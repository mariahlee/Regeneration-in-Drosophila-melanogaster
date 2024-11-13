#!/bin/bash

# Rna-Seq Pipeline

## Time Tracking
SECONDS=0

## Function to print the error and duration
handle_error() {
  echo "An error occurred. Exiting."
  duration=$SECONDS
  hours=$((duration / 3600))
  minutes=$(( (duration % 3600) / 60))
  seconds=$((duration % 60))
  echo "Script completed in $hours hours, $minutes minutes, and $seconds seconds."
  notify-send "Script canceled" "Your script had an error."
  paplay /usr/share/sounds/freedesktop/stereo/complete.oga
  exit 1
}

## Set the trap to catch errors
trap 'handle_error' ERR

## Assign working directory
cd ~/Regeneration_in_Dmelanogaster/

## Define directories
data_dir="./reads/rawdata/"               # where original fastq files are located 
firstqc_dir="./qc_reports/raw/"           # for outputting QC reports of original reads
trimmed_dir="./reads/trimmed_reads/"      # where trimmed reads will be saved
secondqc_dir="./qc_reports/trimmed/"      # for outputting QC reports of trimmed reads
hisat2_index_dir="./genome/hisat2_index/" # directory where Hisat2 index files are stored
hisat2_index_prefix="${hisat2_index_dir}dmel-all-chromosome-r6.57"
aligned_dir="./reads/aligned_bam/"        # where aligned BAM files will be saved
read_counts_dir="./reads/read_counts/"    # where read count files will be saved
multiqc_dir="./qc_reports/multiqc/"       # for compiled multiqc reports
metrics_dir="./qc_reports/metrics/"       # for other metrics

## Ensure the necessary directories exist
mkdir -p $firstqc_dir $secondqc_dir $trimmed_dir $hisat2_index_dir \
         $aligned_dir $read_counts_dir $multiqc_dir $metrics_dir

## Custom adapter sequences
ADAPTER_FWD="AGATCGGAAGAG"
ADAPTER_REV="GATCGTCGGACT"

## Input genome and annotation file paths
genome_file="./genome/dmel-all-chromosome-r6.57.fasta"
annotation_file="./genome/annotation/dmel-all-r6.57.gtf.gz"

## Load environment
echo "Activating drosphila environment..."
source /home/mariah/anaconda3/etc/profile.d/conda.sh
conda activate drosophila                 # An exiting environment with needed software

## Check if genome is indexed, and index if not
if [ ! -f "${hisat2_index_prefix}.1.ht2" ]; then
    echo "Indexing the reference genome..."
    hisat2-build "$genome_file" \
    $hisat2_index_prefix
fi

## Function to calculate mean read depth
calculate_mean_depth() {
  BAM_FILE=$1
  DEPTH_FILE="${metrics_dir}$(basename ${BAM_FILE%.bam}_depth.txt)"
  MEAN_DEPTH_FILE="${metrics_dir}$(basename ${BAM_FILE%.bam}_mean_depth.txt)"

  # Calculate depth at each position
  samtools depth $BAM_FILE > $DEPTH_FILE

  # Calculate mean depth
  awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print 0}' \
  $DEPTH_FILE > $MEAN_DEPTH_FILE

  echo "Mean read depth for $BAM_FILE: $(cat $MEAN_DEPTH_FILE)"
}

## Function to generate coverage histogram and table per chromosome
generate_coverage_stats() {
  BAM_FILE=$1
  COVERAGE_FILE="${metrics_dir}$(basename ${BAM_FILE%.bam}_coverage.txt)"
  COVERAGE_HIST_FILE="${metrics_dir}$(basename ${BAM_FILE%.bam}_coverage_histogram.txt)"

  # Generate coverage per chromosome
  samtools coverage $BAM_FILE > $COVERAGE_FILE

  # Generate histogram of coverage
  samtools coverage -H $BAM_FILE > $COVERAGE_HIST_FILE

  echo "Coverage stats for $BAM_FILE generated."
}

## Process each pair of R1 and R2 files
for i in {1..8}; do
  R1_FILE="${data_dir}RH${i}_S${i}_L007_R1_001.fastq.gz"
  R2_FILE="${data_dir}RH${i}_S${i}_L007_R2_001.fastq.gz"
  TRIMMED_R1_FILE="${trimmed_dir}RH${i}_S${i}_L007_R1_001_val_1.fq.gz"
  TRIMMED_R2_FILE="${trimmed_dir}RH${i}_S${i}_L007_R2_001_val_2.fq.gz"
  BAM_FILE="${aligned_dir}RH${i}_S${i}_aligned.bam"
  SORTED_BAM_FILE="${aligned_dir}RH${i}_S${i}_aligned_sorted.bam"
  COUNT_FILE="${read_counts_dir}RH${i}_S${i}_counts.tsv"
  METRICS_FILE="${metrics_dir}RH${i}_S${i}_metrics.txt"

  echo "Processing pair: $R1_FILE and $R2_FILE"
  if [[ ! -f "$R1_FILE" ]] || [[ ! -f "$R2_FILE" ]]; then
    echo "Missing file: $R1_FILE or $R2_FILE"
    exit 1
  fi

  # Generate FastQC report for the original reads and save it in $firstqc_dir
   echo "Running FastQC for sample RH${i}..."
   fastqc -o $firstqc_dir $R1_FILE $R2_FILE

  # Run Trim Galore to trim adapters and handle FastQC reports
    echo "Trimming reads using Trim Galore and generating final FastQC reports..."
    trim_galore --paired $R1_FILE $R2_FILE \
                --adapter $ADAPTER_FWD --adapter2 $ADAPTER_REV \
                --fastqc --fastqc_args "-o $secondqc_dir" \
                --trim-n --gzip \
                --output_dir $trimmed_dir

  # Align trimmed reads to the reference genome and pipe to Samtools to get BAM
    echo "Running HISAT2 for sample RH${i}..."
    hisat2 -x $hisat2_index_prefix \
       -1 $TRIMMED_R1_FILE \
       -2 $TRIMMED_R2_FILE \
       -p 8 2> ${aligned_dir}RH${i}_S${i}_alignment_summary.txt | \
    samtools view -bS - > $BAM_FILE
  
  # Generate read counts using featureCounts
    echo "Running featureCounts for sample RH${i}..."
    featureCounts -a $annotation_file \
                  -o "${COUNT_FILE%.tsv}.tsv" \
                  -T 8 \
                  -g gene_id \
                  -p \
                  --primary \
                  -s 0 \
                $BAM_FILE

  # Sort BAM file to ensure it is position sorted
  echo "Sorting BAM file for sample RH${i}..."
  samtools sort -o $SORTED_BAM_FILE $BAM_FILE
  
  # Calculate mean read depth and generate coverage stats using the sorted BAM file
  echo "Calculating mean read depth for sample RH${i}..."
  calculate_mean_depth $SORTED_BAM_FILE
  echo "Generating coverage stats for sample RH${i}..."
  generate_coverage_stats $SORTED_BAM_FILE
done

## Run multiqc for all steps
multiqc $firstqc_dir -o $multiqc_dir/raw_multiqc_report
multiqc $secondqc_dir -o $multiqc_dir/trimmed_multiqc_report
multiqc $aligned_dir -o $aligned_dir/multiqc_report
multiqc $read_counts_dir -o $multiqc_dir/counts_multiqc_report
multiqc $metrics_dir -o $multiqc_dir/metrics_multiqc_report

## Save versions of tools used in the pipeline
echo "Saving versions of tools..."
declare -A tools
tools=(
  ["conda"]="conda --version"
  ["python"]="python --version"
  ["fastqc"]="fastqc --version"
  ["awk"]="awk -Wversion 2>/dev/null || awk --version"  
  ["trim_galore"]="trim_galore --version | grep version"
  ["hisat2"]="hisat2 --version"
  ["samtools"]="samtools --version"
  ["featureCounts"]="featureCounts -v 2>&1 | awk 'NR==2 {print $1, $2}'"
)

output_file="${metrics_dir}regeneration_pipeline_versions.tsv"

for tool in "${!tools[@]}"; do
  version_info=$(${tools[$tool]} | head -n 1)
  echo -e "${tool}\t${version_info}"
done > $output_file

## Keep track of how long the script takes to run
duration=$SECONDS
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60))
seconds=$((duration % 60))
echo "Script completed in $hours hours, $minutes minutes, and $seconds seconds."

## Add notification and sound at the end of the script
notify-send "Script finished" "Your script has completed running."
paplay /usr/share/sounds/freedesktop/stereo/complete.oga
