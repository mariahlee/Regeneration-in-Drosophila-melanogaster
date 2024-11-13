# Regeneration in *Drosophila melanogaster*
This repository contains the code and resources for RNA-Seq data processing, differential gene expression (DGE), and gene ontology (GO) analysis in Drosophila melanogaster, specifically focused on genes involved in regeneration.

## RNA-Seq Pipeline

### Overview
The `rnaseq_pipeline.sh` script automates the RNA-Seq data processing workflow. It includes quality control, adapter trimming, alignment, read counting, and report generation, with runtime logging and tool version tracking.

### Requirements
**Environment**: Anaconda  
**Tools**:  
FastQC == v0.12.1  
Trim Galore == 0.6.10  
HISAT2 == 2.2.1  
samtools == 1.20  
featureCounts == v2.0.6  
multiqc == v1.22.3

## Differential Gene Expression and Gene Ontology Analysis

### Overview
This part of the project performs DGE and GO enrichment analysis on the gene expression data. It includes comparisons between early and late stages, as well as damaged and undamaged conditions. Analysis is centered on the JAK/STAT and JNK pathways, with various visualizations to explore expression patterns.

### Versions
**R**: 4.4.1  
**Packages**: Listed in `DGE_GO_session_info.txt`

### Usage
Counts from all samples are combined into a single csv by calling the `combining_counts.R` script.
