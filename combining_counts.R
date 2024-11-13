#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)

# Define the working directory
working_directory <- "/home/mariah/Regeneration_in_Dmelanogaster" 
# there is no / at the end because of the use of file.path
setwd(working_directory)

# Define the directories for differential gene expression and read counts
read_counts_directory <- file.path(working_directory, "reads/read_counts")

# Verify the directory exists
if (!dir.exists(read_counts_directory)) {
  stop("Error: read_counts_directory does not exist. 
    Please check the path: ", read_counts_directory)
}

# List the count files
count_files <- list.files(path = read_counts_directory, 
                          pattern = "RH\\d+_S\\d+_counts\\.tsv$", 
                          full.names = TRUE)

# Check if files are being detected
if (length(count_files) == 0) {
  stop("No files matched the pattern in: ", read_counts_directory)
}

# Initialize a data frame for the combined counts
combined_counts <- NULL

# Read each count file and merge directly into the combined_counts data frame
for (file in count_files) {

  counts <- fread(file, skip = "# Program")  # Skips metadata
  
  # Extract only the Geneid and the count data column (7th column)
  gene_id <- counts[, .(Geneid)]
  sample_name <- gsub(".*\\/|_counts\\.tsv$", "", file)  # Sample name from file name
  count_data <- counts[, .SD, .SDcols = 7]
  
  # Rename the count column to the sample name
  setnames(count_data, old = names(count_data), new = sample_name)
  
  # Combine Geneid with the sample count
  temp_df <- cbind(gene_id, count_data)
  
  # Print for debugging
  print(paste("Contents of temp_df from file:", sample_name))
  print(head(temp_df))
  
  # Ensure the column name for gene_id remains "Geneid"
  setnames(temp_df, old = "Geneid.Geneid", new = "Geneid", skip_absent = TRUE)
  
  # Print for debugging
  print(paste("Contents of temp_df from file:", sample_name))
  print(head(temp_df))
  
  # Merge the data frame into the combined_counts
  if (is.null(combined_counts)) {
    combined_counts <- temp_df  # Initialize with the first file
  } else {
    combined_counts <- merge(combined_counts, temp_df, by = "Geneid", all = TRUE)
  }
}

# Check if the length of gene_ids matches the number of rows in combined_counts
gene_ids <- combined_counts$Geneid
if (length(gene_ids) != nrow(combined_counts)) {
  stop("Length of gene_ids does not match the number of rows in combined_counts.")
}

# Assign the Gene IDs as row names
rownames(combined_counts) <- gene_ids
combined_counts$Geneid <- NULL  # Remove the Geneid column

# Save the combined counts to a new file
write.table(combined_counts, 
            file.path(read_counts_directory, "RH_combined_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.csv(combined_counts, file.path(read_counts_directory, "RH_combined_counts.csv"))

