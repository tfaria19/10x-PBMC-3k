# 00_download_data.R
# Download and prepare unsorted PBMC 3k (v2.0.0) dataset
# Author: Thomas Faria
# Date: 2025-06-26

# Load required libraries
library(Seurat)
library(tidyverse)

# Create data directory if it doesn't exist
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# Option 1: Download unsorted PBMC 3k dataset (TAR.GZ format)
if (!file.exists("data/raw/pbmc_unsorted_3k_filtered_feature_bc_matrix.tar.gz")) {
  cat("Downloading UNSORTED PBMC 3k dataset from 10x Genomics...\n")
  
  download.file(
    url = "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_filtered_feature_bc_matrix.tar.gz",
    destfile = "data/raw/pbmc_unsorted_3k_filtered_feature_bc_matrix.tar.gz",
    method = "curl",
    mode = "wb"
  )
  
  # Extract the files
  cat("Extracting files...\n")
  untar("data/raw/pbmc_unsorted_3k_filtered_feature_bc_matrix.tar.gz", 
        exdir = "data/raw/")
  
  cat("Download complete!\n")
} else {
  cat("Data already exists.\n")
}

# List extracted files
cat("\nExtracted files:\n")
list.files("data/raw/", recursive = TRUE)

