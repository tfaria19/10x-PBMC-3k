# 01_quality_control.R
# Quality control for unsorted PBMC 3k (v2.0.0) dataset
# Author: Thomas Faria
# Date: 2025-06-26

# Load libraries
library(Seurat)
library(tidyverse)
library(patchwork)

# Read the 10X data (matrix format)
# Read10X returns a list with "Gene Expression" and "Peaks"
pbmc.data <- Read10X(data.dir = "data/raw/filtered_feature_bc_matrix/")

# Extract just the Gene Expression data
if (is.list(pbmc.data)) {
  pbmc.data <- pbmc.data[["Gene Expression"]]
  cat("Extracted Gene Expression data from multiome dataset\n")
}

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3,    # gene must be in at least 3 cells
                           min.features = 200) # cell must have at least 200 genes

# Initial look at the data
cat(paste("Initial dataset:", ncol(pbmc), "cells and", nrow(pbmc), "genes\n"))

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)

p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red")

p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_vline(xintercept = c(500, 5000), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(250, 2500), linetype = "dashed", color = "red")

# Save QC plots
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/01_qc_violin_plots.png", plot = p1, width = 12, height = 6)
ggsave("results/figures/02_qc_scatter_mt.png", plot = p2, width = 8, height = 6)
ggsave("results/figures/03_qc_scatter_features.png", plot = p3, width = 8, height = 6)

# Print QC summary statistics
qc_stats <- pbmc@meta.data %>%
  summarise(
    median_genes = median(nFeature_RNA),
    median_counts = median(nCount_RNA),
    median_mt = median(percent.mt),
    cells_high_mt = sum(percent.mt > 5),
    cells_low_genes = sum(nFeature_RNA < 200)
  )

cat("\nQC Summary Statistics:\n")
print(qc_stats)

# Save the unfiltered object for comparison
dir.create("data", showWarnings = FALSE)
saveRDS(pbmc, file = "data/pbmc_unfiltered.rds")

cat("\nAnalysis complete! Check results/figures/ for QC plots.\n")
