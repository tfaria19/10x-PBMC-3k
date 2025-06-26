# 10x-PBMC-3k

-   First time making a single-cell RNA-seq analysis pipeline using Seurat

## Single-Cell RNA-seq Analysis Portfolio

### Overview

This repository demonstrates single-cell RNA-seq analysis using R and Seurat. The pipeline walks through quality control, clustering, and cell type annotation of a 10x Genomics tutorial dataset for peripheral blood mononuclear cells (PBMCs).

### Data ([link](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-3-k-1-standard-2-0-0))

-   **Source**: 10x Genomics PBMC 3k dataset
-   **Estimated number of cells**: 3,009
-   **Platform**: 10x Genomics Chromium Controller

### Analysis Pipeline

1.  Quality control and filtering
2.  Normalization and feature selection
3.  Dimensionality reduction (PCA, UMAP)
4.  Clustering and cell type annotation
5.  Differential expression analysis

### Repository Structure

├── data/ \# Raw and processed data ├── scripts/ \# Analysis scripts ├── results/ \# Output figures and tables └── docs/ \# Additional documentation

### Tools

-   R 4.x

-   Seurat 5.0

-   tidyverse

-   Additional packages in renv.lock

### Reproducibility

This project uses `renv` for package management. To reproduce:

``` r
install.packages("renv")
renv::restore()
```

### Contact

Thomas Faria - fariathomas17\@gmail.com
