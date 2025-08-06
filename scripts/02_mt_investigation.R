# 02_mt_investigation.R
# Thorough investigation of mitochondrial percentages in PBMC data
# This script helps determine the appropriate MT% filtering threshold
# Author: Thomas Faria
# Date: 2025-07-15

# Load libraries
library(Seurat)
library(tidyverse)
library(patchwork)

# Set up output directory
dir.create("results/figures/mt_investigation", 
           recursive = TRUE, 
           showWarnings = FALSE)

# Load the raw object with QC metrics from script 01
cat("Loading PBMC data from 01_quality_control.R output...\n")
pbmc <- readRDS("data/pbmc_unfiltered.rds")

# Begin analysis
cat("\n=== Mitochondrial Percentage Analysis ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n\n")

# 1. Distribution stats for MT%
cat("MT% Distribution Statistics:\n")
cat("--------------------------------\n")
mt_summary <- summary(pbmc$percent.mt)
print(mt_summary)
cat("\nPercentiles:\n")
mt_percentiles <- quantile(pbmc$percent.mt, 
                           probs = c(0.01, 0.05, 0.10, 0.25, 0.50,
                                     0.75, 0.90, 0.95, 0.99))
print(round(mt_percentiles, 2))

# 2. Threshold impact analysis
cat("\n\n2. Impact of Different MT% Thresholds:\n")
cat("---------------------------------------\n")
thresholds <- c(3, 5, 7.5, 10, 12.5, 15, 20)
threshold_analysis <- data.frame(
  MT_threshold = thresholds,
  cells_retained = NA,
  percent_retained = NA,
  median_genes_retained = NA,
  median_UMI_retained = NA
)

for (i in 1:length(thresholds)) {
  mask <- pbmc$percent.mt < thresholds[i]
  threshold_analysis$cells_retained[i] <- sum(mask)
  threshold_analysis$percent_retained[i] <- round(sum(mask) / ncol(pbmc) * 100, 1)
  threshold_analysis$median_genes_retained[i] <- median(pbmc$nFeature_RNA[mask])
  threshold_analysis$median_UMI_retained[i] <- median(pbmc$nCount_RNA[mask])
}

print(threshold_analysis)
write.csv(threshold_analysis, 
          "results/tables/mt_threshold_analysis.csv", 
          row.names = FALSE)


# 3. Visualizations

# 3a. Histogram with possible thresholds
p_hist <- ggplot(pbmc@meta.data, aes(x = percent.mt)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = c(5, 10, 15), 
             linetype = c("dashed", "dashed", "dotted"), 
             color = c("purple", "darkgreen", "black"),
             size = c(1.2, 1, 0.8)) +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 50)) +
  labs(title = "Distribution of MT%",
       x = "Mitochondrial Percentage (%)",
       y = "Number of Cells") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

ggsave("results/figures/mt_investigation/mt_distribution_histogram.png", 
       p_hist, width = 10, height = 6, dpi = 300)

# 3b. Cumulative distribution
pbmc_meta <- pbmc@meta.data %>%
  arrange(percent.mt) %>%
  mutate(cumulative_percent = row_number() / n() * 100) # Calculate cumulative %

p_cumulative <- ggplot(pbmc_meta, aes(x = percent.mt, y = cumulative_percent)) +
  geom_line(size = 1.2, color = "darkblue") +
  geom_vline(xintercept = c(5, 10, 15), 
             linetype = "dashed", 
             color = c("purple", "darkgreen", "black")) +
  geom_hline(yintercept = c(90, 95), 
             linetype = "dotted", 
             color = "gray50") +
  scale_x_continuous(breaks = seq(0, 50, 5), limits = c(0, 30)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(title = "Cumulative Distribution of MT%",
       subtitle = "Shows percentage of cells retained at each MT% threshold",
       x = "Mitochondrial Percentage (%)",
       y = "Cumulative Percentage of Cells (%)") +
  theme_minimal()

ggsave("results/figures/mt_investigation/mt_cumulative_distribution.png", 
       p_cumulative, width = 10, height = 6, dpi = 300)

# 3c. MT% vs other QC metrics
p_mt_vs_genes <- ggplot(pbmc@meta.data, aes(x = nFeature_RNA, y = percent.mt)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_density_2d(color = "red", alpha = 0.5) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", 
             color = c("purple", "darkgreen")) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(title = "Mitochondrial % vs Gene Count",
       subtitle = "Checking if high MT% cells have other quality issues",
       x = "Number of Genes Detected",
       y = "Mitochondrial Gene Percentage (%)") +
  theme_minimal()

p_mt_vs_umi <- ggplot(pbmc@meta.data, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_density_2d(color = "red", alpha = 0.5) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", 
             color = c("purple", "darkgreen")) +
  scale_y_continuous(limits = c(0, 30)) +
  scale_x_log10() +
  labs(title = "Mitochondrial % vs UMI Count",
       x = "Total UMI Count (log scale)",
       y = "Mitochondrial Gene Percentage (%)") +
  theme_minimal()

# Combine the scatter plots
p_combined_scatter <- p_mt_vs_genes + p_mt_vs_umi
ggsave("results/figures/mt_investigation/mt_vs_other_metrics.png", 
       p_combined_scatter, width = 14, height = 6, dpi = 300)

# 4. Statistical comparison of high vs low MT% cells
cat("\n\n4. Comparison of High vs Low MT% Cells:\n")
cat("----------------------------------------\n")

# Define groups based on different thresholds
comparisons <- list(
  "5%" = 5,
  "10%" = 10,
  "15%" = 15
)

comparison_results <- data.frame()

for (name in names(comparisons)) {
  threshold <- comparisons[[name]]
  
  low_mt <- pbmc@meta.data[pbmc$percent.mt <= threshold, ]
  high_mt <- pbmc@meta.data[pbmc$percent.mt > threshold, ]
  
  result <- data.frame(
    threshold = name,
    low_mt_cells = nrow(low_mt),
    high_mt_cells = nrow(high_mt),
    low_mt_median_genes = median(low_mt$nFeature_RNA),
    high_mt_median_genes = median(high_mt$nFeature_RNA),
    low_mt_median_umi = median(low_mt$nCount_RNA),
    high_mt_median_umi = median(high_mt$nCount_RNA),
    gene_difference_pct = round((median(high_mt$nFeature_RNA) - 
                                   median(low_mt$nFeature_RNA)) / 
                                  median(low_mt$nFeature_RNA) * 100, 1)
  )
  
  comparison_results <- rbind(comparison_results, result)
}

print(comparison_results)
write.csv(comparison_results, "results/tables/mt_high_vs_low_comparison.csv", 
          row.names = FALSE)

# 5. MAD-based thresholding
cat("\n\n5. MAD-based Thresholding:\n")
cat("-------------------------------------\n")

# Calculate MAD-based threshold
mt_median <- median(pbmc$percent.mt)
mt_mad <- mad(pbmc$percent.mt)

cat("Median MT%:", round(mt_median, 2), "\n")
cat("MAD of MT%:", round(mt_mad, 2), "\n")
cat("\nSuggested thresholds based on median + n*MAD:\n")
for (n in c(2, 2.5, 3)) {
  threshold <- mt_median + n * mt_mad
  cells_retained <- sum(pbmc$percent.mt < threshold)
  percent_retained <- round(cells_retained / ncol(pbmc) * 100, 1)
  cat(sprintf("  Median + %.1f*MAD = %.1f%% (retains %d cells, %.1f%%)\n", 
              n, threshold, cells_retained, percent_retained))
}
