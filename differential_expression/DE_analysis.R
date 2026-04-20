#!/usr/bin/env Rscript
# VDA-miRNA Differential Expression Framework
# Author: Vidit Zainith / Collaborative Team
# ------------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(apeglm)
library(pheatmap)

# 1. Pipeline Inputs
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript DE_analysis.R <counts.csv> <metadata.txt> <output_prefix>")
}

counts_file <- args[1]
meta_file <- args[2]
out_prefix <- args[3]

# 2. Data Preparation
counts <- read.csv(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
coldata <- read.table(meta_file, header = TRUE, row.names = 1)

# Quality Check: Column-Row ordering
counts <- counts[, rownames(coldata)]
stopifnot(all(colnames(counts) == rownames(coldata)))

# 3. RUN DESeq2
# We use a design focused on Vascular Dementia vs Control
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

# 4. Results & LFC Shrinkage (apeglm)
# Essential for smaller miRNA datasets to stabilize low-count variance
res_name <- resultsNames(dds)[2]
res <- lfcShrink(dds, coef = res_name, type = "apeglm")

# Convert to dataframe and add significance tags
res_df <- as.data.frame(res)
res_df$ID <- rownames(res_df)
res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              ifelse(res_df$log2FoldChange > 1, "UP", "DOWN"), "NS")

# 5. Output
write.csv(res_df, paste0(out_prefix, "_full_results.csv"), row.names = FALSE)
write.csv(subset(res_df, Significance != "NS"), paste0(out_prefix, "_significant.csv"), row.names = FALSE)

# 6. Professional Visualization: Volcano Plot
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("UP" = "#d73027", "DOWN" = "#4575b4", "NS" = "grey80")) +
  theme_minimal() +
  labs(title = paste("Differential Expression:", res_name),
       subtitle = "Processed using LFC Shrinkage (apeglm)",
       x = "Log2 Fold Change", y = "-Log10 Padj")

ggsave(paste0(out_prefix, "_volcano.png"), plot = p, width = 7, height = 5)

message(paste("Differential Expression Analysis Completed. Results saved with prefix:", out_prefix))
