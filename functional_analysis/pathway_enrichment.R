#!/usr/bin/env Rscript
# Functional Enrichment for VDA-miRNA Targets
# ------------------------------------------------------------------

library(clusterProfiler)
library(ReactomePA)
library(org.Rn.eg.db) # Using Rat database as per original project
library(enrichplot)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pathway_enrichment.R <gene_list.txt> <output_prefix>")
}

# 1. Load Gene List
gene_symbols <- read.table(args[1], header = FALSE)$V1

# 2. Convert to Entrez IDs (mapping from SYMBOL)
mapping <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db")
entrez_ids <- mapping$ENTREZID

# 3. GO Enrichment - Biological Process
go_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Rn.eg.db, ont = "BP", readable = TRUE)
p1 <- dotplot(go_bp, showCategory = 15) + ggtitle("GO: Biological Process Enrichment")
ggsave(paste0(args[2], "_GO_BP.png"), plot = p1, width = 8, height = 10)

# 4. KEGG Enrichment
kegg <- enrichKEGG(gene = entrez_ids, organism = 'rno')
kegg <- setReadable(kegg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
p2 <- dotplot(kegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment")
ggsave(paste0(args[2], "_KEGG.png"), plot = p2, width = 8, height = 10)

# 5. Reactome Analysis
reactome <- enrichPathway(gene = entrez_ids, organism = "rat", readable = TRUE)
p3 <- dotplot(reactome, showCategory = 15) + ggtitle("Reactome Pathway Enrichment")
ggsave(paste0(args[2], "_Reactome.png"), plot = p3, width = 8, height = 10)

message("Functional Analysis Completed Successfully.")
