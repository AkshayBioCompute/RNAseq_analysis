# Load necessary libraries
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(reshape2)

# Set the working directory (change this to your directory)
setwd("D:/Biostate")

# Load the gene counts data
count_data <- read.table("gene_counts.txt", header = TRUE, row.names = 1)

# Remove rows with all zero counts and NaN values
count_data <- count_data[rowSums(count_data) > 0, ]
count_data <- count_data[!apply(count_data, 1, function(x) any(is.nan(x))), ]

# Define sample conditions (customize based on your experiment)
col_data <- data.frame(
  tissue = c(rep("Heart", 4), rep("Liver", 4)),
  time = c(rep("ZT0", 2), rep("ZT12", 2), rep("ZT0", 2), rep("ZT12", 2)),
  row.names = colnames(count_data)
)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ tissue + time + tissue:time)

# Normalization of counts using size factors (default in DESeq2)
dds <- estimateSizeFactors(dds)

# Transform normalized counts for visualization (optional)
normalized_counts <- counts(dds, normalized = TRUE)

# Run the differential expression analysis
dds <- DESeq(dds)

# Extract results for various comparisons
res_tissue <- results(dds, name = "tissue_Liver_vs_Heart")
res_time <- results(dds, name = "time_ZT12_vs_ZT0")
res_interaction_Liver_vs_Heart_ZT12 <- results(dds, name = "tissueLiver.timeZT12")
res_Liver_vs_Heart_ZT0 <- results(dds, contrast = c("tissue", "Liver", "Heart"))
res_interaction_Liver_ZT12_vs_ZT0 <- results(dds, contrast = c("time", "ZT12", "ZT0"))
res_interaction_Heart_ZT12_vs_ZT0 <- results(dds, contrast = c("time", "ZT12", "ZT0"))

# Save results to CSV files
write.csv(as.data.frame(res_tissue), file = "All_Tissue_Specific_DEGs.csv")
write.csv(as.data.frame(res_time), file = "All_Time_Specific_DEGs.csv")
write.csv(as.data.frame(res_interaction_Liver_vs_Heart_ZT12), file = "Interaction_Liver_vs_Heart_ZT12.csv")
write.csv(as.data.frame(res_Liver_vs_Heart_ZT0), file = "Interaction_Liver_vs_Heart_ZT0.csv")
write.csv(as.data.frame(res_interaction_Liver_ZT12_vs_ZT0), file = "Interaction_Liver_ZT12_vs_ZT0.csv")
write.csv(as.data.frame(res_interaction_Heart_ZT12_vs_ZT0), file = "Interaction_Heart_ZT12_vs_ZT0.csv")

# Function to get top genes
get_top_genes <- function(results, n = 10) {
  top_genes <- results[order(results$padj), ]
  top_genes <- top_genes[!is.na(top_genes$padj), ]
  return(head(top_genes[order(-top_genes$log2FoldChange), ], n))
}

# Get top genes for different comparisons
top_genes_tissue <- get_top_genes(res_tissue, 10)
top_genes_time <- get_top_genes(res_time, 10)
top_genes_interaction_Liver_vs_Heart_ZT12 <- get_top_genes(res_interaction_Liver_vs_Heart_ZT12, 10)
top_genes_Liver_vs_Heart_ZT0 <- get_top_genes(res_Liver_vs_Heart_ZT0, 10)
top_genes_interaction_Liver_ZT12_vs_ZT0 <- get_top_genes(res_interaction_Liver_ZT12_vs_ZT0, 10)
top_genes_interaction_Heart_ZT12_vs_ZT0 <- get_top_genes(res_interaction_Heart_ZT12_vs_ZT0, 10)

# Save top genes to CSV files
write.csv(as.data.frame(top_genes_tissue), file = "Top_Tissue_Specific_DEGs.csv")
write.csv(as.data.frame(top_genes_time), file = "Top_Time_Specific_DEGs.csv")
write.csv(as.data.frame(top_genes_interaction_Liver_vs_Heart_ZT12), file = "Top_Interaction_Liver_vs_Heart_ZT12.csv")
write.csv(as.data.frame(top_genes_Liver_vs_Heart_ZT0), file = "Top_Interaction_Liver_vs_Heart_ZT0.csv")
write.csv(as.data.frame(top_genes_interaction_Liver_ZT12_vs_ZT0), file = "Top_Interaction_Liver_ZT12_vs_ZT0.csv")
write.csv(as.data.frame(top_genes_interaction_Heart_ZT12_vs_ZT0), file = "Top_Interaction_Heart_ZT12_vs_ZT0.csv")

# Prepare heatmap data
heatmap_data_tissue <- assay(dds)[rownames(top_genes_tissue), ]
heatmap_data_time <- assay(dds)[rownames(top_genes_time), ]
heatmap_data_interaction_Liver_vs_Heart_ZT12 <- assay(dds)[rownames(top_genes_interaction_Liver_vs_Heart_ZT12), ]
heatmap_data_Liver_vs_Heart_ZT0 <- assay(dds)[rownames(top_genes_Liver_vs_Heart_ZT0), ]
heatmap_data_interaction_Liver_ZT12_vs_ZT0 <- assay(dds)[rownames(top_genes_interaction_Liver_ZT12_vs_ZT0), ]
heatmap_data_interaction_Heart_ZT12_vs_ZT0 <- assay(dds)[rownames(top_genes_interaction_Heart_ZT12_vs_ZT0), ]

# Volcano Plots
png("Volcano_Plot_Tissue_Specific_DEGs.png", width = 1600, height = 1200, res = 300)
EnhancedVolcano(res_tissue,
                lab = rownames(res_tissue),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Tissue Specific DEGs",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()

png("Volcano_Plot_Time_Specific_DEGs.png", width = 1600, height = 1200, res = 300)
EnhancedVolcano(res_time,
                lab = rownames(res_time),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Time Specific DEGs",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()

png("Volcano_Plot_Interaction_Effect_Liver_ZT12_vs_Heart_ZT12.png", width = 1600, height = 1200, res = 300)
EnhancedVolcano(res_interaction_Liver_vs_Heart_ZT12,
                lab = rownames(res_interaction_Liver_vs_Heart_ZT12),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Interaction: Liver ZT12 vs Heart ZT12",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()

png("Volcano_Plot_Interaction_Effect_Liver_ZT0_vs_Heart_ZT0.png", width = 1600, height = 1200, res = 300)
EnhancedVolcano(res_Liver_vs_Heart_ZT0,
                lab = rownames(res_Liver_vs_Heart_ZT0),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Interaction: Liver vs Heart ZT0",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()

png("Volcano_Plot_Interaction_Effect_Liver_ZT12_vs_ZT0.png", width = 1600, height = 1200, res = 300)
EnhancedVolcano(res_interaction_Liver_ZT12_vs_ZT0,
                lab = rownames(res_interaction_Liver_ZT12_vs_ZT0),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Interaction: Liver ZT12 vs ZT0",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()

png("Volcano_Plot_Interaction_Effect_Heart_ZT12_vs_ZT0.png", width = 1600, height = 1200, res = 300)
EnhancedVolcano(res_interaction_Heart_ZT12_vs_ZT0,
                lab = rownames(res_interaction_Heart_ZT12_vs_ZT0),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Interaction: Heart ZT12 vs ZT0",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()

# Heatmaps for Top DEGs
png("Heatmap_Top_Tissue_Specific_DEGs.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_tissue, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for Top Tissue-Specific DEGs")
dev.off()

png("Heatmap_Top_Time_Specific_DEGs.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_time, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for Top Time-Specific DEGs")
dev.off()

png("Heatmap_Top_Interaction_Effect_Liver_vs_Heart_ZT12.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_interaction_Liver_vs_Heart_ZT12, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for Top Interaction Effect: Liver vs Heart ZT12")
dev.off()

png("Heatmap_Top_Interaction_Effect_Liver_vs_Heart_ZT0.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_Liver_vs_Heart_ZT0, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for Interaction Effect: Liver vs Heart ZT0")
dev.off()

png("Heatmap_Top_Interaction_Effect_Liver_ZT12_vs_ZT0.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_interaction_Liver_ZT12_vs_ZT0, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for Interaction Effect: Liver ZT12 vs ZT0")
dev.off()

png("Heatmap_Top_Interaction_Effect_Heart_ZT12_vs_ZT0.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_interaction_Heart_ZT12_vs_ZT0, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for Interaction Effect: Heart ZT12 vs ZT0")
dev.off()

# Heatmap for All DEGs
heatmap_data_all_DEGs <- assay(dds)

# Save heatmap data for all DEGs
write.csv(heatmap_data_all_DEGs, file = "All_DEGs_Data.csv")

# Clustering analysis heatmap for all DEGs
png("Heatmap_All_DEGs.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_data_all_DEGs, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Clustering Analysis for All DEGs")
dev.off()

