# Set working directory
setwd("D:/Biostate")

# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # For mouse gene annotation
library(enrichplot)
library(ggplot2)
library(dplyr)

# Get a list of all CSV files in the directory
csv_files <- list.files(pattern = "*.csv")

# Loop through each CSV file
for (file in csv_files) {
  # Read the DEG data
  deg_data <- read.csv(file, header = TRUE, row.names = 1)
  
  # Remove rows with any NA values
  deg_data_clean <- na.omit(deg_data)
  
  # Prepare Gene List
  gene_list <- rownames(deg_data_clean)
  
  # Remove version numbers from Ensembl IDs if necessary
  gene_list <- sub("\\..*", "", gene_list)  # This will remove everything after the dot
  
  # Convert Ensembl Gene IDs to Entrez IDs
  gene_ids <- bitr(gene_list, fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
  
  # Perform GO Enrichment Analysis
  go_results <- enrichGO(gene         = gene_ids$ENTREZID,
                         OrgDb        = org.Mm.eg.db,
                         keyType      = "ENTREZID",
                         ont          = "BP",  # Change to "CC" or "MF" if needed
                         pAdjustMethod = "BH",
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  
  # Perform KEGG Enrichment Analysis
  kegg_results <- enrichKEGG(gene         = gene_ids$ENTREZID,
                             organism     = 'mmu',  # Mouse organism code
                             pAdjustMethod = "BH",
                             qvalueCutoff  = 0.05)
  
  # Visualize Results
  # Bubble Plot for GO Terms
  go_plot <- dotplot(go_results, showCategory=10) + ggtitle(paste("Top 10 GO Enriched Terms for", file))
  print(go_plot)
  ggsave(paste0("Top_10_GO_Enriched_Terms_", tools::file_path_sans_ext(file), ".png"), plot = go_plot, width = 8, height = 6)
  
  # Bar Plot for KEGG Pathways
  kegg_plot <- barplot(kegg_results, showCategory=10) + ggtitle(paste("Top 10 KEGG Enriched Pathways for", file))
  print(kegg_plot)
  ggsave(paste0("Top_10_KEGG_Enriched_Pathways_", tools::file_path_sans_ext(file), ".png"), plot = kegg_plot, width = 8, height = 6)
  
  # Summarize Results
  go_summary <- as.data.frame(go_results)
  kegg_summary <- as.data.frame(kegg_results)
  
  # Write summaries to CSV
  write.csv(go_summary, paste0("Mouse_GO_Enrichment_Summary_", tools::file_path_sans_ext(file), ".csv"), row.names = FALSE)
  write.csv(kegg_summary, paste0("Mouse_KEGG_Enrichment_Summary_", tools::file_path_sans_ext(file), ".csv"), row.names = FALSE)
}
