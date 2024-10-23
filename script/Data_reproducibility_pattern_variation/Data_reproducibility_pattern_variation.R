# Set the working directory (if needed)
setwd("D:/Biostate")

# Load necessary libraries
library(pheatmap)
library(factoextra)
library(ggplot2)

# Read the gene count data
data <- read.table("gene_counts.txt", header = TRUE, sep = "\t", row.names = 1)

# View the first few rows of the data
head(data)

# Remove rows where all counts are zero
data_filtered <- data[rowSums(data) > 0, ]

# Log2 transformation (if needed)
data_log <- log2(data_filtered + 1)

# Scatter plot matrix to check reproducibility
pairs(data_log, main = "Scatterplot Matrix for Reproducibility")

# Save scatter plot matrix with higher resolution (300 dpi)
png("scatterplot_matrix_high_res.png", width = 3000, height = 3000, res = 300)
pairs(data_log, main = "Scatterplot Matrix for Reproducibility")
dev.off()

# Compute correlation matrix
cor_matrix <- cor(data_log)

# Generate heatmap
pheatmap(cor_matrix, 
         main = "Correlation Heatmap",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         display_numbers = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100))

# Save heatmap with higher resolution (300 dpi)
png("correlation_heatmap_high_res.png", width = 3000, height = 2400, res = 300)
pheatmap(cor_matrix, 
         main = "Correlation Heatmap",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         display_numbers = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Perform PCA on the log-transformed data (transpose the matrix)
pca_result <- prcomp(t(data_log), scale. = TRUE)

# Visualize PCA result
fviz_pca_ind(pca_result, 
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, 
             title = "PCA of Samples")

# Save PCA plot with higher resolution
ggsave("PCA_plot_high_res.png", plot = last_plot(), width = 10, height = 8, dpi = 300)

# Scree plot to show the explained variance
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# Save Scree plot with higher resolution
ggsave("scree_plot_high_res.png", width = 10, height = 8, dpi = 300)
