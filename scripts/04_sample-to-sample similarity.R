# 1. Set working directory
setwd("C:/Users/Ibrah/Downloads/Transcriptomics")

# 2. Load required libraries
library(DESeq2)
library(readxl)
library(pheatmap)
library(dplyr)
library(ggplot2)

# 3. Read metadata
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Extract only the necessary columns
metadata <- metadata %>%
  dplyr::select(LibraryName = 2, everything())  # Ensure LibraryName is accessible

# 4. Read gene count data
counts <- read_excel("sample names.xlsx")

# Set gene_id as row names and drop gene_name column
rownames(counts) <- counts$gene_id
counts <- counts[, -c(1,2)]  # Remove gene_id and gene_name columns

# 5. Subset columns based on metadata$LibraryName
# Ensure column names in counts exactly match metadata$LibraryName
matched_samples <- metadata$LibraryName[metadata$LibraryName %in% colnames(counts)]
counts <- counts[, matched_samples]
metadata <- metadata[metadata$LibraryName %in% matched_samples, ]

# Reorder count columns to match metadata order
counts <- counts[, match(metadata$LibraryName, colnames(counts))]

# 6. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Condition)

# 7. Variance Stabilizing Transformation
vsd <- vst(dds, blind = TRUE)

# 8. Sample Distance Matrix
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- metadata$LibraryName
colnames(sample_dist_matrix) <- metadata$LibraryName

# 9. Sample-to-Sample Distance Heatmap
# ✅ Correct call with closing parenthesis
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample-to-Sample Distance Heatmap")

# 10. PCA Plot
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = Condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.1, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()
___________________________________________________________________

