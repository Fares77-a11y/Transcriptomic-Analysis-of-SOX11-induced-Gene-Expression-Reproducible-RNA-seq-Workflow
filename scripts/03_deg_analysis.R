#Specifying a working directory
setwd("C:/Users/Ibrah/Downloads/Transcriptomics")
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
#Load Gene count matrix
counts <- read.delim("kallisto/kallisto.merged.gene_counts.tsv", row.names = 1, check.names = FALSE)
head(counts)  # Check structure
# Identify sample columns (all columns except the first column, which is gene names)
sample_columns <- 2:ncol(counts)
#Load metadata
metadata <- read.csv("metadata.csv", row.names = 1)
head(metadata)# Check structure
#Grouping of samples
metadata$Condition <- factor(metadata$Condition, levels = c("Negative", "SOX11-induced"))
# Reorder metadata rows to match the countData column order
metadata <- metadata[colnames(counts)[sample_columns], ]
# Round counts to nearest integer
counts_rounded <- round(counts[, sample_columns])

# Then pass the rounded counts to DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(
  countData = counts_rounded,
  colData = metadata,
  design = ~ Condition
)
# 🔬 Normalize counts (DESeq2 handles this internally)
dds <- DESeq(dds)
# 📊 Extract results: SOX11-induced vs Negative
res <- results(dds, contrast = c("Condition", "SOX11-induced", "Negative"))
# 🏷 Add gene names
res <- as.data.frame(res) %>% arrange(padj)
res$gene_id <- rownames(res)
#Export DESeq2 results
write.csv(res, file = "DESeq2_DEGs_results.csv", row.names = FALSE)
# 📌 Extract significant DEGs (padj < 0.05)
sig_DEGs <- res %>% filter(padj < 0.05)
# 📜 Save significant DEGs separately
write.csv(sig_DEGs, file = "DESeq2_DEGs_significant.csv", row.names = FALSE)
#PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="Condition")
#PCA to Classify samples according to specific gene
plotCounts(dds, gene="A1BG", intgroup="Condition")
#PCA
res <- results(dds)
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup="Condition")


# 🎨 **Volcano Plot**
#UP, Down, and NS
# Prepare your result dataframe
res_df <- as.data.frame(res)

# Classify genes based on thresholds
logFC_cutoff <- 1
pval_cutoff <- 0.05

res_df <- res_df %>%
  mutate(Change = case_when(
    log2FoldChange >= logFC_cutoff & pvalue < pval_cutoff ~ "UP",
    log2FoldChange <= -logFC_cutoff & pvalue < pval_cutoff ~ "DOWN",
    TRUE ~ "NOT"
  ))

# Count the number of genes in each group
gene_counts <- table(res_df$Change)

# Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Change)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NOT" = "grey30")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right") +
  labs(
    title = "Volcano Plot (SOX11-induced vs Negative)",
    x = expression(log[2]~"(FC)"),
    y = expression(-log[10]~"(P-value)"),
    color = "Change",
    caption = paste0(
      "Cutoff for logFC is ", logFC_cutoff, "\n",
      "UP genes: ", gene_counts["UP"], "\n",
      "DOWN genes: ", gene_counts["DOWN"], "\n",
      "Non-significant: ", gene_counts["NOT"]
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
#Significant Gene names
p_cutoff <- 0.05
fc_cutoff <- 1
selectLab <- rownames(res)[which(res$padj < p_cutoff & abs(res$log2FoldChange) > fc_cutoff)]

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                selectLab = rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)],
)
_____________________________________________________________________________
# Select genes with padj < 0.01 and abs(log2FC) > 2
#Gene names
sigGenes <- rownames(res)[which(res$padj < 0.01 & abs(res$log2FoldChange) > 2)]

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 2,
                selectLab = sigGenes,  # <-- Labels shown only for selected significant genes
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                labSize = 4,
                colAlpha = 0.85,
                col = c("grey30", "blue3", "red3", "purple"),
                title = "SOX11-induced vs Negative",
                subtitle = "Volcano plot with gene labels",
                caption = "Source: DESeq2",
                legendPosition = "right"
)

__________________________________________________________________________


# 🔥 **Heatmap (Top 50 DEGs)**
#  Extract top 50 significant genes
top100 <- head(order(res$padj), 100)
mat <- assay(vst(dds))[top100, ]  # Use variance-stabilized data for smoother heatmap
# Essential annotations only
anno <- metadata[, c("Condition", "SampleName")]
rownames(anno) <- rownames(metadata)

# Custom color palette for condition
ann_colors <- list(
  Condition = c(Negative = "#1f77b4", `SOX11-induced` = "#ff7f0e")
)
# Create annotation dataframe (make sure rownames match column names of 'mat')
anno <- metadata[, c("Condition"), drop = FALSE]# Keep only the 'Condition' column for simplicity
rownames(anno) <- rownames(metadata)    # Ensure rownames are SRX IDs
show_rownames = TRUE


# Double-check that column names of 'mat' and rownames of 'anno' align
all(colnames(mat) == rownames(anno))  # Should return TRUE
# Generate the heatmap
pheatmap(mat,
         scale = "row",  # Scale per gene
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,      # ✅ Show gene names
         fontsize_row = 6,          # ✅ Adjust gene label size for clarity
         fontsize_col = 8,          # ✅ Adjust sample label size
         annotation_col = anno,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         main = "Heatmap of Top 100 DEGs (Scaled)")

___________________________________________-
  
#Heatmap with TOP 50 variably expressed genes mapped to library names
  
  # 1. Extract top 50 DEGs
top50 <- head(order(res$padj), 50)
mat <- assay(vst(dds))[top50, ]  # variance-stabilized counts

# 2. Map SRX IDs (colnames in mat) to LibraryName
# Assuming 'metadata' has rownames as SRX IDs and a column named "LibraryName"
colnames(mat) <- metadata[colnames(mat), "LibraryName"]

# 3. Create condition annotation dataframe using SRX IDs
anno <- metadata[, c("Condition"), drop = FALSE]
rownames(anno) <- metadata$Experiment  # this ensures it matches the SRX IDs
anno <- anno[colnames(vst(dds)), , drop = FALSE]  # reorder to match expression data
rownames(anno) <- colnames(mat)  # assign new rownames matching new column names in mat

# 4. Define annotation colors
ann_colors <- list(
  Condition = c(Negative = "#1f77b4", `SOX11-induced` = "#ff7f0e")
)

# 5. Plot the heatmap with LibraryNames as column labels
pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         annotation_col = anno,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         main = "Heatmap of Top 50 DEGs")
