# ======================================================================================
# Title:  Differential Gene Expression Analysis Pipeline
# Description: RNA-Seq data analysis from TCGA-BRCA, identifying differentially 
#              expressed genes (DEGs) between Luminal A subtype and Normal tissue,
#              with a specific focus on the ITGB1 target gene.

# ======================================================================================

# ======================================================================================
# Phase 0: Environment Setup & Package Management
# ======================================================================================

# Define required packages
required_packages <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
  "BiocParallel", "apeglm", "AnnotationDbi", "org.Hs.eg.db",
  "EnhancedVolcano", "ggplot2", "pheatmap"
)

# Install missing packages via BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
} 
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages, update = FALSE, ask = FALSE)

# Load libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(BiocParallel)
  library(apeglm)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(EnhancedVolcano)
  library(ggplot2)
  library(pheatmap)
})

# Create an output directory for saving results and plots
dir.create("results", showWarnings = FALSE)

# ======================================================================================
# Phase 1: Data Acquisition & Preprocessing (TCGA-BRCA)
# ======================================================================================

# 1. Query GDC for BRCA Transcriptome Profiling
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 2. Identify Luminal A subtype patients
subtypes <- TCGAquery_subtype(tumor = "BRCA")
lumA_patients <- subtypes$patient[subtypes$BRCA_Subtype_PAM50 == "LumA"]

# 3. Filter query for Normal Tissue and Luminal A Primary Tumors
all_cases <- query$results[[1]]$cases
is_normal <- query$results[[1]]$sample_type == "Solid Tissue Normal"
is_lumA <- query$results[[1]]$sample_type == "Primary Tumor" & substr(all_cases, 1, 12) %in% lumA_patients

query$results[[1]] <- query$results[[1]][is_normal | is_lumA, ]

# 4. Download and prepare SummarizedExperiment (SE) object
options(timeout = 300) # Increase timeout for large downloads
se <- GDCprepare(query)

# 5. Define conditions for DESeq2 and save raw SE object
se$condition <- factor(se$shortLetterCode)
saveRDS(se, "results/se_BRCA_LumA_Normal_Filtered.rds")

# Clean up memory
rm(query, subtypes, lumA_patients)
gc()

# ======================================================================================
# Phase 2: DESeqDataSet Construction & Pre-filtering
# ======================================================================================

dds <- DESeqDataSet(se, design = ~ condition)

# Pre-filtering: keep rows with at least 10 counts across minimum group size
smallestGroupSize <- min(table(dds$condition))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Relevel to set Normal Tissue (NT) as reference
dds$condition <- relevel(dds$condition, ref = "NT")

# ======================================================================================
# Phase 3: Differential Gene Expression Analysis (DESeq2)
# ======================================================================================

# Run parallelized DESeq analysis
register(MulticoreParam(workers = 4)) # Adjust workers based on CPU
dds <- DESeq(dds)

res <- results(dds, alpha = 0.05)
# ======================================================================================
# Phase 4: Alternative Shrinkage Estimators (LFC Shrinkage)
# ======================================================================================

# Apply apeglm shrinkage for robust fold-change estimates
res_shrunk <- lfcShrink(dds, coef = "condition_TP_vs_NT", type = "apeglm")

# ======================================================================================
# Phase 5: Gene Annotation
# ======================================================================================

res_df <- as.data.frame(res_shrunk)

# Clean ENSEMBL IDs and map to HGNC Symbols
clean_ensembl <- gsub("\\..*", "", rownames(res_df))
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = clean_ensembl,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Annotate Expression Status (UP/DOWN/Stable)
res_df$status <- "Stable"
res_df$status[res_df$log2FoldChange > 0.5 & res_df$padj < 0.05] <- "UP"
res_df$status[res_df$log2FoldChange < -0.5 & res_df$padj < 0.05] <- "DOWN"

# Order by adjusted p-value and save results
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file = "results/DESeq2_LumA_vs_Normal_Annotated.csv")
saveRDS(dds, file = "results/dds_calculated_backup.rds")

# ======================================================================================
# Phase 6: Visualization (PCA, Volcano, Heatmap)
# ======================================================================================

# ---------------------------------------------------------
# 6.1 PCA Plot
# ---------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

custom_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, type = "norm", level = 0.95) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Luminal A vs Normal Breast Tissue") +
  theme_gray() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_blank(), 
    legend.position = "right",
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

ggsave(filename = "results/PCA_Plot_Ellipses.pdf", plot = custom_pca, width = 8, height = 6)

# ---------------------------------------------------------
# 6.2 Volcano Plot (Targeting ITGB1)
# ---------------------------------------------------------
keyvals <- rep('grey', nrow(res_df))
names(keyvals) <- rep('NS', nrow(res_df))
keyvals[which(abs(res_df$log2FoldChange) > 0.5 & res_df$padj >= 0.05)] <- 'green'
names(keyvals)[which(keyvals == 'green')] <- 'Log2 FC'
keyvals[which(abs(res_df$log2FoldChange) <= 0.5 & res_df$padj < 0.05)] <- 'royalblue'
names(keyvals)[which(keyvals == 'royalblue')] <- 'padj'
keyvals[which(abs(res_df$log2FoldChange) > 0.5 & res_df$padj < 0.05)] <- 'red'
names(keyvals)[which(keyvals == 'red')] <- 'padj & log2FC'

# Highlight target gene
keyvals[which(res_df$symbol == "ITGB1")] <- 'black'
names(keyvals)[which(res_df$symbol == "ITGB1")] <- 'ITGB1 (Target)'

my_volcano <- EnhancedVolcano(res_df,
                              lab = res_df$symbol,
                              selectLab = c("ITGB1"),
                              x = 'log2FoldChange', y = 'padj',
                              ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                              pCutoff = 0.05, FCcutoff = 0.5,
                              colCustom = keyvals,
                              pointSize = 2.5, labSize = 6.0, labFace = 'bold',
                              title = 'Luminal A vs Normal Breast Tissue',
                              subtitle = 'Volcano plot (FC Cutoff = 0.5, padj < 0.05)',
                              legendPosition = 'right',
                              drawConnectors = TRUE, colConnectors = 'black') + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(filename = "results/Volcano_Plot_ITGB1.pdf", plot = my_volcano, width = 10, height = 8)

# ---------------------------------------------------------
# 6.3 Minimal Heatmap
# ---------------------------------------------------------
ens_col <- names(res_df)[sapply(res_df, function(x) any(grepl("^ENSG", as.character(x))))]
if(length(ens_col) > 0 && !any(grepl("^ENSG", rownames(res_df)))) {
  rownames(res_df) <- make.unique(as.character(res_df[[ens_col[1]]]))
}

# Filter significant genes with robust baseline expression
sig_genes <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5 & res_df$baseMean > 50), ]
sig_ordered <- sig_genes[order(sig_genes$padj), ]
valid_ids <- head(rownames(sig_ordered), 40)

# Ensure target gene inclusion
target_gene_id <- rownames(res_df[which(res_df$symbol == "ITGB1"), ])
if (length(target_gene_id) > 0 && !(target_gene_id %in% valid_ids)) {
  valid_ids <- c(valid_ids, target_gene_id)
}

valid_ids <- intersect(valid_ids, rownames(vsd))
mat <- assay(vsd)[valid_ids, , drop=FALSE]
mat <- mat[complete.cases(mat), , drop=FALSE]
mat <- mat[apply(mat, 1, var) > 0, , drop=FALSE]
valid_ids <- rownames(mat)

# Assign gene symbols
gene_symbols <- res_df[valid_ids, "symbol"]
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- valid_ids[is.na(gene_symbols) | gene_symbols == ""]
rownames(mat) <- make.unique(as.character(gene_symbols))

# Annotations
anno_col <- as.data.frame(colData(vsd)[, "condition", drop=FALSE])
colnames(anno_col) <- c("Tissue_Status")
anno_row <- data.frame(Log10_BaseMean = log10(res_df[valid_ids, "baseMean"] + 1))
rownames(anno_row) <- rownames(mat)

my_anno_colors <- list(Log10_BaseMean = colorRampPalette(c("red", "white", "forestgreen"))(50))
patient_order <- order(anno_col$Tissue_Status)
mat_clean <- mat[, patient_order]
anno_col_clean <- anno_col[patient_order, , drop=FALSE]

# Generate Heatmap
pheatmap(mat_clean,
         annotation_col = anno_col_clean,
         annotation_row = anno_row,
         annotation_colors = my_anno_colors,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = NA,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename = "results/Heatmap_Top40_ITGB1.pdf",
         width = 9, height = 10)