# TCGA Differential Gene Expression (DGE) Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21188431.svg)](https://doi.org/10.5281/zenodo.21188431)

## Overview
A complete, automated, and reproducible R pipeline for performing Differential Gene Expression (DGE) analysis using RNA-Seq data from The Cancer Genome Atlas (TCGA). This script handles everything from raw data acquisition to publication-ready data visualization.
While this pipeline is highly adaptable for various cancer types, it is currently configured for a case study comparing Luminal A Breast Cancer (BRCA) primary tumors with normal breast tissues, with a specific focus on highlighting the ITGB1 gene signature.

While this pipeline is highly adaptable for various cancer types, it is currently configured for a case study comparing Luminal A Breast Cancer (BRCA) primary tumors with normal breast tissues, with a specific focus on highlighting the *ITGB1* gene signature.

## Key Features
*   Direct TCGA Integration: Automated downloading and structuring of STAR-Counts data using TCGAbiolinks.
*   Robust DGE Analysis: Utilizes DESeq2 with internal pre-filtering and apeglm LFC shrinkage for accurate transcriptomic profiling.
*   Gene Annotation: Seamless conversion of ENSEMBL IDs to HGNC symbols via org.Hs.eg.db.
*   Publication-Ready Visualizations: Automatically generates highly customized, high-resolution plots:
    *   PCA Plots with 95% confidence ellipses (ggplot2).
    *   Target-highlighted Volcano Plots (EnhancedVolcano).
    *   Clean, cluster-free, baseMean-annotated Heatmaps (pheatmap).

## Prerequisites
To run this pipeline, ensure you have R installed along with the following packages:
`R
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "BiocParallel", "apeglm", "AnnotationDbi", "org.Hs.eg.db", "EnhancedVolcano", "pheatmap"))
install.packages("ggplot2")

## Usage
Simply clone the repository and run the TCGA_LumA_DESeq2_Pipeline.R script. The pipeline will automatically create a results/ directory in your working environment to store all downloaded datasets, intermediate .rds files, CSV outputs, and PDF visualizations.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
