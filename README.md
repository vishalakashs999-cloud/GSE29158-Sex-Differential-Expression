# GSE29158: Sex-Based Differential Gene Expression Analysis

## Project Overview

Computational analysis of sex-associated gene expression differences in human B-cell transcriptomic data using R and DESeq2.

**Dataset:** GSE29158 / SRP001563  
**Samples:** 20 (10 female, 10 male)  
**Comparison:** Male vs Female

## Bioinformatics Workflow

1. Retrieved gene-level expression data using recount3.
2. Matched SRA runs with GEO sample metadata.
3. Aggregated sequencing runs by biological sample.
4. Performed differential expression analysis using DESeq2.
5. Annotated genes using Ensembl and Bioconductor resources.
6. Generated volcano plots and heatmaps using ggplot2 and pheatmap.

## Differential Expression Criteria

- Adjusted p-value (FDR) < 0.05
- Absolute log2 fold change > 0.5
- Positive log2FC: higher expression in males
- Negative log2FC: higher expression in females

## Results

The analysis identified sex-associated expression patterns, including differences involving Y-chromosome genes such as RPS4Y1, DDX3Y and KDM5D, and X-linked genes including XIST.

The volcano plot visualizes the magnitude and statistical significance of differential expression. The heatmap displays expression patterns across male and female samples using row-scaled, variance-stabilized expression values.

## Research Focus

Exploration of sex-associated transcriptional differences in human B cells, including candidate genes involved in glycosylation-related biological pathways.

## Tools

R, recount3, GEOquery, DESeq2, SummarizedExperiment, AnnotationDbi, biomaRt, ggplot2, pheatmap and writexl.

## Data Source

Publicly available GEO dataset GSE29158 and associated SRA study SRP001563.
