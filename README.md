# GSE29158: Sex-Based Differential Gene Expression Analysis

## Project Overview

An exploratory computational transcriptomics project investigating sex-associated gene expression differences in human B-cell samples using R, recount3, and DESeq2.

**Dataset:** GSE29158 / SRP001563  
**Samples:** 20 biological samples (10 female, 10 male)  
**Comparison:** Male vs Female  
**Analysis:** Differential gene expression using DESeq2

## Research Objective

To investigate transcriptional differences between male and female human B-cell samples, identify sex-associated differentially expressed genes, and explore candidate genes potentially involved in glycosylation-related biological processes.

## Bioinformatics Workflow

1. Retrieved RNA-seq gene-level count data using recount3.
2. Matched SRA sequencing runs to GEO biological samples using metadata.
3. Aggregated sequencing runs belonging to the same biological sample.
4. Prepared sample metadata and performed count-based differential expression analysis using DESeq2.
5. Annotated gene identifiers using available annotation resources.
6. Generated volcano plots and heatmaps for visualization.
7. Exported differential expression results for downstream interpretation.

## Differential Expression Criteria

- Adjusted p-value (Benjamini–Hochberg FDR) < 0.05
- Absolute log2 fold change > 0.5
- Positive log2 fold change: higher expression in males
- Negative log2 fold change: higher expression in females

## Results and Visualization

The analysis identified sex-associated transcriptional patterns, including higher expression of Y-linked genes such as **RPS4Y1, DDX3Y, and KDM5D** in male samples and **XIST** in female samples.

**Volcano plot:** Visualizes gene-level log2 fold changes and adjusted p-values, distinguishing male-enriched and female-enriched expression.

**Heatmap:** Displays selected differentially expressed genes across the 20 samples using row-scaled, variance-stabilized expression values.

The project also explores glycosylation-associated candidate genes as a hypothesis-generating extension of the differential expression analysis.

## Tools and Technologies

- **Programming:** R
- **Data retrieval:** recount3, GEOquery
- **Differential expression:** DESeq2
- **Data handling and annotation:** SummarizedExperiment and Bioconductor annotation resources
- **Visualization:** ggplot2, pheatmap
- **Results export:** writexl

## Reproducibility

The R analysis script documents the computational workflow.

The analysis requires the corresponding public sequencing data and correctly matched sample metadata, including the local metadata and SRA run-mapping CSV files used in the original workflow.

The uploaded script is a reconstructed and organized version of the analysis code. It has not been independently rerun and validated end-to-end in this repository.

## Limitations

- The analysis is exploratory and based on 20 biological samples.
- Sex-associated transcriptional differences do not independently establish causal biological mechanisms.
- Glycosylation-related candidate genes require additional statistical assessment and experimental validation before conclusions about altered glycosylation can be drawn.

## Data Source

Publicly available NCBI Gene Expression Omnibus dataset **GSE29158** and associated sequencing study **SRP001563**.

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29158
- SRA: https://www.ncbi.nlm.nih.gov/sra/?term=SRP001563

## Project Status

Exploratory analysis completed; code organized for portfolio presentation. Independent end-to-end reproducibility verification remains pending.
