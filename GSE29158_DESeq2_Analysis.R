dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData = sample_metadata,
  design = ~ SEX
)
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("SEX", "M", "F"))
res_df <- as.data.frame(res)
res_df$Gene_ID <- rownames(res_df)

# ---- 7. Gene annotations from recount3 ----
gene_data <- as.data.frame(SummarizedExperiment::rowData(rse))
gene_data$Gene_ID <- rownames(gene_data)
idx <- match(res_df$Gene_ID, gene_data$Gene_ID)
if ("gene_name" %in% names(gene_data)) {
  res_df$Gene_Name <- as.character(gene_data$gene_name[idx])
} else {
  res_df$Gene_Name <- res_df$Gene_ID
}
missing_names <- is.na(res_df$Gene_Name) | !nzchar(res_df$Gene_Name)
res_df$Gene_Name[missing_names] <- res_df$Gene_ID[missing_names]
res_df$Gene_Type <- if ("gene_type" %in% names(gene_data)) {
  as.character(gene_data$gene_type[idx])
} else {
  NA_character_
}

# ---- 8. Consistent significance classification ----
res_df$Regulation <- dplyr::case_when(
  !is.na(res_df$padj) & res_df$padj < padj_cutoff &
    res_df$log2FoldChange > lfc_cutoff ~ "Male Upregulated",
  !is.na(res_df$padj) & res_df$padj < padj_cutoff &
    res_df$log2FoldChange < -lfc_cutoff ~ "Female Upregulated",
  TRUE ~ "Not Significant"
)
male_up <- res_df[res_df$Regulation == "Male Upregulated", , drop = FALSE]
female_up <- res_df[res_df$Regulation == "Female Upregulated", , drop = FALSE]
male_up <- male_up[order(male_up$padj), , drop = FALSE]
female_up <- female_up[order(female_up$padj), , drop = FALSE]
message("Male upregulated: ", nrow(male_up),
        "; female upregulated: ", nrow(female_up))

# ---- 9. Volcano plot ----
volcano_df <- res_df
volcano_df$negLog10Padj <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))
label_genes <- dplyr::bind_rows(
  head(male_up, 10), head(female_up, 10)
)
label_genes$negLog10Padj <- -log10(pmax(label_genes$padj, .Machine$double.xmin))
volcano_plot <- ggplot2::ggplot(
  volcano_df, ggplot2::aes(x = log2FoldChange, y = negLog10Padj)
) +
  ggplot2::geom_point(ggplot2::aes(color = Regulation), alpha = 0.7, size = 1.7, na.rm = TRUE) +
  ggrepel::geom_text_repel(
    data = label_genes, ggplot2::aes(label = Gene_Name),
    size = 3, max.overlaps = Inf, show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c(
    "Male Upregulated" = "darkgreen",
    "Female Upregulated" = "red",
    "Not Significant" = "grey70"
  )) +
  ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  ggplot2::labs(
    title = "GSE29158: Male vs Female",
    subtitle = "DESeq2 | FDR < 0.05 and |log2FC| > 0.5",
    x = "Log2 Fold Change (Male / Female)",
    y = "-Log10 Adjusted P-value", color = "Expression"
  ) +
  ggplot2::theme_bw(base_size = 12)
ggplot2::ggsave(
  file.path(output_dir, "GSE29158_Volcano_Plot.png"),
  volcano_plot, width = 11, height = 8, dpi = 300
)

# ---- 10. Heatmap: top 20 in each direction ----
top_gene_ids <- unique(c(head(male_up$Gene_ID, 20), head(female_up$Gene_ID, 20)))
if (length(top_gene_ids) > 1L) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  vst_matrix <- SummarizedExperiment::assay(vsd)
  top_gene_ids <- intersect(top_gene_ids, rownames(vst_matrix))
  heatmap_matrix <- vst_matrix[top_gene_ids, , drop = FALSE]
  labels <- res_df$Gene_Name[match(rownames(heatmap_matrix), res_df$Gene_ID)]
  rownames(heatmap_matrix) <- make.unique(ifelse(
    is.na(labels) | !nzchar(labels), rownames(heatmap_matrix), labels
  ))
  sample_order <- rownames(sample_metadata)[order(sample_metadata$SEX)]
  heatmap_matrix <- heatmap_matrix[, sample_order, drop = FALSE]
  annotation_col <- data.frame(SEX = sample_metadata[sample_order, "SEX"])
  rownames(annotation_col) <- sample_order
  pheatmap::pheatmap(
    heatmap_matrix, scale = "row", cluster_rows = TRUE,
    cluster_cols = FALSE, annotation_col = annotation_col,
    gaps_col = sum(sample_metadata$SEX == "F"),
    fontsize_col = 7, fontsize_row = 8, angle_col = 90,
    color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = "GSE29158: Male vs Female | FDR < 0.05, |log2FC| > 0.5",
    filename = file.path(output_dir, "GSE29158_Heatmap_log2FC_0.5.png"),
    width = 12, height = 11
  )
} else {
  warning("Too few significant genes to create the top-gene heatmap")
}

# ---- 11. Export DESeq2 results and glycosylation candidates ----
glycan_genes <- c(
  "FUT8", "ST6GAL1", "B4GALT1", "MGAT1", "MGAT2", "MGAT3",
  "MGAT5", "STT3A", "STT3B", "B3GNT2", "MAN2A1"
)
res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
writexl::write_xlsx(
  list(
    "All Filtered Genes" = res_df,
    "All Significant" = subset(res_df, !is.na(padj) & padj < padj_cutoff),
    "Male Upregulated" = male_up,
    "Female Upregulated" = female_up,
    "N-Glycan Candidates" = subset(res_df, Gene_Name %in% glycan_genes)
  ),
  file.path(output_dir, "GSE29158_Updated_log2FC_0.5.xlsx")
)
message("Saved outputs to: ", normalizePath(output_dir))
