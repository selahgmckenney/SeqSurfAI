plot_expression_heatmap <- function(
  vsd,
  metadata,
  genes,
  annotation_cols,
  sample_col = "sample",
  silent = FALSE
) {
  genes <- intersect(genes, rownames(vsd))
  validate(need(length(genes) > 1, "At least two requested genes must be present in the expression matrix."))

  mat <- vsd[genes, , drop = FALSE]
  mat <- mat[, colnames(mat) %in% metadata[[sample_col]], drop = FALSE]
  mat <- mat[, match(metadata[[sample_col]][metadata[[sample_col]] %in% colnames(mat)], colnames(mat)), drop = FALSE]
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0

  annotation_cols <- intersect(annotation_cols, names(metadata))
  annotation <- NULL
  if (length(annotation_cols) > 0) {
    annotation <- metadata[, annotation_cols, drop = FALSE]
    rownames(annotation) <- metadata[[sample_col]]
    annotation <- annotation[colnames(mat_scaled), , drop = FALSE]
  }

  pheatmap::pheatmap(
    mat_scaled,
    annotation_col = annotation,
    show_colnames = FALSE,
    fontsize_row = 8,
    border_color = NA,
    main = "Scaled expression",
    silent = silent
  )
}

plot_gene_expression <- function(expr_df, group_col) {
  validate(need(group_col %in% names(expr_df), "Selected grouping column is not available."))

  ggplot2::ggplot(expr_df, ggplot2::aes(x = .data[[group_col]], y = expression, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, show.legend = FALSE) +
    ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.75, show.legend = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(gene), scales = "free_y") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(x = group_col, y = "Normalized expression")
}
