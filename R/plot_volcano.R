plot_volcano_ggplot <- function(deg, padj_cutoff, lfc_cutoff) {
  deg$neg_log10_padj <- safe_neg_log10(deg$padj)
  deg$direction <- dplyr::case_when(
    !is.na(deg$padj) & deg$padj <= padj_cutoff & deg$log2FoldChange >= lfc_cutoff ~ "Up",
    !is.na(deg$padj) & deg$padj <= padj_cutoff & deg$log2FoldChange <= -lfc_cutoff ~ "Down",
    TRUE ~ "Not significant"
  )

  p <- ggplot2::ggplot(
    deg,
    ggplot2::aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      color = direction,
      text = paste0(
        "Gene: ", gene,
        "<br>log2FC: ", signif(log2FoldChange, 3),
        "<br>padj: ", signif(padj, 3)
      )
    )
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.7) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey45") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey45") +
    ggplot2::scale_color_manual(
      values = c("Up" = "#2c7fb8", "Down" = "#d95f0e", "Not significant" = "grey70")
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    )

  p
}

plot_volcano <- function(deg, padj_cutoff, lfc_cutoff) {
  p <- plot_volcano_ggplot(deg, padj_cutoff, lfc_cutoff)
  if (requireNamespace("plotly", quietly = TRUE)) {
    plotly::ggplotly(p, tooltip = "text")
  } else {
    p
  }
}

pathway_plot_data <- function(gsea, padj_cutoff, direction = "Both", top_n = 20) {
  out <- gsea |>
    add_pathway_overlap_stats() |>
    dplyr::filter(!is.na(.data$padj), .data$padj <= padj_cutoff)

  if (direction == "Positive") {
    out <- dplyr::filter(out, .data$NES > 0)
  } else if (direction == "Negative") {
    out <- dplyr::filter(out, .data$NES < 0)
  }

  out |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES))) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::arrange(.data$NES)
}

plot_pathway_nes_ggplot <- function(gsea, padj_cutoff, direction = "Both", top_n = 20) {
  df <- pathway_plot_data(gsea, padj_cutoff, direction, top_n)

  validate(need(nrow(df) > 0, "No pathways pass the selected adjusted p-value cutoff."))

  df$pathway_label <- gsub("^(HALLMARK|REACTOME|GOBP|WP)_", "", df$pathway)
  df$pathway_label <- gsub("_", " ", df$pathway_label)
  df$pathway_label <- stats::reorder(df$pathway_label, df$NES)
  df$direction <- ifelse(df$NES > 0, "Positive NES", "Negative NES")
  df$hover <- paste0(
    "Pathway: ", df$pathway,
    "<br>NES: ", signif(df$NES, 3),
    "<br>padj: ", signif(df$padj, 3),
    "<br>Leading-edge genes: ", df$leading_edge_genes,
    if ("size" %in% names(df)) paste0(" of ", df$size, " pathway genes") else ""
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = pathway_label, y = NES, fill = direction, text = hover)
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::scale_fill_manual(values = c("Positive NES" = "#2c7fb8", "Negative NES" = "#d95f0e")) +
    ggplot2::labs(x = NULL, y = "Normalized enrichment score")
}

plot_pathway_nes <- function(gsea, padj_cutoff, direction = "Both", top_n = 20) {
  p <- plot_pathway_nes_ggplot(gsea, padj_cutoff, direction, top_n)
  if (requireNamespace("plotly", quietly = TRUE)) {
    plotly::ggplotly(p, tooltip = "text")
  } else {
    p
  }
}
