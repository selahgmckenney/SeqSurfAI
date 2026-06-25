pc_axis_label <- function(pc, variance) {
  if (!is.null(variance) && pc %in% names(variance)) {
    paste0(pc, " (", round(variance[[pc]], 1), "%)")
  } else {
    pc
  }
}

plot_pca <- function(
  pca_df,
  x_pc,
  y_pc,
  color_by,
  shape_by = NULL,
  label_by = NULL,
  variance = NULL
) {
  pca_df <- pca_df[, !duplicated(names(pca_df)), drop = FALSE]
  if (is.null(shape_by) || !nzchar(shape_by) || !shape_by %in% names(pca_df)) shape_by <- NULL
  if (is.null(label_by) || !nzchar(label_by) || !label_by %in% names(pca_df)) label_by <- NULL
  validate(
    need(x_pc %in% names(pca_df), "Selected x-axis PC is not available."),
    need(y_pc %in% names(pca_df), "Selected y-axis PC is not available."),
    need(color_by %in% names(pca_df), "Selected color column is not available.")
  )

  aes_args <- ggplot2::aes(
    x = .data[[x_pc]],
    y = .data[[y_pc]],
    color = .data[[color_by]]
  )

  if (!is.null(shape_by) && nzchar(shape_by) && shape_by %in% names(pca_df)) {
    aes_args$shape <- rlang::expr(.data[[shape_by]])
  }

  p <- ggplot2::ggplot(pca_df, aes_args) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      x = pc_axis_label(x_pc, variance),
      y = pc_axis_label(y_pc, variance),
      color = color_by,
      shape = shape_by
    )

  if (!is.null(label_by) && nzchar(label_by) && label_by %in% names(pca_df)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[[label_by]]),
      size = 3,
      vjust = -0.7,
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }

  p
}

top_loading_genes <- function(loadings, pc, n, direction) {
  df <- loadings |>
    dplyr::filter(.data$PC == pc)

  if (direction == "positive") {
    df <- df |> dplyr::filter(.data$loading > 0) |> dplyr::arrange(dplyr::desc(.data$loading))
  } else if (direction == "negative") {
    df <- df |> dplyr::filter(.data$loading < 0) |> dplyr::arrange(.data$loading)
  } else {
    df <- df |> dplyr::arrange(dplyr::desc(abs(.data$loading)))
  }

  df |> dplyr::slice_head(n = n)
}

plot_top_loadings <- function(loadings_df) {
  validate(need(nrow(loadings_df) > 0, "No loading genes matched this selection."))

  loadings_df$gene <- stats::reorder(loadings_df$gene, loadings_df$loading)
  ggplot2::ggplot(loadings_df, ggplot2::aes(x = gene, y = loading, fill = loading > 0)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "#d95f0e")) +
    ggplot2::labs(x = NULL, y = "PCA loading")
}
