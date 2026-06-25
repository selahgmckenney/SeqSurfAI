select_signature_genes <- function(
  deg,
  contrast,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  direction = "Up",
  max_genes = 500
) {
  out <- deg
  if ("contrast" %in% names(out)) {
    out <- out[out$contrast == contrast, , drop = FALSE]
  }
  out <- out[
    !is.na(out$padj) &
      out$padj <= padj_cutoff &
      !is.na(out$log2FoldChange) &
      abs(out$log2FoldChange) >= lfc_cutoff,
    ,
    drop = FALSE
  ]

  if (direction == "Up") {
    out <- out[out$log2FoldChange > 0, , drop = FALSE]
  } else if (direction == "Down") {
    out <- out[out$log2FoldChange < 0, , drop = FALSE]
  }

  out <- out[order(out$padj, -abs(out$log2FoldChange)), , drop = FALSE]
  out <- out[!duplicated(toupper(out$gene)), , drop = FALSE]
  utils::head(out, max_genes)
}

enrichment_barplot <- function(
  results,
  term_col,
  padj_col,
  score_col = NULL,
  title = NULL,
  top_n = 20
) {
  validate(need(nrow(results) > 0, "No enriched terms passed the selected cutoff."))
  results <- results[!is.na(results[[padj_col]]), , drop = FALSE]
  results <- results[order(results[[padj_col]]), , drop = FALSE]
  results <- utils::head(results, top_n)
  results$display_term <- factor(
    results[[term_col]],
    levels = rev(results[[term_col]])
  )
  results$display_score <- if (is.null(score_col)) {
    -log10(pmax(results[[padj_col]], .Machine$double.xmin))
  } else {
    results[[score_col]]
  }

  ggplot2::ggplot(
    results,
    ggplot2::aes(x = display_term, y = display_score, fill = -log10(.data[[padj_col]]))
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "C", end = 0.9, name = "-log10 adj. p") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = NULL,
      y = if (is.null(score_col)) "-log10 adjusted p-value" else score_col,
      title = title
    )
}

run_clusterprofiler_go <- function(
  signature,
  universe_genes,
  ontology = "BP",
  padj_cutoff = 0.05
) {
  validate(need(nrow(signature) >= 5, "At least five selected genes are required."))

  selected_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(signature$gene),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  universe_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(universe_genes),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  selected_ids <- unique(stats::na.omit(unname(selected_map)))
  universe_ids <- unique(stats::na.omit(unname(universe_map)))
  validate(need(length(selected_ids) >= 5, "Fewer than five selected genes mapped to Entrez IDs."))

  result <- clusterProfiler::enrichGO(
    gene = selected_ids,
    universe = universe_ids,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  out <- as.data.frame(result)
  out <- out[!is.na(out$p.adjust) & out$p.adjust <= padj_cutoff, , drop = FALSE]
  out
}

