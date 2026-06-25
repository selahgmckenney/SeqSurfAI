test_that("DEG filtering assigns direction and preserves column order", {
  dat <- load_test_dataset("synthetic")
  out <- filter_deg(
    dat$deg,
    "treated vs control",
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )

  expect_equal(names(out)[1:3], c("contrast", "gene", "direction"))
  expect_true(all(out$contrast == "treated vs control"))
  expect_equal(out$direction[out$gene == "MYCN"], "Up")
  expect_equal(out$direction[out$gene == "NTRK1"], "Down")
})

test_that("PCA labels and loading selectors are internally consistent", {
  dat <- load_test_dataset("synthetic")
  expect_match(pc_axis_label("PC1", dat$pca_variance), "^PC1 \\([0-9.]+%\\)$")

  top <- top_loading_genes(dat$loadings, "PC1", 4, "absolute")
  expect_equal(nrow(top), 4)
  expect_true(all(diff(abs(top$loading)) <= 0))

  plot <- plot_pca(
    dat$pca, "PC1", "PC2", "condition",
    "sex", "", dat$pca_variance
  )
  expect_s3_class(plot, "ggplot")
})

test_that("GSEA plot filtering obeys collection, direction, cutoff, and top N", {
  dat <- load_test_dataset("synthetic")
  gsea <- dat$gsea[
    dat$gsea$contrast == "treated vs control" &
      dat$gsea$collection == "Hallmark",
    ,
    drop = FALSE
  ]

  positive <- pathway_plot_data(gsea, 0.1, "Positive", 10)
  negative <- pathway_plot_data(gsea, 0.1, "Negative", 10)

  expect_lte(nrow(positive), 10)
  expect_lte(nrow(negative), 10)
  expect_true(all(positive$NES > 0))
  expect_true(all(negative$NES < 0))
})

test_that("gene expression joins preserve all samples", {
  dat <- load_test_dataset("synthetic")
  out <- expression_long(
    dat$vsd, dat$metadata, "MYCN",
    sample_col = dat$info$sample_id_col
  )

  expect_equal(nrow(out), nrow(dat$metadata))
  expect_false(anyNA(out$expression))
  expect_setequal(out[[dat$info$sample_id_col]], dat$metadata[[dat$info$sample_id_col]])
})
