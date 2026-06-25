test_that("clean SeqSurf app can start without bundled datasets", {
  ids <- available_datasets(data_root)
  expect_equal(ids, character())
  expect_equal(dataset_choice_list(ids, data_root), c("No active dataset" = ""))
})

test_that("app-ready datasets satisfy the core data contract", {
  dat <- load_test_dataset("synthetic")
  sample_col <- dat$info$sample_id_col

  expect_true(sample_col %in% names(dat$metadata))
  expect_false(anyDuplicated(dat$metadata[[sample_col]]) > 0)
  expect_setequal(dat$metadata[[sample_col]], colnames(dat$vsd))
  expect_true(all(c("PC1", "PC2") %in% names(dat$pca)))
  expect_true(all(c("gene", "log2FoldChange", "padj", "contrast") %in% names(dat$deg)))
  expect_true(all(c("pathway", "NES", "padj", "contrast", "collection") %in% names(dat$gsea)))
  expect_false(anyDuplicated(dat$deg[, c("contrast", "gene")]) > 0)
  expect_false(anyDuplicated(dat$gsea[, c("contrast", "collection", "pathway")]) > 0)
})
