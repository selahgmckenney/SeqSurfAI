test_that("uploaded expression and metadata write an app-ready dataset", {
  expr <- data.frame(
    gene = c("MYCN", "NTRK1", "GENEUP", "GENEDOWN", "STABLE1", "STABLE2"),
    S1 = c(5, 9, 2, 7, 4, 6),
    S2 = c(5, 9, 3, 8, 4, 5),
    S3 = c(5, 9, 2, 7, 5, 6),
    S4 = c(9, 4, 7, 2, 4, 6),
    S5 = c(9, 4, 8, 3, 5, 5),
    S6 = c(9, 4, 7, 2, 4, 6),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample_id = paste0("S", 1:6),
    condition = rep(c("control", "treated"), each = 3),
    patient = paste0("P", 1:6),
    stringsAsFactors = FALSE
  )
  settings <- list(
    dataset_id = "uploaded_unit_test",
    display_name = "Uploaded unit test",
    gene_col = "gene",
    sample_col = "sample_id",
    group_col = "condition",
    positive_group = "treated",
    negative_group = "control",
    patient_id_col = "patient",
    display_sample_col = "sample_id",
    contrast_label = "treated vs control"
  )
  root <- file.path(tempdir(), paste0("seqsurf_upload_test_", Sys.getpid()))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  result <- write_uploaded_app_dataset(expr, metadata, settings, root)
  dat <- load_rnaseq_dataset(result$dataset_id, root)

  expect_equal(result$rows, 6)
  expect_equal(result$genes, 6)
  expect_equal(dat$info$default_contrast, "treated vs control")
  expect_true(all(c("metadata.rds", "vsd_matrix.rds", "deg_results.csv", "dataset_info.rds") %in% basename(list.files(result$dataset_dir))))
  expect_true(all(c("gene", "log2FoldChange", "padj", "contrast") %in% names(dat$deg)))
  expect_equal(dat$deg$log2FoldChange[dat$deg$gene == "MYCN"], 4)
})

test_that("uploaded count-like matrices use count-aware app metadata", {
  expr <- data.frame(
    gene = c("MYCN", "NTRK1", "GENEUP", "GENEDOWN", "STABLE1", "STABLE2"),
    S1 = c(100, 800, 50, 900, 400, 500),
    S2 = c(120, 820, 55, 880, 410, 520),
    S3 = c(110, 790, 60, 910, 405, 510),
    S4 = c(900, 100, 700, 90, 420, 505),
    S5 = c(920, 110, 740, 95, 430, 515),
    S6 = c(880, 120, 710, 85, 415, 500),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample_id = paste0("S", 1:6),
    condition = rep(c("control", "treated"), each = 3),
    stringsAsFactors = FALSE
  )
  settings <- list(
    dataset_id = "uploaded_count_unit_test",
    display_name = "Uploaded count unit test",
    gene_col = "gene",
    sample_col = "sample_id",
    group_col = "condition",
    positive_group = "treated",
    negative_group = "control",
    patient_id_col = "",
    display_sample_col = "sample_id",
    contrast_label = "treated vs control"
  )
  root <- file.path(tempdir(), paste0("seqsurf_upload_count_test_", Sys.getpid()))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  result <- write_uploaded_app_dataset(expr, metadata, settings, root)
  dat <- load_rnaseq_dataset(result$dataset_id, root)

  expect_equal(detect_uploaded_matrix_type(prepare_uploaded_expression(expr, "gene")), "raw_counts")
  expect_true(grepl("count", dat$info$expression_unit, ignore.case = TRUE))
  expect_true(grepl("Matrix type: raw_counts", dat$info$notes, fixed = TRUE))
  expect_true(all(is.finite(dat$vsd)))
})
