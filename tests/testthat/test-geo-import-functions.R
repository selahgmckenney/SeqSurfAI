test_that("GEO import design helpers produce app-safe defaults", {
  design <- parse_design_from_inputs(
    accession = "gse123",
    platform_type = "normalized_expression",
    expression_data_type = "normalized_expression",
    group_col = "condition",
    positive_group = "treated",
    negative_group = "control",
    covariates = "sex, batch",
    sample_id_col = "geo_accession",
    patient_id_col = "patient",
    display_sample_col = "title",
    default_timepoint_col = "condition",
    contrast_label = "treated vs control"
  )

  expect_equal(design$accession, "GSE123")
  expect_equal(geo_dataset_id("GSE123"), "geo_gse123")
  expect_equal(design$covariates, c("sex", "batch"))
  expect_true(all(c("Field", "Value") %in% names(design_to_table(design))))
})

test_that("GEO report writer creates standalone HTML", {
  dat <- load_test_dataset("synthetic")
  kb <- list(
    accession = "GSE123",
    geo = list(title = "Synthetic study", pubmed_id = "0", sample_count = nrow(dat$metadata), platform_accessions = "GPL"),
    pubmed = list(title = "Synthetic paper"),
    extracted = list(
      study_objective = "Study treatment.",
      cohort_characteristics = "Synthetic samples.",
      sample_groups = "Treated and control groups.",
      statistical_methods = "limma",
      reported_pathways = "MYCN",
      reported_biomarkers = "MYCN",
      study_limitations = "Processed expression."
    ),
    claims = data.frame(
      claim_type = "biomarker",
      entity = "MYCN",
      comparison = "amplified vs non-amplified",
      reported_direction = "up",
      reported_significance = "significant",
      evidence_sentence = "MYCN was higher in amplified tumors.",
      source = "test",
      confidence = "high"
    )
  )
  out_dir <- file.path(tempdir(), "geo_report_test")
  dir.create(out_dir, showWarnings = FALSE)
  path <- write_geo_reproducibility_report(dat, kb, out_dir)
  expect_true(file.exists(path))
  report_lines <- readLines(path, warn = FALSE)
  expect_true(any(grepl("SeqSurf AI reproducibility report", report_lines, fixed = TRUE)))
  expect_true(any(grepl("Study Memory", report_lines, fixed = TRUE)))
  expect_true(any(grepl("Analyzability Route", report_lines, fixed = TRUE)))
  expect_true(any(grepl("Extracted Claims", report_lines, fixed = TRUE)))
  expect_true(any(grepl("Agreement And Disagreement", report_lines, fixed = TRUE)))
  expect_true(any(grepl("Discovery Plan", report_lines, fixed = TRUE)))
})

test_that("demo walkthrough exists for the curated GEO example", {
  walkthrough <- file.path("docs", "demo_walkthrough_GSE19804.md")
  if (!file.exists(walkthrough)) {
    walkthrough <- file.path("..", "..", "docs", "demo_walkthrough_GSE19804.md")
  }
  expect_true(file.exists(walkthrough))
  text <- readLines(walkthrough, warn = FALSE)
  expect_true(any(grepl("GSE19804", text, fixed = TRUE)))
  expect_true(any(grepl("In-app GEO processed-matrix reanalysis", text, fixed = TRUE)))
})

test_that("GEO supplementary count helpers classify and prepare count matrices", {
  kb <- list(
    accession = "GSE999",
    geo = list(
      supplementary_files = c(
        "https://example.org/GSE999_gene_counts.txt.gz",
        "https://example.org/GSE999_RAW.tar",
        "https://example.org/SRR999.fastq.gz"
      )
    )
  )
  candidates <- geo_supplementary_count_candidates(kb)
  expect_true(candidates$File[[1]] == "GSE999_gene_counts.txt.gz")
  expect_equal(candidates$Status[[1]], "Likely count matrix")

  count_path <- file.path(tempdir(), "GSE999_gene_counts.txt")
  writeLines(
    c(
      "gene\tGSM1\tGSM2\tGSM3\tGSM4",
      "MYCN\t100\t120\t700\t750",
      "NTRK1\t800\t820\t100\t90",
      "GENE3\t50\t60\t55\t65",
      "GENE4\t900\t880\t85\t95"
    ),
    count_path
  )
  count_df <- read_geo_supplementary_table(count_path)
  metadata <- data.frame(
    geo_accession = paste0("GSM", 1:4),
    condition = c("control", "control", "treated", "treated"),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$geo_accession
  inputs <- prepare_geo_count_upload_inputs(count_df, metadata)
  expect_equal(inputs$gene_col, "gene")
  expect_equal(ncol(inputs$expr_df), 5)
  expect_true(all(c("app_sample_id", "condition") %in% names(inputs$metadata)))
  expect_true(select_count_de_method("DESeq2") %in% c("deseq2", "edger", "limma_voom", "fallback"))
  expect_equal(infer_count_de_method_from_text("The authors used DESeq2 for differential expression."), "deseq2")

  design <- rule_based_geo_metadata_design("GSE999", metadata)
  expect_equal(design$expression_data_type, "raw_counts")
  expect_equal(design$sample_id_col, "geo_accession")

  out_dir <- file.path(tempdir(), "count_matrix_pipeline_dataset")
  settings <- list(
    dataset_id = "count_matrix_pipeline_dataset",
    display_name = "Count matrix pipeline dataset",
    gene_col = "gene",
    sample_col = "app_sample_id",
    group_col = "condition",
    positive_group = "treated",
    negative_group = "control",
    patient_id_col = "app_sample_id",
    display_sample_col = "geo_accession",
    contrast_label = "treated vs control",
    count_method = "limma_voom"
  )
  result <- write_uploaded_app_dataset(inputs$expr_df, inputs$metadata, settings, tempdir())
  expect_true(file.exists(file.path(result$dataset_dir, "count_matrix_input.csv")))
  expect_true(file.exists(file.path(result$dataset_dir, "sample_metadata_input.csv")))
  expect_true(file.exists(file.path(result$dataset_dir, "count_matrix_pipeline.R")))
  script <- readLines(file.path(result$dataset_dir, "count_matrix_pipeline.R"), warn = FALSE)
  expect_true(any(grepl("SeqSurf count-matrix pipeline handoff", script, fixed = TRUE)))
  expect_true(file.exists(file.path(out_dir, "deg_results.csv")))
})
