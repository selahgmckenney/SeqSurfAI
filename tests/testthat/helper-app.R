suppressPackageStartupMessages({
  library(testthat)
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(bslib)
})

app_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
data_root <- file.path(app_root, "data")

source(file.path(app_root, "R", "helper_functions.R"))
source(file.path(app_root, "R", "plot_pca.R"))
source(file.path(app_root, "R", "plot_volcano.R"))
source(file.path(app_root, "R", "plot_heatmap.R"))
source(file.path(app_root, "R", "enrichment_functions.R"))
source(file.path(app_root, "R", "ai_assistant_functions.R"))
source(file.path(app_root, "R", "geo_import_functions.R"))
source(file.path(app_root, "R", "upload_functions.R"))

load_test_dataset <- function(id) {
  if (identical(id, "synthetic")) return(make_synthetic_dataset())
  load_rnaseq_dataset(id, data_root)
}

lookup_deg <- function(deg, contrast, gene) {
  deg[deg$contrast == contrast & toupper(deg$gene) == toupper(gene), , drop = FALSE]
}

make_synthetic_dataset <- function() {
  samples <- paste0("S", 1:6)
  metadata <- data.frame(
    app_sample_id = samples,
    title = samples,
    patient = samples,
    condition = rep(c("control", "treated"), each = 3),
    sex = rep(c("Female", "Male"), 3),
    stringsAsFactors = FALSE
  )
  expr <- matrix(
    c(
      5, 5, 5, 9, 9, 9,
      9, 9, 9, 4, 4, 4,
      2, 3, 2, 7, 8, 7,
      7, 8, 7, 2, 3, 2,
      4, 4, 5, 4, 5, 4,
      6, 5, 6, 6, 5, 6,
      8, 8, 7, 4, 4, 5,
      7, 7, 8, 4, 5, 4,
      8, 7, 8, 5, 4, 5,
      3, 3, 4, 8, 8, 7,
      4, 3, 4, 7, 8, 8,
      3, 4, 3, 8, 7, 8
    ),
    nrow = 12,
    byrow = TRUE,
    dimnames = list(c("MYCN", "NTRK1", "GENEUP", "GENEDOWN", "STABLE1", "STABLE2", "EPCAM", "KRT8", "CDH1", "VIM", "FN1", "CDH2"), samples)
  )
  pca_df <- data.frame(
    app_sample_id = samples,
    PC1 = c(-2, -1.5, -1, 1, 1.5, 2),
    PC2 = c(-1, 0, 1, -1, 0, 1),
    metadata,
    check.names = FALSE
  )
  deg <- data.frame(
    contrast = "treated vs control",
    gene = rownames(expr),
    log2FoldChange = c(3, -3, 2, -2, 0.1, -0.1, -3, -3, -3, 3, 3, 3),
    baseMean = rowMeans(expr),
    stat = c(8, -8, 5, -5, 0.2, -0.2, -4, -4, -4, 4, 4, 4),
    pvalue = c(1e-5, 1e-5, 0.001, 0.001, 0.8, 0.8, rep(0.002, 6)),
    padj = c(1e-4, 1e-4, 0.005, 0.005, 0.9, 0.9, rep(0.01, 6)),
    stringsAsFactors = FALSE
  )
  gsea <- data.frame(
    collection = "Hallmark",
    msigdb_version = "test",
    contrast = "treated vs control",
    pathway = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_APOPTOSIS"),
    pval = c(0.001, 0.002),
    padj = c(0.01, 0.02),
    log2err = 0,
    ES = c(0.7, -0.6),
    NES = c(2.5, -2.1),
    size = c(30, 25),
    leadingEdge = c("MYCN;GENEUP", "NTRK1;GENEDOWN"),
    stringsAsFactors = FALSE
  )
  list(
    id = "synthetic",
    info = list(
      dataset_id = "synthetic",
      display_name = "Synthetic test dataset",
      sample_id_col = "app_sample_id",
      display_sample_col = "title",
      patient_id_col = "patient",
      default_group_col = "condition",
      default_timepoint_col = "condition",
      default_contrast = "treated vs control",
      notes = "Tiny synthetic dataset for SeqSurf AI tests."
    ),
    metadata = metadata,
    vsd = expr,
    pca = pca_df,
    pca_variance = c(PC1 = 70, PC2 = 20),
    loadings = data.frame(
      PC = rep("PC1", nrow(expr)),
      gene = rownames(expr),
      loading = c(0.6, -0.6, 0.4, -0.4, 0.1, -0.1),
      stringsAsFactors = FALSE
    ),
    deg = deg,
    gsea = gsea,
    pathway_scores = matrix(1, nrow = 1, ncol = length(samples), dimnames = list("HALLMARK_TEST", samples))
  )
}
