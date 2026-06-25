set.seed(19804)

data_root <- "data"
base_id <- "geo_gse19804"
base_dir <- file.path(data_root, base_id)

if (!dir.exists(base_dir)) {
  stop("Base demo dataset is missing: ", base_dir, call. = FALSE)
}

read_base <- function(file) readRDS(file.path(base_dir, file))
base_metadata <- read_base("metadata.rds")
base_vsd <- read_base("vsd_matrix.rds")
base_pca <- read_base("pca_df.rds")
base_var <- read_base("pca_variance.rds")
base_loadings <- read_base("pca_loadings.rds")
base_scores <- read.csv(file.path(base_dir, "pathway_scores.csv"), row.names = 1, check.names = FALSE)

make_deg <- function(contrast, positive_group, negative_group, effects) {
  genes <- unique(c(
    names(effects),
    "ACTB", "GAPDH", "RPLP0", "HSP90AB1", "B2M", "COL1A1", "KRT18", "CD44",
    "MKI67", "TOP2A", "FOXA1", "ESR1", "ERBB2", "SOX2", "JUN", "FOS"
  ))
  n <- length(genes)
  logfc <- rnorm(n, mean = 0, sd = 0.22)
  names(logfc) <- genes
  logfc[names(effects)] <- effects
  pval <- pmin(0.95, exp(-abs(logfc) * 3.2) * runif(n, 0.0005, 0.12))
  names(pval) <- genes
  pval[names(effects)] <- pmin(pval[names(effects)], runif(length(effects), 1e-7, 0.012))
  padj <- p.adjust(pval, method = "BH")
  names(padj) <- genes
  padj[names(effects)] <- pmin(padj[names(effects)], runif(length(effects), 1e-6, 0.03))
  data.frame(
    contrast = contrast,
    gene = genes,
    log2FoldChange = unname(logfc),
    baseMean = round(runif(n, 5, 800), 3),
    stat = round(unname(logfc) / runif(n, 0.08, 0.28), 3),
    pvalue = pval,
    padj = padj,
    positive_group = positive_group,
    negative_group = negative_group,
    positive_n = 6,
    negative_n = 6,
    check.names = FALSE
  )
}

make_gsea <- function(contrast, nes) {
  pathways <- unique(c(
    names(nes),
    "HALLMARK_APOPTOSIS",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
  ))
  n <- length(pathways)
  nes_values <- rnorm(n, 0, 0.65)
  names(nes_values) <- pathways
  nes_values[names(nes)] <- nes
  pval <- pmin(0.99, exp(-abs(nes_values) * 1.7) * runif(n, 0.001, 0.18))
  names(pval) <- pathways
  pval[names(nes)] <- pmin(pval[names(nes)], runif(length(nes), 1e-5, 0.025))
  padj <- p.adjust(pval, method = "BH")
  names(padj) <- pathways
  padj[names(nes)] <- pmin(padj[names(nes)], runif(length(nes), 1e-4, 0.05))
  data.frame(
    collection = "Hallmark",
    msigdb_version = "2026.1.Hs",
    contrast = contrast,
    pathway = pathways,
    pval = pval,
    padj = padj,
    log2err = runif(n, 0.15, 0.65),
    ES = nes_values / 3,
    NES = unname(nes_values),
    size = sample(35:220, n, replace = TRUE),
    leadingEdge = c(
      "MYCN;MKI67;TOP2A;CCND1",
      "EPCAM;KRT8;CDH1;DDR1",
      "VIM;FN1;COL1A1;CD44",
      rep("JUN;FOS;B2M;STAT5A", max(0, n - 3))
    )[seq_len(n)],
    check.names = FALSE
  )
}

write_dataset <- function(dataset_id, display_name, contrast, positive_group, negative_group, deg, gsea) {
  out_dir <- file.path(data_root, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  info <- read_base("dataset_info.rds")
  info$dataset_id <- dataset_id
  info$display_name <- display_name
  info$default_contrast <- contrast
  info$positive_group <- positive_group
  info$negative_group <- negative_group
  info$notes <- paste(
    "Precomputed comparison-demo dataset derived from the SeqSurf app contract.",
    "Use for fast Compare Studies demonstrations; not a claim of exact paper reproduction."
  )

  saveRDS(base_metadata, file.path(out_dir, "metadata.rds"))
  saveRDS(base_vsd, file.path(out_dir, "vsd_matrix.rds"))
  saveRDS(base_pca, file.path(out_dir, "pca_df.rds"))
  saveRDS(base_var, file.path(out_dir, "pca_variance.rds"))
  saveRDS(base_loadings, file.path(out_dir, "pca_loadings.rds"))
  write.csv(deg, file.path(out_dir, "deg_results.csv"), row.names = FALSE)
  write.csv(gsea, file.path(out_dir, "gsea_hallmark.csv"), row.names = FALSE)
  write.csv(base_scores, file.path(out_dir, "pathway_scores.csv"))
  saveRDS(info, file.path(out_dir, "dataset_info.rds"))
}

claim_table <- function(genes, pathways, comparison, source) {
  gene_claims <- data.frame(
    claim_type = "biomarker",
    entity = genes,
    comparison = comparison,
    reported_direction = "up",
    reported_significance = "reported demo claim",
    evidence_sentence = paste(genes, "is included as a clean SeqSurf comparison-demo biomarker."),
    source = source,
    confidence = "demo",
    check.names = FALSE
  )
  pathway_claims <- data.frame(
    claim_type = "pathway",
    entity = pathways,
    comparison = comparison,
    reported_direction = "up",
    reported_significance = "reported demo claim",
    evidence_sentence = paste(pathways, "is included as a clean SeqSurf comparison-demo pathway."),
    source = source,
    confidence = "demo",
    check.names = FALSE
  )
  rbind(gene_claims, pathway_claims)
}

update_kb_claims <- function(accession, genes, pathways, comparison) {
  path <- file.path(data_root, "ai_knowledge_base", paste0(accession, ".rds"))
  if (!file.exists(path)) return(invisible(FALSE))
  kb <- readRDS(path)
  demo_claims <- claim_table(genes, pathways, comparison, "SeqSurf comparison demo curation")
  kb$claims <- demo_claims
  kb$extracted$reported_biomarkers <- paste(genes, collapse = "; ")
  kb$extracted$reported_pathways <- paste(pathways, collapse = "; ")
  saveRDS(kb, path)
  invisible(TRUE)
}

gse60450_contrast <- "Demo contrast: lactation-associated epithelial survival signal"
gse16476_contrast <- "Demo contrast: neuroblastoma MYCN/proliferation signal"

deg_60450 <- make_deg(
  gse60450_contrast,
  "lactation survival-like signal",
  "baseline mammary cell state",
  c(MCL1 = 1.7, EGF = 1.25, STAT5A = 1.05, EPCAM = 0.9, KRT8 = 0.85, CDH1 = 0.72,
    DDR1 = 0.65, SEMA5A = -0.95, MYCN = 0.58, ALK = 0.45, VIM = -0.62, FN1 = -0.7)
)

gsea_60450 <- make_gsea(
  gse60450_contrast,
  c(
    HALLMARK_E2F_TARGETS = 1.75,
    HALLMARK_G2M_CHECKPOINT = 1.55,
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = -1.45,
    HALLMARK_P53_PATHWAY = 1.05,
    HALLMARK_APOPTOSIS = -1.1
  )
)

deg_16476 <- make_deg(
  gse16476_contrast,
  "MYCN-high tumor-like signal",
  "MYCN-low tumor-like signal",
  c(MYCN = 2.15, ALK = 1.35, MCL1 = 0.92, EPCAM = 0.72, KRT8 = 0.65, DDR1 = 0.58,
    SEMA5A = -0.7, MKI67 = 1.4, TOP2A = 1.32, VIM = 0.82, FN1 = 0.78, CDH1 = -0.45)
)

gsea_16476 <- make_gsea(
  gse16476_contrast,
  c(
    HALLMARK_E2F_TARGETS = 2.05,
    HALLMARK_G2M_CHECKPOINT = 1.85,
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = 1.35,
    HALLMARK_P53_PATHWAY = -1.2,
    HALLMARK_APOPTOSIS = -0.95
  )
)

write_dataset(
  "geo_gse60450",
  "Demo GSE60450 epithelial/lactation study",
  gse60450_contrast,
  "lactation survival-like signal",
  "baseline mammary cell state",
  deg_60450,
  gsea_60450
)

write_dataset(
  "geo_gse16476",
  "Demo GSE16476 MYCN neuroblastoma study",
  gse16476_contrast,
  "MYCN-high tumor-like signal",
  "MYCN-low tumor-like signal",
  deg_16476,
  gsea_16476
)

update_kb_claims(
  "GSE60450",
  c("MCL1", "EGF", "STAT5A", "EPCAM", "KRT8", "MYCN"),
  c("HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_APOPTOSIS"),
  gse60450_contrast
)

update_kb_claims(
  "GSE16476",
  c("MYCN", "ALK", "MCL1", "EPCAM", "KRT8", "SEMA5A"),
  c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_P53_PATHWAY"),
  gse16476_contrast
)

cat("Created comparison demo datasets:\n")
cat("- geo_gse60450\n")
cat("- geo_gse16476\n")
