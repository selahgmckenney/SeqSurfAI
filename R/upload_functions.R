read_uploaded_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

prepare_uploaded_expression <- function(expr_df, gene_col = NULL) {
  if (is.null(gene_col) || !nzchar(gene_col) || !gene_col %in% names(expr_df)) {
    gene_col <- names(expr_df)[[1]]
  }
  genes <- as.character(expr_df[[gene_col]])
  expr <- expr_df[, setdiff(names(expr_df), gene_col), drop = FALSE]
  expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(x)))
  keep_samples <- vapply(expr, function(x) sum(is.finite(x)) > 0, logical(1))
  expr <- expr[, keep_samples, drop = FALSE]
  expr <- as.matrix(expr)
  keep_genes <- !is.na(genes) & nzchar(genes) & apply(expr, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  expr <- expr[keep_genes, , drop = FALSE]
  genes <- toupper(make.names(genes[keep_genes], unique = FALSE))
  rownames(expr) <- genes
  collapsed <- rowsum(expr, group = rownames(expr), reorder = FALSE)
  counts <- as.numeric(table(rownames(expr))[rownames(collapsed)])
  sweep(collapsed, 1, counts, "/")
}

detect_uploaded_matrix_type <- function(expr) {
  finite <- as.numeric(expr[is.finite(expr)])
  if (length(finite) == 0) return("normalized_expression")
  integer_like <- mean(abs(finite - round(finite)) < 1e-6) > 0.95
  nonnegative <- min(finite, na.rm = TRUE) >= 0
  high_dynamic_range <- stats::quantile(finite, 0.95, na.rm = TRUE) > 50
  if (integer_like && nonnegative && high_dynamic_range) "raw_counts" else "normalized_expression"
}

transform_uploaded_for_visualization <- function(expr, matrix_type) {
  if (identical(matrix_type, "raw_counts") && requireNamespace("edgeR", quietly = TRUE)) {
    return(edgeR::cpm(round(pmax(expr, 0)), log = TRUE, prior.count = 1))
  }
  expr
}

select_count_de_method <- function(preferred = "auto") {
  preferred <- tolower(first_nonempty(preferred, "auto"))
  preferred <- gsub("-", "_", preferred, fixed = TRUE)
  if (preferred %in% c("deseq2", "edger", "limma_voom")) return(preferred)
  if (requireNamespace("DESeq2", quietly = TRUE)) return("deseq2")
  if (requireNamespace("edgeR", quietly = TRUE) && requireNamespace("limma", quietly = TRUE)) return("limma_voom")
  if (requireNamespace("edgeR", quietly = TRUE)) return("edger")
  "fallback"
}

infer_count_de_method_from_text <- function(text = "") {
  text <- tolower(paste(text, collapse = " "))
  if (grepl("deseq2|deseq", text)) return("deseq2")
  if (grepl("edger|edge r|edge-r", text)) return("edger")
  if (grepl("limma|voom", text)) return("limma_voom")
  "auto"
}

run_uploaded_deg <- function(expr, metadata, sample_col, group_col, positive_group, negative_group, contrast_label, de_method = "auto") {
  keep <- metadata[[group_col]] %in% c(positive_group, negative_group)
  meta <- metadata[keep, , drop = FALSE]
  samples <- intersect(as.character(meta[[sample_col]]), colnames(expr))
  meta <- meta[match(samples, as.character(meta[[sample_col]])), , drop = FALSE]
  expr <- expr[, samples, drop = FALSE]
  group <- as.character(meta[[group_col]])
  pos <- group == positive_group
  neg <- group == negative_group

  matrix_type <- detect_uploaded_matrix_type(expr)
  selected_method <- if (identical(matrix_type, "raw_counts")) select_count_de_method(de_method) else "normalized_expression"

  if (identical(matrix_type, "raw_counts") && identical(selected_method, "deseq2") && requireNamespace("DESeq2", quietly = TRUE)) {
    deseq_out <- tryCatch({
      group_factor <- factor(group, levels = c(negative_group, positive_group))
      dds_meta <- data.frame(group_factor = group_factor, row.names = samples, check.names = TRUE)
      dds <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(countData = round(pmax(expr, 0)), colData = dds_meta, design = ~ group_factor))
      dds <- suppressWarnings(DESeq2::DESeq(dds, quiet = TRUE))
      res <- as.data.frame(DESeq2::results(dds, contrast = c("group_factor", positive_group, negative_group)))
      out <- data.frame(
        contrast = contrast_label,
        gene = rownames(res),
        log2FoldChange = res$log2FoldChange,
        baseMean = res$baseMean,
        stat = res$stat,
        pvalue = res$pvalue,
        padj = res$padj,
        positive_group = positive_group,
        negative_group = negative_group,
        positive_n = sum(pos),
        negative_n = sum(neg),
        check.names = FALSE
      )
      attr(out, "analysis_method") <- "DESeq2 count-matrix automation"
      out
    }, error = function(e) NULL)
    if (!is.null(deseq_out)) return(deseq_out)
  }

  if (identical(matrix_type, "raw_counts") && identical(selected_method, "edger") && requireNamespace("edgeR", quietly = TRUE)) {
    edger_out <- tryCatch({
      group_factor <- factor(group, levels = c(negative_group, positive_group))
      dge <- edgeR::DGEList(counts = round(pmax(expr, 0)), group = group_factor)
      dge <- edgeR::calcNormFactors(dge)
      dge <- edgeR::estimateDisp(dge)
      et <- edgeR::exactTest(dge, pair = c(negative_group, positive_group))
      tab <- edgeR::topTags(et, n = Inf, sort.by = "none")$table
      out <- data.frame(
        contrast = contrast_label,
        gene = rownames(tab),
        log2FoldChange = tab$logFC,
        baseMean = rowMeans(expr, na.rm = TRUE),
        stat = tab$logCPM,
        pvalue = tab$PValue,
        padj = tab$FDR,
        positive_group = positive_group,
        negative_group = negative_group,
        positive_n = sum(pos),
        negative_n = sum(neg),
        check.names = FALSE
      )
      attr(out, "analysis_method") <- "edgeR exactTest count-matrix automation"
      out
    }, error = function(e) NULL)
    if (!is.null(edger_out)) return(edger_out)
  }

  if (identical(matrix_type, "raw_counts") && requireNamespace("edgeR", quietly = TRUE) && requireNamespace("limma", quietly = TRUE)) {
    group_factor <- factor(group, levels = c(negative_group, positive_group))
    model <- stats::model.matrix(~ 0 + group_factor)
    colnames(model) <- make.names(levels(group_factor))
    dge <- edgeR::DGEList(counts = round(pmax(expr, 0)))
    dge <- edgeR::calcNormFactors(dge)
    fit_expr <- limma::voom(dge, model, plot = FALSE)
    fit <- limma::lmFit(fit_expr, model)
    contrast <- limma::makeContrasts(
      contrasts = paste0(make.names(positive_group), "-", make.names(negative_group)),
      levels = model
    )
    fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
    out <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
    out$contrast <- contrast_label
    out$gene <- rownames(out)
    out$log2FoldChange <- out$logFC
    out$baseMean <- rowMeans(expr, na.rm = TRUE)
    out$stat <- out$t
    out$pvalue <- out$P.Value
    out$padj <- out$adj.P.Val
    out$positive_group <- positive_group
    out$negative_group <- negative_group
    out$positive_n <- sum(pos)
    out$negative_n <- sum(neg)
    out <- out[, c("contrast", "gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj", "positive_group", "negative_group", "positive_n", "negative_n")]
    attr(out, "analysis_method") <- "limma-voom count-matrix automation"
    return(out)
  }

  rows <- lapply(seq_len(nrow(expr)), function(i) {
    x_pos <- as.numeric(expr[i, pos])
    x_neg <- as.numeric(expr[i, neg])
    pvalue <- tryCatch(stats::t.test(x_pos, x_neg)$p.value, error = function(e) NA_real_)
    data.frame(
      contrast = contrast_label,
      gene = rownames(expr)[[i]],
      log2FoldChange = mean(x_pos, na.rm = TRUE) - mean(x_neg, na.rm = TRUE),
      baseMean = mean(expr[i, ], na.rm = TRUE),
      stat = NA_real_,
      pvalue = pvalue,
      positive_group = positive_group,
      negative_group = negative_group,
      positive_n = sum(pos),
      negative_n = sum(neg),
      check.names = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out <- out[, c("contrast", "gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj", "positive_group", "negative_group", "positive_n", "negative_n")]
  attr(out, "analysis_method") <- if (identical(matrix_type, "raw_counts")) "Welch t-test fallback on count-like matrix; install edgeR and limma for voom." else "Welch t-test on uploaded expression values"
  out
}

count_matrix_pipeline_script <- function(settings, method = "auto") {
  method <- select_count_de_method(method)
  group_col <- first_nonempty(settings$group_col, "condition")
  positive_group <- first_nonempty(settings$positive_group, "treated")
  negative_group <- first_nonempty(settings$negative_group, "control")
  contrast_label <- first_nonempty(settings$contrast_label, paste(positive_group, "vs", negative_group))
  paste(
    "# SeqSurf count-matrix pipeline handoff",
    "# Run from the generated dataset directory.",
    "",
    "counts <- read.csv('count_matrix_input.csv', check.names = FALSE)",
    "metadata <- read.csv('sample_metadata_input.csv', check.names = FALSE)",
    "rownames(counts) <- counts$gene",
    "counts$gene <- NULL",
    "counts <- as.matrix(counts)",
    "rownames(metadata) <- metadata$app_sample_id",
    "metadata <- metadata[colnames(counts), , drop = FALSE]",
    paste0("metadata$", group_col, " <- factor(metadata$", group_col, ", levels = c('", negative_group, "', '", positive_group, "'))"),
    "",
    if (identical(method, "deseq2")) paste(
      "suppressPackageStartupMessages(library(DESeq2))",
      paste0("dds <- DESeqDataSetFromMatrix(round(pmax(counts, 0)), metadata, design = ~ ", group_col, ")"),
      "dds <- DESeq(dds)",
      paste0("res <- as.data.frame(results(dds, contrast = c('", group_col, "', '", positive_group, "', '", negative_group, "')))"),
      "deg <- data.frame(",
      paste0("  contrast = '", contrast_label, "',"),
      "  gene = rownames(res),",
      "  log2FoldChange = res$log2FoldChange,",
      "  baseMean = res$baseMean,",
      "  stat = res$stat,",
      "  pvalue = res$pvalue,",
      "  padj = res$padj,",
      paste0("  positive_group = '", positive_group, "',"),
      paste0("  negative_group = '", negative_group, "',"),
      paste0("  positive_n = sum(metadata$", group_col, " == '", positive_group, "'),"),
      paste0("  negative_n = sum(metadata$", group_col, " == '", negative_group, "')"),
      ")",
      sep = "\n"
    ) else if (identical(method, "edger")) paste(
      "suppressPackageStartupMessages(library(edgeR))",
      paste0("group <- metadata$", group_col),
      "dge <- DGEList(counts = round(pmax(counts, 0)), group = group)",
      "dge <- calcNormFactors(dge)",
      "dge <- estimateDisp(dge)",
      paste0("et <- exactTest(dge, pair = c('", negative_group, "', '", positive_group, "'))"),
      "tab <- topTags(et, n = Inf, sort.by = 'none')$table",
      "deg <- data.frame(",
      paste0("  contrast = '", contrast_label, "',"),
      "  gene = rownames(tab), log2FoldChange = tab$logFC, baseMean = rowMeans(counts),",
      "  stat = tab$logCPM, pvalue = tab$PValue, padj = tab$FDR,",
      paste0("  positive_group = '", positive_group, "', negative_group = '", negative_group, "',"),
      paste0("  positive_n = sum(group == '", positive_group, "'), negative_n = sum(group == '", negative_group, "')"),
      ")",
      sep = "\n"
    ) else paste(
      "suppressPackageStartupMessages(library(edgeR))",
      "suppressPackageStartupMessages(library(limma))",
      paste0("group <- factor(metadata$", group_col, ", levels = c('", negative_group, "', '", positive_group, "'))"),
      "model <- model.matrix(~ 0 + group)",
      "colnames(model) <- make.names(levels(group))",
      "dge <- DGEList(counts = round(pmax(counts, 0)))",
      "dge <- calcNormFactors(dge)",
      "v <- voom(dge, model, plot = FALSE)",
      "fit <- lmFit(v, model)",
      paste0("contr <- makeContrasts(", make.names(positive_group), "-", make.names(negative_group), ", levels = model)"),
      "fit <- eBayes(contrasts.fit(fit, contr))",
      "tab <- topTable(fit, coef = 1, number = Inf, sort.by = 'none')",
      "deg <- data.frame(",
      paste0("  contrast = '", contrast_label, "',"),
      "  gene = rownames(tab), log2FoldChange = tab$logFC, baseMean = rowMeans(counts),",
      "  stat = tab$t, pvalue = tab$P.Value, padj = tab$adj.P.Val,",
      paste0("  positive_group = '", positive_group, "', negative_group = '", negative_group, "',"),
      paste0("  positive_n = sum(group == '", positive_group, "'), negative_n = sum(group == '", negative_group, "')"),
      ")",
      sep = "\n"
    ),
    "",
    "write.csv(deg, 'deg_results_recomputed.csv', row.names = FALSE)",
    sep = "\n"
  )
}

write_uploaded_app_dataset <- function(expr_df, metadata, settings, data_root = "data") {
  expr <- prepare_uploaded_expression(expr_df, settings$gene_col)
  metadata <- clean_column_names(metadata)
  validate(need(settings$sample_col %in% names(metadata), "Selected sample ID column is not in metadata."))
  validate(need(settings$group_col %in% names(metadata), "Selected group column is not in metadata."))
  metadata$app_sample_id <- make.names(as.character(metadata[[settings$sample_col]]), unique = TRUE)
  colnames(expr) <- make.names(colnames(expr), unique = TRUE)
  common <- intersect(metadata$app_sample_id, colnames(expr))
  validate(need(length(common) >= 4, "At least four matched expression/metadata samples are required."))
  metadata <- metadata[match(common, metadata$app_sample_id), , drop = FALSE]
  expr <- expr[, common, drop = FALSE]

  matrix_type <- detect_uploaded_matrix_type(expr)
  expr_for_app <- transform_uploaded_for_visualization(expr, matrix_type)
  pca_genes <- names(sort(matrixStats::rowVars(expr_for_app), decreasing = TRUE))[seq_len(min(2000, nrow(expr_for_app)))]
  pca <- stats::prcomp(t(expr_for_app[pca_genes, , drop = FALSE]), scale. = FALSE)
  pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_variance <- setNames(pca_variance, paste0("PC", seq_along(pca_variance)))
  pca_df <- data.frame(app_sample_id = rownames(pca$x), sample = rownames(pca$x), pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE], metadata[rownames(pca$x), , drop = FALSE], row.names = NULL, check.names = FALSE)
  pca_loadings <- dplyr::bind_rows(lapply(seq_len(min(10, ncol(pca$rotation))), function(i) {
    data.frame(PC = paste0("PC", i), gene = rownames(pca$rotation), loading = pca$rotation[, i], stringsAsFactors = FALSE)
  }))

  contrast_label <- first_nonempty(settings$contrast_label, paste(settings$positive_group, "vs", settings$negative_group))
  deg <- run_uploaded_deg(expr, metadata, "app_sample_id", settings$group_col, settings$positive_group, settings$negative_group, contrast_label, first_nonempty(settings$count_method, "auto"))
  analysis_method <- first_nonempty(attr(deg, "analysis_method"), "Uploaded expression workflow")
  pathway_scores <- matrix(nrow = 0, ncol = ncol(expr_for_app), dimnames = list(character(), colnames(expr_for_app)))
  gsea <- empty_gsea_table()

  dataset_id <- sanitize_filename(tolower(first_nonempty(settings$dataset_id, paste0("uploaded_", format(Sys.time(), "%Y%m%d_%H%M%S")))))
  dataset_dir <- file.path(data_root, dataset_id)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  dataset_info <- list(
    dataset_id = dataset_id,
    display_name = first_nonempty(settings$display_name, paste("Uploaded", dataset_id)),
    sample_id_col = "app_sample_id",
    display_sample_col = first_nonempty(settings$display_sample_col, "app_sample_id"),
    patient_id_col = first_nonempty(settings$patient_id_col, "app_sample_id"),
    default_group_col = settings$group_col,
    default_timepoint_col = settings$group_col,
    default_contrast = contrast_label,
    positive_group = settings$positive_group,
    negative_group = settings$negative_group,
    expression_unit = if (identical(matrix_type, "raw_counts")) "logCPM from uploaded count matrix" else "User-uploaded expression values",
    notes = paste("Uploaded expression matrix and metadata.", analysis_method, "Matrix type:", matrix_type)
  )

  saveRDS(metadata, file.path(dataset_dir, "metadata.rds"))
  saveRDS(expr_for_app, file.path(dataset_dir, "vsd_matrix.rds"))
  saveRDS(pca_df, file.path(dataset_dir, "pca_df.rds"))
  saveRDS(pca_variance, file.path(dataset_dir, "pca_variance.rds"))
  saveRDS(pca_loadings, file.path(dataset_dir, "pca_loadings.rds"))
  saveRDS(dataset_info, file.path(dataset_dir, "dataset_info.rds"))
  if (identical(matrix_type, "raw_counts")) {
    write.csv(data.frame(gene = rownames(expr), expr, check.names = FALSE), file.path(dataset_dir, "count_matrix_input.csv"), row.names = FALSE)
    write.csv(metadata, file.path(dataset_dir, "sample_metadata_input.csv"), row.names = FALSE)
    writeLines(count_matrix_pipeline_script(settings, first_nonempty(settings$count_method, "auto")), file.path(dataset_dir, "count_matrix_pipeline.R"))
  }
  write.csv(deg, file.path(dataset_dir, "deg_results.csv"), row.names = FALSE)
  write.csv(gsea, file.path(dataset_dir, "gsea_hallmark.csv"), row.names = FALSE)
  write.csv(pathway_scores, file.path(dataset_dir, "pathway_scores.csv"))
  list(dataset_id = dataset_id, dataset_dir = dataset_dir, rows = nrow(metadata), genes = nrow(expr))
}
