geo_import_root <- function(data_root = "data") {
  file.path(data_root, "_geo_imports")
}

geo_import_cache_dir <- function(accession, data_root = "data") {
  file.path(geo_import_root(data_root), normalize_geo_accession(accession))
}

geo_dataset_id <- function(accession) {
  tolower(paste0("geo_", normalize_geo_accession(accession)))
}

pick_geo_expression_set <- function(gsets) {
  if (inherits(gsets, "ExpressionSet")) return(gsets)
  if (!is.list(gsets) || length(gsets) == 0) {
    stop("GEOquery did not return an ExpressionSet.", call. = FALSE)
  }
  sizes <- vapply(gsets, function(x) {
    if (inherits(x, "ExpressionSet")) {
      expr <- Biobase::exprs(x)
      nrow(expr) * ncol(expr)
    } else {
      0L
    }
  }, numeric(1))
  if (max(sizes) == 0) {
    stop(
      "GEO did not provide an analyzable series-matrix expression table for this accession. ",
      "It may require raw/supplementary-file processing before it can be imported.",
      call. = FALSE
    )
  }
  gsets[[which.max(sizes)]]
}

fetch_geo_expression_set <- function(accession, data_root = "data") {
  accession <- normalize_geo_accession(accession)
  if (!requireNamespace("GEOquery", quietly = TRUE) || !requireNamespace("Biobase", quietly = TRUE)) {
    stop("GEO import requires optional packages GEOquery and Biobase.", call. = FALSE)
  }
  cache_dir <- geo_import_cache_dir(accession, data_root)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(cache_dir, paste0(accession, "_expression_set.rds"))
  if (file.exists(cache_path)) return(readRDS(cache_path))
  gsets <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = TRUE)
  eset <- pick_geo_expression_set(gsets)
  saveRDS(eset, cache_path)
  eset
}

geo_pheno_summary <- function(eset, max_cols = 40) {
  pheno <- Biobase::pData(eset)
  cols <- utils::head(names(pheno), max_cols)
  summaries <- vapply(cols, function(col) {
    x <- pheno[[col]]
    vals <- utils::head(sort(table(x, useNA = "ifany"), decreasing = TRUE), 8)
    paste(paste(names(vals), vals, sep = " n="), collapse = "; ")
  }, character(1))
  paste(paste(cols, summaries, sep = ": "), collapse = "\n")
}

detect_expression_data_type <- function(expr, eset = NULL) {
  finite <- expr[is.finite(expr)]
  if (length(finite) == 0) return("normalized_expression")
  integer_like <- mean(abs(finite - round(finite)) < 1e-6) > 0.95
  if (integer_like && min(finite) >= 0 && stats::quantile(finite, 0.95, na.rm = TRUE) > 50) {
    return("raw_counts")
  }
  platform_text <- ""
  if (!is.null(eset)) {
    platform_text <- paste(
      Biobase::annotation(eset),
      Biobase::experimentData(eset)@title,
      collapse = " "
    )
  }
  if (grepl("affymetrix|illumina|agilent|array|microarray", platform_text, ignore.case = TRUE)) {
    return("microarray_or_processed_expression")
  }
  "normalized_expression"
}

candidate_group_columns <- function(pheno) {
  cols <- names(pheno)[vapply(pheno, function(x) {
    x <- as.character(x)
    n <- length(unique(stats::na.omit(x)))
    n >= 2 && n <= 8 && n < length(x) * 0.75
  }, logical(1))]
  preferred <- grep("group|condition|treat|status|phenotype|response|disease|mycn|stage|class", cols, ignore.case = TRUE, value = TRUE)
  unique(c(preferred, cols))
}

rule_based_geo_design <- function(accession, eset) {
  pheno <- Biobase::pData(eset)
  expr <- Biobase::exprs(eset)
  candidates <- candidate_group_columns(pheno)
  group_col <- if (length(candidates) > 0) candidates[[1]] else names(pheno)[[1]]
  levels <- names(sort(table(pheno[[group_col]], useNA = "no"), decreasing = TRUE))
  levels <- levels[nzchar(levels)]
  negative <- if (length(levels) >= 1) levels[[1]] else NA_character_
  positive <- if (length(levels) >= 2) levels[[2]] else NA_character_
  list(
    accession = normalize_geo_accession(accession),
    platform_type = detect_expression_data_type(expr, eset),
    expression_data_type = detect_expression_data_type(expr, eset),
    group_col = group_col,
    positive_group = positive,
    negative_group = negative,
    covariates = character(),
    sample_id_col = "geo_accession",
    patient_id_col = "geo_accession",
    display_sample_col = "title",
    default_timepoint_col = group_col,
    contrast_label = paste(positive, "vs", negative),
    rationale = "Rule-based proposal from phenotype columns. Review before running analysis."
  )
}

rule_based_geo_metadata_design <- function(accession, metadata) {
  pheno <- clean_column_names(metadata)
  candidates <- candidate_group_columns(pheno)
  group_col <- if (length(candidates) > 0) candidates[[1]] else names(pheno)[[1]]
  levels <- names(sort(table(pheno[[group_col]], useNA = "no"), decreasing = TRUE))
  levels <- levels[nzchar(levels)]
  negative <- if (length(levels) >= 1) levels[[1]] else NA_character_
  positive <- if (length(levels) >= 2) levels[[2]] else NA_character_
  list(
    accession = normalize_geo_accession(accession),
    platform_type = "supplementary_count_matrix",
    expression_data_type = "raw_counts",
    group_col = group_col,
    positive_group = positive,
    negative_group = negative,
    covariates = character(),
    sample_id_col = "geo_accession",
    patient_id_col = "geo_accession",
    display_sample_col = if ("title" %in% names(pheno)) "title" else "geo_accession",
    default_timepoint_col = group_col,
    contrast_label = paste(positive, "vs", negative),
    rationale = "Rule-based proposal from GEO sample metadata for a supplementary count-matrix route. Review before running analysis."
  )
}

geo_design_schema <- function() {
  list(
    type = "object",
    additionalProperties = FALSE,
    required = c(
      "platform_type", "expression_data_type", "group_col", "positive_group",
      "negative_group", "covariates", "sample_id_col", "patient_id_col",
      "display_sample_col", "default_timepoint_col", "contrast_label", "rationale"
    ),
    properties = list(
      platform_type = list(type = "string"),
      expression_data_type = list(type = "string", enum = c("raw_counts", "normalized_expression", "microarray_or_processed_expression")),
      group_col = list(type = "string"),
      positive_group = list(type = "string"),
      negative_group = list(type = "string"),
      covariates = list(type = "array", items = list(type = "string")),
      sample_id_col = list(type = "string"),
      patient_id_col = list(type = "string"),
      display_sample_col = list(type = "string"),
      default_timepoint_col = list(type = "string"),
      contrast_label = list(type = "string"),
      rationale = list(type = "string")
    )
  )
}

propose_geo_analysis_design <- function(accession, kb = NULL, api_key = "", model = "gpt-4o-mini", data_root = "data") {
  eset <- fetch_geo_expression_set(accession, data_root)
  fallback <- rule_based_geo_design(accession, eset)
  if (!openai_available(api_key) || is.null(kb)) return(fallback)

  prompt <- paste(
    "You are designing a reproducible GEO transcriptomics reanalysis for a Shiny RNA-seq Explorer.",
    "Use the paper methods and GEO phenotype columns to propose one primary contrast.",
    "Only choose columns that exist in the phenotype metadata. Prefer the original paper's main biological comparison.",
    "",
    "GEO accession:", normalize_geo_accession(accession),
    "",
    "Paper/method context:",
    kb_evidence_text(kb, max_chars = 35000),
    "",
    "Phenotype columns and values:",
    geo_pheno_summary(eset),
    "",
    "Rule-based fallback:",
    paste(names(fallback), unlist(fallback), sep = "=", collapse = "\n"),
    sep = "\n"
  )
  result <- tryCatch(openai_structured_json(prompt, geo_design_schema(), api_key, model), error = function(e) NULL)
  if (is.null(result)) return(fallback)
  result$accession <- normalize_geo_accession(accession)
  result$covariates <- unlist(result$covariates)
  result
}

design_to_table <- function(design) {
  if (is.null(design)) return(data.frame(Field = character(), Value = character()))
  data.frame(Field = names(design), Value = vapply(design, function(x) paste(x, collapse = ", "), character(1)), check.names = FALSE)
}

parse_design_from_inputs <- function(accession, platform_type, expression_data_type, group_col, positive_group, negative_group, covariates, sample_id_col, patient_id_col, display_sample_col, default_timepoint_col, contrast_label, rationale = "User-confirmed app design.") {
  covariate_values <- unique(trimws(unlist(strsplit(first_nonempty(covariates, ""), "[,;\\n\\r]+"))))
  covariate_values <- covariate_values[nzchar(covariate_values)]
  list(
    accession = normalize_geo_accession(accession),
    platform_type = platform_type,
    expression_data_type = expression_data_type,
    group_col = group_col,
    positive_group = positive_group,
    negative_group = negative_group,
    covariates = covariate_values,
    sample_id_col = sample_id_col,
    patient_id_col = patient_id_col,
    display_sample_col = display_sample_col,
    default_timepoint_col = default_timepoint_col,
    contrast_label = contrast_label,
    rationale = rationale
  )
}

parse_gene_symbol_from_feature <- function(row) {
  values <- as.character(row)
  names(values) <- names(row)
  symbol_cols <- grep("gene.?symbol|symbol|gene_assignment|gene.?name", names(values), ignore.case = TRUE, value = TRUE)
  for (col in symbol_cols) {
    x <- values[[col]]
    if (is.na(x) || !nzchar(x) || x == "---") next
    if (grepl("///|//", x)) {
      pieces <- trimws(unlist(strsplit(x, "///|//|;|,|\\s+")))
      pieces <- pieces[grepl("^[A-Za-z][A-Za-z0-9.-]{1,15}$", pieces)]
      if (length(pieces) > 0) return(toupper(pieces[[1]]))
    }
    pieces <- trimws(unlist(strsplit(x, "///|//|;|,|\\s+")))
    pieces <- pieces[grepl("^[A-Za-z][A-Za-z0-9.-]{1,15}$", pieces)]
    if (length(pieces) > 0) return(toupper(pieces[[1]]))
  }
  NA_character_
}

expression_to_gene_matrix <- function(eset) {
  expr <- Biobase::exprs(eset)
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    stop(
      "This GEO ExpressionSet has no expression features to analyze. ",
      "Choose a GEO series with a processed expression matrix, or add a custom raw-file preparation script.",
      call. = FALSE
    )
  }
  expr <- suppressWarnings(matrix(as.numeric(expr), nrow = nrow(expr), ncol = ncol(expr), dimnames = dimnames(expr)))
  keep_numeric <- rowSums(is.finite(expr)) >= 2
  expr <- expr[keep_numeric, , drop = FALSE]
  if (nrow(expr) == 0) {
    stop(
      "The GEO expression table could not be converted to numeric expression values.",
      call. = FALSE
    )
  }
  feature <- Biobase::fData(eset)
  genes <- if (nrow(feature) == nrow(expr) && ncol(feature) > 0) {
    apply(feature, 1, parse_gene_symbol_from_feature)
  } else {
    rep(NA_character_, nrow(expr))
  }
  genes[is.na(genes) | !nzchar(genes)] <- rownames(expr)[is.na(genes) | !nzchar(genes)]
  keep <- !is.na(genes) & nzchar(genes)
  expr <- expr[keep, , drop = FALSE]
  genes <- make.names(genes[keep], unique = FALSE)
  genes <- toupper(gsub("\\.", "-", genes))
  rownames(expr) <- genes
  expr <- expr[apply(expr, 1, function(x) stats::sd(x, na.rm = TRUE) > 0), , drop = FALSE]
  if (nrow(expr) == 0) {
    stop("No variable expression features remain after filtering.", call. = FALSE)
  }
  collapsed <- rowsum(expr, group = rownames(expr), reorder = FALSE)
  counts <- as.numeric(table(rownames(expr))[rownames(collapsed)])
  sweep(collapsed, 1, counts, "/")
}

prepare_geo_metadata <- function(eset, design) {
  pheno <- clean_column_names(Biobase::pData(eset))
  sample_ids <- if (design$sample_id_col %in% names(pheno)) as.character(pheno[[design$sample_id_col]]) else rownames(pheno)
  sample_ids[is.na(sample_ids) | !nzchar(sample_ids)] <- rownames(pheno)[is.na(sample_ids) | !nzchar(sample_ids)]
  pheno$app_sample_id <- make.names(sample_ids, unique = TRUE)
  if (!"geo_accession" %in% names(pheno)) pheno$geo_accession <- rownames(pheno)
  if (!design$patient_id_col %in% names(pheno)) pheno[[design$patient_id_col]] <- pheno$app_sample_id
  if (!design$display_sample_col %in% names(pheno)) pheno[[design$display_sample_col]] <- pheno$app_sample_id
  rownames(pheno) <- pheno$app_sample_id
  pheno
}

run_geo_deg <- function(expr_gene, metadata, design) {
  keep <- !is.na(metadata[[design$group_col]]) &
    metadata[[design$group_col]] %in% c(design$negative_group, design$positive_group)
  meta <- metadata[keep, , drop = FALSE]
  expr <- expr_gene[, meta$app_sample_id, drop = FALSE]
  group <- factor(meta[[design$group_col]], levels = c(design$negative_group, design$positive_group))
  design_df <- data.frame(group = group, meta[, intersect(design$covariates, names(meta)), drop = FALSE], check.names = TRUE)
  model <- stats::model.matrix(~ 0 + group, data = design_df)
  colnames(model) <- make.names(sub("^group", "", colnames(model)))

  if (design$expression_data_type == "raw_counts") {
    dge <- edgeR::DGEList(counts = round(pmax(expr, 0)))
    dge <- edgeR::calcNormFactors(dge)
    fit_expr <- limma::voom(dge, model, plot = FALSE)
  } else {
    fit_expr <- expr
  }

  fit <- limma::lmFit(fit_expr, model)
  contrast <- limma::makeContrasts(
    contrasts = paste0(make.names(design$positive_group), "-", make.names(design$negative_group)),
    levels = model
  )
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
  out <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
  out$gene <- rownames(out)
  out$contrast <- design$contrast_label
  out$log2FoldChange <- out$logFC
  out$pvalue <- out$P.Value
  out$padj <- out$adj.P.Val
  out$stat <- out$t
  out$baseMean <- rowMeans(expr, na.rm = TRUE)
  out$positive_group <- design$positive_group
  out$negative_group <- design$negative_group
  out$positive_n <- sum(group == design$positive_group)
  out$negative_n <- sum(group == design$negative_group)
  out[, c("contrast", "gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj", "positive_group", "negative_group", "positive_n", "negative_n")]
}

run_geo_gsea <- function(deg, expr_gene) {
  collections <- list(
    Hallmark = list(collection = "H", subcollection = NULL),
    Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
    `GO Biological Process` = list(collection = "C5", subcollection = "GO:BP"),
    WikiPathways = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS")
  )
  gsea_rows <- list()
  idx <- 1L
  for (collection_label in names(collections)) {
    genesets <- msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = collections[[collection_label]]$collection,
      subcollection = collections[[collection_label]]$subcollection
    )
    pathways <- split(genesets$gene_symbol, genesets$gs_name)
    ranks <- deg$stat
    names(ranks) <- deg$gene
    ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
    result <- fgsea::fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500, eps = 0)
    result <- as.data.frame(result)
    result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", character(1))
    result$contrast <- unique(deg$contrast)[[1]]
    result$collection <- collection_label
    result$msigdb_version <- unique(genesets$db_version)[[1]]
    gsea_rows[[idx]] <- result[, c("collection", "msigdb_version", "contrast", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")]
    idx <- idx + 1L
  }
  dplyr::bind_rows(gsea_rows)
}

compute_geo_pathway_scores <- function(expr_gene) {
  genesets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  pathways <- split(genesets$gene_symbol, genesets$gs_name)
  gene_z <- t(scale(t(expr_gene)))
  gene_z[!is.finite(gene_z)] <- 0
  scores <- vapply(pathways, function(genes) {
    present <- intersect(genes, rownames(gene_z))
    if (length(present) < 5) return(rep(NA_real_, ncol(gene_z)))
    colMeans(gene_z[present, , drop = FALSE])
  }, numeric(ncol(gene_z)))
  out <- t(scores)
  colnames(out) <- colnames(expr_gene)
  out
}

geo_supplementary_count_candidates <- function(kb = NULL, accession = NULL) {
  supp <- if (!is.null(kb) && !is.null(kb$geo$supplementary_files)) {
    kb$geo$supplementary_files
  } else {
    character()
  }
  if (length(supp) == 0) {
    return(data.frame(
      File = character(), URL = character(), Score = numeric(),
      Candidate_type = character(), Status = character(), check.names = FALSE
    ))
  }
  urls <- unique(as.character(supp))
  file <- basename(gsub("\\?.*$", "", urls))
  low <- tolower(file)
  data_like <- grepl("\\.(txt|tsv|csv|tab|dat)(\\.gz)?$", low)
  count_like <- grepl("count|counts|featurecounts|htseq|readcount|gene_count|raw", low)
  raw_like <- grepl("fastq|fq|sra|bam|sam|cel|idat", low)
  archive_like <- grepl("\\.(tar|tgz|zip|tar\\.gz)$", low)
  score <- as.numeric(data_like) * 40 + as.numeric(count_like) * 45 - as.numeric(raw_like) * 60 - as.numeric(archive_like) * 40
  status <- ifelse(score >= 60, "Likely count matrix", ifelse(data_like & !raw_like & !archive_like, "Possible table", "Not count-matrix candidate"))
  out <- data.frame(
    File = file,
    URL = urls,
    Score = score,
    Candidate_type = ifelse(count_like, "count-like supplementary table", "supplementary table"),
    Status = status,
    check.names = FALSE
  )
  out[order(out$Score, decreasing = TRUE), , drop = FALSE]
}

download_geo_supplementary_file <- function(url, accession, data_root = "data", timeout_sec = 120) {
  accession <- normalize_geo_accession(accession)
  if (!grepl("^https?://|^ftp://", url)) {
    stop("Supplementary candidate URL is not downloadable: ", url, call. = FALSE)
  }
  dir <- file.path(geo_import_cache_dir(accession, data_root), "supplementary")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(dir, safe_basename_from_url(url, paste0(accession, "_supplementary_count_candidate.txt.gz")))
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout_sec, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

read_geo_supplementary_table <- function(path) {
  open <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(open), add = TRUE)
  ext_path <- sub("\\.gz$", "", path, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(ext_path))
  if (identical(ext, "csv")) {
    utils::read.csv(open, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  } else {
    utils::read.delim(open, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  }
}

guess_count_gene_column <- function(df) {
  preferred <- grep("^(gene|gene_id|geneid|symbol|gene_symbol|ensembl|ensembl_id|feature|id)$", names(df), ignore.case = TRUE, value = TRUE)
  if (length(preferred) > 0) return(preferred[[1]])
  names(df)[[1]]
}

sample_name_matches <- function(table_names, metadata_samples) {
  table_clean <- make.names(table_names)
  meta_clean <- make.names(metadata_samples)
  table_names[table_clean %in% meta_clean]
}

prepare_geo_count_upload_inputs <- function(count_df, metadata) {
  gene_col <- guess_count_gene_column(count_df)
  sample_ids <- unique(c(
    rownames(metadata),
    if ("geo_accession" %in% names(metadata)) as.character(metadata$geo_accession),
    if ("app_sample_id" %in% names(metadata)) as.character(metadata$app_sample_id)
  ))
  matched_cols <- sample_name_matches(setdiff(names(count_df), gene_col), sample_ids)
  numeric_cols <- names(count_df)[vapply(count_df, function(x) sum(is.finite(suppressWarnings(as.numeric(x)))) >= 4, logical(1))]
  numeric_cols <- setdiff(numeric_cols, gene_col)
  sample_cols <- if (length(matched_cols) >= 4) matched_cols else numeric_cols
  if (length(sample_cols) < 4) {
    stop("No count-like supplementary table with at least four numeric sample columns could be parsed.", call. = FALSE)
  }
  expr_df <- count_df[, c(gene_col, sample_cols), drop = FALSE]
  names(expr_df)[[1]] <- "gene"

  metadata <- clean_column_names(metadata)
  metadata$app_sample_id <- make.names(if ("geo_accession" %in% names(metadata)) as.character(metadata$geo_accession) else rownames(metadata), unique = TRUE)
  sample_map <- data.frame(original = sample_cols, app_sample_id = make.names(sample_cols, unique = TRUE), check.names = FALSE)
  common <- intersect(metadata$app_sample_id, sample_map$app_sample_id)
  if (length(common) < 4 && "geo_accession" %in% names(metadata)) {
    common_original <- intersect(as.character(metadata$geo_accession), sample_cols)
    common <- make.names(common_original, unique = TRUE)
  }
  if (length(common) < 4) {
    metadata <- data.frame(app_sample_id = make.names(sample_cols, unique = TRUE), geo_accession = sample_cols, stringsAsFactors = FALSE)
  }
  colnames(expr_df)[-1] <- make.names(colnames(expr_df)[-1], unique = TRUE)
  list(expr_df = expr_df, metadata = metadata, gene_col = "gene")
}

infer_geo_count_method <- function(kb = NULL, preferred = "auto") {
  preferred <- first_nonempty(preferred, "auto")
  if (!identical(tolower(preferred), "auto")) return(select_count_de_method(preferred))
  method_text <- if (is.null(kb)) "" else paste(
    kb$extracted$statistical_methods,
    kb$extra_papers,
    kb_evidence_text(kb, max_chars = 20000),
    collapse = " "
  )
  select_count_de_method(infer_count_de_method_from_text(method_text))
}

fetch_geo_sample_metadata <- function(accession, data_root = "data") {
  accession <- normalize_geo_accession(accession)
  eset <- tryCatch(fetch_geo_expression_set(accession, data_root), error = function(e) NULL)
  if (!is.null(eset)) {
    pheno <- clean_column_names(Biobase::pData(eset))
    if (!"geo_accession" %in% names(pheno)) pheno$geo_accession <- rownames(pheno)
    pheno$app_sample_id <- make.names(as.character(pheno$geo_accession), unique = TRUE)
    return(pheno)
  }
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("GEO sample metadata retrieval requires GEOquery.", call. = FALSE)
  }
  cache_dir <- geo_import_cache_dir(accession, data_root)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(cache_dir, paste0(accession, "_sample_metadata.rds"))
  if (file.exists(cache_path)) return(readRDS(cache_path))
  gse <- GEOquery::getGEO(accession, GSEMatrix = FALSE)
  samples <- GEOquery::GSMList(gse)
  rows <- lapply(names(samples), function(id) {
    meta <- GEOquery::Meta(samples[[id]])
    vals <- vapply(meta, function(x) paste(unique(as.character(x)), collapse = "; "), character(1))
    data.frame(as.list(vals), check.names = FALSE, stringsAsFactors = FALSE) |>
      transform(geo_accession = id, app_sample_id = make.names(id))
  })
  pheno <- clean_column_names(dplyr::bind_rows(rows))
  rownames(pheno) <- pheno$app_sample_id
  saveRDS(pheno, cache_path)
  pheno
}

write_geo_supplementary_count_dataset <- function(accession, design, data_root = "data", kb = NULL, candidate_url = NULL) {
  accession <- normalize_geo_accession(accession)
  kb <- if (is.null(kb)) read_study_kb(accession, data_root) else kb
  candidates <- geo_supplementary_count_candidates(kb, accession)
  candidates <- candidates[candidates$Status %in% c("Likely count matrix", "Possible table"), , drop = FALSE]
  if (!is.null(candidate_url) && nzchar(candidate_url)) {
    candidates <- data.frame(File = basename(candidate_url), URL = candidate_url, Score = 100, Candidate_type = "user-selected", Status = "Likely count matrix", check.names = FALSE)
  }
  if (nrow(candidates) == 0) {
    stop("No GEO supplementary count-matrix candidates were detected for ", accession, ".", call. = FALSE)
  }
  metadata <- fetch_geo_sample_metadata(accession, data_root)
  errors <- character()
  for (i in seq_len(nrow(candidates))) {
    url <- candidates$URL[[i]]
    result <- tryCatch({
      local_path <- download_geo_supplementary_file(url, accession, data_root)
      count_df <- read_geo_supplementary_table(local_path)
      inputs <- prepare_geo_count_upload_inputs(count_df, metadata)
      settings <- list(
        dataset_id = paste0(geo_dataset_id(accession), "_counts"),
        display_name = paste("Imported", accession, "supplementary counts"),
        gene_col = inputs$gene_col,
        sample_col = "app_sample_id",
        group_col = design$group_col,
        positive_group = design$positive_group,
        negative_group = design$negative_group,
        patient_id_col = design$patient_id_col,
        display_sample_col = design$display_sample_col,
        contrast_label = design$contrast_label,
        count_method = infer_geo_count_method(kb)
      )
      result <- write_uploaded_app_dataset(inputs$expr_df, inputs$metadata, settings, data_root)
      dat <- load_rnaseq_dataset(result$dataset_id, data_root)
      report_path <- write_geo_reproducibility_report(dat, kb, result$dataset_dir)
      result$report_path <- report_path
      result$candidate_url <- url
      result$candidate_file <- candidates$File[[i]]
      result$count_method <- settings$count_method
      result$pipeline_script <- file.path(result$dataset_dir, "count_matrix_pipeline.R")
      result
    }, error = function(e) e)
    if (!inherits(result, "error")) return(result)
    errors <- c(errors, paste(candidates$File[[i]], conditionMessage(result), sep = ": "))
  }
  stop("Supplementary count import failed for all candidates. ", paste(errors, collapse = " | "), call. = FALSE)
}

write_geo_app_dataset <- function(accession, design, data_root = "data", kb = NULL) {
  eset <- fetch_geo_expression_set(accession, data_root)
  expr_gene <- expression_to_gene_matrix(eset)
  metadata <- prepare_geo_metadata(eset, design)
  common <- intersect(metadata$app_sample_id, colnames(expr_gene))
  if (length(common) == 0) {
    colnames(expr_gene) <- metadata$app_sample_id[seq_len(ncol(expr_gene))]
    common <- metadata$app_sample_id
  }
  metadata <- metadata[common, , drop = FALSE]
  expr_gene <- expr_gene[, common, drop = FALSE]

  pca_genes <- names(sort(matrixStats::rowVars(expr_gene), decreasing = TRUE))[seq_len(min(2000, nrow(expr_gene)))]
  pca <- stats::prcomp(t(expr_gene[pca_genes, , drop = FALSE]), scale. = FALSE)
  pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_variance <- setNames(pca_variance, paste0("PC", seq_along(pca_variance)))
  pca_df <- data.frame(app_sample_id = rownames(pca$x), sample = rownames(pca$x), pca$x[, seq_len(min(10, ncol(pca$x))), drop = FALSE], metadata[rownames(pca$x), , drop = FALSE], row.names = NULL, check.names = FALSE)
  pca_loadings <- dplyr::bind_rows(lapply(seq_len(min(10, ncol(pca$rotation))), function(i) {
    data.frame(PC = paste0("PC", i), gene = rownames(pca$rotation), loading = pca$rotation[, i], stringsAsFactors = FALSE)
  }))

  deg <- run_geo_deg(expr_gene, metadata, design)
  gsea <- run_geo_gsea(deg, expr_gene)
  pathway_scores <- compute_geo_pathway_scores(expr_gene)

  dataset_id <- geo_dataset_id(accession)
  dataset_dir <- file.path(data_root, dataset_id)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  dataset_info <- list(
    dataset_id = dataset_id,
    display_name = paste("Imported", normalize_geo_accession(accession)),
    sample_id_col = "app_sample_id",
    display_sample_col = design$display_sample_col,
    patient_id_col = design$patient_id_col,
    default_group_col = design$group_col,
    default_timepoint_col = design$default_timepoint_col,
    default_contrast = design$contrast_label,
    positive_group = design$positive_group,
    negative_group = design$negative_group,
    expression_unit = if (design$expression_data_type == "raw_counts") "voom logCPM for visualization" else "GEO processed expression",
    notes = paste("Imported from", normalize_geo_accession(accession), "using", design$expression_data_type, "pipeline.", design$rationale)
  )
  saveRDS(metadata, file.path(dataset_dir, "metadata.rds"))
  saveRDS(expr_gene, file.path(dataset_dir, "vsd_matrix.rds"))
  saveRDS(pca_df, file.path(dataset_dir, "pca_df.rds"))
  saveRDS(pca_variance, file.path(dataset_dir, "pca_variance.rds"))
  saveRDS(pca_loadings, file.path(dataset_dir, "pca_loadings.rds"))
  saveRDS(dataset_info, file.path(dataset_dir, "dataset_info.rds"))
  write.csv(deg, file.path(dataset_dir, "deg_results.csv"), row.names = FALSE)
  write.csv(gsea, file.path(dataset_dir, "gsea_hallmark.csv"), row.names = FALSE)
  write.csv(pathway_scores, file.path(dataset_dir, "pathway_scores.csv"))

  dat <- load_rnaseq_dataset(dataset_id, data_root)
  report_path <- write_geo_reproducibility_report(dat, kb, dataset_dir)
  list(dataset_id = dataset_id, dataset_dir = dataset_dir, report_path = report_path, rows = nrow(metadata), genes = nrow(expr_gene))
}

write_geo_reproducibility_report <- function(dat, kb, dataset_dir) {
  comparison <- compare_findings_against_publication(dat, kb)
  score <- score_reproducibility_from_comparison(comparison)
  triage <- summarize_relevant_different_results(dat, kb, comparison)
  gaps <- gap_analysis_report(dat, kb)
  suggestions <- suggest_ai_analyses(dat, kb)
  new_comparisons <- suggest_unexplored_comparisons(dat, kb)
  study_memory <- study_kb_summary_table(kb)
  prior_work <- prior_work_report_table(kb)
  analyzability <- geo_analyzability_route(kb, expression_set_status = "Report generated after app import")
  feasibility <- exact_reproduction_feasibility_table(kb)
  assessment <- reanalysis_assessment_table(kb)
  claims <- normalize_claim_table(kb$claims)
  virtual_lab <- virtual_lab_summary_table(kb)
  virtual_recs <- virtual_lab_recommendation_table(kb)
  citation_evidence <- citation_evidence_table(kb)
  report_path <- file.path(dataset_dir, "ai_reproducibility_report.html")
  escape <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x
  }
  html_table <- function(df, class = "data-table") {
    if (is.null(df) || nrow(df) == 0) return("<p class='empty'>No rows available.</p>")
    df[] <- lapply(df, function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      x
    })
    header <- paste0("<tr>", paste0("<th>", escape(names(df)), "</th>", collapse = ""), "</tr>")
    rows <- apply(df, 1, function(row) paste0("<tr>", paste0("<td>", escape(row), "</td>", collapse = ""), "</tr>"))
    paste0("<table class='", class, "'>", header, paste(rows, collapse = "\n"), "</table>")
  }
  metric_value <- function(df, metric) {
    if (is.null(df) || !all(c("Metric", "Value") %in% names(df))) return("")
    hit <- df$Value[df$Metric == metric]
    if (length(hit) == 0) "" else as.character(hit[[1]])
  }
  dataset_title <- first_nonempty(dat$info$display_name, kb$geo$title, kb$pubmed$title, kb$accession, "SeqSurf dataset")
  accession <- first_nonempty(kb$accession, dat$info$geo_accession, "")
  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M %Z")
  default_contrast <- first_nonempty(dat$info$default_contrast, "")
  sample_count <- if (!is.null(dat$metadata)) nrow(dat$metadata) else first_nonempty(kb$geo$sample_count, "")
  top_summary <- data.frame(
    Item = c("Dataset", "GEO accession", "Default contrast", "Samples", "Generated"),
    Value = c(dataset_title, accession, default_contrast, sample_count, generated_at),
    check.names = FALSE
  )
  route_label <- first_nonempty(
    if ("Route" %in% names(analyzability)) analyzability$Route[[1]] else "",
    if ("route" %in% names(analyzability)) analyzability$route[[1]] else "",
    ""
  )
  summary_cards <- paste0(
    "<div class='summary-grid'>",
    "<div class='summary-card'><span>Reproducibility score</span><strong>", escape(metric_value(score, "Post-analysis reproducibility score")), "</strong></div>",
    "<div class='summary-card'><span>Supported findings</span><strong>", escape(metric_value(score, "Supported findings")), "</strong></div>",
    "<div class='summary-card'><span>Direction mismatches</span><strong>", escape(metric_value(score, "Direction mismatches")), "</strong></div>",
    "<div class='summary-card'><span>Analyzability route</span><strong>", escape(route_label), "</strong></div>",
    "</div>"
  )
  section <- function(title, subtitle = NULL, body) {
    paste0(
      "<section class='report-section'>",
      "<h2>", escape(title), "</h2>",
      if (!is.null(subtitle)) paste0("<p class='section-note'>", escape(subtitle), "</p>") else "",
      body,
      "</section>"
    )
  }
  discovery_plan <- dplyr::bind_rows(
    data.frame(Source = "Discovery suggestions", suggestions, check.names = FALSE),
    data.frame(Source = "New comparisons", new_comparisons, check.names = FALSE)
  )
  if (nrow(discovery_plan) == 0) {
    discovery_plan <- data.frame(Source = "Discovery", Recommendation = "Run Validate first, then return to Discovery.", check.names = FALSE)
  }
  css <- "
    :root { --ink:#17253a; --muted:#5d6878; --line:#d9e2ea; --soft:#f5f8fb; --blue:#1268a8; --teal:#0b7c86; }
    * { box-sizing: border-box; }
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 0; color: var(--ink); background: #fff; line-height: 1.45; }
    main { max-width: 1180px; margin: 0 auto; padding: 36px; }
    header { border-bottom: 4px solid var(--blue); padding-bottom: 20px; margin-bottom: 24px; }
    h1 { font-size: 34px; margin: 0 0 8px; letter-spacing: 0; }
    h2 { font-size: 22px; margin: 0 0 8px; }
    .subtitle, .section-note { color: var(--muted); margin: 0 0 14px; }
    .summary-grid { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin: 18px 0 4px; }
    .summary-card { border: 1px solid var(--line); border-radius: 8px; padding: 14px 16px; background: var(--soft); min-height: 92px; }
    .summary-card span { display: block; color: var(--muted); font-weight: 700; font-size: 12px; text-transform: uppercase; margin-bottom: 8px; }
    .summary-card strong { display: block; font-size: 20px; }
    .report-section { border: 1px solid var(--line); border-left: 5px solid var(--teal); border-radius: 8px; padding: 20px; margin: 18px 0; break-inside: avoid; }
    table.data-table { width: 100%; border-collapse: collapse; margin: 12px 0; font-size: 13px; }
    .data-table th { text-align: left; background: #edf4fa; color: var(--ink); border-bottom: 1px solid var(--line); padding: 9px; vertical-align: top; }
    .data-table td { border-bottom: 1px solid #e8edf2; padding: 9px; vertical-align: top; }
    .data-table tr:nth-child(even) td { background: #fbfcfd; }
    .empty { color: var(--muted); font-style: italic; }
    .print-note { margin-top: 28px; padding: 12px 14px; background: #f8fafc; border: 1px solid var(--line); color: var(--muted); border-radius: 8px; }
    @media print { main { max-width: none; padding: 18mm; } .report-section { page-break-inside: avoid; } .summary-grid { grid-template-columns: repeat(2, 1fr); } }
  "
  html <- paste(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'>",
    paste0("<title>SeqSurf Reproducibility Report - ", escape(accession), "</title>"),
    paste0("<style>", css, "</style>"),
    "</head><body><main>",
    "<header>",
    "<p class='subtitle'>SeqSurf AI reproducibility report</p>",
    paste0("<h1>", escape(dataset_title), "</h1>"),
    paste0("<p class='subtitle'>", escape(accession), " | ", escape(default_contrast), "</p>"),
    summary_cards,
    "</header>",
    section("Study Memory", "Dataset identity, original paper context, related literature, extracted methods, biomarkers, pathways, and limitations.", paste0(html_table(top_summary), html_table(study_memory), html_table(prior_work))),
    section("Analyzability Route", "Whether the dataset can be analyzed inside SeqSurf, needs a count-matrix pipeline, or requires raw FASTQ/SRA handoff.", paste0(html_table(analyzability), html_table(feasibility), html_table(assessment))),
    section("Study Memory 2.0 Evidence", "Virtual-lab agent summaries and citation/evidence snippets used to support the generated memory.", paste0(html_table(virtual_lab), html_table(virtual_recs), html_table(citation_evidence))),
    section("Extracted Claims", "Structured gene, pathway, method, and result claims extracted from paper text, related papers, and user notes.", html_table(claims)),
    section("Reanalysis Results", "Computed app results available for comparison after importing or reanalyzing the dataset.", paste0(html_table(score), html_table(triage))),
    section("Agreement And Disagreement", "Where app results support, weaken, miss, or conflict with the original publication claims.", html_table(comparison)),
    section("Discovery Plan", "Suggested modern follow-up analyses, metadata-derived contrasts, and gaps worth addressing next.", paste0(html_table(gaps), html_table(discovery_plan))),
    "<p class='print-note'>This HTML report is print/PDF-ready. Use the browser print dialog and choose Save as PDF for submission or sharing.</p>",
    "</main></body></html>",
    sep = "\n"
  )
  writeLines(html, report_path)
  report_path
}
