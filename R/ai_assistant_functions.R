normalize_geo_accession <- function(accession) {
  accession <- toupper(trimws(accession))
  accession <- gsub("[^A-Z0-9]", "", accession)
  accession
}

is_geo_accession <- function(accession) {
  grepl("^GS[EM][0-9]+$", normalize_geo_accession(accession))
}

assistant_kb_dir <- function(data_root = "data") {
  file.path(data_root, "ai_knowledge_base")
}

assistant_kb_path <- function(accession, data_root = "data") {
  file.path(assistant_kb_dir(data_root), paste0(normalize_geo_accession(accession), ".rds"))
}

assistant_document_dir <- function(accession, data_root = "data") {
  file.path(assistant_kb_dir(data_root), normalize_geo_accession(accession), "documents")
}

assistant_paperqa_dir <- function(accession, data_root = "data") {
  file.path(assistant_kb_dir(data_root), normalize_geo_accession(accession), "paperqa")
}

available_study_kbs <- function(data_root = "data") {
  root <- assistant_kb_dir(data_root)
  if (!dir.exists(root)) return(character())
  sub("\\.rds$", "", basename(list.files(root, pattern = "\\.rds$", full.names = FALSE)))
}

read_study_kb <- function(accession, data_root = "data") {
  path <- assistant_kb_path(accession, data_root)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

first_nonempty <- function(...) {
  values <- list(...)
  for (value in values) {
    if (!is.null(value) && length(value) > 0 && !is.na(value[[1]]) && nzchar(trimws(value[[1]]))) {
      return(as.character(value[[1]]))
    }
  }
  NA_character_
}

shorten_text <- function(x, max_chars = 420) {
  x <- trimws(gsub("\\s+", " ", first_nonempty(x, "")))
  if (is.na(x)) x <- ""
  if (nchar(x) <= max_chars) return(x)
  paste0(substr(x, 1, max_chars - 3), "...")
}

extract_first_sentence <- function(x) {
  x <- trimws(gsub("\\s+", " ", first_nonempty(x, "")))
  if (!nzchar(x)) return(NA_character_)
  sentence <- sub("^(.{20,}?[.!?])\\s+.*$", "\\1", x)
  shorten_text(sentence, 360)
}

extract_geo_metadata <- function(accession) {
  accession <- normalize_geo_accession(accession)
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Install optional package GEOquery to retrieve GEO metadata from inside the app.", call. = FALSE)
  }

  gse <- GEOquery::getGEO(accession, GSEMatrix = FALSE)
  meta <- GEOquery::Meta(gse)
  samples <- GEOquery::GSMList(gse)
  sample_ids <- names(samples)
  sample_organisms <- unique(unlist(lapply(samples, function(sample) GEOquery::Meta(sample)$organism_ch1), use.names = FALSE))
  sample_organisms <- sample_organisms[!is.na(sample_organisms) & nzchar(sample_organisms)]
  platform_ids <- names(GEOquery::GPLList(gse))
  supplementary_files <- unique(unlist(meta$supplementary_file, use.names = FALSE))
  supplementary_files <- supplementary_files[!is.na(supplementary_files) & nzchar(supplementary_files)]

  list(
    accession = accession,
    title = first_nonempty(meta$title),
    summary = first_nonempty(meta$summary),
    overall_design = first_nonempty(meta$overall_design),
    pubmed_id = first_nonempty(meta$pubmed_id),
    submission_date = first_nonempty(meta$submission_date),
    last_update_date = first_nonempty(meta$last_update_date),
    organism = first_nonempty(paste(sample_organisms, collapse = "; "), paste(unique(unlist(meta$sample_organism_ch1)), collapse = "; ")),
    sample_count = length(sample_ids),
    sample_accessions = sample_ids,
    platform_accessions = platform_ids,
    supplementary_files = supplementary_files,
    source = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession)
  )
}

extract_pubmed_metadata <- function(pubmed_id) {
  if (is.null(pubmed_id) || is.na(pubmed_id) || !nzchar(pubmed_id)) {
    return(list(pubmed_id = NA_character_, title = NA_character_, abstract = NA_character_, journal = NA_character_))
  }
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    return(list(
      pubmed_id = pubmed_id,
      title = NA_character_,
      abstract = NA_character_,
      journal = NA_character_,
      note = "Install optional package rentrez to retrieve PubMed title and abstract."
    ))
  }

  summary <- rentrez::entrez_summary(db = "pubmed", id = pubmed_id)
  abstract <- tryCatch(
    rentrez::entrez_fetch(db = "pubmed", id = pubmed_id, rettype = "abstract", retmode = "text"),
    error = function(e) NA_character_
  )
  links <- tryCatch(rentrez::entrez_link(dbfrom = "pubmed", db = "pmc", id = pubmed_id), error = function(e) NULL)
  pmc_id <- if (!is.null(links$links$pubmed_pmc) && length(links$links$pubmed_pmc) > 0) {
    paste0("PMC", links$links$pubmed_pmc[[1]])
  } else {
    NA_character_
  }

  list(
    pubmed_id = pubmed_id,
    title = first_nonempty(summary$title),
    abstract = first_nonempty(abstract),
    journal = first_nonempty(summary$fulljournalname, summary$source),
    publication_date = first_nonempty(summary$pubdate),
    pmc_id = pmc_id,
    source = paste0("https://pubmed.ncbi.nlm.nih.gov/", pubmed_id, "/")
  )
}

empty_related_publications <- function(status = "No related publications searched yet.") {
  data.frame(
    PubMed_ID = "",
    Title = "",
    Journal = "",
    Publication_date = "",
    Abstract = "",
    Source = "",
    Status = status,
    check.names = FALSE
  )
}

discover_related_pubmed_publications <- function(accession, max_papers = 12) {
  accession <- normalize_geo_accession(accession)
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    out <- empty_related_publications("Install optional package rentrez to search PubMed for papers mentioning this GEO accession.")
    return(out)
  }
  search <- tryCatch(
    rentrez::entrez_search(
      db = "pubmed",
      term = paste0("\"", accession, "\"[All Fields]"),
      retmax = max_papers,
      sort = "relevance"
    ),
    error = function(e) e
  )
  if (inherits(search, "error") || length(search$ids) == 0) {
    out <- empty_related_publications(paste("No PubMed papers found mentioning", accession))
    return(out)
  }
  rows <- lapply(search$ids, function(id) {
    meta <- extract_pubmed_metadata(id)
    data.frame(
      PubMed_ID = first_nonempty(meta$pubmed_id, id),
      Title = first_nonempty(meta$title, ""),
      Journal = first_nonempty(meta$journal, ""),
      Publication_date = first_nonempty(meta$publication_date, ""),
      Abstract = first_nonempty(meta$abstract, ""),
      Source = first_nonempty(meta$source, paste0("https://pubmed.ncbi.nlm.nih.gov/", id, "/")),
      Status = "PubMed paper mentioning GEO accession",
      check.names = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

related_publications_evidence_text <- function(related_publications, max_chars = 30000) {
  if (is.null(related_publications) || nrow(related_publications) == 0 || !"Title" %in% names(related_publications)) return("")
  related_publications <- related_publications[nzchar(related_publications$Title) | nzchar(related_publications$Abstract), , drop = FALSE]
  if (nrow(related_publications) == 0) return("")
  text <- paste(
    paste0(
      "PMID ", related_publications$PubMed_ID,
      "\nTitle: ", related_publications$Title,
      "\nJournal: ", related_publications$Journal,
      "\nDate: ", related_publications$Publication_date,
      "\nAbstract: ", related_publications$Abstract
    ),
    collapse = "\n\n"
  )
  substr(text, 1, max_chars)
}

build_study_knowledge_base <- function(accession, extra_papers = "", data_root = "data") {
  accession <- normalize_geo_accession(accession)
  if (!is_geo_accession(accession)) {
    stop("Enter a GEO Series or Sample accession such as GSE49711.", call. = FALSE)
  }

  geo <- extract_geo_metadata(accession)
  pubmed <- extract_pubmed_metadata(geo$pubmed_id)
  related_publications <- discover_related_pubmed_publications(accession)
  related_text <- related_publications_evidence_text(related_publications)
  evidence_text <- paste(geo$title, geo$summary, geo$overall_design, pubmed$title, pubmed$abstract, related_text, extra_papers, sep = "\n")
  claims <- extract_structured_claims(evidence_text, source = "GEO/PubMed/related papers/additional papers")

  kb <- list(
    accession = accession,
    created_at = Sys.time(),
    geo = geo,
    pubmed = pubmed,
    related_publications = related_publications,
    extra_papers = extra_papers,
    extracted = list(
      study_objective = first_nonempty(extract_first_sentence(pubmed$abstract), extract_first_sentence(geo$summary)),
      cohort_characteristics = shorten_text(first_nonempty(geo$overall_design, geo$summary), 650),
      sample_groups = shorten_text(first_nonempty(geo$overall_design), 500),
      statistical_methods = detect_methods(evidence_text),
      key_findings = shorten_text(first_nonempty(pubmed$abstract, geo$summary), 700),
      reported_pathways = detect_terms(evidence_text, pathway_dictionary()),
      reported_biomarkers = detect_biomarkers(evidence_text),
      study_limitations = detect_limitations(evidence_text)
    ),
    claims = claims,
    documents = data.frame(Source = character(), URL = character(), Local_path = character(), Text_path = character(), Status = character())
  )

  dir.create(assistant_kb_dir(data_root), recursive = TRUE, showWarnings = FALSE)
  saveRDS(kb, assistant_kb_path(accession, data_root))
  kb
}

get_or_build_study_knowledge_base <- function(accession, extra_papers = "", data_root = "data", refresh = FALSE) {
  accession <- normalize_geo_accession(accession)
  extra_papers <- trimws(first_nonempty(extra_papers, ""))
  if (!refresh && !nzchar(extra_papers)) {
    cached <- read_study_kb(accession, data_root)
    if (!is.null(cached)) {
      cached$loaded_from_cache <- TRUE
      return(cached)
    }
  }
  kb <- build_study_knowledge_base(accession, extra_papers, data_root)
  kb$loaded_from_cache <- FALSE
  kb
}

structured_claim_columns <- function() {
  c(
    "claim_type", "entity", "comparison", "reported_direction",
    "reported_significance", "evidence_sentence", "source", "confidence", "source_section"
  )
}

detect_source_section <- function(text) {
  if (grepl("abstract", text, ignore.case = TRUE)) return("Abstract")
  if (grepl("method", text, ignore.case = TRUE)) return("Methods")
  if (grepl("result", text, ignore.case = TRUE)) return("Results")
  if (grepl("discussion", text, ignore.case = TRUE)) return("Discussion")
  "Unclassified"
}

empty_claim_table <- function() {
  as.data.frame(stats::setNames(
    rep(list(character()), length(structured_claim_columns())),
    structured_claim_columns()
  ))
}

normalize_claim_table <- function(claims) {
  if (is.null(claims) || length(claims) == 0) return(empty_claim_table())
  if (is.list(claims) && is.null(dim(claims)) && length(claims) > 0 && is.list(claims[[1]])) {
    claims <- dplyr::bind_rows(claims)
  } else {
    claims <- as.data.frame(claims, stringsAsFactors = FALSE)
  }
  missing <- setdiff(structured_claim_columns(), names(claims))
  for (col in missing) claims[[col]] <- NA_character_
  claims <- claims[, structured_claim_columns(), drop = FALSE]
  claims[] <- lapply(claims, as.character)
  claims <- claims[!is.na(claims$entity) & nzchar(trimws(claims$entity)), , drop = FALSE]
  rownames(claims) <- NULL
  claims
}

extract_structured_claims <- function(text, source = "rule-based extraction") {
  biomarkers <- split_detected_terms(detect_biomarkers(text))
  pathways <- split_detected_terms(detect_terms(text, pathway_dictionary()))
  rows <- list()
  sec <- detect_source_section(source)

  for (gene in biomarkers) {
    rows[[length(rows) + 1]] <- data.frame(
      claim_type = "biomarker",
      entity = gene,
      comparison = NA_character_,
      reported_direction = "not specified",
      reported_significance = "reported or mentioned",
      evidence_sentence = extract_first_sentence(text),
      source = source,
      confidence = "low",
      source_section = sec,
      check.names = FALSE
    )
  }

  for (pathway in pathways) {
    rows[[length(rows) + 1]] <- data.frame(
      claim_type = "pathway",
      entity = pathway,
      comparison = NA_character_,
      reported_direction = "not specified",
      reported_significance = "reported or mentioned",
      evidence_sentence = extract_first_sentence(text),
      source = source,
      confidence = "low",
      source_section = sec,
      check.names = FALSE
    )
  }

  if (length(rows) == 0) return(empty_claim_table())
  normalize_claim_table(dplyr::bind_rows(rows))
}

claim_terms <- function(kb, claim_type = NULL) {
  claims <- normalize_claim_table(kb$claims)
  if (!is.null(claim_type)) claims <- claims[claims$claim_type == claim_type, , drop = FALSE]
  unique(claims$entity[nzchar(claims$entity)])
}

assistant_setup_status <- function(api_key = "") {
  optional_packages <- c("GEOquery", "rentrez", "pdftools", "httr2", "jsonlite", "survival", "immunedeconv")
  data.frame(
    Capability = c(
      "GEO metadata retrieval",
      "PubMed/PMC retrieval",
      "PDF text extraction",
      "Jina Reader URL/PDF ingestion",
      "LLM API transport",
      "LLM JSON parsing",
      "Survival analysis",
      "Immune deconvolution",
      "PaperQA Python package",
      "OpenAI API key"
    ),
    Requirement = c(optional_packages[1:3], "httr2; optional JINA_API_KEY for higher rate limits", optional_packages[4:7], "Python package paper-qa or paperqa", "OPENAI_API_KEY environment variable"),
    Status = c(
      ifelse(vapply(optional_packages[1:3], requireNamespace, logical(1), quietly = TRUE), "Ready", "Missing optional dependency"),
      ifelse(requireNamespace("httr2", quietly = TRUE), ifelse(nzchar(jina_api_key()), "Ready with API key", "Ready without API key"), "Missing optional dependency"),
      ifelse(vapply(optional_packages[4:7], requireNamespace, logical(1), quietly = TRUE), "Ready", "Missing optional dependency"),
      paperqa_status()$Status[[1]],
      ifelse(nzchar(openai_api_key(api_key)), "Ready", "Missing API key")
    ),
    check.names = FALSE
  )
}

jina_api_key <- function(api_key = "") {
  api_key <- trimws(first_nonempty(api_key, ""))
  if (is.na(api_key)) api_key <- ""
  if (nzchar(api_key)) return(api_key)
  env_key <- trimws(Sys.getenv("JINA_API_KEY", unset = ""))
  if (is.na(env_key)) "" else env_key
}

jina_reader_url <- function(url) {
  url <- trimws(first_nonempty(url, ""))
  if (!nzchar(url)) return("")
  paste0("https://r.jina.ai/", url)
}

jina_reader_available <- function() {
  requireNamespace("httr2", quietly = TRUE)
}

fetch_url_with_jina_reader <- function(url, api_key = "", timeout_sec = 90) {
  if (!jina_reader_available()) {
    stop("Jina Reader ingestion requires optional package httr2.", call. = FALSE)
  }
  req <- httr2::request(jina_reader_url(url)) |>
    httr2::req_headers(
      Accept = "text/markdown",
      `X-No-Cache` = "true",
      `X-Timeout` = as.character(timeout_sec),
      `X-Remove-Selector` = "nav,footer,header,.sidebar,#references"
    ) |>
    httr2::req_timeout(timeout_sec) |>
    httr2::req_error(is_error = function(resp) FALSE)
  key <- jina_api_key(api_key)
  if (nzchar(key)) {
    req <- httr2::req_auth_bearer_token(req, key)
  }
  response <- httr2::req_perform(req)
  body <- httr2::resp_body_string(response)
  if (httr2::resp_status(response) >= 300) {
    stop(paste("Jina Reader failed with status", httr2::resp_status(response), shorten_text(body, 280)), call. = FALSE)
  }
  body
}

paperqa_status <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(data.frame(Status = "Optional: install reticulate + paper-qa for citation-grounded paper chat.", Detail = "reticulate not installed", check.names = FALSE))
  }
  available <- FALSE
  detail <- "Python package not found"
  for (module in c("paperqa", "paper_qa")) {
    found <- tryCatch(reticulate::py_module_available(module), error = function(e) FALSE)
    if (isTRUE(found)) {
      available <- TRUE
      detail <- paste("Python module available:", module)
      break
    }
  }
  data.frame(
    Status = if (available) "Ready" else "Optional dependency missing",
    Detail = detail,
    check.names = FALSE
  )
}

discover_document_sources <- function(kb) {
  if (is.null(kb)) {
    return(data.frame(Source = character(), URL = character(), Status = character()))
  }
  urls <- character()
  labels <- character()

  if (!is.null(kb$geo$supplementary_files) && length(kb$geo$supplementary_files) > 0) {
    urls <- c(urls, kb$geo$supplementary_files)
    labels <- c(labels, rep("GEO supplementary file", length(kb$geo$supplementary_files)))
  }
  if (!is.null(kb$pubmed$source) && !is.na(kb$pubmed$source)) {
    urls <- c(urls, kb$pubmed$source)
    labels <- c(labels, "PubMed record")
  }
  if (!is.null(kb$pubmed$pmc_id) && !is.na(kb$pubmed$pmc_id) && nzchar(kb$pubmed$pmc_id)) {
    urls <- c(urls, paste0("https://pmc.ncbi.nlm.nih.gov/articles/", kb$pubmed$pmc_id, "/pdf/"))
    labels <- c(labels, "PMC full-text PDF")
  }
  if (!is.null(kb$related_publications) && nrow(kb$related_publications) > 0 && "Source" %in% names(kb$related_publications)) {
    related <- kb$related_publications
    related <- related[nzchar(related$Source) & grepl("^https?://", related$Source), , drop = FALSE]
    if (nrow(related) > 0) {
      urls <- c(urls, related$Source)
      labels <- c(labels, paste("Related paper:", vapply(related$Title, shorten_text, character(1), max_chars = 90)))
    }
  }

  if (length(urls) == 0) {
    return(data.frame(Source = "No automatic document sources found", URL = "", Status = "Paste paper text or add PDFs manually."))
  }

  data.frame(Source = labels, URL = urls, Status = "Discovered", check.names = FALSE)
}

is_paper_document_url <- function(url) {
  if (!nzchar(url)) return(FALSE)
  url_low <- tolower(url)
  paper_like <- grepl("\\.pdf($|[?])|\\.txt($|[?])|\\.html?($|[?])|\\.xml($|[?])|pmc\\.ncbi|pubmed\\.ncbi", url_low)
  archive_like <- grepl("\\.tar($|[?])|\\.tar\\.gz($|[?])|\\.tgz($|[?])|\\.gz($|[?])|\\.zip($|[?])|\\.cel($|[?])|sra|fastq|counts", url_low)
  paper_like && !archive_like
}

safe_basename_from_url <- function(url, fallback) {
  name <- basename(strsplit(url, "\\?", fixed = FALSE)[[1]][[1]])
  name <- sanitize_filename(first_nonempty(name, fallback))
  if (!nzchar(name)) fallback else name
}

extract_document_text <- function(local_path) {
  ext <- tolower(tools::file_ext(local_path))
  if (ext == "pdf") {
    if (!requireNamespace("pdftools", quietly = TRUE)) {
      return("PDF text extraction requires optional package pdftools.")
    }
    return(paste(pdftools::pdf_text(local_path), collapse = "\n"))
  }
  if (ext %in% c("txt", "csv", "tsv", "soft", "html", "htm", "xml")) {
    return(paste(readLines(local_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  }
  paste("Text extraction is not configured for file type:", ext)
}

write_text_document <- function(text, text_dir, base_name) {
  base_name <- sanitize_filename(first_nonempty(base_name, "document"))
  text_path <- file.path(text_dir, paste0(tools::file_path_sans_ext(base_name), ".txt"))
  writeLines(text, text_path, useBytes = TRUE)
  text_path
}

ingest_document_source <- function(source, url, doc_dir, text_dir, timeout_sec = 45) {
  if (!nzchar(url) || !grepl("^https?://|^ftp://", url)) {
    return(data.frame(Source = source, URL = url, Local_path = "", Text_path = "", Ingestion = "", Status = "No retrievable URL", check.names = FALSE))
  }

  if (grepl("^https?://", url) && jina_reader_available()) {
    jina_result <- tryCatch({
      text <- fetch_url_with_jina_reader(url, timeout_sec = max(timeout_sec, 90))
      text_path <- write_text_document(text, text_dir, paste0(safe_basename_from_url(url, "jina_document"), "_jina"))
      data.frame(Source = source, URL = url, Local_path = "", Text_path = text_path, Ingestion = "Jina Reader", Status = "Jina Reader indexed", check.names = FALSE)
    }, error = function(e) e)
    if (!inherits(jina_result, "error")) return(jina_result)
    jina_status <- paste("Jina failed:", conditionMessage(jina_result))
  } else {
    jina_status <- "Jina Reader unavailable"
  }

  local_path <- file.path(doc_dir, safe_basename_from_url(url, "document"))
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout_sec, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  direct_result <- tryCatch({
    utils::download.file(url, local_path, mode = "wb", quiet = TRUE)
    text <- extract_document_text(local_path)
    text_path <- write_text_document(text, text_dir, basename(local_path))
    data.frame(Source = source, URL = url, Local_path = if (file.exists(local_path)) local_path else "", Text_path = text_path, Ingestion = "Direct download", Status = paste("Downloaded and indexed;", jina_status), check.names = FALSE)
  }, error = function(e) e)
  if (!inherits(direct_result, "error")) return(direct_result)

  data.frame(
    Source = source,
    URL = url,
    Local_path = if (file.exists(local_path)) local_path else "",
    Text_path = "",
    Ingestion = "Failed",
    Status = paste(jina_status, "Direct download failed:", conditionMessage(direct_result)),
    check.names = FALSE
  )
}

prepare_paperqa_corpus <- function(kb, data_root = "data") {
  if (is.null(kb)) stop("Train or select a study knowledge base first.", call. = FALSE)
  accession <- normalize_geo_accession(kb$accession)
  paperqa_dir <- assistant_paperqa_dir(accession, data_root)
  corpus_dir <- file.path(paperqa_dir, "corpus")
  dir.create(assistant_kb_dir(data_root), recursive = TRUE, showWarnings = FALSE)
  dir.create(corpus_dir, recursive = TRUE, showWarnings = FALSE)

  docs <- if (!is.null(kb$documents) && nrow(kb$documents) > 0 && "Text_path" %in% names(kb$documents)) {
    kb$documents
  } else {
    data.frame(Source = character(), URL = character(), Text_path = character(), check.names = FALSE)
  }
  docs <- docs[nzchar(docs$Text_path) & file.exists(docs$Text_path), , drop = FALSE]
  copied <- lapply(seq_len(nrow(docs)), function(i) {
    dest <- file.path(corpus_dir, sanitize_filename(paste0(i, "_", tools::file_path_sans_ext(basename(docs$Text_path[[i]])), ".txt")))
    file.copy(docs$Text_path[[i]], dest, overwrite = TRUE)
    data.frame(
      Source = first_nonempty(docs$Source[[i]], paste("Document", i)),
      URL = first_nonempty(docs$URL[[i]], ""),
      Text_file = dest,
      Status = "Ready for PaperQA indexing",
      check.names = FALSE
    )
  })

  if (length(copied) == 0) {
    manifest <- data.frame(Source = "No ingested paper text yet", URL = "", Text_file = "", Status = "Fetch evidence with Jina Reader first", check.names = FALSE)
  } else {
    manifest <- dplyr::bind_rows(copied)
  }
  manifest_path <- file.path(paperqa_dir, "paperqa_manifest.csv")
  write.csv(manifest, manifest_path, row.names = FALSE)
  kb$paperqa <- list(
    prepared_at = Sys.time(),
    corpus_dir = corpus_dir,
    manifest_path = manifest_path,
    status = paperqa_status()
  )
  saveRDS(kb, assistant_kb_path(accession, data_root))
  kb
}

detect_document_sections <- function(text, max_chars = 1400) {
  text <- gsub("\\r", "\n", first_nonempty(text, ""))
  text <- gsub("-\\s*\n\\s*", "", text)
  text_flat <- gsub("\\s+", " ", text)
  lines <- trimws(unlist(strsplit(text, "\n+")))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    return(data.frame(Section = character(), Excerpt = character(), Status = character()))
  }
  normalized <- tolower(gsub("^\\s*([0-9]+\\.?\\s*)+", "", lines))
  normalized <- gsub("[[:punct:]]+$", "", normalized)
  section_patterns <- list(
    Abstract = "^(abstract|summary)$",
    Methods = "^(materials and methods|methods|methodology|experimental procedures|statistical analysis|bioinformatics analysis|data processing|differential expression analysis|rna-seq analysis|microarray analysis)$",
    Results = "^(results|findings|main results)$",
    `Tables/Figures` = "^(tables?|figures?|figure legends?|table legends?)$",
    Discussion = "^(discussion|conclusions?|interpretation)$",
    Limitations = "^(limitations|study limitations|limitations of the study)$",
    `Data availability` = "^(data availability|availability of data|data accessibility|accession numbers?)",
    Supplement = "^(supplementary|supplemental|supporting information)"
  )
  fallback_patterns <- list(
    Abstract = "objective|purpose|background|we investigated|we examined|this study",
    Methods = "rna-seq|microarray|sequenc|affymetrix|illumina|normaliz|differential expression|deseq2|edger|limma|voom|gsea|statistical analysis|false discovery",
    Results = "we found|identified|significant|differentially expressed|upregulated|downregulated|enriched|associated with|correlated",
    `Tables/Figures` = "figure [0-9]|table [0-9]|supplementary figure|supplementary table",
    Discussion = "suggest|indicate|consistent with|these results|in conclusion|we conclude",
    Limitations = "limitation|limited by|small sample|retrospective|validation|further study|not available|missing",
    `Data availability` = "geo|gene expression omnibus|accession|deposited|available at|sra|arrayexpress",
    Supplement = "supplementary|supplemental|supporting information|additional file"
  )
  sentence_excerpt <- function(pattern) {
    sentences <- unlist(strsplit(text_flat, "(?<=[.!?])\\s+", perl = TRUE))
    keep <- grepl(pattern, sentences, ignore.case = TRUE)
    hits <- unique(sentences[keep])
    hits <- hits[nchar(hits) >= 35]
    shorten_text(paste(utils::head(hits, 5), collapse = " "), max_chars)
  }
  hits <- lapply(names(section_patterns), function(section) {
    idx <- grep(section_patterns[[section]], normalized, perl = TRUE)
    if (length(idx) == 0) {
      fallback <- sentence_excerpt(fallback_patterns[[section]])
      return(data.frame(
        Section = section,
        Excerpt = fallback,
        Status = if (nzchar(fallback)) "Inferred from sentences" else "Not detected",
        check.names = FALSE
      ))
    }
    start <- idx[[1]]
    next_heading <- length(lines) + 1
    for (pattern in unlist(section_patterns)) {
      later <- grep(pattern, normalized, perl = TRUE)
      later <- later[later > start]
      if (length(later) > 0) next_heading <- min(next_heading, later[[1]])
    }
    end <- min(next_heading - 1, start + 110, length(lines))
    excerpt <- if (start + 1 <= end) paste(lines[(start + 1):end], collapse = " ") else ""
    excerpt <- shorten_text(excerpt, max_chars)
    if (!nzchar(excerpt)) excerpt <- sentence_excerpt(fallback_patterns[[section]])
    data.frame(Section = section, Excerpt = excerpt, Status = if (nzchar(excerpt)) "Detected" else "Heading only", check.names = FALSE)
  })
  dplyr::bind_rows(hits)
}

section_claim_candidate_table <- function(text, max_sentences = 18) {
  text <- gsub("\\s+", " ", first_nonempty(text, ""))
  if (!nzchar(text)) {
    return(data.frame(Claim_candidate = character(), Evidence_type = character(), Status = character(), check.names = FALSE))
  }
  sentences <- unlist(strsplit(text, "(?<=[.!?])\\s+", perl = TRUE))
  patterns <- c(
    gene = "\\b[A-Z0-9]{2,}[A-Z0-9-]*\\b",
    method = "deseq2|edger|limma|voom|microarray|rna-seq|sequenc|normaliz|model|covariate",
    significance = "significant|p\\s*[<=>]|adjusted|fdr|q-value|fold change|differential|log2",
    pathway = "pathway|gsea|gene set|enrichment|ontology|reactome|kegg|hallmark",
    outcome = "survival|prognosis|response|risk|recurrence|metastasis",
    limitation = "limitation|limited|missing|unavailable|not available|small sample|further study"
  )
  keep <- vapply(sentences, function(sentence) {
    sentence_low <- tolower(sentence)
    nchar(sentence) >= 40 && nchar(sentence) <= 450 &&
      (grepl(patterns[["significance"]], sentence_low) ||
         grepl(patterns[["pathway"]], sentence_low) ||
         grepl(patterns[["outcome"]], sentence_low))
  }, logical(1))
  candidates <- unique(sentences[keep])
  candidates <- utils::head(candidates, max_sentences)
  if (length(candidates) == 0) {
    return(data.frame(Claim_candidate = "No claim-like sentences detected.", Evidence_type = "Missing", Status = "Needs LLM/user extraction", check.names = FALSE))
  }
  evidence_type <- vapply(candidates, function(sentence) {
    sentence_low <- tolower(sentence)
    paste(names(patterns)[vapply(patterns, grepl, logical(1), x = sentence_low)], collapse = "; ")
  }, character(1))
  data.frame(Claim_candidate = candidates, Evidence_type = evidence_type, Status = "Candidate", check.names = FALSE)
}

summarize_document_sections <- function(kb) {
  if (is.null(kb)) {
    return(data.frame(Section = "Paper sections", Excerpt = "Train or select a study knowledge base first.", Status = "Not ready", check.names = FALSE))
  }
  if (!is.null(kb$document_sections) && nrow(kb$document_sections) > 0) {
    return(kb$document_sections)
  }
  fallback <- kb_evidence_text(kb, max_chars = 30000)
  sections <- detect_document_sections(fallback)
  if (nrow(sections) == 0 || all(sections$Status == "Not detected")) {
    return(data.frame(
      Section = c("Methods", "Results", "Limitations", "Supplement"),
      Excerpt = c(
        first_nonempty(kb$extracted$statistical_methods, "No methods section detected."),
        first_nonempty(kb$extracted$key_findings, "No results section detected."),
        first_nonempty(kb$extracted$study_limitations, "No limitations section detected."),
        paste(length(if (is.null(kb$geo$supplementary_files)) character() else kb$geo$supplementary_files), "GEO supplementary links parsed.")
      ),
      Status = c("Fallback", "Fallback", "Fallback", "Fallback"),
      check.names = FALSE
    ))
  }
  sections
}

fetch_and_ingest_documents <- function(kb, data_root = "data", max_docs = 4, timeout_sec = 45) {
  if (is.null(kb)) stop("Train or select a study knowledge base first.", call. = FALSE)
  sources <- discover_document_sources(kb)
  sources$paper_like <- vapply(sources$URL, is_paper_document_url, logical(1))
  if (any(sources$paper_like)) {
    sources <- sources[sources$paper_like, , drop = FALSE]
  }
  sources <- utils::head(sources, max_docs)
  accession <- normalize_geo_accession(kb$accession)
  doc_dir <- assistant_document_dir(accession, data_root)
  text_dir <- file.path(doc_dir, "text")
  dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)

  ingested <- lapply(seq_len(nrow(sources)), function(i) {
    if (!isTRUE(sources$paper_like[[i]])) {
      return(data.frame(Source = sources$Source[[i]], URL = sources$URL[[i]], Local_path = "", Text_path = "", Ingestion = "", Status = "Skipped non-paper archive or data file", check.names = FALSE))
    }
    ingest_document_source(
      source = sources$Source[[i]],
      url = sources$URL[[i]],
      doc_dir = doc_dir,
      text_dir = text_dir,
      timeout_sec = timeout_sec
    )
  })

  docs <- dplyr::bind_rows(ingested)
  kb$documents <- docs
  doc_text <- paste(vapply(docs$Text_path, function(path) {
    if (nzchar(path) && file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n") else ""
  }, character(1)), collapse = "\n")
  if (nzchar(doc_text)) {
    sections <- detect_document_sections(doc_text)
    kb$document_sections <- summarize_document_sections(list(
      document_sections = sections,
      geo = kb$geo,
      pubmed = kb$pubmed,
      extracted = kb$extracted,
      documents = docs
    ))
    kb$claim_candidates <- section_claim_candidate_table(paste(sections$Excerpt, collapse = " "))
    kb$claims <- dplyr::bind_rows(kb$claims, extract_structured_claims(doc_text, source = "downloaded documents"))
    kb$claims <- unique(kb$claims)
  }
  kb <- prepare_paperqa_corpus(kb, data_root)
  saveRDS(kb, assistant_kb_path(accession, data_root))
  kb
}

openai_api_key <- function(api_key = "") {
  api_key <- trimws(first_nonempty(api_key, ""))
  if (is.na(api_key)) api_key <- ""
  if (nzchar(api_key)) return(api_key)
  env_key <- trimws(Sys.getenv("OPENAI_API_KEY", unset = ""))
  if (is.na(env_key)) "" else env_key
}

openai_available <- function(api_key = "") {
  nzchar(openai_api_key(api_key)) &&
    requireNamespace("httr2", quietly = TRUE) &&
    requireNamespace("jsonlite", quietly = TRUE)
}

extract_response_text <- function(body) {
  if (!is.null(body$output_text) && nzchar(body$output_text)) return(body$output_text)
  pieces <- character()
  if (!is.null(body$output)) {
    for (item in body$output) {
      if (!is.null(item$content)) {
        for (content in item$content) {
          if (!is.null(content$text)) pieces <- c(pieces, content$text)
        }
      }
    }
  }
  paste(pieces, collapse = "\n")
}

publication_extraction_schema <- function() {
  list(
    type = "object",
    additionalProperties = FALSE,
    required = c(
      "study_objective", "cohort_characteristics", "sample_groups", "statistical_methods",
      "key_findings", "reported_pathways", "reported_biomarkers", "study_limitations", "claims"
    ),
    properties = list(
      study_objective = list(type = "string"),
      cohort_characteristics = list(type = "string"),
      sample_groups = list(type = "string"),
      statistical_methods = list(type = "string"),
      key_findings = list(type = "string"),
      reported_pathways = list(type = "array", items = list(type = "string")),
      reported_biomarkers = list(type = "array", items = list(type = "string")),
      study_limitations = list(type = "string"),
      claims = list(
        type = "array",
        items = list(
          type = "object",
          additionalProperties = FALSE,
          required = structured_claim_columns(),
          properties = list(
            claim_type = list(type = "string", enum = c("gene", "biomarker", "pathway", "method", "cohort", "limitation", "other")),
            entity = list(type = "string"),
            comparison = list(type = "string"),
            reported_direction = list(type = "string", enum = c("up", "down", "positive", "negative", "enriched", "depleted", "higher", "lower", "not specified")),
            reported_significance = list(type = "string"),
            evidence_sentence = list(type = "string"),
            source = list(type = "string"),
            confidence = list(type = "string", enum = c("high", "medium", "low"))
          )
        )
      )
    )
  )
}

seqsurf_virtual_lab_schema <- function() {
  agent_schema <- list(
    type = "object",
    additionalProperties = FALSE,
    required = c("role", "summary", "evidence", "recommendations", "uncertainties"),
    properties = list(
      role = list(type = "string"),
      summary = list(type = "string"),
      evidence = list(type = "array", items = list(type = "string")),
      recommendations = list(type = "array", items = list(type = "string")),
      uncertainties = list(type = "array", items = list(type = "string"))
    )
  )
  list(
    type = "object",
    additionalProperties = FALSE,
    required = c(
      "literature_curator", "methods_extractor", "claims_extractor",
      "reproducibility_judge", "discovery_scientist",
      "study_objective", "cohort_characteristics", "sample_groups",
      "statistical_methods", "key_findings", "reported_pathways",
      "reported_biomarkers", "study_limitations", "claims"
    ),
    properties = list(
      literature_curator = agent_schema,
      methods_extractor = agent_schema,
      claims_extractor = agent_schema,
      reproducibility_judge = agent_schema,
      discovery_scientist = agent_schema,
      study_objective = list(type = "string"),
      cohort_characteristics = list(type = "string"),
      sample_groups = list(type = "string"),
      statistical_methods = list(type = "string"),
      key_findings = list(type = "string"),
      reported_pathways = list(type = "array", items = list(type = "string")),
      reported_biomarkers = list(type = "array", items = list(type = "string")),
      study_limitations = list(type = "string"),
      claims = publication_extraction_schema()$properties$claims
    )
  )
}

openai_structured_json <- function(prompt, schema, api_key = "", model = "gpt-4o-mini") {
  if (!openai_available(api_key)) {
    stop("OpenAI extraction requires OPENAI_API_KEY or a pasted API key plus optional packages httr2 and jsonlite.", call. = FALSE)
  }
  request_body <- list(
    model = model,
    input = prompt,
    text = list(
      format = list(
        type = "json_schema",
        name = "seqsurf_structured_extraction",
        strict = TRUE,
        schema = schema
      )
    )
  )
  response <- httr2::request("https://api.openai.com/v1/responses") |>
    httr2::req_auth_bearer_token(openai_api_key(api_key)) |>
    httr2::req_headers(`Content-Type` = "application/json") |>
    httr2::req_body_json(request_body, auto_unbox = TRUE) |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_timeout(120) |>
    httr2::req_perform()
  body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  if (httr2::resp_status(response) >= 300) {
    message <- first_nonempty(body$error$message, paste("OpenAI API request failed with status", httr2::resp_status(response)))
    stop(message, call. = FALSE)
  }
  text <- extract_response_text(body)
  if (!nzchar(text)) stop("OpenAI API returned no text output.", call. = FALSE)
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

kb_evidence_text <- function(kb, max_chars = 60000) {
  doc_text <- ""
  if (!is.null(kb$documents) && "Text_path" %in% names(kb$documents)) {
    doc_text <- paste(vapply(kb$documents$Text_path, function(path) {
      if (nzchar(path) && file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n") else ""
    }, character(1)), collapse = "\n")
  }
  text <- paste(
    kb$geo$title, kb$geo$summary, kb$geo$overall_design,
    kb$pubmed$title, kb$pubmed$abstract,
    related_publications_evidence_text(kb$related_publications),
    kb$extra_papers,
    doc_text,
    sep = "\n\n"
  )
  substr(text, 1, max_chars)
}

llm_extract_publication_knowledge <- function(kb, api_key = "", model = "gpt-4o-mini", data_root = "data") {
  if (is.null(kb)) stop("Train or select a study knowledge base first.", call. = FALSE)
  prompt <- paste(
    "You are a careful scientific curator for an RNA-seq reproducibility app.",
    "Extract publication-grounded facts from the supplied GEO, PubMed, supplementary, and user-provided paper text.",
    "Only extract claims supported by the text. Preserve uncertainty. For gene/pathway claims, include direction and comparison only when stated.",
    "Return schema-valid JSON only.",
    "",
    "TEXT:",
    kb_evidence_text(kb),
    sep = "\n"
  )
  result <- openai_structured_json(prompt, publication_extraction_schema(), api_key, model)
  kb$extracted <- list(
    study_objective = first_nonempty(result$study_objective),
    cohort_characteristics = first_nonempty(result$cohort_characteristics),
    sample_groups = first_nonempty(result$sample_groups),
    statistical_methods = first_nonempty(result$statistical_methods),
    key_findings = first_nonempty(result$key_findings),
    reported_pathways = paste(unlist(result$reported_pathways), collapse = "; "),
    reported_biomarkers = paste(unlist(result$reported_biomarkers), collapse = "; "),
    study_limitations = first_nonempty(result$study_limitations)
  )
  kb$claims <- normalize_claim_table(result$claims)
  kb$llm <- list(model = model, extracted_at = Sys.time(), source = "OpenAI Responses API structured output")
  saveRDS(kb, assistant_kb_path(kb$accession, data_root))
  kb
}

virtual_lab_agent_names <- function() {
  c(
    literature_curator = "Literature Curator",
    methods_extractor = "Methods Extractor",
    claims_extractor = "Claims Extractor",
    reproducibility_judge = "Reproducibility Judge",
    discovery_scientist = "Discovery Scientist"
  )
}

build_virtual_lab_prompt <- function(kb) {
  paste(
    "You are SeqSurf Study Memory 2.0: a small virtual lab of specialist scientific agents for GEO transcriptomic reanalysis.",
    "Use only the supplied evidence. Do not invent papers, methods, datasets, claims, organisms, sample groups, or results.",
    "Each agent must produce a concise scientific review:",
    "1. Literature Curator: what has been done before, original paper/related-paper context, and citation-quality evidence.",
    "2. Methods Extractor: assay type, preprocessing, normalization, statistical methods, covariates, contrasts, and exact-reproduction feasibility.",
    "3. Claims Extractor: structured genes, pathways, comparisons, directions, and evidence sentences that can be validated.",
    "4. Reproducibility Judge: whether exact reproduction, paper-style approximation, or external workflow is appropriate, and why.",
    "5. Discovery Scientist: better modern follow-up analyses that are actually supported by the public data and metadata.",
    "Return schema-valid JSON only. Keep evidence passages short and attributable to the provided text.",
    "",
    "GEO ACCESSION:",
    kb$accession,
    "",
    "EVIDENCE:",
    kb_evidence_text(kb, max_chars = 100000),
    sep = "\n"
  )
}

run_study_memory_virtual_lab <- function(kb, api_key = "", model = "gpt-4o-mini", data_root = "data") {
  if (is.null(kb)) stop("Train or select a study knowledge base first.", call. = FALSE)
  result <- openai_structured_json(
    build_virtual_lab_prompt(kb),
    seqsurf_virtual_lab_schema(),
    api_key,
    model
  )
  agents <- virtual_lab_agent_names()
  kb$virtual_lab <- lapply(names(agents), function(agent_id) {
    agent <- result[[agent_id]]
    list(
      agent_id = agent_id,
      role = agents[[agent_id]],
      summary = first_nonempty(agent$summary, ""),
      evidence = unlist(agent$evidence),
      recommendations = unlist(agent$recommendations),
      uncertainties = unlist(agent$uncertainties)
    )
  })
  names(kb$virtual_lab) <- names(agents)
  kb$extracted <- list(
    study_objective = first_nonempty(result$study_objective),
    cohort_characteristics = first_nonempty(result$cohort_characteristics),
    sample_groups = first_nonempty(result$sample_groups),
    statistical_methods = first_nonempty(result$statistical_methods),
    key_findings = first_nonempty(result$key_findings),
    reported_pathways = paste(unlist(result$reported_pathways), collapse = "; "),
    reported_biomarkers = paste(unlist(result$reported_biomarkers), collapse = "; "),
    study_limitations = first_nonempty(result$study_limitations)
  )
  kb$claims <- normalize_claim_table(result$claims)
  kb$llm <- list(model = model, extracted_at = Sys.time(), source = "SeqSurf Study Memory 2.0 virtual lab structured output")
  saveRDS(kb, assistant_kb_path(kb$accession, data_root))
  kb
}

virtual_lab_summary_table <- function(kb) {
  if (is.null(kb) || is.null(kb$virtual_lab) || length(kb$virtual_lab) == 0) {
    return(data.frame(
      Agent = names(virtual_lab_agent_names()),
      Role = unname(virtual_lab_agent_names()),
      Summary = "Run Study Memory 2.0 after fetching evidence text.",
      Evidence = "",
      Uncertainty = "Not built yet",
      check.names = FALSE
    ))
  }
  rows <- lapply(kb$virtual_lab, function(agent) {
    data.frame(
      Agent = first_nonempty(agent$agent_id, ""),
      Role = first_nonempty(agent$role, ""),
      Summary = first_nonempty(agent$summary, ""),
      Evidence = paste(utils::head(agent$evidence, 3), collapse = " | "),
      Uncertainty = paste(utils::head(agent$uncertainties, 3), collapse = " | "),
      check.names = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

virtual_lab_recommendation_table <- function(kb) {
  if (is.null(kb) || is.null(kb$virtual_lab) || length(kb$virtual_lab) == 0) {
    return(data.frame(Role = "Study Memory 2.0", Recommendation = "Run the virtual lab extraction after evidence text is available.", Status = "Not ready", check.names = FALSE))
  }
  rows <- lapply(kb$virtual_lab, function(agent) {
    recs <- unlist(agent$recommendations)
    if (length(recs) == 0) recs <- "No recommendation returned."
    data.frame(
      Role = first_nonempty(agent$role, ""),
      Recommendation = recs,
      Status = "Agent recommendation",
      check.names = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

citation_evidence_table <- function(kb, max_rows = 40) {
  empty <- data.frame(
    Source = "Study Memory",
    Evidence = "Run evidence ingestion and Study Memory 2.0 to create citation-grounded evidence.",
    Used_for = "Not ready",
    Confidence = "",
    check.names = FALSE
  )
  if (is.null(kb)) return(empty)

  claim_rows <- data.frame()
  claims <- normalize_claim_table(kb$claims)
  if (nrow(claims) > 0) {
    claim_rows <- data.frame(
      Source = ifelse(nzchar(claims$source), claims$source, "Extracted claim"),
      Evidence = claims$evidence_sentence,
      Used_for = paste(
        claims$claim_type,
        ifelse(nzchar(claims$entity), paste0(":", claims$entity), ""),
        ifelse(nzchar(claims$reported_direction), paste0(" (", claims$reported_direction, ")"), ""),
        sep = ""
      ),
      Confidence = claims$confidence,
      check.names = FALSE
    )
  }

  agent_rows <- data.frame()
  if (!is.null(kb$virtual_lab) && length(kb$virtual_lab) > 0) {
    agent_rows <- dplyr::bind_rows(lapply(kb$virtual_lab, function(agent) {
      evidence <- unlist(agent$evidence)
      evidence <- evidence[nzchar(evidence)]
      if (length(evidence) == 0) return(NULL)
      data.frame(
        Source = first_nonempty(agent$role, agent$agent_id, "Virtual lab agent"),
        Evidence = evidence,
        Used_for = paste(utils::head(unlist(agent$recommendations), 2), collapse = " | "),
        Confidence = "Agent-cited evidence",
        check.names = FALSE
      )
    }))
  }

  document_rows <- data.frame()
  if (!is.null(kb$documents) && nrow(kb$documents) > 0) {
    document_rows <- data.frame(
      Source = kb$documents$Source,
      Evidence = if ("Text_path" %in% names(kb$documents)) kb$documents$Text_path else kb$documents$URL,
      Used_for = "PaperQA / source corpus",
      Confidence = kb$documents$Status,
      check.names = FALSE
    )
  }

  rows <- dplyr::bind_rows(claim_rows, agent_rows, document_rows)
  if (nrow(rows) == 0) return(empty)
  rows <- rows[nzchar(rows$Evidence), , drop = FALSE]
  utils::head(rows, max_rows)
}

retrieve_context_chunks <- function(kb, question, n = 6) {
  text <- kb_evidence_text(kb, max_chars = 120000)
  chunks <- unlist(strsplit(text, "(?<=\\.)\\s+|\\n\\n+", perl = TRUE))
  chunks <- chunks[nchar(chunks) > 80]
  if (length(chunks) == 0) return("")
  terms <- unique(tolower(unlist(strsplit(question, "[^A-Za-z0-9]+"))))
  terms <- terms[nchar(terms) >= 4]
  scores <- vapply(chunks, function(chunk) {
    sum(terms %in% unique(tolower(unlist(strsplit(chunk, "[^A-Za-z0-9]+")))))
  }, numeric(1))
  paste(utils::head(chunks[order(scores, decreasing = TRUE)], n), collapse = "\n\n")
}

chat_with_study_assistant <- function(dat, kb, question, api_key = "", model = "gpt-4o-mini") {
  if (is.null(kb)) stop("Train or select a study knowledge base first.", call. = FALSE)
  if (!openai_available(api_key)) {
    stop("Study chat requires OPENAI_API_KEY or a pasted API key plus optional packages httr2 and jsonlite.", call. = FALSE)
  }
  context <- retrieve_context_chunks(kb, question)
  prompt <- paste(
    "You are a study-specific scientific AI assistant embedded in an RNA-seq Explorer app.",
    "Answer using the study knowledge base and active app dataset summaries. Be explicit when evidence is missing.",
    "For Discovery questions, avoid generic wish-list methods unless the current metadata/results can support them.",
    "Prioritize analyses that can be run in the current app tabs or by the generated external workflow.",
    "Format the answer as concise markdown sections with headings:",
    "## Best next analysis",
    "## Why it matters",
    "## Tabs or files to use",
    "## Checks before trusting it",
    "",
    "ACTIVE DATASET:",
    paste("Dataset:", dat$info$display_name),
    paste("Default contrast:", dat$info$default_contrast),
    paste("Samples:", nrow(dat$metadata)),
    paste("Top DEG genes:", paste(utils::head(dat$deg$gene[order(dat$deg$padj, na.last = TRUE)], 20), collapse = ", ")),
    "",
    "RETRIEVED STUDY CONTEXT:",
    context,
    "",
    "QUESTION:",
    question,
    sep = "\n"
  )
  response <- httr2::request("https://api.openai.com/v1/responses") |>
    httr2::req_auth_bearer_token(openai_api_key(api_key)) |>
    httr2::req_body_json(list(model = model, input = prompt), auto_unbox = TRUE) |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_timeout(120) |>
    httr2::req_perform()
  body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  if (httr2::resp_status(response) >= 300) {
    message <- first_nonempty(body$error$message, paste("OpenAI API request failed with status", httr2::resp_status(response)))
    stop(message, call. = FALSE)
  }
  extract_response_text(body)
}

detect_terms <- function(text, dictionary) {
  text_low <- tolower(first_nonempty(text, ""))
  found <- dictionary[vapply(dictionary, function(term) grepl(tolower(term), text_low, fixed = TRUE), logical(1))]
  if (length(found) == 0) "Not automatically detected; add paper notes to enrich the knowledge base." else paste(found, collapse = "; ")
}

pathway_dictionary <- function() {
  c(
    "MYCN", "ALK", "p53", "apoptosis", "cell cycle", "immune response", "interferon",
    "inflammation", "hypoxia", "angiogenesis", "epithelial mesenchymal transition",
    "oxidative phosphorylation", "glycolysis", "DNA repair", "PI3K", "MAPK", "WNT",
    "TGF-beta", "Notch", "neural crest"
  )
}

detect_biomarkers <- function(text) {
  candidates <- unique(unlist(regmatches(text, gregexpr("\\b[A-Z0-9]{2,}[A-Z0-9-]*\\b", text))))
  candidates <- setdiff(candidates, c("RNA", "DNA", "GEO", "NCBI", "PMID", "GSE", "GPL", "GSM"))
  candidates <- candidates[nchar(candidates) <= 12]
  if (length(candidates) == 0) return("Not automatically detected; add paper notes to enrich the knowledge base.")
  paste(utils::head(candidates, 24), collapse = "; ")
}

split_detected_terms <- function(x) {
  x <- first_nonempty(x, "")
  if (!nzchar(x) || grepl("^Not automatically detected", x)) return(character())
  terms <- unique(trimws(unlist(strsplit(x, ";|,|\\n"))))
  terms[nzchar(terms)]
}

detect_methods <- function(text) {
  methods <- detect_terms(text, c(
    "DESeq2", "edgeR", "limma", "voom", "GSEA", "Gene Set Enrichment Analysis",
    "WGCNA", "Cox regression", "Kaplan-Meier", "PCA", "RMA", "batch correction",
    "negative binomial", "false discovery rate", "Benjamini"
  ))
  methods
}

detect_limitations <- function(text) {
  text_low <- tolower(first_nonempty(text, ""))
  hits <- c()
  if (grepl("retrospective", text_low)) hits <- c(hits, "Retrospective cohort")
  if (grepl("batch", text_low)) hits <- c(hits, "Batch effects require review")
  if (grepl("microarray", text_low)) hits <- c(hits, "Platform is not RNA-seq")
  if (grepl("processed", text_low)) hits <- c(hits, "Uses processed public expression")
  if (length(hits) == 0) hits <- "Not automatically detected; review the original methods."
  paste(unique(hits), collapse = "; ")
}

study_kb_summary_table <- function(kb) {
  if (is.null(kb)) return(data.frame(Field = character(), Value = character()))
  data.frame(
    Field = c(
      "GEO accession", "GEO title", "PubMed ID", "PubMed title", "Study objective",
      "Cohort characteristics", "Sample groups", "Statistical methods",
      "Reported pathways", "Reported biomarkers", "Study limitations", "Related papers found"
    ),
    Value = c(
      kb$accession,
      first_nonempty(kb$geo$title),
      first_nonempty(kb$geo$pubmed_id),
      first_nonempty(kb$pubmed$title),
      first_nonempty(kb$extracted$study_objective),
      first_nonempty(kb$extracted$cohort_characteristics),
      first_nonempty(kb$extracted$sample_groups),
      first_nonempty(kb$extracted$statistical_methods),
      first_nonempty(kb$extracted$reported_pathways),
      first_nonempty(kb$extracted$reported_biomarkers),
      first_nonempty(kb$extracted$study_limitations),
      as.character(if (!is.null(kb$related_publications) && nrow(kb$related_publications) > 0 && "PubMed_ID" %in% names(kb$related_publications)) {
        sum(nzchar(kb$related_publications$PubMed_ID))
      } else {
        0
      })
    ),
    check.names = FALSE
  )
}

prior_work_report_table <- function(kb) {
  if (is.null(kb)) {
    return(data.frame(Section = "Prior work", Finding = "Train study context first.", Evidence = "No GEO/PubMed knowledge base selected.", Status = "Not ready", check.names = FALSE))
  }
  claims <- normalize_claim_table(kb$claims)
  docs <- if (is.null(kb$documents)) 0 else nrow(kb$documents)
  related <- if (is.null(kb$related_publications)) empty_related_publications("Not searched.") else kb$related_publications
  related_titles <- if (nrow(related) > 0 && !"Not searched." %in% related$Status) {
    paste(utils::head(related$Title, 3), collapse = " | ")
  } else {
    "No related PubMed reuse papers detected yet."
  }
  data.frame(
    Section = c(
      "Original study",
      "Methods detected",
      "Reported results",
      "Related papers",
      "Dataset reuse clues",
      "Limitations",
      "Evidence depth"
    ),
    Finding = c(
      first_nonempty(kb$extracted$study_objective, kb$geo$title),
      first_nonempty(kb$extracted$statistical_methods, "No statistical methods automatically detected."),
      paste(
        "Biomarkers:", first_nonempty(kb$extracted$reported_biomarkers, "none detected"),
        "| Pathways:", first_nonempty(kb$extracted$reported_pathways, "none detected")
      ),
      paste(nrow(related), "PubMed papers mentioning the GEO accession."),
      if (grepl("validation|replication|cohort|reanalysis|reuse|public", kb_evidence_blob(kb), ignore.case = TRUE)) {
        "Study text suggests reuse, validation, cohort structure, or public reanalysis context."
      } else {
        "No obvious reuse/validation language detected; manual literature search may still be needed."
      },
      first_nonempty(kb$extracted$study_limitations, "No limitations automatically detected."),
      paste(nrow(claims), "structured claims;", docs, "downloaded/indexed documents.")
    ),
    Evidence = c(
      first_nonempty(kb$pubmed$title, kb$geo$title),
      first_nonempty(kb$extracted$statistical_methods, "Rule/LLM extraction did not identify methods."),
      first_nonempty(kb$extracted$key_findings, "Use LLM extraction or add paper notes for richer results."),
      related_titles,
      "Heuristic scan of GEO/PubMed/notes text.",
      first_nonempty(kb$extracted$study_limitations, "Not stated in available context."),
      "Knowledge-base extraction status."
    ),
    Status = c("Review", "Review", "Review", "Review", "Check", "Check", if (nrow(claims) > 0) "Usable" else "Needs richer extraction"),
    check.names = FALSE
  )
}

exact_reproduction_feasibility_table <- function(kb) {
  if (is.null(kb)) {
    return(data.frame(Check = "Feasibility", Assessment = "Train study context first.", Status = "Not ready", Rationale = "No study evidence available.", check.names = FALSE))
  }
  text <- tolower(kb_evidence_blob(kb))
  methods <- split_detected_terms(kb$extracted$statistical_methods)
  has_methods <- length(methods) > 0
  has_counts <- grepl("count matrix|counts matrix|raw counts|htseq|featurecounts|read counts|gene counts", text)
  has_raw <- grepl("fastq|sra|bam|sam", text)
  has_processed <- grepl("processed|normalized|series matrix|rma|microarray|expression profiling by array", text)
  has_groups <- nzchar(first_nonempty(kb$extracted$sample_groups, ""))
  has_claims <- nrow(normalize_claim_table(kb$claims)) > 0
  has_covariate_detail <- grepl("covariate|batch|adjust|model|linear model|design matrix|paired|blocking", text)
  has_thresholds <- grepl("fdr|adjusted p|q-value|fold change|log2|p\\s*<", text)
  file_mode <- if (has_raw) {
    "Raw FASTQ/SRA/BAM signal"
  } else if (has_counts) {
    "Count matrix signal"
  } else if (has_processed) {
    "Processed expression/series matrix signal"
  } else {
    "Public expression format unclear"
  }
  status <- if (has_raw && has_methods && has_groups && has_covariate_detail) {
    "Exact reproduction possible only with external raw-data workflow"
  } else if (has_counts && has_methods && has_groups) {
    "Near-exact count-level reproduction likely"
  } else if (has_processed && has_methods && has_groups) {
    "Paper-style approximation likely"
  } else {
    "Not directly reproducible from public evidence yet"
  }
  expected_agreement <- if (grepl("Exact|Near-exact", status)) {
    "High agreement is expected if software/reference versions and filtering choices can be matched; mismatches should be investigated as reproducibility gaps."
  } else if (grepl("approximation", status)) {
    "Moderate agreement is the fair target; direction and broad pathway support matter more than identical p-values."
  } else {
    "Low or mixed agreement should not be overinterpreted until missing files, sample groups, and methods are curated."
  }
  data.frame(
    Check = c("Overall feasibility", "Public expression files", "Methods detail", "Group definitions", "Model/covariates", "Thresholds/effect criteria", "Claim extraction", "Expected agreement", "Recommended mode"),
    Assessment = c(
      status,
      file_mode,
      if (has_methods) paste(methods, collapse = "; ") else "Methods not detected.",
      if (has_groups) first_nonempty(kb$extracted$sample_groups) else "Sample groups not detected.",
      if (has_covariate_detail) "Model/covariate language detected." else "Model/covariate details not clearly detected.",
      if (has_thresholds) "Statistical thresholds/effect criteria detected." else "Thresholds not clearly detected.",
      if (has_claims) paste(nrow(normalize_claim_table(kb$claims)), "claims available.") else "No structured claims yet.",
      expected_agreement,
      if (grepl("external raw-data", status)) "Generate FASTQ/SRA workflow handoff, then import count/expression outputs." else if (grepl("Near-exact", status)) "Run count-aware DESeq2/edgeR/limma-voom automation and document deviations." else if (grepl("approximation", status)) "Run paper-style approximation and disclose public-data limits." else "Do not present as exact reproduction until files/methods are curated."
    ),
    Status = c(status, if (has_raw || has_counts || has_processed) "Usable" else "Manual check", if (has_methods) "Usable" else "Missing", if (has_groups) "Usable" else "Missing", if (has_covariate_detail) "Usable" else "Approximate", if (has_thresholds) "Usable" else "Approximate", if (has_claims) "Usable" else "Needs extraction", "Interpretation guide", "Decision"),
    Rationale = c(
      "Direct reproducibility depends on public files, clear methods, sample grouping, covariates, thresholds, and claim extraction.",
      "Determines whether GEO import can automate analysis or needs custom preprocessing.",
      "Needed to match normalization, model, covariates, and thresholds.",
      "Needed to recreate contrasts.",
      "Missing covariates make exact p-value reproduction unlikely.",
      "Missing thresholds make agreement scoring more qualitative.",
      "Needed for verification scoring.",
      "Defines how strictly SeqSurf should compare reproduced results to the paper.",
      "This wording should be used in demos, abstracts, and reports."
    ),
    check.names = FALSE
  )
}

raw_supplementary_triage_table <- function(kb) {
  if (is.null(kb)) {
    return(data.frame(Resource = "Supplementary/raw files", Signal = "Train study context first.", Action = "Not ready", Status = "Not ready", check.names = FALSE))
  }
  supp <- if (is.null(kb$geo$supplementary_files)) character() else kb$geo$supplementary_files
  text <- tolower(paste(kb_evidence_blob(kb), supp, collapse = " "))
  data.frame(
    Resource = c("Series matrix", "Count matrix", "Raw sequencing", "Supplementary archives", "Manual pipeline handoff"),
    Signal = c(
      if (grepl("series matrix|processed|normalized", text)) "Processed/series-matrix language detected." else "Not clearly detected.",
      if (grepl("count|counts|htseq|featurecounts", text)) "Count-matrix language detected." else "Not clearly detected.",
      if (grepl("fastq|sra|bam|sam", text)) "Raw sequencing language detected." else "Not clearly detected.",
      paste(length(supp), "supplementary file links parsed from GEO."),
      "If raw/count files are present but not series-matrix compatible, generate a paper-to-code brief."
    ),
    Action = c(
      "Use automated GEO import first.",
      "Prefer count-aware DESeq2/edgeR script outside Shiny, then write app dataset contract.",
      "Do not process FASTQ inside Shiny; hand off to reproducible pipeline.",
      "Fetch/index small paper-like documents; avoid pulling huge archives during interactive demos.",
      "Use exported brief as the coding-agent/spec handoff."
    ),
    Status = c(
      if (grepl("series matrix|processed|normalized", text)) "Automatable" else "Check",
      if (grepl("count|counts|htseq|featurecounts", text)) "Pipeline needed" else "Unknown",
      if (grepl("fastq|sra|bam|sam", text)) "Pipeline needed" else "Unknown",
      if (length(supp) > 0) "Review links" else "No links parsed",
      "Recommended"
    ),
    check.names = FALSE
  )
}

count_matrix_pipeline_plan <- function(kb = NULL) {
  route <- geo_analyzability_route(kb)
  feasibility <- exact_reproduction_feasibility_table(kb)
  accession <- if (is.null(kb)) "GEO_ACCESSION" else first_nonempty(kb$accession, "GEO_ACCESSION")
  method_hint <- if (!is.null(kb)) first_nonempty(kb$extracted$statistical_methods, "DESeq2/edgeR/limma-voom; choose based on paper methods and count format.") else "DESeq2/edgeR/limma-voom"
  groups <- if (!is.null(kb)) first_nonempty(kb$extracted$sample_groups, "Define group column from metadata before running.") else "Define group column from metadata before running."
  paste(
    "# Count-matrix automation plan",
    "",
    paste0("- GEO/study: ", accession),
    paste0("- Route decision: ", route$Route, " (", route$Can_run_in_app, ", ", route$Confidence, "% confidence)"),
    paste0("- Reproduction level: ", feasibility$Status[[1]]),
    paste0("- Method hints from paper/context: ", method_hint),
    paste0("- Grouping evidence: ", groups),
    "",
    "## In-app path",
    "1. Upload a gene-by-sample count matrix and a sample metadata table on Start.",
    "2. Map gene, sample ID, group, reference group, and comparison group columns.",
    "3. SeqSurf detects integer-like nonnegative counts and runs count-aware limma-voom when edgeR/limma are installed.",
    "4. The app writes metadata.rds, vsd_matrix.rds/logCPM matrix, PCA objects, DEG results, GSEA placeholder/results, and dataset_info.rds.",
    "5. Validate paper agreement against extracted claims; interpret agreement using the reproduction level above.",
    "",
    "## External exactness checks",
    "- Record raw file URL, count-generation method, genome/annotation version, filtering thresholds, normalization, model formula, contrast, and exclusions.",
    "- If the paper used DESeq2 specifically, run a DESeq2 script externally for exact-method comparison and import normalized values/results into the app contract.",
    "- If covariates or sample pairing are missing from public metadata, label the result paper-style approximation, not exact reproduction.",
    sep = "\n"
  )
}

raw_sra_handoff_plan <- function(kb = NULL) {
  accession <- if (is.null(kb)) "GEO_ACCESSION" else first_nonempty(kb$accession, "GEO_ACCESSION")
  supp <- if (is.null(kb) || is.null(kb$geo$supplementary_files)) character() else kb$geo$supplementary_files
  supp_lines <- if (length(supp) == 0) "- No supplementary/raw URLs parsed yet." else paste0("- ", supp, collapse = "\n")
  paste(
    "# FASTQ/SRA external workflow handoff",
    "",
    paste0("- GEO/study: ", accession),
    "- SeqSurf does not process arbitrary FASTQ/SRA inside Shiny. Use an external workflow, then import the generated count matrix/app contract.",
    "",
    "## Raw/supplementary signals",
    supp_lines,
    "",
    "## Nextflow/Snakemake-style steps",
    "1. Resolve SRA run accessions and download with fasterq-dump/prefetch or nf-core/fetchngs.",
    "2. Run FastQC/MultiQC and record read counts, quality failures, and excluded samples.",
    "3. Quantify with Salmon/STAR+featureCounts using a documented genome and GTF annotation version.",
    "4. Summarize to gene counts, join GEO sample metadata, and verify sample counts against the paper.",
    "5. Run DESeq2/edgeR/limma-voom with the paper's contrast/covariates where public metadata allow.",
    "6. Export the SeqSurf app contract files and a provenance table.",
    "",
    "## Interpretation rule",
    "- If reference genome, annotation, preprocessing, covariates, or sample exclusions differ from the paper, report the workflow as an approximate reproduction and let the validation score reflect claim agreement.",
    sep = "\n"
  )
}

extract_sra_accessions <- function(kb = NULL) {
  text <- if (is.null(kb)) "" else paste(kb_evidence_blob(kb), kb$geo$supplementary_files, collapse = " ")
  hits <- unlist(regmatches(text, gregexpr("\\b(SRR|ERR|DRR|SRX|ERX|DRX|SRS|ERS|DRS|SRP|ERP|DRP)[0-9]+\\b", text, perl = TRUE)))
  unique(hits[nzchar(hits)])
}

no_code_ai_handoff_guide <- function(kb = NULL) {
  accession <- if (is.null(kb)) "GEO_ACCESSION" else first_nonempty(kb$accession, "GEO_ACCESSION")
  paste(
    "# SeqSurf no-code external analysis guide",
    "",
    paste0("Study: ", accession),
    "",
    "This guide is for users who do not code but need to analyze a GEO dataset that cannot be fully processed inside SeqSurf.",
    "SeqSurf prepares the scientific plan, workflow files, and AI prompt. The user runs the workflow with help from an AI coding assistant, then uploads the resulting count matrix and metadata back into SeqSurf.",
    "",
    "## What to install",
    "",
    "1. Install R from CRAN: https://cran.r-project.org/",
    "2. Install Positron, a free R/Python data-science IDE: https://positron.posit.co/download.html",
    "3. Choose an AI coding assistant:",
    "   - Claude Code web or desktop: https://claude.ai/code",
    "   - Claude Code docs/install: https://code.claude.com/docs/en/overview",
    "   - Or another coding assistant that can read folders and run terminal/R commands.",
    "",
    "## What to download from SeqSurf",
    "",
    "Download these from the Start -> Plan or Train tabs:",
    "",
    "1. Raw/SRA workflow bundle ZIP",
    "2. Paper-to-code brief",
    "3. This no-code AI guide",
    "4. AI assistant prompt",
    "",
    "Put all downloaded files into one project folder, then unzip the workflow bundle.",
    "",
    "## What the AI assistant should do",
    "",
    "Paste `PROMPT_FOR_AI_ASSISTANT.md` into Claude Code, Positron Assistant, Codex, or another AI coding assistant.",
    "The AI should:",
    "",
    "- inspect the workflow folder",
    "- install missing R packages or tell the user exactly what is missing",
    "- resolve GEO/SRA run accessions",
    "- download FASTQ files when raw SRA processing is required",
    "- quantify or count reads using a documented method",
    "- create `outputs/seqsurf_gene_counts.csv`",
    "- create `outputs/seqsurf_metadata.csv`",
    "- document every assumption and any deviation from the original paper",
    "",
    "## What to upload back into SeqSurf",
    "",
    "After the AI-assisted workflow finishes, return to SeqSurf and upload:",
    "",
    "- `outputs/seqsurf_gene_counts.csv` as the expression/count matrix",
    "- `outputs/seqsurf_metadata.csv` as the metadata table",
    "",
    "Then use Validate to run DESeq2/edgeR/limma-voom-style analysis and compare against paper claims.",
    "",
    "## Safety rules",
    "",
    "- Do not paste private API keys, passwords, or patient identifiers into the AI assistant.",
    "- Treat the external workflow as approximate unless genome build, annotation, covariates, sample exclusions, and software versions match the paper.",
    "- Keep the generated logs and scripts. They are part of the reproducibility record.",
    sep = "\n"
  )
}

external_ai_analysis_prompt <- function(kb = NULL) {
  accession <- if (is.null(kb)) "GEO_ACCESSION" else first_nonempty(kb$accession, "GEO_ACCESSION")
  route <- geo_analyzability_route(kb)
  route_text <- paste0(
    "- Route: ", route$Route[[1]],
    "\n- Can run in SeqSurf app: ", route$Can_run_in_app[[1]],
    "\n- Evidence: ", route$Evidence[[1]],
    "\n- Blocking reason: ", route$Blocking_reason[[1]],
    "\n- Recommended next step: ", route$Recommended_next_step[[1]]
  )
  paste(
    "You are helping a non-coder run a reproducible transcriptomics reanalysis prepared by SeqSurf AI.",
    "",
    "Goal:",
    paste0("Analyze GEO study ", accession, " using the files in this project folder, then produce files that can be uploaded back into SeqSurf."),
    "",
    "Important constraints:",
    "- Explain every step in plain language before running commands.",
    "- Ask before installing software or downloading very large files.",
    "- Do not delete user files.",
    "- Do not claim exact reproduction unless the public files, genome/annotation, software, covariates, and sample exclusions match the original paper.",
    "- Save a short provenance note for every command, input file, and assumption.",
    "",
    "SeqSurf route decision:",
    route_text,
    "",
    "Files you may see:",
    "- README.md",
    "- AI_HANDOFF_GUIDE.md",
    "- accessions/sra_seed_accessions.txt",
    "- metadata/sample_sheet_template.csv",
    "- scripts/01_resolve_sra_runs.R",
    "- scripts/02_download_fastq.sh",
    "- scripts/03_quant_salmon.sh",
    "- scripts/04_build_gene_counts.R",
    "- scripts/05_prepare_seqsurf_import.R",
    "- paper-to-code brief from SeqSurf",
    "",
    "Do this:",
    "1. Inspect the folder and summarize what files are present.",
    "2. Check whether R is installed and whether required R packages are installed.",
    "3. Check whether command line tools are available: sra-tools, fastqc, multiqc, salmon.",
    "4. If this is an SRA/raw route, resolve SRR/ERR/DRR run accessions.",
    "5. Help the user review and edit `metadata/sample_sheet.csv` so sample groups are correct.",
    "6. Download FASTQ files only after warning the user about disk/time requirements.",
    "7. Quantify or count reads with a documented method.",
    "8. Create `outputs/seqsurf_gene_counts.csv` with genes as rows and samples as columns.",
    "9. Create `outputs/seqsurf_metadata.csv` with sample IDs and group labels.",
    "10. Write `outputs/provenance_notes.md` summarizing files, commands, versions, assumptions, and deviations from the paper.",
    "11. Tell the user to upload `outputs/seqsurf_gene_counts.csv` and `outputs/seqsurf_metadata.csv` back into SeqSurf.",
    "",
    "If something fails:",
    "- Diagnose the error.",
    "- Explain it in plain language.",
    "- Propose the smallest safe fix.",
    "- If a file or metadata column is missing, ask the user before guessing.",
    sep = "\n"
  )
}

raw_sra_workflow_files <- function(kb = NULL) {
  accession <- if (is.null(kb)) "GEO_ACCESSION" else first_nonempty(kb$accession, "GEO_ACCESSION")
  sra_accessions <- extract_sra_accessions(kb)
  sra_seed <- if (length(sra_accessions) == 0) accession else paste(sra_accessions, collapse = "\n")
  sample_sheet <- data.frame(
    run = if (length(grep("^SRR|^ERR|^DRR", sra_accessions, value = TRUE)) > 0) grep("^SRR|^ERR|^DRR", sra_accessions, value = TRUE) else "RUN_ACCESSION",
    sample_id = "sample_1",
    group = "GROUP_A",
    fastq_1 = "fastq/RUN_ACCESSION_1.fastq.gz",
    fastq_2 = "fastq/RUN_ACCESSION_2.fastq.gz",
    stringsAsFactors = FALSE
  )
  list(
    "README.md" = paste(
      "# SeqSurf external raw/SRA workflow bundle",
      "",
      paste0("Study: ", accession),
      "",
      "This bundle is for datasets that SeqSurf should not process inside Shiny because they require raw FASTQ/SRA/BAM/CEL/IDAT-style preprocessing.",
      "",
      "## What this does",
      "",
      "1. Resolve SRA run accessions from GEO/SRA metadata.",
      "2. Download FASTQ files with SRA Toolkit.",
      "3. Quantify reads with Salmon, or adapt the shell script to STAR/featureCounts.",
      "4. Build a gene-level count matrix.",
      "5. Prepare count matrix + metadata files that can be imported back into SeqSurf.",
      "",
      "## Requirements",
      "",
      "- R packages: GEOquery, readr, dplyr. Optional for Salmon import: tximport.",
      "- Command line tools: sra-tools, salmon, multiqc, fastqc.",
      "- A Salmon transcriptome index and a transcript-to-gene map.",
      "",
      "## Run order",
      "",
      "```bash",
      "Rscript scripts/01_resolve_sra_runs.R",
      "bash scripts/02_download_fastq.sh",
      "bash scripts/03_quant_salmon.sh /path/to/salmon_index",
      "Rscript scripts/04_build_gene_counts.R tx2gene.csv",
      "Rscript scripts/05_prepare_seqsurf_import.R",
      "```",
      "",
      "After step 5, upload `outputs/seqsurf_gene_counts.csv` and `outputs/seqsurf_metadata.csv` in SeqSurf.",
      sep = "\n"
    ),
    "AI_HANDOFF_GUIDE.md" = no_code_ai_handoff_guide(kb),
    "PROMPT_FOR_AI_ASSISTANT.md" = external_ai_analysis_prompt(kb),
    "START_HERE.command" = paste(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "cd \"$(dirname \"$0\")\"",
      "echo 'SeqSurf external workflow workspace'",
      "echo '1) Running preflight checks now.'",
      "Rscript scripts/00_preflight_check.R || true",
      "echo ''",
      "echo '2) Next, open this folder in Positron/Claude Code/Codex and paste PROMPT_FOR_AI_ASSISTANT.md.'",
      "echo '3) Raw FASTQ download and quantification require your confirmation because they may use substantial disk and compute.'",
      "read -r -p 'Press Return to close this window.' _",
      sep = "\n"
    ),
    "accessions/sra_seed_accessions.txt" = sra_seed,
    "metadata/sample_sheet_template.csv" = paste(capture.output(utils::write.csv(sample_sheet, row.names = FALSE, quote = TRUE)), collapse = "\n"),
    "scripts/00_preflight_check.R" = paste(
      "cat('# SeqSurf external workflow preflight\\n\\n')",
      "cat('Working directory:', normalizePath(getwd(), mustWork = FALSE), '\\n\\n')",
      "required_r <- c('GEOquery', 'readr', 'dplyr')",
      "optional_r <- c('tximport', 'BiocManager')",
      "pkg_status <- function(pkgs) data.frame(Package = pkgs, Installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE), check.names = FALSE)",
      "print(pkg_status(required_r))",
      "cat('\\nOptional R packages:\\n')",
      "print(pkg_status(optional_r))",
      "tools <- c('prefetch', 'fasterq-dump', 'fastqc', 'multiqc', 'salmon')",
      "tool_status <- data.frame(Tool = tools, Path = vapply(tools, Sys.which, character(1)), check.names = FALSE)",
      "tool_status$Available <- nzchar(tool_status$Path)",
      "cat('\\nCommand line tools:\\n')",
      "print(tool_status)",
      "dir.create('outputs', showWarnings = FALSE, recursive = TRUE)",
      "utils::write.csv(tool_status, 'outputs/preflight_tool_status.csv', row.names = FALSE)",
      "cat('\\nWrote outputs/preflight_tool_status.csv\\n')",
      sep = "\n"
    ),
    "scripts/01_resolve_sra_runs.R" = paste(
      paste0("accession <- Sys.getenv('SEQSURF_GEO_ACCESSION', unset = '", accession, "')"),
      "dir.create('metadata', showWarnings = FALSE, recursive = TRUE)",
      "dir.create('accessions', showWarnings = FALSE, recursive = TRUE)",
      "seed <- readLines('accessions/sra_seed_accessions.txt', warn = FALSE)",
      "seed <- unique(seed[nzchar(seed)])",
      "message('Resolving SRA runs for ', accession)",
      "runs <- character()",
      "if (requireNamespace('GEOquery', quietly = TRUE)) {",
      "  gse <- tryCatch(GEOquery::getGEO(accession, GSEMatrix = FALSE), error = function(e) NULL)",
      "  if (!is.null(gse)) {",
      "    gsms <- GEOquery::GSMList(gse)",
      "    rel <- unlist(lapply(gsms, function(x) GEOquery::Meta(x)$relation))",
      "    runs <- unique(unlist(regmatches(rel, gregexpr('\\\\b(SRR|ERR|DRR)[0-9]+\\\\b', rel, perl = TRUE))))",
      "  }",
      "}",
      "runs <- unique(c(runs, grep('^(SRR|ERR|DRR)[0-9]+$', seed, value = TRUE)))",
      "if (length(runs) == 0) {",
      "  message('No SRR/ERR/DRR runs were auto-resolved. Add runs to accessions/sra_seed_accessions.txt or metadata/sample_sheet_template.csv.')",
      "} else {",
      "  writeLines(runs, 'accessions/sra_runs.txt')",
      "  sample_sheet <- data.frame(run = runs, sample_id = runs, group = 'GROUP_A', fastq_1 = paste0('fastq/', runs, '_1.fastq.gz'), fastq_2 = paste0('fastq/', runs, '_2.fastq.gz'))",
      "  utils::write.csv(sample_sheet, 'metadata/sample_sheet.csv', row.names = FALSE)",
      "  message('Wrote ', length(runs), ' runs to accessions/sra_runs.txt and metadata/sample_sheet.csv')",
      "}",
      sep = "\n"
    ),
    "scripts/02_download_fastq.sh" = paste(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "mkdir -p sra fastq logs",
      "if [ ! -f accessions/sra_runs.txt ]; then",
      "  echo 'Missing accessions/sra_runs.txt. Run scripts/01_resolve_sra_runs.R first or create it manually.'",
      "  exit 1",
      "fi",
      "while read -r run; do",
      "  [ -z \"$run\" ] && continue",
      "  echo \"Downloading $run\"",
      "  prefetch \"$run\" --output-directory sra 2>&1 | tee \"logs/${run}_prefetch.log\"",
      "  fasterq-dump \"sra/${run}\" --split-files --threads \"${THREADS:-8}\" --outdir fastq 2>&1 | tee \"logs/${run}_fasterq.log\"",
      "  gzip -f fastq/${run}*.fastq",
      "done < accessions/sra_runs.txt",
      "fastqc fastq/*.fastq.gz -o logs || true",
      "multiqc logs -o logs || true",
      sep = "\n"
    ),
    "scripts/03_quant_salmon.sh" = paste(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "INDEX=\"${1:-}\"",
      "if [ -z \"$INDEX\" ]; then echo 'Usage: bash scripts/03_quant_salmon.sh /path/to/salmon_index'; exit 1; fi",
      "mkdir -p salmon logs",
      "SHEET=\"metadata/sample_sheet.csv\"",
      "if [ ! -f \"$SHEET\" ]; then SHEET=\"metadata/sample_sheet_template.csv\"; fi",
      "tail -n +2 \"$SHEET\" | while IFS=, read -r run sample_id group fastq1 fastq2; do",
      "  run=$(echo \"$run\" | tr -d '\"')",
      "  fastq1=$(echo \"$fastq1\" | tr -d '\"')",
      "  fastq2=$(echo \"$fastq2\" | tr -d '\"')",
      "  echo \"Quantifying $run\"",
      "  if [ -f \"$fastq2\" ]; then",
      "    salmon quant -i \"$INDEX\" -l A -1 \"$fastq1\" -2 \"$fastq2\" -p \"${THREADS:-8}\" -o \"salmon/$run\" 2>&1 | tee \"logs/${run}_salmon.log\"",
      "  else",
      "    salmon quant -i \"$INDEX\" -l A -r \"$fastq1\" -p \"${THREADS:-8}\" -o \"salmon/$run\" 2>&1 | tee \"logs/${run}_salmon.log\"",
      "  fi",
      "done",
      sep = "\n"
    ),
    "scripts/04_build_gene_counts.R" = paste(
      "args <- commandArgs(trailingOnly = TRUE)",
      "tx2gene_path <- if (length(args) >= 1) args[[1]] else 'tx2gene.csv'",
      "dir.create('outputs', showWarnings = FALSE, recursive = TRUE)",
      "sample_sheet <- if (file.exists('metadata/sample_sheet.csv')) read.csv('metadata/sample_sheet.csv') else read.csv('metadata/sample_sheet_template.csv')",
      "quant_files <- file.path('salmon', sample_sheet$run, 'quant.sf')",
      "names(quant_files) <- sample_sheet$sample_id",
      "if (requireNamespace('tximport', quietly = TRUE) && file.exists(tx2gene_path) && all(file.exists(quant_files))) {",
      "  tx2gene <- read.csv(tx2gene_path)",
      "  txi <- tximport::tximport(quant_files, type = 'salmon', tx2gene = tx2gene)",
      "  counts <- round(txi$counts)",
      "} else {",
      "  stop('Install tximport, provide tx2gene.csv, and confirm salmon/*/quant.sf files exist.')",
      "}",
      "counts_df <- data.frame(gene = rownames(counts), counts, check.names = FALSE)",
      "write.csv(counts_df, 'outputs/seqsurf_gene_counts.csv', row.names = FALSE)",
      "message('Wrote outputs/seqsurf_gene_counts.csv')",
      sep = "\n"
    ),
    "scripts/05_prepare_seqsurf_import.R" = paste(
      "dir.create('outputs', showWarnings = FALSE, recursive = TRUE)",
      "sample_sheet <- if (file.exists('metadata/sample_sheet.csv')) read.csv('metadata/sample_sheet.csv') else read.csv('metadata/sample_sheet_template.csv')",
      "metadata <- data.frame(app_sample_id = sample_sheet$sample_id, geo_accession = sample_sheet$run, group = sample_sheet$group, stringsAsFactors = FALSE)",
      "write.csv(metadata, 'outputs/seqsurf_metadata.csv', row.names = FALSE)",
      "message('Upload outputs/seqsurf_gene_counts.csv and outputs/seqsurf_metadata.csv in SeqSurf.')",
      sep = "\n"
    )
  )
}

write_raw_sra_workflow_bundle <- function(kb = NULL, zipfile) {
  files <- raw_sra_workflow_files(kb)
  bundle_dir <- tempfile("seqsurf_raw_sra_bundle_")
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  for (path in names(files)) {
    out <- file.path(bundle_dir, path)
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[path]], out)
  }
  Sys.chmod(c(
    file.path(bundle_dir, "START_HERE.command"),
    file.path(bundle_dir, "scripts", c("02_download_fastq.sh", "03_quant_salmon.sh"))
  ), mode = "0755")
  setwd(bundle_dir)
  utils::zip(zipfile, files = names(files), flags = "-r9Xq")
  zipfile
}

external_workflow_root <- function(data_root = "data") {
  file.path(data_root, "external_workflows")
}

external_workflow_dir <- function(kb = NULL, data_root = "data") {
  accession <- if (is.null(kb)) "study" else normalize_geo_accession(first_nonempty(kb$accession, "study"))
  file.path(external_workflow_root(data_root), accession)
}

write_external_workflow_workspace <- function(kb = NULL, data_root = "data") {
  dir <- external_workflow_dir(kb, data_root)
  files <- raw_sra_workflow_files(kb)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (path in names(files)) {
    out <- file.path(dir, path)
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[path]], out, useBytes = TRUE)
  }
  Sys.chmod(c(
    file.path(dir, "START_HERE.command"),
    file.path(dir, "scripts", c("02_download_fastq.sh", "03_quant_salmon.sh"))
  ), mode = "0755")
  list(
    directory = normalizePath(dir, mustWork = FALSE),
    prompt = normalizePath(file.path(dir, "PROMPT_FOR_AI_ASSISTANT.md"), mustWork = FALSE),
    starter = normalizePath(file.path(dir, "START_HERE.command"), mustWork = FALSE),
    preflight = normalizePath(file.path(dir, "scripts", "00_preflight_check.R"), mustWork = FALSE),
    command = paste0(
      "cd ", shQuote(normalizePath(dir, mustWork = FALSE)),
      "\nRscript scripts/00_preflight_check.R",
      "\nopen -a Positron ", shQuote(normalizePath(dir, mustWork = FALSE))
    )
  )
}

learning_root <- function(data_root = "data") {
  file.path(data_root, "ai_learning")
}

learning_log_path <- function(data_root = "data") {
  file.path(learning_root(data_root), "seqsurf_learning_log.csv")
}

shared_learning_pool_path <- function(data_root = "data") {
  file.path(learning_root(data_root), "seqsurf_shared_learning_pool.csv")
}

training_export_root <- function(data_root = "data") {
  file.path(learning_root(data_root), "training_exports")
}

empty_learning_log <- function() {
  data.frame(
    timestamp = character(),
    accession = character(),
    dataset_id = character(),
    route = character(),
    assay_signal = character(),
    sample_bin = character(),
    method_terms = character(),
    learned_tags = character(),
    event_type = character(),
    outcome = character(),
    usefulness_rating = numeric(),
    reproducibility_score = numeric(),
    proceed_score = numeric(),
    claims_checked = numeric(),
    supported_claims = numeric(),
    missing_claims = numeric(),
    direction_mismatches = numeric(),
    notes = character(),
    check.names = FALSE
  )
}

read_learning_log <- function(data_root = "data") {
  path <- learning_log_path(data_root)
  if (!file.exists(path)) return(empty_learning_log())
  out <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) empty_learning_log())
  missing <- setdiff(names(empty_learning_log()), names(out))
  for (col in missing) out[[col]] <- empty_learning_log()[[col]]
  out[, names(empty_learning_log()), drop = FALSE]
}

privacy_safe_learning_log <- function(log) {
  if (is.null(log) || nrow(log) == 0) return(empty_learning_log())
  safe <- log[, names(empty_learning_log()), drop = FALSE]
  safe$timestamp <- as.character(as.Date(Sys.time()))
  safe$accession <- ""
  safe$dataset_id <- ""
  safe$notes <- ""
  safe$event_type <- paste0("shared_", first_nonempty(safe$event_type, "learning"))
  safe
}

read_shared_learning_pool <- function(data_root = "data") {
  path <- shared_learning_pool_path(data_root)
  if (!file.exists(path)) return(empty_learning_log())
  out <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) empty_learning_log())
  missing <- setdiff(names(empty_learning_log()), names(out))
  for (col in missing) out[[col]] <- empty_learning_log()[[col]]
  out <- out[, names(empty_learning_log()), drop = FALSE]
  out$accession <- ""
  out$dataset_id <- ""
  out$notes <- ""
  out
}

combined_learning_log <- function(data_root = "data", include_shared = TRUE) {
  local <- read_learning_log(data_root)
  if (!include_shared) return(local)
  dplyr::bind_rows(local, read_shared_learning_pool(data_root))
}

write_shared_learning_export <- function(file, data_root = "data") {
  log <- read_learning_log(data_root)
  safe <- privacy_safe_learning_log(log)
  if (nrow(safe) == 0) {
    safe <- empty_learning_log()
  }
  metadata <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    field = c("format", "privacy", "rows"),
    value = c(
      "SeqSurf shared learning packet v1",
      "No accession, dataset_id, or free-text notes included.",
      as.character(nrow(safe))
    ),
    check.names = FALSE
  )
  tmp_dir <- tempfile("seqsurf_shared_learning_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(metadata, file.path(tmp_dir, "README_metadata.csv"), row.names = FALSE)
  write.csv(safe, file.path(tmp_dir, "shared_learning_events.csv"), row.names = FALSE)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(tmp_dir)
  utils::zip(file, files = c("README_metadata.csv", "shared_learning_events.csv"), flags = "-r9Xq")
  file
}

read_shared_learning_import_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "zip")) {
    tmp <- tempfile("seqsurf_shared_import_")
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(path, exdir = tmp)
    csv <- file.path(tmp, "shared_learning_events.csv")
    if (!file.exists(csv)) stop("Shared learning ZIP must contain shared_learning_events.csv.", call. = FALSE)
    out <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  missing <- setdiff(names(empty_learning_log()), names(out))
  for (col in missing) out[[col]] <- empty_learning_log()[[col]]
  privacy_safe_learning_log(out[, names(empty_learning_log()), drop = FALSE])
}

import_shared_learning_packet <- function(path, data_root = "data") {
  imported <- read_shared_learning_import_file(path)
  if (nrow(imported) == 0) {
    return(data.frame(Imported_rows = 0, Total_shared_rows = nrow(read_shared_learning_pool(data_root)), Status = "No events in packet."))
  }
  dir.create(learning_root(data_root), recursive = TRUE, showWarnings = FALSE)
  pool <- read_shared_learning_pool(data_root)
  combined <- unique(dplyr::bind_rows(pool, imported))
  write.csv(combined, shared_learning_pool_path(data_root), row.names = FALSE)
  data.frame(
    Imported_rows = nrow(imported),
    Total_shared_rows = nrow(combined),
    Status = "Shared learning imported. Future recommendations can use these aggregate patterns.",
    check.names = FALSE
  )
}

json_escape <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) return("")
  x <- as.character(x[[1]])
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("\n", "\\\\n", x)
  x <- gsub("\r", "\\\\r", x)
  x <- gsub("\t", "\\\\t", x)
  x
}

json_value <- function(x) {
  if (is.null(x) || length(x) == 0) return("null")
  if (length(x) > 1) return(paste0("[", paste(vapply(x, json_value, character(1)), collapse = ","), "]"))
  if (is.na(x)) return("null")
  if (is.logical(x)) return(if (isTRUE(x)) "true" else "false")
  if (is.numeric(x)) return(if (is.finite(x)) as.character(x) else "null")
  paste0("\"", json_escape(x), "\"")
}

json_object <- function(x) {
  fields <- paste0("\"", names(x), "\":", vapply(x, json_value, character(1)))
  paste0("{", paste(fields, collapse = ","), "}")
}

sanitize_training_text <- function(x, max_chars = 1200) {
  shorten_text(gsub("\\s+", " ", first_nonempty(x, "")), max_chars)
}

training_example_from_learning_event <- function(event, kb = NULL, dat = NULL, data_root = "data", include_private = FALSE) {
  route <- first_nonempty(event$route, if (!is.null(kb)) geo_analyzability_route(kb)$Route else "Unknown")
  features <- if (is.null(kb)) {
    list(
      assay_signal = first_nonempty(event$assay_signal, "unknown assay"),
      sample_bin = first_nonempty(event$sample_bin, "unknown sample size"),
      method_terms = split_detected_terms(first_nonempty(event$method_terms, "")),
      learned_tags = split_detected_terms(first_nonempty(event$learned_tags, ""))
    )
  } else {
    study_learning_features(kb)
  }
  assessment <- reanalysis_assessment_table(kb)
  proceed_values <- assessment$Value[assessment$Score == "Proceed recommendation"]
  proceed <- if (length(proceed_values) == 0) suppressWarnings(as.numeric(event$proceed_score)) else suppressWarnings(as.numeric(proceed_values[[1]]))
  claims <- if (is.null(kb)) empty_claim_table() else normalize_claim_table(kb$claims)
  comparison <- if (!is.null(dat)) tryCatch(compare_findings_against_publication(dat, kb), error = function(e) NULL) else NULL
  score <- if (!is.null(dat)) tryCatch(post_analysis_reproducibility_score(dat, kb), error = function(e) data.frame()) else data.frame()
  pitfalls <- tryCatch(sequencing_pitfall_table(dat, kb, data_root), error = function(e) data.frame())
  checklist <- tryCatch(standardized_pipeline_checklist(dat, kb), error = function(e) data.frame())
  signatures <- tryCatch(signature_screen_table(dat, kb), error = function(e) data.frame())
  adaptive <- tryCatch(adaptive_reanalysis_recommendation_table(kb, data_root), error = function(e) data.frame())
  next_plan <- tryCatch(adaptive_ranked_discovery_plan(dat, kb, comparison, data_root), error = function(e) data.frame())

  top_values <- function(df, col, n = 5) {
    if (is.null(df) || nrow(df) == 0 || !col %in% names(df)) return("")
    paste(utils::head(df[[col]], n), collapse = " | ")
  }
  label_outcome <- first_nonempty(event$outcome, "unlabeled")
  label_useful <- suppressWarnings(as.numeric(event$usefulness_rating))
  label_reproduced <- suppressWarnings(as.numeric(event$reproducibility_score))
  target <- paste(
    "Recommendation:",
    top_values(adaptive, "Recommendation", 4),
    "Next analyses:",
    top_values(next_plan, "Analysis", 4),
    "Warnings:",
    top_values(pitfalls, "Pitfall", 4),
    "Signature plan:",
    top_values(signatures, "Recommendation", 3)
  )
  input_summary <- paste(
    "Route:", route,
    "Assay:", features$assay_signal,
    "Sample bin:", features$sample_bin,
    "Methods:", paste(features$method_terms, collapse = "; "),
    "Tags:", paste(features$learned_tags, collapse = "; "),
    "Proceed score:", first_nonempty(as.character(proceed), "NA"),
    "Claims:", nrow(claims),
    "Checklist statuses:", top_values(checklist, "Status", 8),
    if (include_private && !is.null(kb)) paste("Study context:", sanitize_training_text(kb_evidence_blob(kb))) else ""
  )
  list(
    example_id = paste0("seqsurf_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999, 1))),
    privacy_mode = if (include_private) "private_full_context" else "privacy_safe",
    task = "seqsurf_reanalysis_recommendation",
    input = input_summary,
    label_outcome = label_outcome,
    label_usefulness_rating = label_useful,
    label_reproducibility_score = label_reproduced,
    label_claims_checked = suppressWarnings(as.numeric(event$claims_checked)),
    label_supported_claims = suppressWarnings(as.numeric(event$supported_claims)),
    label_missing_claims = suppressWarnings(as.numeric(event$missing_claims)),
    label_direction_mismatches = suppressWarnings(as.numeric(event$direction_mismatches)),
    target = target,
    route = route,
    assay_signal = features$assay_signal,
    sample_bin = features$sample_bin,
    learned_tags = paste(features$learned_tags, collapse = "; "),
    method_terms = paste(features$method_terms, collapse = "; ")
  )
}

training_examples_from_learning_log <- function(kb = NULL, dat = NULL, data_root = "data", include_shared = TRUE, include_private = FALSE) {
  log <- combined_learning_log(data_root, include_shared)
  if (nrow(log) == 0) return(list())
  lapply(seq_len(nrow(log)), function(i) {
    training_example_from_learning_event(log[i, , drop = FALSE], kb, dat, data_root, include_private)
  })
}

training_examples_summary_table <- function(kb = NULL, dat = NULL, data_root = "data", include_shared = TRUE, include_private = FALSE) {
  examples <- training_examples_from_learning_log(kb, dat, data_root, include_shared, include_private)
  if (length(examples) == 0) {
    return(data.frame(
      Metric = c("Training examples", "Privacy mode", "What is needed"),
      Value = c(0, if (include_private) "Private full context" else "Privacy safe", "Save at least one Teach SeqSurf learning event."),
      check.names = FALSE
    ))
  }
  data.frame(
    Metric = c("Training examples", "Privacy mode", "Tasks", "Routes represented", "Labels included"),
    Value = c(
      length(examples),
      if (include_private) "Private full context" else "Privacy safe",
      paste(unique(vapply(examples, function(x) x$task, character(1))), collapse = "; "),
      paste(unique(vapply(examples, function(x) x$route, character(1))), collapse = "; "),
      "outcome, usefulness rating, reproducibility score, claim counts, target recommendation"
    ),
    check.names = FALSE
  )
}

write_training_examples_jsonl <- function(file, kb = NULL, dat = NULL, data_root = "data", include_shared = TRUE, include_private = FALSE) {
  examples <- training_examples_from_learning_log(kb, dat, data_root, include_shared, include_private)
  lines <- vapply(examples, json_object, character(1))
  writeLines(lines, file, useBytes = TRUE)
  file
}

write_training_data_export <- function(file, kb = NULL, dat = NULL, data_root = "data", include_shared = TRUE, include_private = FALSE) {
  tmp_dir <- tempfile("seqsurf_training_export_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  jsonl <- file.path(tmp_dir, "seqsurf_training_examples.jsonl")
  write_training_examples_jsonl(jsonl, kb, dat, data_root, include_shared, include_private)
  manifest <- data.frame(
    Field = c("format", "privacy_mode", "example_count", "intended_use", "excluded_fields"),
    Value = c(
      "SeqSurf training examples JSONL v1",
      if (include_private) "Private full context; may include study context excerpts." else "Privacy safe; no accession, dataset ID, notes, expression matrix, metadata table, or paper text.",
      length(readLines(jsonl, warn = FALSE)),
      "Fine-tuning/ranking data for reanalysis recommendation, warning, reproducibility, and signature-planning models.",
      if (include_private) "Raw expression matrices and full metadata tables." else "GEO accession, dataset_id, free-text notes, raw data, metadata tables, paper text."
    ),
    check.names = FALSE
  )
  write.csv(manifest, file.path(tmp_dir, "TRAINING_EXPORT_MANIFEST.csv"), row.names = FALSE)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(tmp_dir)
  utils::zip(file, files = c("seqsurf_training_examples.jsonl", "TRAINING_EXPORT_MANIFEST.csv"), flags = "-r9Xq")
  file
}

record_learning_event <- function(
    accession = "",
    dataset_id = "",
    event_type = "user_feedback",
    outcome = "",
    usefulness_rating = NA_real_,
    notes = "",
    kb = NULL,
    dat = NULL,
    comparison = NULL,
    data_root = "data") {
  dir.create(learning_root(data_root), recursive = TRUE, showWarnings = FALSE)
  route <- if (is.null(kb)) "Unknown" else first_nonempty(geo_analyzability_route(kb)$Route, "Unknown")
  features <- study_learning_features(kb)
  learned_tags <- unique(c(features$learned_tags, learning_tags_from_text(outcome, notes)))
  assessment <- reanalysis_assessment_table(kb)
  proceed_values <- assessment$Value[assessment$Score == "Proceed recommendation"]
  proceed <- if (length(proceed_values) == 0) NA_real_ else suppressWarnings(as.numeric(proceed_values[[1]]))
  score <- if (!is.null(dat)) {
    suppressWarnings(post_analysis_reproducibility_score(dat, kb)$Value[1])
  } else {
    NA_real_
  }
  comparison <- if (is.null(comparison) && !is.null(dat)) {
    tryCatch(compare_findings_against_publication(dat, kb), error = function(e) NULL)
  } else {
    comparison
  }
  categories <- classify_publication_comparison(comparison)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    accession = normalize_geo_accession(first_nonempty(accession, if (!is.null(kb)) kb$accession else "")),
    dataset_id = first_nonempty(dataset_id, if (!is.null(dat)) dat$id else ""),
    route = route,
    assay_signal = features$assay_signal,
    sample_bin = features$sample_bin,
    method_terms = paste(features$method_terms, collapse = "; "),
    learned_tags = paste(learned_tags, collapse = "; "),
    event_type = event_type,
    outcome = first_nonempty(outcome, ""),
    usefulness_rating = suppressWarnings(as.numeric(usefulness_rating)),
    reproducibility_score = suppressWarnings(as.numeric(score)),
    proceed_score = proceed,
    claims_checked = length(categories),
    supported_claims = sum(categories == "Supported"),
    missing_claims = sum(categories == "Missing"),
    direction_mismatches = sum(categories == "Direction mismatch"),
    notes = first_nonempty(notes, ""),
    check.names = FALSE
  )
  log <- dplyr::bind_rows(read_learning_log(data_root), row)
  write.csv(log, learning_log_path(data_root), row.names = FALSE)
  row
}

learning_summary_table <- function(data_root = "data", include_shared = TRUE) {
  log <- combined_learning_log(data_root, include_shared)
  if (nrow(log) == 0) {
    return(data.frame(
      Metric = c("Learning events", "Local events", "Shared events", "Datasets seen", "Average usefulness", "Best current lesson"),
      Value = c(0, 0, 0, 0, "No ratings yet", "SeqSurf will start adapting once users save feedback after datasets."),
      check.names = FALSE
    ))
  }
  avg_rating <- suppressWarnings(mean(as.numeric(log$usefulness_rating), na.rm = TRUE))
  avg_text <- if (is.nan(avg_rating)) "No ratings yet" else round(avg_rating, 2)
  common_outcome <- names(sort(table(log$outcome[nzchar(log$outcome)]), decreasing = TRUE))[1]
  if (is.na(common_outcome) || !nzchar(common_outcome)) common_outcome <- "Keep collecting feedback."
  local_n <- nrow(read_learning_log(data_root))
  shared_n <- nrow(read_shared_learning_pool(data_root))
  data.frame(
    Metric = c("Learning events", "Local events", "Shared events", "Datasets seen", "Average usefulness", "Most common outcome"),
    Value = c(nrow(log), local_n, shared_n, length(unique(log$accession[nzchar(log$accession)])), avg_text, common_outcome),
    check.names = FALSE
  )
}

learning_lessons_table <- function(data_root = "data", include_shared = TRUE) {
  log <- combined_learning_log(data_root, include_shared)
  if (nrow(log) == 0) {
    return(data.frame(
      Lesson = "No learned lessons yet",
      Evidence = "Save feedback after a dataset to teach SeqSurf what helped, failed, or was scientifically useful.",
      Action = "Use Discovery -> Teach SeqSurf after validation.",
      check.names = FALSE
    ))
  }
  route_stats <- log |>
    dplyr::group_by(.data$route) |>
    dplyr::summarise(
      Events = dplyr::n(),
      Mean_usefulness = suppressWarnings(mean(as.numeric(.data$usefulness_rating), na.rm = TRUE)),
      Mean_reproducibility = suppressWarnings(mean(as.numeric(.data$reproducibility_score), na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$Events))
  route_stats$Mean_usefulness[is.nan(route_stats$Mean_usefulness)] <- NA_real_
  route_stats$Mean_reproducibility[is.nan(route_stats$Mean_reproducibility)] <- NA_real_
  data.frame(
    Lesson = paste("Route pattern:", route_stats$route),
    Evidence = paste0(
      route_stats$Events, " events; mean usefulness ",
      ifelse(is.na(route_stats$Mean_usefulness), "unrated", round(route_stats$Mean_usefulness, 2)),
      "; mean reproducibility ",
      ifelse(is.na(route_stats$Mean_reproducibility), "not scored", round(route_stats$Mean_reproducibility, 1))
    ),
    Action = ifelse(
      !is.na(route_stats$Mean_usefulness) & route_stats$Mean_usefulness >= 4,
      "Prefer this route for similar datasets when data availability checks pass.",
      "Keep collecting examples before changing recommendations."
    ),
    check.names = FALSE
  )
}

study_learning_features <- function(kb = NULL) {
  text <- tolower(if (is.null(kb)) "" else kb_evidence_blob(kb))
  sample_count <- suppressWarnings(as.integer(first_nonempty(if (!is.null(kb)) kb$geo$sample_count else NA_integer_, NA_integer_)))
  sample_bin <- if (is.na(sample_count)) {
    "unknown sample size"
  } else if (sample_count >= 100) {
    "large cohort"
  } else if (sample_count >= 30) {
    "medium cohort"
  } else {
    "small cohort"
  }
  assay_signal <- if (grepl("rna-seq|rnaseq|sequenc|hiseq|novaseq|nextseq|count", text)) {
    "rna-seq/count"
  } else if (grepl("microarray|affymetrix|agilent|beadchip|array|probe", text)) {
    "array/processed"
  } else {
    "unknown assay"
  }
  method_terms <- if (is.null(kb)) character() else split_detected_terms(first_nonempty(kb$extracted$statistical_methods, ""))
  tags <- c(
    if (grepl("survival|progression|relapse|response|death|prognos|risk", text)) "outcome",
    if (grepl("batch|center|platform|cohort|validation|replication", text)) "batch",
    if (grepl("immune|microenvironment|cell type|deconvolution", text)) "immune",
    if (grepl("epithelial|mesenchymal|emt", text)) "epithelial_signature",
    if (grepl("count|counts|featurecounts|htseq|read count", text)) "count_matrix",
    if (grepl("fastq|sra|bam|sam", text)) "raw_pipeline"
  )
  list(
    assay_signal = assay_signal,
    sample_bin = sample_bin,
    method_terms = method_terms,
    learned_tags = unique(tags)
  )
}

learning_tags_from_text <- function(...) {
  text <- tolower(paste(..., collapse = " "))
  tags <- c(
    if (grepl("metadata|group|label|phenotype|sample", text)) "metadata_curation",
    if (grepl("count|counts|featurecounts|htseq|voom|deseq2|edger", text)) "count_matrix",
    if (grepl("raw|fastq|sra|bam|salmon|star", text)) "raw_pipeline",
    if (grepl("survival|outcome|response|relapse|progression|risk", text)) "outcome",
    if (grepl("immune|cell type|microenvironment|deconvolution", text)) "immune",
    if (grepl("epithelial|mesenchymal|emt", text)) "epithelial_signature",
    if (grepl("batch|platform|center|confound", text)) "batch",
    if (grepl("disagree|mismatch|different|failed|not reproduce", text)) "disagreement",
    if (grepl("reproduced|matched|supported", text)) "reproduced"
  )
  unique(tags)
}

learning_prior_for_study <- function(kb = NULL, data_root = "data", include_shared = TRUE) {
  log <- combined_learning_log(data_root, include_shared)
  if (nrow(log) == 0) {
    return(data.frame(
      Signal = "No local learning yet",
      Evidence = "No feedback events have been saved.",
      Recommendation_bias = "Use deterministic scorecard only.",
      Confidence = "None",
      check.names = FALSE
    ))
  }
  route <- if (is.null(kb)) "Unknown" else first_nonempty(geo_analyzability_route(kb)$Route, "Unknown")
  features <- study_learning_features(kb)
  same_route <- log[log$route == route, , drop = FALSE]
  assay_mask <- log$assay_signal == features$assay_signal & !is.na(log$assay_signal) & nzchar(log$assay_signal)
  assay_mask[is.na(assay_mask)] <- FALSE
  same_assay <- log[assay_mask, , drop = FALSE]
  wanted_tags <- features$learned_tags
  tag_pattern <- paste(wanted_tags, collapse = "|")
  same_tags <- if (length(wanted_tags) == 0) {
    log[FALSE, , drop = FALSE]
  } else {
    tag_mask <- grepl(tag_pattern, log$learned_tags, ignore.case = TRUE)
    tag_mask[is.na(tag_mask)] <- FALSE
    log[tag_mask, , drop = FALSE]
  }
  relevant <- unique(dplyr::bind_rows(same_route, same_assay, same_tags))
  if (nrow(relevant) == 0) {
    return(data.frame(
      Signal = paste("No close local precedent for", route),
      Evidence = paste(nrow(log), "learning events exist, but none match this route/assay/tags closely."),
      Recommendation_bias = "Use deterministic scorecard and ask user to save feedback after this dataset.",
      Confidence = "Low",
      check.names = FALSE
    ))
  }
  avg_use <- suppressWarnings(mean(as.numeric(relevant$usefulness_rating), na.rm = TRUE))
  avg_rep <- suppressWarnings(mean(as.numeric(relevant$reproducibility_score), na.rm = TRUE))
  avg_use_text <- if (is.nan(avg_use)) "unrated" else round(avg_use, 2)
  avg_rep_text <- if (is.nan(avg_rep)) "not scored" else round(avg_rep, 1)
  outcomes <- relevant$outcome[nzchar(relevant$outcome)]
  top_outcome <- if (length(outcomes) == 0) "No outcome labels yet" else names(sort(table(outcomes), decreasing = TRUE))[[1]]
  failures <- sum(grepl("failed|manual|low-value", relevant$outcome, ignore.case = TRUE))
  successes <- sum(grepl("useful|reproduced|interesting", relevant$outcome, ignore.case = TRUE))
  bias <- if (!is.nan(avg_use) && avg_use >= 4 && successes >= failures) {
    "Upgrade confidence: similar local datasets were useful."
  } else if (failures > successes) {
    "Add caution: similar local datasets often needed manual review or failed."
  } else if (!is.nan(avg_rep) && avg_rep < 45) {
    "Prioritize discrepancy analysis: similar datasets had weak reproduction."
  } else {
    "Neutral: local precedents are mixed or sparse."
  }
  data.frame(
    Signal = paste("Matched local precedent:", route),
    Evidence = paste0(
      nrow(relevant), " relevant events; average usefulness ", avg_use_text,
      "; average reproducibility ", avg_rep_text,
      "; most common outcome: ", top_outcome
    ),
    Recommendation_bias = bias,
    Confidence = if (nrow(relevant) >= 5) "Moderate" else "Early",
    check.names = FALSE
  )
}

adaptive_reanalysis_recommendation_table <- function(kb, data_root = "data") {
  base <- reanalysis_recommendation_table(kb)
  base$Source <- "Scorecard"
  prior <- learning_prior_for_study(kb, data_root)
  learned_rows <- list()
  bias <- first_nonempty(prior$Recommendation_bias[[1]], "")
  if (grepl("Upgrade confidence", bias, fixed = TRUE)) {
    learned_rows[[length(learned_rows) + 1]] <- data.frame(
      Recommendation = "Proceed sooner for similar high-performing datasets",
      Why = prior$Evidence[[1]],
      Status = "Learned boost",
      Source = "Local learning",
      check.names = FALSE
    )
  } else if (grepl("Add caution", bias, fixed = TRUE)) {
    learned_rows[[length(learned_rows) + 1]] <- data.frame(
      Recommendation = "Add a manual metadata/file audit before full reanalysis",
      Why = prior$Evidence[[1]],
      Status = "Learned caution",
      Source = "Local learning",
      check.names = FALSE
    )
  } else if (grepl("discrepancy", bias, ignore.case = TRUE)) {
    learned_rows[[length(learned_rows) + 1]] <- data.frame(
      Recommendation = "Plan for paper-disagreement analysis from the start",
      Why = prior$Evidence[[1]],
      Status = "Learned priority",
      Source = "Local learning",
      check.names = FALSE
    )
  } else {
    learned_rows[[length(learned_rows) + 1]] <- data.frame(
      Recommendation = "Collect feedback after this dataset",
      Why = prior$Evidence[[1]],
      Status = "Learning",
      Source = "Local learning",
      check.names = FALSE
    )
  }
  tags <- study_learning_features(kb)$learned_tags
  if ("epithelial_signature" %in% tags) {
    learned_rows[[length(learned_rows) + 1]] <- data.frame(
      Recommendation = "Track epithelial/EMT signature behavior explicitly",
      Why = "Study context or local notes mention epithelial/mesenchymal biology; compare signature scores, DEG direction, and pathway enrichment across groups.",
      Status = "Adaptive hypothesis",
      Source = "Feature-aware learning",
      check.names = FALSE
    )
  }
  if ("outcome" %in% tags) {
    learned_rows[[length(learned_rows) + 1]] <- data.frame(
      Recommendation = "Prioritize outcome-aware validation if metadata supports it",
      Why = "Outcome-like language is present and prior feedback can make these datasets especially valuable when labels are curated.",
      Status = "Adaptive hypothesis",
      Source = "Feature-aware learning",
      check.names = FALSE
    )
  }
  learned <- dplyr::bind_rows(learned_rows)
  learned$Priority <- seq_len(nrow(learned))
  out <- dplyr::bind_rows(
    learned[, c("Priority", "Recommendation", "Why", "Status", "Source"), drop = FALSE],
    base[, c("Priority", "Recommendation", "Why", "Status", "Source"), drop = FALSE]
  )
  out$Priority <- seq_len(nrow(out))
  out
}

geo_analyzability_route <- function(kb, expression_set_status = NULL) {
  if (is.null(kb)) {
    return(data.frame(
      Route = "Not ready",
      Can_run_in_app = "No",
      Confidence = 0,
      Evidence = "Train study context first.",
      Recommended_next_step = "Enter a GEO accession and click Get study info.",
      Blocking_reason = "No GEO/PubMed metadata available.",
      check.names = FALSE
    ))
  }
  supp <- if (is.null(kb$geo$supplementary_files)) character() else kb$geo$supplementary_files
  text <- tolower(paste(kb_evidence_blob(kb), supp, collapse = " "))
  has_series <- grepl("series matrix|processed|normalized|rma|expression profiling by array|expression profiling by high throughput sequencing", text)
  has_count <- grepl("count matrix|counts matrix|raw counts|htseq|featurecounts|read counts|gene counts", text)
  has_raw <- grepl("fastq|sra|bam|sam|cel.gz|idat", text)
  has_supp <- length(supp) > 0
  has_methods <- length(split_detected_terms(kb$extracted$statistical_methods)) > 0
  has_groups <- nzchar(first_nonempty(kb$extracted$sample_groups, ""))
  eset_ok <- identical(expression_set_status, "available")
  eset_failed <- identical(expression_set_status, "unavailable")

  route <- "Manual curation required"
  can_run <- "No"
  confidence <- 25
  evidence <- "Public expression format is unclear from GEO/PubMed context."
  next_step <- "Review GEO files manually or paste methods/supplement notes."
  block <- "No processed expression, count matrix, or raw-file signal was confidently detected."

  if (eset_ok || (has_series && !has_count && !has_raw)) {
    route <- "In-app GEO processed-matrix reanalysis"
    can_run <- "Yes"
    confidence <- if (eset_ok) 90 else 70
    evidence <- if (eset_ok) "GEOquery returned an analyzable ExpressionSet/series matrix." else "Study context suggests processed or normalized expression suitable for GEO series-matrix import."
    next_step <- "Use Analyze, review the proposed design, then run reanalysis."
    block <- "None for processed-matrix route; still review sample grouping and paper-method match."
  } else if (has_count) {
    route <- "Count-matrix pipeline then app import"
    can_run <- "Partial"
    confidence <- 75
    evidence <- "Study context or supplementary links suggest gene/read count files."
    next_step <- "Run a count-aware DESeq2/edgeR/limma-voom script outside Shiny, then import the app dataset contract."
    block <- "Interactive Shiny should not silently choose normalization/modeling for arbitrary count files."
  } else if (has_raw) {
    route <- "Raw sequencing pipeline handoff"
    can_run <- "No"
    confidence <- 80
    evidence <- "Study context or supplementary links suggest FASTQ/SRA/BAM/CEL/IDAT-style raw data."
    next_step <- "Generate a paper-to-code brief and run a reproducible external workflow before loading outputs into SeqSurf AI."
    block <- "Raw processing is compute-heavy and requires external tools, references, and workflow provenance."
  } else if (has_supp) {
    route <- "Supplementary-file review"
    can_run <- "Maybe"
    confidence <- 55
    evidence <- paste(length(supp), "supplementary file links exist, but file type is not clearly classified.")
    next_step <- "Review supplementary links for processed matrices or count tables before running automated import."
    block <- "File type and matrix shape are unknown."
  }

  if (eset_failed && can_run == "Yes") {
    route <- "Processed-matrix route failed"
    can_run <- "No"
    confidence <- 85
    evidence <- "GEOquery did not return an analyzable ExpressionSet despite processed-data language."
    next_step <- "Use supplementary-file triage or an external preprocessing pipeline."
    block <- "Automated GEO series-matrix import failed."
  }

  if (!has_methods || !has_groups) {
    block <- paste(
      block,
      if (!has_methods) "Methods are not clearly extracted." else "",
      if (!has_groups) "Sample groups are not clearly extracted." else ""
    )
    confidence <- max(15, confidence - 10)
  }

  data.frame(
    Route = route,
    Can_run_in_app = can_run,
    Confidence = confidence,
    Evidence = evidence,
    Recommended_next_step = next_step,
    Blocking_reason = trimws(block),
    check.names = FALSE
  )
}

geo_analyzability_checks <- function(kb) {
  if (is.null(kb)) {
    return(data.frame(Check = "Study context", Signal = "Missing", Status = "Not ready", Impact = "Cannot classify analyzability yet.", check.names = FALSE))
  }
  supp <- if (is.null(kb$geo$supplementary_files)) character() else kb$geo$supplementary_files
  text <- tolower(paste(kb_evidence_blob(kb), supp, collapse = " "))
  checks <- data.frame(
    Check = c("Processed/series-matrix signal", "Count-matrix signal", "Raw-file signal", "Supplementary links", "Methods extracted", "Sample groups extracted"),
    Signal = c(
      if (grepl("series matrix|processed|normalized|rma", text)) "Detected" else "Not detected",
      if (grepl("count matrix|counts matrix|raw counts|htseq|featurecounts|read counts|gene counts", text)) "Detected" else "Not detected",
      if (grepl("fastq|sra|bam|sam|cel.gz|idat", text)) "Detected" else "Not detected",
      paste(length(supp), "links"),
      if (length(split_detected_terms(kb$extracted$statistical_methods)) > 0) "Detected" else "Not detected",
      if (nzchar(first_nonempty(kb$extracted$sample_groups, ""))) "Detected" else "Not detected"
    ),
    Status = c(
      if (grepl("series matrix|processed|normalized|rma", text)) "Supports in-app route" else "No support",
      if (grepl("count matrix|counts matrix|raw counts|htseq|featurecounts|read counts|gene counts", text)) "External count pipeline" else "No count route",
      if (grepl("fastq|sra|bam|sam|cel.gz|idat", text)) "External raw pipeline" else "No raw route",
      if (length(supp) > 0) "Review" else "None parsed",
      if (length(split_detected_terms(kb$extracted$statistical_methods)) > 0) "Usable" else "Needs extraction",
      if (nzchar(first_nonempty(kb$extracted$sample_groups, ""))) "Usable" else "Needs extraction"
    ),
    Impact = c(
      "Best case for automated Shiny reanalysis.",
      "Requires count-aware modeling before app import.",
      "Requires external workflow before app import.",
      "May contain matrices or raw archives.",
      "Needed to decide exact vs approximate reproduction.",
      "Needed to recreate paper contrasts."
    ),
    check.names = FALSE
  )
  checks
}

ai_status_badge <- function(status = "Ready", detail = "Enter a GEO accession or upload a dataset to begin.") {
  tags$div(
    class = "seq-avatar-panel",
    tags$div(class = "seq-avatar", tags$i(class = "fa-solid fa-water")),
    tags$div(
      class = "seq-avatar-copy",
      tags$strong("SeqSurfer"),
      tags$span(status),
      tags$p(detail)
    )
  )
}

score_band <- function(score) {
  if (is.na(score)) return("Not available")
  if (score >= 80) return("High")
  if (score >= 60) return("Moderate-high")
  if (score >= 40) return("Moderate")
  if (score >= 20) return("Low-moderate")
  "Low"
}

clip_score <- function(x) {
  max(0, min(100, round(x)))
}

kb_evidence_blob <- function(kb) {
  if (is.null(kb)) return("")
  paste(
    first_nonempty(kb$geo$title, ""),
    first_nonempty(kb$geo$summary, ""),
    first_nonempty(kb$geo$overall_design, ""),
    first_nonempty(kb$pubmed$title, ""),
    first_nonempty(kb$pubmed$abstract, ""),
    related_publications_evidence_text(if (!is.null(kb$related_publications)) kb$related_publications else NULL),
    first_nonempty(kb$extra_papers, ""),
    first_nonempty(kb$extracted$study_objective, ""),
    first_nonempty(kb$extracted$cohort_characteristics, ""),
    first_nonempty(kb$extracted$sample_groups, ""),
    first_nonempty(kb$extracted$statistical_methods, ""),
    first_nonempty(kb$extracted$key_findings, ""),
    first_nonempty(kb$extracted$reported_pathways, ""),
    first_nonempty(kb$extracted$reported_biomarkers, ""),
    first_nonempty(kb$extracted$study_limitations, ""),
    collapse = "\n"
  )
}

reanalysis_score_inputs <- function(kb) {
  text <- kb_evidence_blob(kb)
  text_low <- tolower(text)
  claims <- if (is.null(kb)) empty_claim_table() else normalize_claim_table(kb$claims)
  methods <- split_detected_terms(if (is.null(kb)) "" else kb$extracted$statistical_methods)
  biomarkers <- split_detected_terms(if (is.null(kb)) "" else kb$extracted$reported_biomarkers)
  pathways <- split_detected_terms(if (is.null(kb)) "" else kb$extracted$reported_pathways)
  sample_count <- suppressWarnings(as.integer(first_nonempty(if (!is.null(kb)) kb$geo$sample_count else NA_integer_, NA_integer_)))
  supp_files <- if (is.null(kb$geo$supplementary_files)) character() else kb$geo$supplementary_files
  platform_text <- paste(if (!is.null(kb$geo$platform_accessions)) kb$geo$platform_accessions else "", text, collapse = " ")

  list(
    text = text,
    text_low = text_low,
    claims = claims,
    methods = methods,
    biomarkers = biomarkers,
    pathways = pathways,
    sample_count = sample_count,
    supplementary_file_count = length(supp_files),
    has_pubmed = !is.null(kb$pubmed$pubmed_id) && !is.na(kb$pubmed$pubmed_id) && nzchar(as.character(kb$pubmed$pubmed_id)),
    has_abstract = !is.null(kb$pubmed$abstract) && !is.na(kb$pubmed$abstract) && nzchar(kb$pubmed$abstract),
    has_extra_papers = !is.null(kb$extra_papers) && !is.na(kb$extra_papers) && nzchar(trimws(kb$extra_papers)),
    mentions_counts = grepl("count|counts|rna-seq|rnaseq|sequenc", text_low),
    mentions_processed = grepl("processed|normalized|series matrix|microarray|array|rma", text_low),
    mentions_outcome = grepl("survival|progression|relapse|response|death|prognos|risk", text_low),
    mentions_batch = grepl("batch|center|platform|cohort|validation|replication", text_low),
    platform_is_array = grepl("microarray|affymetrix|agilent|illumina array|array", platform_text, ignore.case = TRUE),
    platform_is_rnaseq = grepl("rna-seq|rnaseq|sequenc|hiseq|novaseq|nextseq", platform_text, ignore.case = TRUE)
  )
}

assess_reanalysis_opportunity <- function(kb) {
  if (is.null(kb)) {
    return(list(
      scores = data.frame(
        Score = "Readiness",
        Value = 0,
        Band = "Not available",
        Rationale = "Enter a GEO accession or train a study knowledge base first.",
        check.names = FALSE
      ),
      recommendations = data.frame(
        Priority = 1,
        Recommendation = "Train study context first",
        Why = "The app needs GEO/PubMed context before it can score reanalysis value.",
        Status = "Not ready",
        check.names = FALSE
      )
    ))
  }

  x <- reanalysis_score_inputs(kb)
  sample_score <- if (is.na(x$sample_count)) 5 else if (x$sample_count >= 100) 25 else if (x$sample_count >= 30) 18 else if (x$sample_count >= 10) 10 else 4
  evidence_score <- 0
  evidence_score <- evidence_score + if (x$has_pubmed) 12 else 0
  evidence_score <- evidence_score + if (x$has_abstract) 10 else 0
  evidence_score <- evidence_score + min(10, x$supplementary_file_count * 2)
  evidence_score <- evidence_score + if (x$has_extra_papers) 8 else 0
  claim_score <- min(20, nrow(x$claims) * 4 + length(x$biomarkers) * 2 + length(x$pathways) * 2)
  method_score <- min(15, length(x$methods) * 3)

  data_availability <- clip_score(20 + sample_score + evidence_score + if (x$mentions_counts || x$mentions_processed) 15 else 0)
  reproducibility <- clip_score(25 + evidence_score + claim_score + method_score - if (x$platform_is_array) 5 else 0)
  usefulness <- clip_score(25 + sample_score + claim_score + if (x$mentions_outcome) 12 else 0 + if (x$mentions_batch) 6 else 0)
  novelty <- clip_score(35 + if (x$mentions_outcome) 15 else 0 + if (x$mentions_batch) 10 else 0 + if (x$platform_is_rnaseq) 8 else 0 - if (length(x$methods) >= 5) 5 else 0)
  technical_risk <- clip_score(35 + if (x$platform_is_array) 15 else 0 + if (!x$mentions_counts && !x$mentions_processed) 20 else 0 - sample_score / 2 - evidence_score / 3)
  proceed_score <- clip_score((data_availability + reproducibility + usefulness + novelty + (100 - technical_risk)) / 5)

  scores <- data.frame(
    Score = c(
      "Data availability",
      "Original-method reproducibility",
      "Reanalysis usefulness",
      "Novelty opportunity",
      "Technical risk",
      "Proceed recommendation"
    ),
    Value = c(data_availability, reproducibility, usefulness, novelty, technical_risk, proceed_score),
    Band = vapply(c(data_availability, reproducibility, usefulness, novelty, technical_risk, proceed_score), score_band, character(1)),
    Rationale = c(
      paste0("GEO reports ", first_nonempty(as.character(x$sample_count), "unknown"), " samples; PubMed/supplement/paper context contributes ", evidence_score, " evidence points."),
      paste0("Detected ", nrow(x$claims), " structured claims and ", length(x$methods), " method terms for publication matching."),
      paste0("Reanalysis value is driven by sample count, reported biomarkers/pathways, and outcome/batch terms in the study context."),
      paste0("Opportunity increases when outcome, validation, batch, or RNA-seq context suggests modern extensions beyond the original analysis."),
      paste0("Risk increases when public data availability is unclear or platform/method details imply extra preprocessing."),
      if (proceed_score >= 70) "Proceed: strong candidate for automated reanalysis." else if (proceed_score >= 50) "Proceed with review: likely useful, but inspect design and data availability first." else "Caution: use manual review before spending time on full reanalysis."
    ),
    check.names = FALSE
  )

  recommendations <- reanalysis_recommendations_from_inputs(x)
  list(scores = scores, recommendations = recommendations)
}

reanalysis_assessment_table <- function(kb) {
  assess_reanalysis_opportunity(kb)$scores
}

reanalysis_recommendations_from_inputs <- function(x) {
  rows <- list()
  add <- function(recommendation, why, status = "Recommended") {
    rows[[length(rows) + 1]] <<- data.frame(
      Recommendation = recommendation,
      Why = why,
      Status = status,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (length(x$biomarkers) > 0) {
    add(
      "Reproduce reported biomarker directions",
      paste("Check DEG direction and significance for:", paste(utils::head(x$biomarkers, 8), collapse = ", "))
    )
  }
  if (length(x$pathways) > 0) {
    add(
      "Reproduce reported pathway signals",
      paste("Run preranked GSEA and compare NES direction for:", paste(utils::head(x$pathways, 6), collapse = ", "))
    )
  }
  if (x$mentions_outcome) {
    add(
      "Add outcome-aware reanalysis",
      "Study context mentions survival, progression, response, risk, or prognosis; test whether expression/pathway scores track outcome metadata."
    )
  }
  if (x$mentions_batch) {
    add(
      "Run sensitivity and batch checks",
      "Study context mentions platform, batch, validation, center, or cohort structure; test whether PCA/DEG results are robust to these variables."
    )
  }
  if (x$platform_is_rnaseq) {
    add(
      "Use modern RNA-seq pathway and gene-set layers",
      "RNA-seq context supports updated GSEA, gene-set scoring, regulon/network screens, and sample-level pathway analysis."
    )
  }
  if (x$platform_is_array) {
    add(
      "Review probe-to-gene mapping before interpreting genes",
      "Array platforms often need careful probe annotation and duplicate-gene handling before DEG/GSEA interpretation.",
      "Caution"
    )
  }
  if (!x$mentions_counts && !x$mentions_processed) {
    add(
      "Confirm expression data availability",
      "The study context does not clearly mention counts, processed expression, or a series matrix; GEO import may require manual raw-file preparation.",
      "Check first"
    )
  }
  add(
    "Compare reproduced results against extracted claims",
    "Use publication claims as a verification layer, then flag supported, missing, weak, or direction-mismatched findings."
  )
  add(
    "Generate a modern follow-up analysis plan",
    "After paper-style reproduction, prioritize unexplored metadata contrasts, immune/cell-composition screens, survival screens, and cross-cohort validation where data permit."
  )

  out <- dplyr::bind_rows(rows)
  out$Priority <- seq_len(nrow(out))
  out[, c("Priority", "Recommendation", "Why", "Status"), drop = FALSE]
}

reanalysis_recommendation_table <- function(kb) {
  assess_reanalysis_opportunity(kb)$recommendations
}

build_reanalysis_narrative_prompt <- function(kb) {
  if (is.null(kb)) stop("Train study context before generating an AI narrative.", call. = FALSE)
  assessment <- reanalysis_assessment_table(kb)
  recommendations <- reanalysis_recommendation_table(kb)
  score_lines <- paste0(
    "- ", assessment$Score, ": ", assessment$Value, "/100 (", assessment$Band, "). ",
    assessment$Rationale,
    collapse = "\n"
  )
  recommendation_lines <- paste0(
    "- ", recommendations$Priority, ". ", recommendations$Recommendation,
    " [", recommendations$Status, "]: ", recommendations$Why,
    collapse = "\n"
  )
  summary_lines <- paste0(
    "- ", study_kb_summary_table(kb)$Field, ": ", study_kb_summary_table(kb)$Value,
    collapse = "\n"
  )

  paste(
    "You are a careful transcriptomics reanalysis advisor embedded in SeqSurf AI.",
    "Write a concise, publication-aware reanalysis decision memo for a scientist deciding whether this GEO dataset is worth reanalyzing.",
    "Use only the supplied study context, deterministic scorecard, and recommendation table.",
    "Do not invent paper findings, available files, methods, or biological conclusions.",
    "If evidence is weak, say exactly what needs manual checking before proceeding.",
    "Be specific about what new value reanalysis could add compared with the original study.",
    "Name concrete analyses or hypotheses, not generic phrases like 'updated techniques' unless you say which technique and why.",
    "For array studies, frame RNA-seq as a cross-platform validation idea only if RNA-seq data are not actually available in this GEO record.",
    "Use Markdown with exactly five level-4 headings. Do not include a title heading.",
    "Keep each section to 1-3 short sentences. Use bullets in the Questions to test and Risks/checks sections.",
    "",
    "Use these headings exactly:",
    "#### Why reanalyze this dataset?",
    "#### What has already been done?",
    "#### Questions to test next",
    "#### Best first reanalysis",
    "#### Risks and checks",
    "",
    "For 'Why reanalyze this dataset?', explicitly mention the strongest drivers from the scorecard: sample size, public matrix availability, dated platform/methods, reported biomarkers/pathways, outcome/metadata opportunity, or reproducibility risk.",
    "For 'Questions to test next', provide 2-4 bullet questions that are directly testable from the GEO metadata/expression or clearly label when a question requires added clinical metadata.",
    "For 'Best first reanalysis', give one concrete first analysis and the expected decision it would support.",
    "",
    "STUDY SUMMARY:",
    summary_lines,
    "",
    "DETERMINISTIC SCORECARD:",
    score_lines,
    "",
    "RECOMMENDED REANALYSIS TARGETS:",
    recommendation_lines,
    "",
    "SOURCE CONTEXT EXCERPT:",
    substr(kb_evidence_text(kb, max_chars = 18000), 1, 18000),
    sep = "\n"
  )
}

generate_reanalysis_narrative <- function(kb, api_key = "", model = "gpt-4o-mini") {
  if (is.null(kb)) stop("Train study context before generating an AI narrative.", call. = FALSE)
  if (!openai_available(api_key)) {
    stop("AI narrative requires OPENAI_API_KEY or a pasted API key plus optional packages httr2 and jsonlite.", call. = FALSE)
  }
  response <- httr2::request("https://api.openai.com/v1/responses") |>
    httr2::req_auth_bearer_token(openai_api_key(api_key)) |>
    httr2::req_body_json(
      list(
        model = model,
        input = build_reanalysis_narrative_prompt(kb)
      ),
      auto_unbox = TRUE
    ) |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_timeout(120) |>
    httr2::req_perform()
  body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  if (httr2::resp_status(response) >= 300) {
    message <- first_nonempty(body$error$message, paste("OpenAI API request failed with status", httr2::resp_status(response)))
    stop(message, call. = FALSE)
  }
  text <- extract_response_text(body)
  if (!nzchar(text)) stop("OpenAI API returned no narrative text.", call. = FALSE)
  text
}

verify_dataset_against_kb <- function(dat, kb) {
  if (is.null(kb)) {
    return(data.frame(Check = "Knowledge base", Result = "No study knowledge base selected.", Status = "Needs training"))
  }

  meta <- dat$metadata
  info <- dat$info
  geo_mentions <- unique(na.omit(c(
    info$dataset_id,
    info$display_name,
    info$notes,
    if ("geo_accession" %in% names(meta)) meta$geo_accession else NA_character_
  )))
  accession_match <- any(grepl(tolower(kb$accession), tolower(geo_mentions), fixed = TRUE))
  sample_delta <- if (is.null(kb$geo$sample_count) || is.na(kb$geo$sample_count)) NA_integer_ else nrow(meta) - kb$geo$sample_count
  platform_note <- if (length(kb$geo$platform_accessions) == 0) "No GEO platform parsed" else paste(kb$geo$platform_accessions, collapse = "; ")

  data.frame(
    Check = c("Accession trace", "Sample count", "Default comparison", "Platform record", "Expression unit"),
    Result = c(
      if (accession_match) paste("Current dataset mentions", kb$accession) else paste("Current dataset does not explicitly mention", kb$accession),
      if (is.na(sample_delta)) "GEO sample count unavailable" else paste0("App has ", nrow(meta), " samples; GEO reports ", kb$geo$sample_count, " samples; delta ", sample_delta),
      paste(info$default_contrast, "using", info$default_group_col),
      platform_note,
      first_nonempty(info$expression_unit, "Not specified")
    ),
    Status = c(
      if (accession_match) "Pass" else "Review",
      if (is.na(sample_delta) || sample_delta == 0) "Pass" else "Review",
      "Review",
      "Review",
      "Review"
    ),
    check.names = FALSE
  )
}

compare_findings_against_publication <- function(dat, kb, padj_cutoff = 0.1, lfc_cutoff = 0.5) {
  if (is.null(kb)) {
    return(data.frame(
      Finding = "Publication knowledge",
      Evidence_in_app = "Train or select a study knowledge base first.",
      Interpretation = "Not ready",
      check.names = FALSE
    ))
  }

  claims <- normalize_claim_table(kb$claims)
  if (nrow(claims) == 0) {
    claims <- dplyr::bind_rows(
      extract_structured_claims(kb$extracted$reported_biomarkers, "extracted biomarker summary"),
      extract_structured_claims(kb$extracted$reported_pathways, "extracted pathway summary")
    )
  }
  deg <- dat$deg
  gsea <- dat$gsea
  rows <- list()

  gene_claims <- claims[claims$claim_type %in% c("gene", "biomarker"), , drop = FALSE]
  pathway_claims <- claims[claims$claim_type == "pathway", , drop = FALSE]

  for (i in seq_len(nrow(gene_claims))) {
    gene <- gene_claims$entity[[i]]
    matches <- deg[toupper(deg$gene) == toupper(gene), , drop = FALSE]
    if (nrow(matches) == 0) {
      rows[[length(rows) + 1]] <- data.frame(
        Finding = paste("Reported biomarker:", gene),
        Evidence_in_app = "Gene not found in DEG table.",
        Interpretation = "Missing from app result",
        check.names = FALSE
      )
    } else {
      best <- matches[order(matches$padj, -abs(matches$log2FoldChange), na.last = TRUE), , drop = FALSE][1, ]
      app_direction <- ifelse(best$log2FoldChange > 0, "up", ifelse(best$log2FoldChange < 0, "down", "flat"))
      paper_direction <- tolower(first_nonempty(gene_claims$reported_direction[[i]], "not specified"))
      direction_status <- if (paper_direction %in% c("up", "higher", "positive") && app_direction == "up") {
        "direction agrees"
      } else if (paper_direction %in% c("down", "lower", "negative") && app_direction == "down") {
        "direction agrees"
      } else if (paper_direction == "not specified") {
        "direction not specified in paper claim"
      } else {
        "direction mismatch"
      }
      rows[[length(rows) + 1]] <- data.frame(
        Finding = paste("Reported biomarker:", gene),
        Evidence_in_app = paste0(
          best$contrast,
          "; log2FC=", signif(best$log2FoldChange, 3),
          "; padj=", signif(best$padj, 3),
          "; paper direction=", paper_direction,
          "; app direction=", app_direction
        ),
        Interpretation = if (direction_status == "direction mismatch") {
          "Direction mismatch"
        } else if (!is.na(best$padj) && best$padj <= padj_cutoff && abs(best$log2FoldChange) >= lfc_cutoff) {
          paste("Supported by app DEG result;", direction_status)
        } else {
          paste("Detected but not significant at selected thresholds;", direction_status)
        },
        check.names = FALSE
      )
    }
  }

  for (i in seq_len(nrow(pathway_claims))) {
    term <- pathway_claims$entity[[i]]
    pattern <- tolower(term)
    matches <- gsea[grepl(pattern, tolower(gsea$pathway), fixed = TRUE), , drop = FALSE]
    if (nrow(matches) == 0) {
      rows[[length(rows) + 1]] <- data.frame(
        Finding = paste("Reported pathway:", term),
        Evidence_in_app = "No matching pathway name found in GSEA table.",
        Interpretation = "Missing from app pathway result",
        check.names = FALSE
      )
    } else {
      best <- matches[order(matches$padj, -abs(matches$NES), na.last = TRUE), , drop = FALSE][1, ]
      app_direction <- ifelse(best$NES > 0, "enriched", ifelse(best$NES < 0, "depleted", "flat"))
      paper_direction <- tolower(first_nonempty(pathway_claims$reported_direction[[i]], "not specified"))
      direction_status <- if (paper_direction %in% c("enriched", "up", "higher", "positive") && best$NES > 0) {
        "direction agrees"
      } else if (paper_direction %in% c("depleted", "down", "lower", "negative") && best$NES < 0) {
        "direction agrees"
      } else if (paper_direction == "not specified") {
        "direction not specified in paper claim"
      } else {
        "direction mismatch"
      }
      rows[[length(rows) + 1]] <- data.frame(
        Finding = paste("Reported pathway:", term),
        Evidence_in_app = paste0(
          best$contrast,
          "; ", best$collection,
          "; NES=", signif(best$NES, 3),
          "; padj=", signif(best$padj, 3),
          "; paper direction=", paper_direction,
          "; app direction=", app_direction
        ),
        Interpretation = if (direction_status == "direction mismatch") {
          "Direction mismatch"
        } else if (!is.na(best$padj) && best$padj <= padj_cutoff) {
          paste("Supported by app GSEA result;", direction_status)
        } else {
          paste("Detected but not significant at selected threshold;", direction_status)
        },
        check.names = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      Finding = "Extracted paper findings",
      Evidence_in_app = "No reported biomarkers or pathways were automatically extracted. Add paper notes with specific genes/pathways, then train again.",
      Interpretation = "Needs richer paper extraction",
      check.names = FALSE
    ))
  }

  dplyr::bind_rows(rows)
}

classify_publication_comparison <- function(comparison) {
  if (is.null(comparison) || nrow(comparison) == 0 || !"Interpretation" %in% names(comparison)) {
    return(character())
  }
  text <- tolower(comparison$Interpretation)
  dplyr::case_when(
    grepl("supported by app", text, fixed = TRUE) ~ "Supported",
    grepl("direction mismatch", text, fixed = TRUE) ~ "Direction mismatch",
    grepl("missing from app", text, fixed = TRUE) ~ "Missing",
    grepl("not significant", text, fixed = TRUE) ~ "Detected but weak",
    grepl("needs richer", text, fixed = TRUE) ~ "Needs more paper extraction",
    grepl("not ready", text, fixed = TRUE) ~ "Not ready",
    TRUE ~ "Review"
  )
}

score_reproducibility_from_comparison <- function(comparison) {
  if (is.null(comparison) || nrow(comparison) == 0) {
    return(data.frame(
      Metric = c("Post-analysis reproducibility score", "Claims checked", "Supported claims", "Direction mismatches", "Missing claims", "Weak detections"),
      Value = c(0, 0, 0, 0, 0, 0),
      Status = c("Not ready", rep("No comparison rows", 5)),
      Rationale = c(
        "Run publication claim extraction and compare against the selected dataset.",
        rep("No paper-to-result comparison has been generated.", 5)
      ),
      check.names = FALSE
    ))
  }

  categories <- classify_publication_comparison(comparison)
  n <- length(categories)
  supported <- sum(categories == "Supported")
  mismatch <- sum(categories == "Direction mismatch")
  missing <- sum(categories == "Missing")
  weak <- sum(categories == "Detected but weak")
  not_ready <- sum(categories %in% c("Needs more paper extraction", "Not ready"))
  review <- sum(categories == "Review")

  raw <- if (n == 0) {
    0
  } else {
    (supported * 1 + weak * 0.45 + review * 0.25 - mismatch * 0.35 - missing * 0.2 - not_ready * 0.1) / n * 100
  }
  score <- clip_score(raw)
  status <- if (not_ready == n) {
    "Needs claims"
  } else if (score >= 75) {
    "Strong reproduction"
  } else if (score >= 55) {
    "Partial reproduction"
  } else if (score >= 35) {
    "Weak reproduction"
  } else {
    "Poor or inconclusive reproduction"
  }

  data.frame(
    Metric = c(
      "Post-analysis reproducibility score",
      "Claims checked",
      "Supported claims",
      "Direction mismatches",
      "Missing claims",
      "Weak detections"
    ),
    Value = c(score, n, supported, mismatch, missing, weak),
    Status = c(
      status,
      paste(score_band(score), "overall"),
      if (supported > 0) "Reproduced" else "None",
      if (mismatch > 0) "Review urgently" else "None",
      if (missing > 0) "Review data/annotation" else "None",
      if (weak > 0) "Below selected thresholds" else "None"
    ),
    Rationale = c(
      "Weighted score from supported, weak, missing, and direction-mismatched paper claims.",
      "Number of extracted publication claims checked against DEG/GSEA outputs.",
      "Claims with matching gene/pathway signal and selected significance/effect threshold.",
      "Claims where app direction conflicts with reported publication direction.",
      "Claims whose gene/pathway entity was not found in the current app result tables.",
      "Claims detected in results but not significant at selected thresholds."
    ),
    check.names = FALSE
  )
}

post_analysis_reproducibility_score <- function(dat, kb, padj_cutoff = 0.1, lfc_cutoff = 0.5) {
  comparison <- compare_findings_against_publication(dat, kb, padj_cutoff, lfc_cutoff)
  score_reproducibility_from_comparison(comparison)
}

claim_entity_set <- function(kb) {
  if (is.null(kb)) return(character())
  claims <- normalize_claim_table(kb$claims)
  extracted <- c(
    split_detected_terms(first_nonempty(kb$extracted$reported_biomarkers, "")),
    split_detected_terms(first_nonempty(kb$extracted$reported_pathways, ""))
  )
  unique(toupper(trimws(c(claims$entity, extracted))))
}

summarize_relevant_different_results <- function(dat, kb, comparison = NULL, padj_cutoff = 0.1, lfc_cutoff = 0.5, n = 10) {
  if (is.null(comparison)) {
    comparison <- compare_findings_against_publication(dat, kb, padj_cutoff, lfc_cutoff)
  }
  categories <- classify_publication_comparison(comparison)
  rows <- list()
  add_rows <- function(df, result_type) {
    if (is.null(df) || nrow(df) == 0) return()
    df$Result_type <- result_type
    rows[[length(rows) + 1]] <<- df
  }

  if (nrow(comparison) > 0 && length(categories) == nrow(comparison)) {
    comparison_subset <- function(mask, priority) {
      sub <- comparison[mask, , drop = FALSE]
      if (nrow(sub) == 0) {
        return(data.frame(Finding = character(), Evidence = character(), Interpretation = character(), Priority = character(), check.names = FALSE))
      }
      data.frame(
        Finding = sub$Finding,
        Evidence = sub$Evidence_in_app,
        Interpretation = sub$Interpretation,
        Priority = rep(priority, nrow(sub)),
        check.names = FALSE
      )
    }
    add_rows(comparison_subset(categories == "Supported", "High"), "Reproduced paper finding")
    add_rows(comparison_subset(categories == "Direction mismatch", "Highest"), "Different from paper")
    add_rows(comparison_subset(categories %in% c("Missing", "Detected but weak"), "Medium"), "Weak or missing paper finding")
  }

  claimed <- claim_entity_set(kb)
  deg <- dat$deg
  if (nrow(deg) > 0 && all(c("gene", "padj", "log2FoldChange") %in% names(deg))) {
    novel_deg <- deg |>
      dplyr::filter(!is.na(.data$padj), .data$padj <= padj_cutoff, abs(.data$log2FoldChange) >= lfc_cutoff) |>
      dplyr::filter(!toupper(.data$gene) %in% claimed) |>
      dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange))) |>
      dplyr::slice_head(n = n)
    if (nrow(novel_deg) > 0) {
      add_rows(
        data.frame(
          Finding = paste("High-signal DEG:", novel_deg$gene),
          Evidence = paste0(
            novel_deg$contrast,
            "; log2FC=", signif(novel_deg$log2FoldChange, 3),
            "; padj=", signif(novel_deg$padj, 3)
          ),
          Interpretation = "Strong app result not currently represented in extracted publication claims.",
          Priority = "Medium",
          check.names = FALSE
        ),
        "Potentially new DEG result"
      )
    }
  }

  gsea <- dat$gsea
  if (nrow(gsea) > 0 && all(c("pathway", "padj", "NES") %in% names(gsea))) {
    novel_gsea <- gsea |>
      dplyr::filter(!is.na(.data$padj), .data$padj <= padj_cutoff) |>
      dplyr::filter(!toupper(.data$pathway) %in% claimed) |>
      dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES))) |>
      dplyr::slice_head(n = n)
    if (nrow(novel_gsea) > 0) {
      add_rows(
        data.frame(
          Finding = paste("High-signal pathway:", novel_gsea$pathway),
          Evidence = paste0(
            novel_gsea$contrast,
            "; ", novel_gsea$collection,
            "; NES=", signif(novel_gsea$NES, 3),
            "; padj=", signif(novel_gsea$padj, 3)
          ),
          Interpretation = "Strong pathway result not currently represented in extracted publication claims.",
          Priority = "Medium",
          check.names = FALSE
        ),
        "Potentially new pathway result"
      )
    }
  }

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(data.frame(
      Result_type = "No triage rows",
      Priority = "Needs input",
      Finding = "No supported, different, weak, missing, or novel high-signal results could be summarized.",
      Evidence = "Run claim extraction and compare to paper, or inspect whether DEG/GSEA tables contain significant rows.",
      Interpretation = "Not ready",
      check.names = FALSE
    ))
  }
  out <- out[, c("Result_type", "Priority", "Finding", "Evidence", "Interpretation"), drop = FALSE]
  priority_order <- c("Highest", "High", "Medium", "Low", "Needs input")
  out$Priority <- factor(out$Priority, levels = priority_order, ordered = TRUE)
  out <- out[order(out$Priority, out$Result_type), , drop = FALSE]
  out$Priority <- as.character(out$Priority)
  rownames(out) <- NULL
  out
}

survival_candidate_columns <- function(metadata) {
  nms <- names(metadata)
  time_cols <- nms[grepl("survival|follow|time|days|month", nms, ignore.case = TRUE)]
  event_cols <- nms[grepl("event|death|progression|relapse|status", nms, ignore.case = TRUE)]
  list(time = time_cols, event = event_cols)
}

run_survival_summary <- function(dat, time_col = NULL, event_col = NULL, group_col = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(data.frame(Status = "Missing optional package survival."))
  }
  meta <- dat$metadata
  candidates <- survival_candidate_columns(meta)
  time_col <- first_nonempty(time_col, candidates$time[1])
  event_col <- first_nonempty(event_col, candidates$event[1])
  group_col <- first_nonempty(group_col, dat$info$default_group_col)
  if (!all(c(time_col, event_col, group_col) %in% names(meta))) {
    return(data.frame(Status = "No usable survival time, event, and group columns were detected."))
  }
  df <- meta[, c(time_col, event_col, group_col), drop = FALSE]
  names(df) <- c("time", "event", "group")
  df$time <- suppressWarnings(as.numeric(df$time))
  df$event <- suppressWarnings(as.integer(as.character(df$event)))
  df <- df[!is.na(df$time) & !is.na(df$event) & !is.na(df$group), , drop = FALSE]
  if (nrow(df) < 5 || length(unique(df$group)) < 2) {
    return(data.frame(Status = "Not enough complete survival records across at least two groups."))
  }
  fit <- survival::survdiff(survival::Surv(time, event) ~ group, data = df)
  pvalue <- stats::pchisq(fit$chisq, df = length(fit$n) - 1, lower.tail = FALSE)
  data.frame(
    Time_column = time_col,
    Event_column = event_col,
    Group_column = group_col,
    Complete_records = nrow(df),
    Groups = paste(names(fit$n), fit$n, sep = " n=", collapse = "; "),
    Logrank_pvalue = pvalue,
    Status = "Survival screen complete",
    check.names = FALSE
  )
}

run_immune_deconvolution_summary <- function(dat, method = "xcell") {
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    return(data.frame(Status = "Missing optional package immunedeconv. Install it to run immune cell deconvolution."))
  }
  expr <- dat$vsd
  expr_df <- data.frame(gene_symbol = rownames(expr), expr, check.names = FALSE)
  result <- tryCatch(
    immunedeconv::deconvolute(expr_df, method),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(data.frame(Status = paste("Immune deconvolution failed:", conditionMessage(result))))
  }
  data.frame(
    Status = "Immune deconvolution complete",
    Method = method,
    Cell_types = max(0, nrow(result) - 1),
    Samples = max(0, ncol(result) - 1),
    check.names = FALSE
  )
}

gap_analysis_report <- function(dat, kb = NULL) {
  available <- c("PCA", "DEG", "Volcano", "Heatmap", "GSEA", "Enrichr", "clusterProfiler", "Gene-level lookup")
  recommended <- c(
    "Immune cell deconvolution",
    "Survival analysis",
    "Kinase activity inference",
    "iLINCS perturbation analysis",
    "Network or regulon analysis",
    "Batch/sensitivity analysis",
    "Cross-cohort validation"
  )
  text <- paste(
    first_nonempty(if (!is.null(kb)) kb$extracted$statistical_methods else ""),
    first_nonempty(if (!is.null(kb)) kb$extracted$key_findings else ""),
    collapse = " "
  )
  missing <- recommended[!vapply(recommended, function(x) grepl(tolower(sub(" analysis$", "", x)), tolower(text), fixed = TRUE), logical(1))]
  data.frame(
    Category = c(rep("Already supported in app", length(available)), rep("Potential gap / next analysis", length(missing))),
    Analysis = c(available, missing),
    Rationale = c(
      rep("Available in the current RNA-seq Explorer workflow.", length(available)),
      rep("Not clearly detected in the trained publication knowledge; useful for hypothesis generation or reproducibility extension.", length(missing))
    ),
    check.names = FALSE
  )
}

standardized_pipeline_checklist <- function(dat = NULL, kb = NULL) {
  route <- if (is.null(kb)) "Not classified" else first_nonempty(geo_analyzability_route(kb)$Route, "Not classified")
  methods <- if (is.null(kb)) character() else split_detected_terms(first_nonempty(kb$extracted$statistical_methods, ""))
  claims <- if (is.null(kb)) empty_claim_table() else normalize_claim_table(kb$claims)
  has_dat <- !is.null(dat)
  has_deg <- has_dat && !is.null(dat$deg) && nrow(dat$deg) > 0
  has_gsea <- has_dat && !is.null(dat$gsea) && nrow(dat$gsea) > 0
  metadata_cols <- if (has_dat && !is.null(dat$metadata)) names(dat$metadata) else character()
  group_col <- if (has_dat) first_nonempty(dat$info$default_group_col, "") else ""
  rows <- data.frame(
    Stage = c(
      "1. Study intake",
      "2. Public-data route",
      "3. Metadata contract",
      "4. Original-method reproduction",
      "5. Result verification",
      "6. Signature/pathway layer",
      "7. Pitfall/provenance record",
      "8. Learning feedback"
    ),
    Required_standard = c(
      "GEO/PubMed context, linked/related papers, added notes, sample count, organism, and assay type are captured.",
      "Dataset is classified as processed matrix, count matrix, raw pipeline, supplementary review, or manual curation.",
      "Sample IDs, group column, contrast, patient/display columns, and exclusions are explicit.",
      "Primary paper-style contrast is run with the closest public-data method and deviations are named.",
      "Paper biomarkers/pathways/claims are compared with DEG/GSEA direction and significance.",
      "Common pathway/signature scores are checked using fixed gene sets and reported with coverage.",
      "Known pitfalls are flagged: batch, probe mapping, missing covariates, raw/reference mismatch, metadata ambiguity.",
      "Outcome/usefulness/reproducibility feedback is saved locally and can be shared anonymously with consent."
    ),
    Current_evidence = c(
      if (is.null(kb)) "No study memory selected." else paste("Study memory:", first_nonempty(kb$accession, "trained study")),
      route,
      if (!has_dat) "No active app dataset." else paste0("Dataset ", first_nonempty(dat$id, dat$info$dataset_id, "active"), "; group column ", group_col, "; ", length(metadata_cols), " metadata columns."),
      if (length(methods) == 0) "No method terms extracted yet." else paste("Detected methods:", paste(utils::head(methods, 6), collapse = ", ")),
      if (nrow(claims) == 0) "No structured claims yet." else paste(nrow(claims), "structured claims available."),
      if (has_gsea) paste(nrow(dat$gsea), "GSEA rows available.") else "No pathway/signature result layer yet.",
      if (is.null(kb)) "Pitfall scan needs study context." else paste("Route/paper context available for pitfall scan."),
      "Use Discovery -> Teach SeqSurf after each dataset."
    ),
    Status = c(
      if (is.null(kb)) "Not ready" else "Ready",
      if (is.null(kb)) "Not ready" else "Ready",
      if (!has_dat) "Needs dataset" else if (nzchar(group_col)) "Ready" else "Needs group mapping",
      if (!has_deg) "Needs reanalysis" else "Ready",
      if (nrow(claims) == 0) "Needs claim extraction" else if (has_deg || has_gsea) "Ready" else "Needs results",
      if (has_gsea) "Ready" else "Needs pathway/signature scoring",
      if (is.null(kb)) "Needs study context" else "Review",
      "Recommended"
    ),
    check.names = FALSE
  )
  rows
}

sequencing_pitfall_table <- function(dat = NULL, kb = NULL, data_root = "data") {
  route <- if (is.null(kb)) "Unknown" else first_nonempty(geo_analyzability_route(kb)$Route, "Unknown")
  text <- tolower(paste(
    if (is.null(kb)) "" else kb_evidence_blob(kb),
    if (is.null(kb) || is.null(kb$geo$supplementary_files)) "" else paste(kb$geo$supplementary_files, collapse = " "),
    collapse = " "
  ))
  meta_names <- if (!is.null(dat) && !is.null(dat$metadata)) names(dat$metadata) else character()
  learned <- learning_prior_for_study(kb, data_root)
  add <- function(pitfall, signal, impact, action, severity = "Review") {
    data.frame(Pitfall = pitfall, Signal = signal, Impact = impact, Action = action, Severity = severity, check.names = FALSE)
  }
  rows <- list(
    add(
      "Metadata/group ambiguity",
      if (length(meta_names) == 0) "No active metadata table." else paste(length(meta_names), "metadata columns available."),
      "Wrong group labels can invert DEG direction and make paper agreement meaningless.",
      "Confirm sample IDs, group levels, reference/comparison direction, exclusions, and duplicated patients before interpreting.",
      if (length(meta_names) == 0) "High" else "Review"
    ),
    add(
      "Original-method mismatch",
      route,
      "Public reanalysis may approximate the paper rather than exactly reproduce it.",
      "Record model formula, normalization, covariates, genome/annotation, and any missing paper-specific settings.",
      if (grepl("raw|count|manual|supplementary", route, ignore.case = TRUE)) "High" else "Review"
    ),
    add(
      "Batch/platform confounding",
      if (grepl("batch|platform|center|cohort|validation", text) || any(grepl("batch|platform|center|cohort|dataset", meta_names, ignore.case = TRUE))) "Detected" else "Not obvious",
      "Technical structure can masquerade as biology.",
      "Color PCA by technical columns and rerun/sensitivity-check models where metadata allow.",
      if (grepl("batch|platform|center|cohort|validation", text) || any(grepl("batch|platform|center|cohort|dataset", meta_names, ignore.case = TRUE))) "High" else "Review"
    ),
    add(
      "Probe/gene identifier mapping",
      if (grepl("microarray|array|affymetrix|agilent|probe", text)) "Array/probe signal detected" else "No array signal detected",
      "Probe collapsing can change apparent biomarker direction and pathway membership.",
      "Document platform annotation, duplicate-gene handling, unmapped probes, and symbol updates.",
      if (grepl("microarray|array|affymetrix|agilent|probe", text)) "High" else "Low"
    ),
    add(
      "Raw/reference reconstruction",
      if (grepl("fastq|sra|bam|sam|salmon|star|featurecounts", text)) "Raw/count workflow signal detected" else "No raw/count signal detected",
      "Genome build, GTF, quantifier, and QC choices affect every downstream conclusion.",
      "Use the external workflow folder, run preflight, and keep provenance logs before import.",
      if (grepl("fastq|sra|bam|sam|salmon|star|featurecounts", text)) "High" else "Low"
    ),
    add(
      "Learned local/shared risk",
      first_nonempty(learned$Evidence[[1]], "No learning prior yet."),
      "Past datasets with similar route/assay/tags can warn about recurring pitfalls.",
      first_nonempty(learned$Recommendation_bias[[1]], "Collect feedback after this dataset."),
      if (grepl("caution|failed|manual", first_nonempty(learned$Recommendation_bias[[1]], ""), ignore.case = TRUE)) "High" else "Review"
    )
  )
  dplyr::bind_rows(rows)
}

signature_library <- function() {
  list(
    epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN3", "CLDN4", "MUC1", "TACSTD2", "DSP"),
    mesenchymal = c("VIM", "FN1", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "COL1A1", "ITGA5")
  )
}

score_gene_signature <- function(expr, genes) {
  present <- intersect(toupper(genes), toupper(rownames(expr)))
  row_lookup <- match(present, toupper(rownames(expr)))
  if (length(row_lookup) < 3) return(rep(NA_real_, ncol(expr)))
  gene_z <- t(scale(t(expr[row_lookup, , drop = FALSE])))
  gene_z[!is.finite(gene_z)] <- 0
  colMeans(gene_z, na.rm = TRUE)
}

signature_screen_table <- function(dat = NULL, kb = NULL) {
  if (is.null(dat) || is.null(dat$vsd) || nrow(dat$vsd) == 0) {
    return(data.frame(
      Signature = "Epithelial/EMT",
      Coverage = "No active expression matrix.",
      Direction = "Not ready",
      Evidence = "Validate or select an app dataset first.",
      Recommendation = "Run reanalysis, then use fixed signature scores for standardized cross-dataset interpretation.",
      check.names = FALSE
    ))
  }
  lib <- signature_library()
  expr <- dat$vsd
  epi_present <- intersect(toupper(lib$epithelial), toupper(rownames(expr)))
  mes_present <- intersect(toupper(lib$mesenchymal), toupper(rownames(expr)))
  epi <- score_gene_signature(expr, lib$epithelial)
  mes <- score_gene_signature(expr, lib$mesenchymal)
  emt_delta <- epi - mes
  meta <- dat$metadata
  group_col <- first_nonempty(dat$info$default_group_col, "")
  group_signal <- "Group comparison unavailable."
  direction <- "Coverage only"
  if (nzchar(group_col) && group_col %in% names(meta) && all(names(emt_delta) %in% rownames(meta))) {
    group <- as.character(meta[names(emt_delta), group_col])
    levels <- names(sort(table(group, useNA = "no"), decreasing = TRUE))
    if (length(levels) >= 2) {
      positive <- first_nonempty(dat$info$positive_group, levels[[2]])
      negative <- first_nonempty(dat$info$negative_group, levels[[1]])
      if (!positive %in% levels) positive <- levels[[2]]
      if (!negative %in% levels) negative <- levels[[1]]
      pos_mean <- mean(emt_delta[group == positive], na.rm = TRUE)
      neg_mean <- mean(emt_delta[group == negative], na.rm = TRUE)
      diff <- pos_mean - neg_mean
      direction <- if (is.na(diff)) "Not enough complete samples" else if (diff > 0.15) "More epithelial in positive group" else if (diff < -0.15) "More mesenchymal in positive group" else "No strong epithelial/mesenchymal shift"
      group_signal <- paste0(positive, " mean delta ", signif(pos_mean, 3), "; ", negative, " mean delta ", signif(neg_mean, 3), "; difference ", signif(diff, 3))
    }
  }
  context <- tolower(if (is.null(kb)) "" else kb_evidence_blob(kb))
  data.frame(
    Signature = c("Epithelial genes", "Mesenchymal genes", "Epithelial-minus-mesenchymal delta"),
    Coverage = c(
      paste(length(epi_present), "/", length(lib$epithelial), "genes present:", paste(utils::head(epi_present, 8), collapse = ", ")),
      paste(length(mes_present), "/", length(lib$mesenchymal), "genes present:", paste(utils::head(mes_present, 8), collapse = ", ")),
      paste("Computed across", ncol(expr), "samples.")
    ),
    Direction = c("Fixed positive epithelial signature", "Fixed mesenchymal/EMT comparator", direction),
    Evidence = c(
      "Standardized gene set used consistently across datasets.",
      "Standardized gene set used consistently across datasets.",
      group_signal
    ),
    Recommendation = c(
      "Report coverage before interpreting epithelial scores.",
      "Report coverage before interpreting mesenchymal/EMT scores.",
      if (grepl("epithelial|mesenchymal|emt", context)) "This matches study context; include as a pre-specified follow-up check." else "Use as a hypothesis-generating standardized screen."
    ),
    check.names = FALSE
  )
}

suggest_unexplored_comparisons <- function(dat, kb = NULL, max_levels = 6) {
  metadata <- dat$metadata
  existing_contrasts <- unique(dat$deg$contrast)
  existing_text <- tolower(paste(existing_contrasts, collapse = " | "))
  preferred_patterns <- c(
    "mycn", "stage", "risk", "response", "progress", "death", "survival",
    "sex", "age", "class", "treatment", "timepoint", "dataset"
  )
  cols <- names(metadata)[vapply(metadata, function(x) {
    values <- unique(stats::na.omit(as.character(x)))
    length(values) >= 2 && length(values) <= max_levels
  }, logical(1))]
  cols <- cols[!cols %in% c(dat$info$sample_id_col, dat$info$display_sample_col, dat$info$patient_id_col)]
  if (length(cols) == 0) {
    return(data.frame(
      Priority = 1,
      Suggested_comparison = "No categorical metadata columns found",
      Group_column = "",
      Rationale = "The selected dataset does not expose obvious categorical sample groups for additional DEG contrasts.",
      Status = "Needs metadata curation",
      check.names = FALSE
    ))
  }

  score_column <- function(col) {
    score <- 0
    if (col == dat$info$default_group_col) score <- score + 3
    if (any(grepl(paste(preferred_patterns, collapse = "|"), col, ignore.case = TRUE))) score <- score + 2
    if (grepl(col, existing_text, fixed = TRUE)) score <- score - 3
    score
  }

  ranked_cols <- cols[order(vapply(cols, score_column, numeric(1)), decreasing = TRUE)]
  rows <- lapply(ranked_cols, function(col) {
    counts <- sort(table(as.character(metadata[[col]]), useNA = "no"), decreasing = TRUE)
    levels <- names(counts)
    if (length(levels) < 2) return(NULL)
    negative <- levels[[1]]
    positive <- levels[[2]]
    label <- paste(positive, "vs", negative)
    already_done <- grepl(tolower(positive), existing_text, fixed = TRUE) &&
      grepl(tolower(negative), existing_text, fixed = TRUE)
    data.frame(
      Suggested_comparison = label,
      Group_column = col,
      Positive_group = positive,
      Negative_group = negative,
      Group_sizes = paste(paste(levels, as.integer(counts), sep = " n="), collapse = "; "),
      Rationale = if (already_done) {
        "This comparison appears related to an existing DEG contrast; use as a sanity check or refine groups."
      } else if (grepl("survival|death|progress", col, ignore.case = TRUE)) {
        "Outcome-associated grouping may reveal biology linked to prognosis or progression."
      } else if (grepl("stage|risk|mycn|response", col, ignore.case = TRUE)) {
        "Clinically relevant grouping likely to be interpretable and publication-relevant."
      } else {
        "Metadata-defined group not clearly represented in current DEG contrasts."
      },
      Status = if (already_done) "Possibly already covered" else "Candidate new comparison",
      check.names = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(data.frame(
      Priority = 1,
      Suggested_comparison = "No usable comparisons found",
      Group_column = "",
      Rationale = "Categorical columns were present but did not contain at least two non-empty levels.",
      Status = "Needs metadata curation",
      check.names = FALSE
    ))
  }
  out$Priority <- seq_len(nrow(out))
  out[, c("Priority", "Suggested_comparison", "Group_column", "Positive_group", "Negative_group", "Group_sizes", "Rationale", "Status"), drop = FALSE]
}

suggest_ai_analyses <- function(dat, kb = NULL) {
  deg <- dat$deg
  gsea <- dat$gsea
  top_deg <- deg |>
    dplyr::filter(!is.na(.data$padj)) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange))) |>
    dplyr::slice_head(n = 8) |>
    dplyr::pull(.data$gene)
  top_pathways <- gsea |>
    dplyr::filter(!is.na(.data$padj)) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES))) |>
    dplyr::slice_head(n = 6) |>
    dplyr::pull(.data$pathway)

  literature_terms <- if (is.null(kb)) {
    character()
  } else {
    unique(unlist(strsplit(first_nonempty(kb$extracted$reported_pathways, ""), ";\\s*")))
  }
  pathway_hits <- top_pathways[tolower(top_pathways) %in% tolower(literature_terms)]

  suggestions <- c(
    paste("Start by validating the default contrast:", dat$info$default_contrast),
    paste("Inspect top differential genes:", paste(top_deg, collapse = ", ")),
    paste("Prioritize enriched pathways:", paste(top_pathways, collapse = " | ")),
    paste("Check whether PCA separates by", dat$info$default_group_col, "or by technical/sample variables."),
    "Use Gene Search for paper biomarkers and compare DEG direction, expression distribution, and leading-edge pathway membership."
  )
  if (length(pathway_hits) > 0) {
    suggestions <- c(suggestions, paste("Paper-to-data pathway overlap detected:", paste(pathway_hits, collapse = ", ")))
  }

  data.frame(Priority = seq_along(suggestions), Suggestion = suggestions, check.names = FALSE)
}

ranked_discovery_plan <- function(dat = NULL, kb = NULL, comparison = NULL, n = 6) {
  if (is.null(dat)) {
    return(data.frame(
      Rank = 1,
      Analysis = "Validate a dataset first",
      Why = "Discovery needs an active app dataset, DEG results, GSEA results, and trained study memory to make specific recommendations.",
      Where_to_open = "Validate",
      Evidence = "No active dataset selected.",
      Trust_check = "Run Validate or select saved results before using Discovery.",
      Readiness = "Not ready",
      check.names = FALSE
    ))
  }

  add <- function(rows, analysis, why, where, evidence, trust, readiness = "Ready") {
    rows[[length(rows) + 1]] <- data.frame(
      Analysis = analysis,
      Why = why,
      Where_to_open = where,
      Evidence = evidence,
      Trust_check = trust,
      Readiness = readiness,
      check.names = FALSE
    )
    rows
  }

  rows <- list()
  deg <- dat$deg
  gsea <- dat$gsea
  meta <- dat$metadata
  default_contrast <- first_nonempty(dat$info$default_contrast, "default contrast")
  top_deg <- deg |>
    dplyr::filter(!is.na(.data$padj)) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange))) |>
    dplyr::slice_head(n = 5)
  top_pathway <- gsea |>
    dplyr::filter(!is.na(.data$padj)) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES))) |>
    dplyr::slice_head(n = 3)
  top_gene_text <- if (nrow(top_deg) > 0) paste(top_deg$gene, collapse = ", ") else "No significant DEG rows available."
  top_pathway_text <- if (nrow(top_pathway) > 0) paste(top_pathway$pathway, collapse = " | ") else "No significant GSEA rows available."

  rows <- add(
    rows,
    paste("Audit the primary contrast:", default_contrast),
    "This is the core paper-style reproduction check: do the strongest DEG and pathway signals match the trained study claims?",
    "Validate -> Score, DEG Explorer, Volcano Plot",
    paste("Top DEG signals:", top_gene_text),
    "Confirm group labels, sample counts, adjusted p-value cutoff, and log2FC direction before presenting."
  )

  if (!is.null(comparison) && nrow(comparison) > 0 && "Interpretation" %in% names(comparison)) {
    different <- comparison[grepl("mismatch|not significant|missing|weak", comparison$Interpretation, ignore.case = TRUE), , drop = FALSE]
    if (nrow(different) > 0) {
      rows <- add(
        rows,
        "Investigate paper-disagreement claims",
        "Disagreement is the most interesting scientific output: it tells you where public reanalysis does not cleanly reproduce the original paper.",
        "Validate -> Paper Agreement, Validate -> Score, Gene Search",
        paste(utils::head(different$Finding, 3), collapse = " | "),
        "Check whether disagreement is caused by sample mismatch, probe/gene mapping, platform normalization, or a truly different result.",
        "High value"
      )
    }
  }

  rows <- add(
    rows,
    "Inspect pathway-level biology",
    "Pathways stabilize gene-level noise and are easier to connect to mechanism, especially when individual biomarkers are platform-dependent.",
    "GSEA Explorer, Heatmap, Gene Search",
    paste("Top pathways:", top_pathway_text),
    "Verify pathway direction, leading-edge genes, and whether the pathway was reported in the original/related papers."
  )

  comparisons <- tryCatch(suggest_unexplored_comparisons(dat, kb), error = function(e) data.frame())
  if (nrow(comparisons) > 0 && !"No usable comparisons found" %in% comparisons$Suggested_comparison) {
    rows <- add(
      rows,
      "Run a metadata-derived comparison",
      "A new comparison can turn reanalysis from reproduction into discovery if the metadata contains a group not tested in the paper.",
      "Validate design controls, DEG Explorer, Volcano Plot",
      first_nonempty(comparisons$Suggested_comparison[[1]], "Candidate comparison from metadata."),
      "Check group sizes and whether the comparison is biologically meaningful, not just a metadata artifact.",
      "Candidate"
    )
  }

  pca_note <- if (!is.null(dat$pca_variance) && length(dat$pca_variance) >= 2) {
    paste0("PC1 ", round(dat$pca_variance[[1]], 1), "%; PC2 ", round(dat$pca_variance[[2]], 1), "% variance.")
  } else {
    "PCA variance not available."
  }
  rows <- add(
    rows,
    "Check PCA and possible confounding",
    "If samples separate by batch, tissue source, platform, sex, or center, DEG results may reflect technical structure rather than biology.",
    "PCA Explorer",
    pca_note,
    "Color PCA by every major metadata column before trusting the default contrast."
  )

  if (any(grepl("time|survival|event|days|months|status", names(meta), ignore.case = TRUE))) {
    rows <- add(
      rows,
      "Run outcome/survival screening",
      "Outcome-aware analysis can make reanalysis more clinically meaningful when time/event metadata exist.",
      "Discovery -> Survival Screen",
      "Time/event-like metadata columns detected.",
      "Confirm censoring definitions and avoid overclaiming small or incomplete clinical metadata.",
      "Conditional"
    )
  }

  out <- dplyr::bind_rows(rows)
  out$Rank <- seq_len(nrow(out))
  utils::head(out[, c("Rank", "Analysis", "Why", "Where_to_open", "Evidence", "Trust_check", "Readiness"), drop = FALSE], n)
}

adaptive_ranked_discovery_plan <- function(dat = NULL, kb = NULL, comparison = NULL, data_root = "data", n = 6) {
  base <- ranked_discovery_plan(dat, kb, comparison, n = n)
  prior <- learning_prior_for_study(kb, data_root)
  bias <- first_nonempty(prior$Recommendation_bias[[1]], "")
  learned <- NULL
  if (grepl("Add caution", bias, fixed = TRUE)) {
    learned <- data.frame(
      Rank = 1,
      Analysis = "Start with a learned-risk audit",
      Why = "Local feedback says similar datasets often needed manual curation or failed before producing useful reanalysis.",
      Where_to_open = "Start -> Analyzability Details, Validate -> Metadata",
      Evidence = prior$Evidence[[1]],
      Trust_check = "Confirm file type, sample groups, sample count, gene IDs, and paper-method match before interpreting DEG/GSEA.",
      Readiness = "Learned caution",
      check.names = FALSE
    )
  } else if (grepl("Upgrade confidence", bias, fixed = TRUE)) {
    learned <- data.frame(
      Rank = 1,
      Analysis = "Use the standard reproduction route first",
      Why = "Local feedback says similar datasets tended to be useful, so the highest-value move is a clean paper-style reproduction before adding novelty.",
      Where_to_open = "Validate -> Run reanalysis, Validate -> Score",
      Evidence = prior$Evidence[[1]],
      Trust_check = "Still verify group labels and paper claims before treating agreement as biological evidence.",
      Readiness = "Learned boost",
      check.names = FALSE
    )
  } else if (grepl("discrepancy", bias, ignore.case = TRUE)) {
    learned <- data.frame(
      Rank = 1,
      Analysis = "Make disagreement a first-class result",
      Why = "Similar local datasets had weak reproduction, so mismatches may be the most informative output.",
      Where_to_open = "Validate -> Paper Agreement, Gene Search, GSEA Explorer",
      Evidence = prior$Evidence[[1]],
      Trust_check = "Separate true biological disagreement from sample mismatch, probe mapping, normalization, and missing covariates.",
      Readiness = "Learned priority",
      check.names = FALSE
    )
  }
  tags <- study_learning_features(kb)$learned_tags
  if ("epithelial_signature" %in% tags) {
    epithelial <- data.frame(
      Rank = 1,
      Analysis = "Score epithelial/EMT signature shifts",
      Why = "This study context matches an epithelial-signature pattern; signature-level analysis can standardize interpretation across sequencing studies.",
      Where_to_open = "GSEA Explorer, Heatmap, Gene Search",
      Evidence = "Epithelial/mesenchymal terms detected in study context or feedback.",
      Trust_check = "Use a fixed gene set and report direction, leading-edge genes, and whether the signature was pre-specified.",
      Readiness = "Adaptive hypothesis",
      check.names = FALSE
    )
    learned <- dplyr::bind_rows(learned, epithelial)
  }
  if (is.null(learned) || nrow(learned) == 0) {
    return(base)
  }
  out <- dplyr::bind_rows(learned, base)
  out$Rank <- seq_len(nrow(out))
  utils::head(out[, c("Rank", "Analysis", "Why", "Where_to_open", "Evidence", "Trust_check", "Readiness"), drop = FALSE], n)
}

paper_to_code_brief <- function(dat = NULL, kb = NULL) {
  info <- if (is.null(dat) || is.null(dat$info)) list() else dat$info
  kb_table <- study_kb_summary_table(kb)
  kb_lines <- if (nrow(kb_table) == 0) {
    "- No trained study knowledge base selected yet."
  } else {
    paste0("- ", kb_table$Field, ": ", kb_table$Value, collapse = "\n")
  }
  suggestions <- if (is.null(dat)) {
    reanalysis_recommendation_table(kb)
  } else {
    suggest_ai_analyses(dat, kb)
  }
  suggestion_col <- if ("Suggestion" %in% names(suggestions)) "Suggestion" else if ("Recommendation" %in% names(suggestions)) "Recommendation" else names(suggestions)[[1]]
  suggestion_lines <- paste0("- ", suggestions[[suggestion_col]], collapse = "\n")
  route <- geo_analyzability_route(kb)
  route_lines <- paste0(
    "- Route: ", route$Route,
    "\n- Can run in SeqSurf app: ", route$Can_run_in_app,
    "\n- Confidence: ", route$Confidence,
    "\n- Evidence: ", route$Evidence,
    "\n- Blocking reason: ", route$Blocking_reason,
    "\n- Recommended next step: ", route$Recommended_next_step
  )
  sections <- summarize_document_sections(kb)
  section_lines <- paste0("- ", sections$Section, " [", sections$Status, "]: ", sections$Excerpt, collapse = "\n")
  feasibility <- exact_reproduction_feasibility_table(kb)
  feasibility_lines <- paste0("- ", feasibility$Check, " [", feasibility$Status, "]: ", feasibility$Assessment, collapse = "\n")
  raw_triage <- raw_supplementary_triage_table(kb)
  raw_lines <- paste0("- ", raw_triage$Resource, " [", raw_triage$Status, "]: ", raw_triage$Action, collapse = "\n")
  claim_candidates <- if (!is.null(kb) && !is.null(kb$claim_candidates)) kb$claim_candidates else section_claim_candidate_table(paste(sections$Excerpt, collapse = " "))
  claim_candidate_lines <- paste0("- [", claim_candidates$Evidence_type, "] ", claim_candidates$Claim_candidate, collapse = "\n")

  paste(
    "# Paper-to-code reproducible reanalysis brief",
    "",
    "Use this brief with Claude Code, Codex, or another coding agent to convert the paper-grounded study plan into reproducible app-ready analysis code.",
    "",
    "## Active app dataset",
    if (is.null(dat)) "- No active app dataset selected yet. This is a pre-validation paper-to-code brief from the trained study memory." else paste0("- Dataset ID: ", first_nonempty(info$dataset_id, dat$id)),
    paste0("- Display name: ", first_nonempty(info$display_name)),
    paste0("- Sample ID column: ", first_nonempty(info$sample_id_col)),
    paste0("- Patient ID column: ", first_nonempty(info$patient_id_col)),
    paste0("- Default group column: ", first_nonempty(info$default_group_col)),
    paste0("- Default contrast: ", first_nonempty(info$default_contrast)),
    paste0("- Expression unit: ", first_nonempty(info$expression_unit, "Not specified")),
    "",
    "## Trained paper knowledge",
    kb_lines,
    "",
    "## Section-aware paper extraction",
    section_lines,
    "",
    "## Reanalysis route decision",
    route_lines,
    "",
    "## Exact-method feasibility",
    feasibility_lines,
    "",
    "## Raw/supplementary file triage",
    raw_lines,
    "",
    "## Claim candidates from paper sections",
    claim_candidate_lines,
    "",
    "## Count-matrix automation plan",
    count_matrix_pipeline_plan(kb),
    "",
    "## FASTQ/SRA external handoff",
    raw_sra_handoff_plan(kb),
    "",
    "## Coding-agent instructions",
    "- If route is processed-matrix reanalysis, write an R script that uses GEOquery/Biobase, maps probes/features to gene symbols where possible, builds the app dataset contract, and records every transformation.",
    "- If route is count-matrix pipeline, write an R script using DESeq2, edgeR, or limma-voom as justified by the paper methods and count matrix format; do not pretend this is exact if paper covariates/files are missing.",
    "- If route is raw sequencing pipeline, write a Nextflow/Snakemake-style plan or shell/R handoff with required tools, reference genome, annotation version, QC, quantification, count summarization, and app-contract export.",
    "- Always include a provenance table listing input URLs/files, sample exclusions, gene/probe mapping decisions, normalization, model formula, contrast, and known deviations from the paper.",
    "",
    "## Reanalysis tasks to code",
    "- Retrieve GEO metadata, sample annotations, expression/count matrices, platform annotations, and supplementary files where available.",
    "- Recreate the paper's primary sample groups and document any samples excluded from the app dataset.",
    "- Reproduce the paper's reported statistical model as closely as public data permits.",
    "- Generate app contract files: metadata.rds, vsd_matrix.rds, pca_df.rds, pca_variance.rds, pca_loadings.rds, deg_results.csv, gsea_hallmark.csv, pathway_scores.csv, dataset_info.rds.",
    "- Add provenance comments for every transformation that changes sample IDs, gene IDs, normalization, filtering, or group labels.",
    "- Add tests that compare sample counts, group counts, top biomarkers, reported pathways, and major paper claims against the generated app objects.",
    "",
    "## Suggested analysis checks",
    suggestion_lines,
    "",
    "## Compare-to-publication report requirements",
    "- Parse reported genes, biomarkers, pathways, and statistical claims from the original paper and supplements.",
    "- Compare reported genes against deg_results.csv by gene, contrast, effect direction, log2 fold change, p-value, and adjusted p-value.",
    "- Compare reported pathways against gsea_hallmark.csv by pathway name, collection, NES direction, and adjusted p-value.",
    "- Produce one table with Supported, Direction mismatch, Not significant, Missing from app result, and Not publicly reproducible statuses.",
    "",
    "## Acceptance criteria",
    "- The generated dataset folder passes the app data contract tests.",
    "- The Start, Train, Validate, and Discovery tabs explain mismatches between GEO, the paper evidence, and the app-ready dataset.",
    "- Any non-reproducible paper claim is flagged as unavailable, ambiguous, or not supported by public GEO files.",
    sep = "\n"
  )
}
