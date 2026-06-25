test_that("GEO accession normalization accepts common series ids", {
  expect_equal(normalize_geo_accession(" gse49711 "), "GSE49711")
  expect_true(is_geo_accession("GSE49711"))
  expect_false(is_geo_accession("not a geo id"))
})

test_that("assistant knowledge base tables and verification are deterministic", {
  dat <- load_test_dataset("synthetic")
  kb <- list(
    accession = "GSE123",
    geo = list(
      title = "Synthetic cohort",
      pubmed_id = "00000000",
      sample_count = nrow(dat$metadata),
      platform_accessions = "GPLTEST"
    ),
    pubmed = list(title = "Synthetic RNA-seq study"),
    extracted = list(
      study_objective = "Evaluate treatment expression.",
      cohort_characteristics = "Synthetic samples with RNA-seq counts and outcome response metadata.",
      sample_groups = "Treated and control groups with validation cohort structure.",
      statistical_methods = "limma; GSEA; PCA; batch correction",
      reported_pathways = "MYCN; cell cycle; apoptosis",
      reported_biomarkers = "MYCN; ALK; NTRK1",
      study_limitations = "Uses processed public expression"
    ),
    claims = data.frame(
      claim_type = c("biomarker", "biomarker", "pathway"),
      entity = c("MYCN", "NTRK1", "apoptosis"),
      comparison = "treated vs control",
      reported_direction = c("up", "down", "depleted"),
      reported_significance = "significant",
      evidence_sentence = "Synthetic claims.",
      source = "test",
      confidence = "high",
      stringsAsFactors = FALSE
    ),
    related_publications = data.frame(
      PubMed_ID = "11111111",
      Title = "Synthetic GEO reuse paper",
      Journal = "Synthetic Journal",
      Publication_date = "2026",
      Abstract = "This paper reuses GSE123 and reports a stronger apoptosis signal.",
      Source = "https://pubmed.ncbi.nlm.nih.gov/11111111/",
      Status = "PubMed paper mentioning GEO accession",
      check.names = FALSE
    ),
    extra_papers = "Additional methods mention survival, progression, response, batch, and RNA-seq reanalysis."
  )

  summary <- study_kb_summary_table(kb)
  assessment <- reanalysis_assessment_table(kb)
  reanalysis_targets <- reanalysis_recommendation_table(kb)
  narrative_prompt <- build_reanalysis_narrative_prompt(kb)
  section_table <- detect_document_sections("Methods\nWe used limma and GSEA.\nResults\nMYCN increased.\nDiscussion\nThe cohort is limited.")
  verification <- verify_dataset_against_kb(dat, kb)
  suggestions <- suggest_ai_analyses(dat, kb)
  comparison_suggestions <- suggest_unexplored_comparisons(dat, kb)
  comparison <- compare_findings_against_publication(dat, kb)
  gaps <- gap_analysis_report(dat, kb)
  post_score <- score_reproducibility_from_comparison(comparison)
  triage <- summarize_relevant_different_results(dat, kb, comparison)
  checklist <- standardized_pipeline_checklist(dat, kb)
  pitfalls <- sequencing_pitfall_table(dat, kb, tempdir())
  signatures <- signature_screen_table(dat, kb)
  prior_work <- prior_work_report_table(kb)
  feasibility <- exact_reproduction_feasibility_table(kb)
  raw_triage <- raw_supplementary_triage_table(kb)
  evidence_blob <- kb_evidence_blob(kb)
  kb$virtual_lab <- list(
    literature_curator = list(
      agent_id = "literature_curator",
      role = "Literature Curator",
      summary = "Prior work reports MYCN and apoptosis.",
      evidence = c("Synthetic evidence sentence."),
      recommendations = c("Check original MYCN claim."),
      uncertainties = c("Full text not available.")
    ),
    methods_extractor = list(
      agent_id = "methods_extractor",
      role = "Methods Extractor",
      summary = "limma and GSEA are detected.",
      evidence = c("Methods mention limma."),
      recommendations = c("Reproduce limma-style contrast."),
      uncertainties = c("Covariates unclear.")
    )
  )
  virtual_lab <- virtual_lab_summary_table(kb)
  virtual_lab_recs <- virtual_lab_recommendation_table(kb)
  citation_evidence <- citation_evidence_table(kb)

  expect_true(all(c("Field", "Value") %in% names(summary)))
  expect_true(any(summary$Field == "Related papers found" & summary$Value == "1"))
  expect_true(grepl("Synthetic GEO reuse paper", evidence_blob, fixed = TRUE))
  expect_true(all(c("Score", "Value", "Band", "Rationale") %in% names(assessment)))
  expect_setequal(
    assessment$Score,
    c("Data availability", "Original-method reproducibility", "Reanalysis usefulness", "Novelty opportunity", "Technical risk", "Proceed recommendation")
  )
  expect_true(all(assessment$Value >= 0 & assessment$Value <= 100))
  expect_true(assessment$Value[assessment$Score == "Proceed recommendation"] >= 50)
  expect_true(all(c("Priority", "Recommendation", "Why", "Status") %in% names(reanalysis_targets)))
  expect_true(any(grepl("biomarker", reanalysis_targets$Recommendation, ignore.case = TRUE)))
  expect_true(any(grepl("outcome", reanalysis_targets$Recommendation, ignore.case = TRUE)))
  expect_true(grepl("What has already been done", narrative_prompt, fixed = TRUE))
  expect_true(grepl("Proceed recommendation", narrative_prompt, fixed = TRUE))
  expect_true(grepl("MYCN", narrative_prompt, fixed = TRUE))
  expect_true(any(section_table$Section == "Methods" & section_table$Status == "Detected"))
  expect_true(any(grepl("limma", section_table$Excerpt, fixed = TRUE)))
  rich_sections <- detect_document_sections("1 Methods\nWe used limma.\nData availability\nGSE files are public.\nFigure legends\nFigure 1 shows MYCN.")
  expect_true(any(rich_sections$Section == "Data availability" & rich_sections$Status == "Detected"))
  expect_true(any(rich_sections$Section == "Tables/Figures" & rich_sections$Status == "Detected"))
  candidates <- section_claim_candidate_table("MYCN was significantly increased with adjusted p < 0.05. Survival was associated with pathway scores.")
  expect_true(any(candidates$Status == "Candidate"))
  expect_true(any(verification$Status == "Pass"))
  expect_true(nrow(suggestions) >= 5)
  expect_true(nrow(comparison_suggestions) >= 1)
  expect_true(all(c("Suggested_comparison", "Group_column", "Status") %in% names(comparison_suggestions)))
  expect_true(any(grepl(dat$info$default_contrast, suggestions$Suggestion, fixed = TRUE)))
  expect_true(all(c("Finding", "Evidence_in_app", "Interpretation") %in% names(comparison)))
  expect_true(any(grepl("MYCN", comparison$Finding, fixed = TRUE)))
  expect_true(all(c("Metric", "Value", "Status", "Rationale") %in% names(post_score)))
  expect_true(post_score$Value[post_score$Metric == "Post-analysis reproducibility score"] >= 40)
  expect_true(all(c("Result_type", "Priority", "Finding", "Evidence", "Interpretation") %in% names(triage)))
  expect_true(any(triage$Result_type == "Reproduced paper finding"))
  expect_true(any(triage$Result_type == "Potentially new DEG result"))
  expect_true(all(c("Stage", "Required_standard", "Current_evidence", "Status") %in% names(checklist)))
  expect_true(any(grepl("Result verification", checklist$Stage)))
  expect_true(all(c("Pitfall", "Signal", "Impact", "Action", "Severity") %in% names(pitfalls)))
  expect_true(any(grepl("Metadata", pitfalls$Pitfall)))
  expect_true(all(c("Signature", "Coverage", "Direction", "Evidence", "Recommendation") %in% names(signatures)))
  expect_true(any(grepl("Epithelial", signatures$Signature)))
  expect_true(all(c("Section", "Finding", "Evidence", "Status") %in% names(prior_work)))
  expect_true(any(prior_work$Section == "Original study"))
  expect_true(all(c("Check", "Assessment", "Status", "Rationale") %in% names(feasibility)))
  expect_true(any(grepl("reproduction", feasibility$Check, ignore.case = TRUE) | grepl("Recommended mode", feasibility$Check)))
  expect_true(all(c("Resource", "Signal", "Action", "Status") %in% names(raw_triage)))
  expect_true(any(raw_triage$Resource == "Raw sequencing"))
  expect_true(grepl("Count-matrix automation plan", count_matrix_pipeline_plan(kb), fixed = TRUE))
  expect_true(grepl("FASTQ/SRA external workflow handoff", raw_sra_handoff_plan(kb), fixed = TRUE))
  raw_files <- raw_sra_workflow_files(kb)
  expect_true(all(c(
    "README.md",
    "AI_HANDOFF_GUIDE.md",
    "PROMPT_FOR_AI_ASSISTANT.md",
    "START_HERE.command",
    "scripts/01_resolve_sra_runs.R",
    "scripts/00_preflight_check.R",
    "scripts/02_download_fastq.sh",
    "scripts/03_quant_salmon.sh",
    "scripts/04_build_gene_counts.R",
    "scripts/05_prepare_seqsurf_import.R"
  ) %in% names(raw_files)))
  workspace_root <- tempfile("seqsurf-workspace-")
  workspace <- write_external_workflow_workspace(kb, workspace_root)
  expect_true(file.exists(file.path(workspace$directory, "PROMPT_FOR_AI_ASSISTANT.md")))
  expect_true(file.exists(file.path(workspace$directory, "START_HERE.command")))
  bundle_path <- file.path(tempdir(), "seqsurf_raw_sra_bundle_test.zip")
  write_raw_sra_workflow_bundle(kb, bundle_path)
  expect_true(file.exists(bundle_path))
  no_code_guide <- no_code_ai_handoff_guide(kb)
  ai_prompt <- external_ai_analysis_prompt(kb)
  expect_true(grepl("https://cran.r-project.org/", no_code_guide, fixed = TRUE))
  expect_true(grepl("outputs/seqsurf_gene_counts.csv", no_code_guide, fixed = TRUE))
  expect_true(grepl("non-coder", ai_prompt, fixed = TRUE))
  expect_true(grepl("outputs/seqsurf_metadata.csv", ai_prompt, fixed = TRUE))
  expect_true(all(c("Category", "Analysis", "Rationale") %in% names(gaps)))
  expect_true(any(gaps$Category == "Potential gap / next analysis"))
  expect_true(all(c("Agent", "Role", "Summary", "Evidence", "Uncertainty") %in% names(virtual_lab)))
  expect_true(any(virtual_lab$Role == "Literature Curator"))
  expect_true(all(c("Role", "Recommendation", "Status") %in% names(virtual_lab_recs)))
  expect_true(any(grepl("MYCN", virtual_lab_recs$Recommendation, fixed = TRUE)))
  expect_true(all(c("Source", "Evidence", "Used_for", "Confidence") %in% names(citation_evidence)))
  expect_true(any(grepl("MYCN", citation_evidence$Evidence, fixed = TRUE) | grepl("Synthetic", citation_evidence$Evidence, fixed = TRUE)))
})

test_that("local learning memory records dataset feedback", {
  root <- tempfile("seqsurf-learning-")
  dat <- load_test_dataset("synthetic")
  kb <- list(
    accession = "GSE123",
    geo = list(
      title = "Synthetic study",
      summary = "Processed normalized series matrix expression study.",
      overall_design = "Treated and control groups.",
      sample_count = nrow(dat$metadata),
      supplementary_files = character(),
      platform_accessions = "GPLTEST"
    ),
    pubmed = list(title = "Synthetic paper", abstract = "Expression analysis with limma and GSEA."),
    extracted = list(
      sample_groups = "Treated and control groups.",
      statistical_methods = "limma; GSEA",
      reported_pathways = "apoptosis",
      reported_biomarkers = "MYCN",
      study_objective = "Test treatment.",
      cohort_characteristics = "Synthetic cohort.",
      key_findings = "MYCN changed.",
      study_limitations = "Processed expression."
    ),
    claims = data.frame(
      claim_type = "biomarker",
      entity = "MYCN",
      comparison = "treated vs control",
      reported_direction = "up",
      reported_significance = "significant",
      evidence_sentence = "MYCN increased.",
      source = "test",
      confidence = "high",
      stringsAsFactors = FALSE
    ),
    related_publications = empty_related_publications("none"),
    extra_papers = ""
  )

  saved <- record_learning_event(
    accession = "GSE123",
    dataset_id = dat$id,
    outcome = "Useful reanalysis candidate",
    usefulness_rating = 5,
    notes = "Synthetic feedback.",
    kb = kb,
    dat = dat,
    data_root = root
  )
  log <- read_learning_log(root)
  summary <- learning_summary_table(root)
  lessons <- learning_lessons_table(root)
  prior <- learning_prior_for_study(kb, root)
  adaptive <- adaptive_reanalysis_recommendation_table(kb, root)
  adaptive_plan <- adaptive_ranked_discovery_plan(dat, kb, comparison = NULL, data_root = root)
  packet <- file.path(tempdir(), "seqsurf_shared_learning_test.zip")
  write_shared_learning_export(packet, root)
  imported_root <- tempfile("seqsurf-shared-learning-")
  imported <- import_shared_learning_packet(packet, imported_root)
  shared <- read_shared_learning_pool(imported_root)
  shared_summary <- learning_summary_table(imported_root)
  training_summary <- training_examples_summary_table(kb, dat, root)
  training_zip <- file.path(tempdir(), "seqsurf_training_examples_test.zip")
  write_training_data_export(training_zip, kb, dat, root)
  training_dir <- tempfile("seqsurf-training-unzip-")
  dir.create(training_dir)
  utils::unzip(training_zip, exdir = training_dir)
  training_jsonl <- file.path(training_dir, "seqsurf_training_examples.jsonl")
  training_lines <- readLines(training_jsonl, warn = FALSE)

  expect_equal(saved$accession[[1]], "GSE123")
  expect_equal(nrow(log), 1)
  expect_true(file.exists(learning_log_path(root)))
  expect_true("assay_signal" %in% names(log))
  expect_true(grepl("array|processed|rna-seq|unknown", log$assay_signal[[1]]))
  expect_true(any(summary$Metric == "Learning events" & summary$Value == 1))
  expect_true(any(grepl("Route pattern", lessons$Lesson)))
  expect_true(grepl("Matched local precedent", prior$Signal[[1]]))
  expect_true(any(adaptive$Source == "Local learning"))
  expect_true(any(grepl("standard reproduction|Audit the primary", adaptive_plan$Analysis, ignore.case = TRUE)))
  expect_true(file.exists(packet))
  expect_equal(imported$Imported_rows[[1]], 1)
  expect_equal(shared$accession[[1]], "")
  expect_equal(shared$dataset_id[[1]], "")
  expect_equal(shared$notes[[1]], "")
  expect_true(any(shared_summary$Metric == "Shared events" & shared_summary$Value == 1))
  expect_true(any(training_summary$Metric == "Training examples" & training_summary$Value >= 1))
  expect_true(file.exists(training_zip))
  expect_true(file.exists(training_jsonl))
  expect_true(any(grepl("\"task\":\"seqsurf_reanalysis_recommendation\"", training_lines, fixed = TRUE)))
  expect_true(any(grepl("\"target\":", training_lines, fixed = TRUE)))
})

test_that("related publication helpers return stable evidence text", {
  empty <- empty_related_publications("Not searched in test.")
  expect_setequal(
    names(empty),
    c("PubMed_ID", "Title", "Journal", "Publication_date", "Abstract", "Source", "Status")
  )
  expect_equal(empty$Status, "Not searched in test.")

  related <- data.frame(
    PubMed_ID = c("123", "456"),
    Title = c("First GEO reuse", "Second GEO reuse"),
    Journal = c("Journal A", "Journal B"),
    Publication_date = c("2024", "2025"),
    Abstract = c("Mentions GSE999.", "Also mentions GSE999."),
    Source = c("https://pubmed.ncbi.nlm.nih.gov/123/", "https://pubmed.ncbi.nlm.nih.gov/456/"),
    Status = "PubMed paper mentioning GEO accession",
    check.names = FALSE
  )
  text <- related_publications_evidence_text(related)
  expect_true(grepl("PMID 123", text, fixed = TRUE))
  expect_true(grepl("Second GEO reuse", text, fixed = TRUE))
})

test_that("Jina and PaperQA helpers prepare deterministic corpus metadata", {
  expect_equal(jina_reader_url("https://example.org/paper.pdf"), "https://r.jina.ai/https://example.org/paper.pdf")
  expect_true(all(c("Status", "Detail") %in% names(paperqa_status())))

  root <- tempfile("seqsurf-paperqa-")
  dir.create(root)
  text_dir <- file.path(root, "texts")
  dir.create(text_dir)
  text_path <- file.path(text_dir, "paper.txt")
  writeLines("Methods: limma. Results: MYCN was significant.", text_path)
  kb <- list(
    accession = "GSE123",
    documents = data.frame(
      Source = "Synthetic paper",
      URL = "https://example.org/paper",
      Text_path = text_path,
      Status = "Jina Reader indexed",
      check.names = FALSE
    )
  )

  prepared <- prepare_paperqa_corpus(kb, root)
  expect_true(file.exists(prepared$paperqa$manifest_path))
  manifest <- read.csv(prepared$paperqa$manifest_path, check.names = FALSE)
  expect_equal(manifest$Source[[1]], "Synthetic paper")
  expect_true(file.exists(manifest$Text_file[[1]]))
})

test_that("post-analysis reproducibility score penalizes direction mismatches", {
  comparison <- data.frame(
    Finding = c("Reported biomarker: MYCN", "Reported biomarker: NTRK1", "Reported pathway: apoptosis"),
    Evidence_in_app = c("supported", "direction mismatch", "missing"),
    Interpretation = c(
      "Supported by app DEG result; direction agrees",
      "Direction mismatch",
      "Missing from app pathway result"
    ),
    check.names = FALSE
  )

  score <- score_reproducibility_from_comparison(comparison)
  categories <- classify_publication_comparison(comparison)

  expect_setequal(categories, c("Supported", "Direction mismatch", "Missing"))
  expect_true(score$Value[score$Metric == "Post-analysis reproducibility score"] < 50)
  expect_equal(score$Value[score$Metric == "Direction mismatches"], 1)
  expect_equal(score$Value[score$Metric == "Missing claims"], 1)
})

test_that("GEO analyzability classifier distinguishes in-app, count, raw, and missing routes", {
  base_kb <- list(
    accession = "GSE123",
    geo = list(
      title = "Synthetic study",
      summary = "Processed normalized series matrix expression study.",
      overall_design = "Treated and control groups.",
      sample_count = 12,
      supplementary_files = character(),
      platform_accessions = "GPLTEST"
    ),
    pubmed = list(title = "Synthetic paper", abstract = "Expression analysis with limma and GSEA."),
    extracted = list(
      sample_groups = "Treated and control groups.",
      statistical_methods = "limma; GSEA",
      reported_pathways = "apoptosis",
      reported_biomarkers = "MYCN",
      study_objective = "Test treatment.",
      cohort_characteristics = "Synthetic cohort.",
      key_findings = "MYCN changed.",
      study_limitations = "Processed expression."
    ),
    claims = empty_claim_table(),
    extra_papers = ""
  )

  processed <- geo_analyzability_route(base_kb)
  processed_verified <- geo_analyzability_route(base_kb, expression_set_status = "available")
  count_kb <- base_kb
  count_kb$geo$summary <- "Supplementary raw counts matrix generated by featureCounts."
  count_kb$geo$supplementary_files <- "GSE123_counts_matrix.txt.gz"
  raw_kb <- base_kb
  raw_kb$geo$summary <- "RNA-seq FASTQ files are available from SRA."
  raw_kb$geo$supplementary_files <- "sra://SRP000000"

  expect_equal(processed$Can_run_in_app, "Yes")
  expect_equal(processed_verified$Confidence, 90)
  expect_equal(geo_analyzability_route(count_kb)$Can_run_in_app, "Partial")
  expect_match(geo_analyzability_route(count_kb)$Route, "Count-matrix")
  expect_equal(geo_analyzability_route(raw_kb)$Can_run_in_app, "No")
  expect_match(geo_analyzability_route(raw_kb)$Route, "Raw sequencing")
  expect_equal(geo_analyzability_route(NULL)$Route, "Not ready")
  expect_true(all(c("Check", "Signal", "Status", "Impact") %in% names(geo_analyzability_checks(base_kb))))

  brief <- paper_to_code_brief(load_test_dataset("synthetic"), base_kb)
  expect_true(grepl("Reanalysis route decision", brief, fixed = TRUE))
  expect_true(grepl("Coding-agent instructions", brief, fixed = TRUE))
  expect_true(grepl("Claude Code, Codex", brief, fixed = TRUE))
})
