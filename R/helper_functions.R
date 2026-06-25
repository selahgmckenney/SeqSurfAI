clean_column_names <- function(x) {
  names(x) <- trimws(names(x))
  x
}

available_datasets <- function(data_root = "data") {
  dirs <- list.dirs(data_root, recursive = FALSE, full.names = FALSE)
  dirs[file.exists(file.path(data_root, dirs, "dataset_info.rds"))]
}

dataset_choice_list <- function(datasets, data_root = "data") {
  if (length(datasets) == 0) {
    return(c("No active dataset" = ""))
  }
  stats::setNames(
    datasets,
    vapply(
      datasets,
      function(id) readRDS(file.path(data_root, id, "dataset_info.rds"))$display_name,
      character(1)
    )
  )
}

empty_gsea_table <- function() {
  data.frame(
    collection = character(),
    msigdb_version = character(),
    contrast = character(),
    pathway = character(),
    pval = numeric(),
    padj = numeric(),
    log2err = numeric(),
    ES = numeric(),
    NES = numeric(),
    size = integer(),
    leadingEdge = character(),
    check.names = FALSE
  )
}

load_rnaseq_dataset <- function(dataset_id, data_root = "data") {
  validate(
    need(!is.null(dataset_id) && nzchar(dataset_id), "Import a GEO dataset or add an app-ready dataset to explore results.")
  )
  dataset_dir <- file.path(data_root, dataset_id)

  required <- c(
    "metadata.rds",
    "vsd_matrix.rds",
    "pca_df.rds",
    "pca_variance.rds",
    "pca_loadings.rds",
    "deg_results.csv",
    "gsea_hallmark.csv",
    "pathway_scores.csv",
    "dataset_info.rds"
  )
  missing <- required[!file.exists(file.path(dataset_dir, required))]
  validate(
    need(length(missing) == 0, paste("Missing dataset files:", paste(missing, collapse = ", ")))
  )

  list(
    id = dataset_id,
    info = readRDS(file.path(dataset_dir, "dataset_info.rds")),
    metadata = clean_column_names(readRDS(file.path(dataset_dir, "metadata.rds"))),
    vsd = readRDS(file.path(dataset_dir, "vsd_matrix.rds")),
    pca = clean_column_names(readRDS(file.path(dataset_dir, "pca_df.rds"))),
    pca_variance = readRDS(file.path(dataset_dir, "pca_variance.rds")),
    loadings = clean_column_names(readRDS(file.path(dataset_dir, "pca_loadings.rds"))),
    deg = clean_column_names(read.csv(file.path(dataset_dir, "deg_results.csv"), check.names = FALSE)),
    gsea = clean_column_names(read.csv(file.path(dataset_dir, "gsea_hallmark.csv"), check.names = FALSE)),
    pathway_scores = read.csv(file.path(dataset_dir, "pathway_scores.csv"), row.names = 1, check.names = FALSE)
  )
}

info_title <- function(label, description) {
  tags$span(
    label,
    tags$sup(" [?]", class = "text-secondary"),
    title = description,
    style = "cursor: help;"
  )
}

external_link <- function(label, url) {
  tags$a(
    label,
    icon("arrow-up-right-from-square"),
    href = url,
    target = "_blank",
    rel = "noopener noreferrer",
    class = "btn btn-outline-primary btn-sm me-2 mb-2"
  )
}

study_stat <- function(value, label) {
  tags$div(
    class = "study-stat",
    tags$strong(value),
    tags$span(label)
  )
}

definition_table <- function(items) {
  tags$dl(
    class = "study-definitions",
    lapply(items, function(item) {
      tagList(tags$dt(item[[1]]), tags$dd(item[[2]]))
    })
  )
}

dataset_study_content <- function(dataset_id, metadata, info) {
  if (dataset_id == "anbl17p1_paxgene") {
    sex_counts <- table(metadata$inferred_sex, useNA = "no")
    mycn_patient_counts <- metadata |>
      dplyr::distinct(.data$patient_id, .data$mycn_expression_group) |>
      dplyr::count(.data$mycn_expression_group)
    response_patient_counts <- metadata |>
      dplyr::distinct(.data$patient_id, .data$inferred_tumor_response_group) |>
      dplyr::count(.data$inferred_tumor_response_group)

    mycn_high <- mycn_patient_counts$n[
      mycn_patient_counts$mycn_expression_group == "MYCN expression-high"
    ]
    mycn_low <- mycn_patient_counts$n[
      mycn_patient_counts$mycn_expression_group == "MYCN not expression-high"
    ]
    response_strong <- response_patient_counts$n[
      response_patient_counts$inferred_tumor_response_group == "Strong tumor-score decrease"
    ]
    response_weak <- response_patient_counts$n[
      response_patient_counts$inferred_tumor_response_group == "Weak/no tumor-score decrease"
    ]

    return(tagList(
      tags$div(
        class = "study-heading",
        tags$span("CLINICAL TRIAL PAXGENE BLOOD COHORT", class = "study-kicker"),
        tags$h2("Neuroblastoma: ANBL17P1 PAXgene Whole Blood"),
        tags$p(
          "Longitudinal PAXgene whole-blood RNA-seq from the ANBL17P1 neuroblastoma ",
          "clinical-trial dataset. This is blood/host transcriptomics with an exploratory ",
          "tumor-RNA signature layer, not bulk tumor tissue."
        )
      ),
      tags$div(
        class = "study-stats",
        study_stat(nrow(metadata), "standard ANBL blood RNA-seq samples"),
        study_stat(length(unique(metadata$patient_id)), "standard patients"),
        study_stat(paste0(unname(mycn_high), " / ", unname(mycn_low)), "MYCN expression-high / not-high patients"),
        study_stat(paste0(unname(sex_counts["Male"]), " / ", unname(sex_counts["Female"])), "male / female samples"),
        study_stat(paste0(unname(response_strong), " / ", unname(response_weak)), "strong / weak-no tumor-score decrease patients")
      ),
      tags$div(
        class = "study-note",
        tags$strong("Interpretation note: "),
        "the response grouping in this app is inferred from the exploratory MRD-core ",
        "tumor-RNA score. It should not be presented as clinical response until linked ",
        "to curated response, imaging, marrow, or outcome metadata."
      ),
      tags$div(
        class = "study-section-grid",
        tags$section(
          tags$h3("Clinical-Trial Context"),
          tags$p(
            "ANBL17P1 includes samples across pretreatment, induction, recovery, and ",
            "post-consolidation immunotherapy windows. DIN-containing induction cycles ",
            "and GM-CSF timing are important for interpreting blood-composition changes."
          ),
          tags$p(
            "The app keeps nonstandard/control records out of the main ANBL dataset and ",
            "uses the 227 standard count-matrix samples from 30 patients."
          )
        ),
        tags$section(
          tags$h3("Comparisons In This App"),
          tags$p(
            "Precomputed contrasts include MYCN expression-high versus not-high, male ",
            "versus female, and strong tumor-score decrease versus weak/no decrease."
          ),
          tags$p(
            "DE contrasts are exploratory. MYCN and sex contrasts use the earliest ",
            "available sample per patient; the tumor-score-decrease contrast uses ",
            "pretreatment samples from evaluable patients."
          )
        )
      ),
      tags$section(
        class = "study-section",
        tags$h3("Inferred Variables"),
        definition_table(list(
          list("MYCN expression group", "Expression-based screen from the earliest available RNA sample. This is not clinical MYCN amplification."),
          list("Inferred sex", "Expression-based sex inference from XIST and Y-marker expression, checked against available metadata where possible."),
          list("Tumor-score decrease", "Patient-level inferred response based on MRD-core neuroblastoma tumor-RNA score decrease from pretreatment. Strong decrease is defined as minimum post-pretreatment delta <= -0.5."),
          list("Expression unit", "DESeq2 VST values for visualization; differential-expression contrasts use log2 CPM with limma.")
        ))
      ),
      tags$section(
        class = "study-section",
        tags$h3("Local Resources"),
        external_link("ANBL17P1 exploratory report", "../neuroblastoma_ANBL17P1/exploratory_analysis/ANBL17P1_story_exploratory_report.html"),
        external_link("Metadata documentation", "../neuroblastoma_ANBL17P1/metadata_master/documentation/neuroblastoma_ANBL17P1_Master_Metadata_Documentation.html")
      )
    ))
  }

  if (dataset_id == "neuroblastoma_gse85047") {
    mycn_counts <- table(metadata$mycn_status, useNA = "no")
    stage_counts <- table(metadata$stage4_status, useNA = "no")
    sex_counts <- table(metadata$inferred_sex, useNA = "no")

    return(tagList(
      tags$div(
        class = "study-heading",
        tags$span("EXTERNAL PUBLIC TUMOR COHORT", class = "study-kicker"),
        tags$h2("Neuroblastoma: GSE85047 / NRC-283"),
        tags$p(
          "A public external neuroblastoma cohort used here to test whether ",
          "signals discovered in GSE49711/GSE62564 replicate outside the original ",
          "SEQC-498 training and validation split."
        )
      ),
      tags$div(
        class = "study-stats",
        study_stat(nrow(metadata), "primary untreated tumors"),
        study_stat(paste0(unname(mycn_counts["Amplified"]), " / ", unname(mycn_counts["Non-amplified"])), "MYCN amplified / non-amplified"),
        study_stat(paste0(unname(stage_counts["Stage 4"]), " / ", unname(stage_counts["Other stages"])), "stage 4 / other stages"),
        study_stat(paste0(unname(sex_counts["Male"]), " / ", unname(sex_counts["Female"])), "inferred male / female"),
        study_stat("Affymetrix HuEx-1_0-st", "microarray platform")
      ),
      tags$div(
        class = "study-note",
        tags$strong("Replication note: "),
        "this cohort is independent of GSE49711/GSE62564. The processed matrix contains ",
        "core transcript-cluster expression summarized by the original submitters. ",
        "Clean MAGEA-family validation is limited because the public processed matrix ",
        "does not retain unambiguous MAGEA transcript clusters."
      ),
      tags$div(
        class = "study-section-grid",
        tags$section(
          tags$h3("Study And Sample Origin"),
          tags$p(
            "GSE85047 contains 283 primary untreated neuroblastoma tumors from the ",
            "Neuroblastoma Research Consortium series. GEO annotations include age at ",
            "diagnosis, INSS stage, MYCN amplification, overall survival, and ",
            "progression-free survival."
          ),
          tags$p(
            "The app uses the public processed series-matrix values rather than the raw ",
            "CEL archive. The submitted data were RMA-sketch normalized and batch corrected."
          )
        ),
        tags$section(
          tags$h3("Data Used In This App"),
          tags$p(
            "GPL5175 transcript-cluster IDs were mapped to gene symbols using GEO's ",
            "`gene_assignment` annotation field. Duplicate gene symbols were collapsed ",
            "by averaging expression across mapped transcript clusters."
          ),
          tags$p(
            "Precomputed contrasts include MYCN amplification, INSS stage 4 status, ",
            "progression, disease death, and age at diagnosis above or below 18 months."
          )
        )
      ),
      tags$section(
        class = "study-section",
        tags$h3("Clinical And Inferred Variables"),
        definition_table(list(
          list("MYCN status", "Submitted MYCN amplification status from GEO."),
          list("INSS stage", "Submitted International Neuroblastoma Staging System stage."),
          list("Age at diagnosis", "Submitted age at diagnosis in days, also grouped above or below 18 months."),
          list("Overall survival", "Submitted event_overall and overall_survival_time fields."),
          list("Progression-free survival", "Submitted event_progression_free and progression-free survival time fields; GEO misspells the time field as `profression_free_survival_time`."),
          list("Inferred sex", "Not submitted in GSE85047. Inferred from Y-linked expression markers DDX3Y and EIF1AY in the processed matrix.")
        ))
      ),
      tags$section(
        class = "study-section",
        tags$h3("Official Resources"),
        external_link("GEO GSE85047", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85047"),
        external_link("GEO GPL5175", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL5175"),
        external_link("BioProject PRJNA336020", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA336020")
      )
    ))
  }

  if (dataset_id == "neuroblastoma_gse49711") {
    sex_counts <- table(metadata$sex, useNA = "no")
    training_counts <- table(metadata$dataset_assignment, useNA = "no")

    return(tagList(
      tags$div(
        class = "study-heading",
        tags$span("PUBLIC TUMOR COHORT", class = "study-kicker"),
        tags$h2("Neuroblastoma: GSE49711"),
        tags$p(
          "A public cohort of primary neuroblastoma tumors assembled to evaluate ",
          "RNA-seq and microarray expression models for clinical endpoint prediction."
        )
      ),
      tags$div(
        class = "study-stats",
        study_stat(nrow(metadata), "primary tumors"),
        study_stat(paste0(unname(sex_counts["Male"]), " / ", unname(sex_counts["Female"])), "male / female"),
        study_stat(
          paste0(unname(training_counts["Training"]), " / ", unname(training_counts["Validation"])),
          "training / validation"
        ),
        study_stat("Illumina HiSeq 2000", "RNA-seq platform")
      ),
      tags$div(
        class = "study-section-grid",
        tags$section(
          tags$h3("Study And Sample Origin"),
          tags$p(
            "GSE49711 contains 498 primary neuroblastoma tumor samples from a German ",
            "clinical cohort. The original study profiled the same tumors by RNA sequencing ",
            "and microarray and divided the cohort evenly into training and validation sets."
          ),
          tags$p(
            "Raw sequencing data were not submitted because of patient-privacy restrictions. ",
            "The app therefore uses public processed expression rather than raw read counts."
          )
        ),
        tags$section(
          tags$h3("Data Used In This App"),
          tags$p(
            "Clinical annotations come from GSE49711. Gene-level expression comes from ",
            "GSE62564, the official reanalysis of the same 498 tumors, and is represented ",
            "as processed log2 RPM."
          ),
          tags$p(
            "RefSeq features were mapped to gene symbols, duplicate mappings were collapsed, ",
            "and genes with expression above zero in at least 10% of tumors were retained. ",
            "Differential expression was performed with limma because the public matrix was ",
            "already normalized."
          )
        )
      ),
      tags$section(
        class = "study-section",
        tags$h3("Clinical Variables"),
        definition_table(list(
          list("Dataset assignment", "Training and validation subsets from the original prediction study; this is a study-design variable, not a biological phenotype."),
          list("Sex", "Submitted sex: 287 male and 211 female tumors."),
          list("Age at diagnosis", "Age in days at initial diagnosis; the app also includes an 18-month age grouping."),
          list("MYCN status", "MYCN amplified, non-amplified, or unavailable."),
          list("High risk", "The submitted clinical high-risk classification."),
          list("INSS stage", "International Neuroblastoma Staging System stage 1, 2, 3, 4, or 4S."),
          list("Class label", "A study-specific favorable or unfavorable label defined only for maximally divergent clinical courses."),
          list("Progression and disease death", "Indicators of tumor progression and death attributed to neuroblastoma.")
        ))
      ),
      tags$section(
        class = "study-section",
        tags$h3("Official Resources"),
        external_link("GEO GSE49711", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711"),
        external_link("GEO GSE62564", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62564"),
        external_link("BioProject PRJNA214797", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA214797")
      )
    ))
  }

  if (dataset_id == "osteosarcoma") {
    sex_counts <- table(metadata$inferred_sex, useNA = "no")
    return(tagList(
      tags$div(
        class = "study-heading",
        tags$span("PHASE 2 CLINICAL TRIAL COHORT", class = "study-kicker"),
        tags$h2("Osteosarcoma: SARC038"),
        tags$p(
          "Whole-blood RNA-seq samples from a trial of regorafenib plus nivolumab ",
          "in patients with refractory or recurrent osteosarcoma."
        )
      ),
      tags$div(
        class = "study-stats",
        study_stat(nrow(metadata), "RNA-seq samples"),
        study_stat(length(unique(metadata$patient_id)), "patients in app dataset"),
        study_stat(paste0(unname(sex_counts["Male"]), " / ", unname(sex_counts["Female"])), "male / female samples"),
        study_stat(length(unique(metadata$timepoint_group)), "sample timepoint groups")
      ),
      tags$div(
        class = "study-note",
        tags$strong("Important distinction: "),
        "the app contains 67 longitudinal blood specimens from 20 patients. ",
        "These sample and patient counts are not the same as the trial's estimated enrollment of 48 participants."
      ),
      tags$div(
        class = "study-section-grid",
        tags$section(
          tags$h3("Trial Overview"),
          tags$p(
            "SARC038 (NCT04803877) is a single-group, open-label phase 2 study sponsored ",
            "by the Sarcoma Alliance for Research through Collaboration. The regimen combines ",
            "regorafenib, a multikinase inhibitor, with nivolumab, an immune checkpoint inhibitor."
          ),
          tags$p(
            "The primary comparison is the 4-month progression-free survival rate against ",
            "historical outcomes from regorafenib alone. Secondary outcomes include response, ",
            "progression-free survival, overall survival, and treatment toxicity."
          )
        ),
        tags$section(
          tags$h3("Samples In This App"),
          tags$p(
            "The molecular dataset consists of whole-blood RNA-seq collected across pre-therapy, ",
            "post-therapy, and progression/end-of-treatment timepoints. Because patients can ",
            "contribute multiple samples, analyses must account for the difference between ",
            "specimen-level and patient-level counts."
          ),
          tags$p(
            "The default expression comparison groups stable disease or partial response ",
            "(SD_PR) against progressive disease (PD). Inferred sex was included as a covariate ",
            "in the differential-expression model."
          )
        )
      ),
      tags$section(
        class = "study-section",
        tags$h3("Trial Context"),
        definition_table(list(
          list("Population", "Patients age 5 years or older with relapsed or refractory high-grade osteosarcoma after at least one prior systemic therapy."),
          list("Design", "Single-arm, Simon two-stage, historically controlled phase 2 trial."),
          list("Interventions", "Oral regorafenib plus intravenous nivolumab, with age-specific dosing."),
          list("Primary outcome", "Four-month progression-free survival assessed against historical regorafenib-alone controls."),
          list("Locations", "Eight participating centers in the United States."),
          list("Record status", "Active, not recruiting; ClinicalTrials.gov record last verified October 2025.")
        ))
      ),
      tags$section(
        class = "study-section",
        tags$h3("Official Resources"),
        external_link("ClinicalTrials.gov NCT04803877", "https://clinicaltrials.gov/study/NCT04803877"),
        external_link("SARC", "https://sarctrials.org/"),
        external_link("NCI Osteosarcoma Information", "https://www.cancer.gov/types/bone/patient/osteosarcoma-treatment-pdq")
      )
    ))
  }

  tagList(tags$h2(info$display_name), tags$p(info$notes))
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

cutoff_label <- function(x) {
  sub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
}

round_for_display <- function(df, digits = 3) {
  numeric <- vapply(df, is.numeric, logical(1))
  df[numeric] <- lapply(df[numeric], function(x) signif(x, digits))
  df
}

pc_columns <- function(pca_df) {
  grep("^PC[0-9]+$", names(pca_df), value = TRUE)
}

numeric_columns <- function(df) {
  names(df)[vapply(df, is.numeric, logical(1))]
}

categorical_columns <- function(df) {
  names(df)[vapply(df, function(x) is.character(x) || is.factor(x) || is.logical(x), logical(1))]
}

metadata_summary <- function(metadata, info) {
  sex_col <- if ("inferred_sex" %in% names(metadata)) {
    "inferred_sex"
  } else if ("sex" %in% names(metadata)) {
    "sex"
  } else {
    NA_character_
  }
  cols <- c(
    Samples = info$sample_id_col,
    Patients = info$patient_id_col,
    Response_groups = info$default_group_col,
    Timepoints = info$default_timepoint_col,
    Sex_categories = sex_col
  )
  data.frame(
    Metric = names(cols),
    Value = vapply(cols, function(col) {
      if (!col %in% names(metadata)) return(NA_character_)
      if (col == info$sample_id_col) return(as.character(nrow(metadata)))
      as.character(length(unique(metadata[[col]][!is.na(metadata[[col]])])))
    }, character(1)),
    check.names = FALSE
  )
}

filter_deg <- function(deg, contrast, padj_cutoff, lfc_cutoff) {
  out <- deg
  if ("contrast" %in% names(out) && !is.null(contrast)) {
    out <- out[out$contrast == contrast, , drop = FALSE]
  }
  out$direction <- dplyr::case_when(
    !is.na(out$padj) & out$padj <= padj_cutoff & out$log2FoldChange >= lfc_cutoff ~ "Up",
    !is.na(out$padj) & out$padj <= padj_cutoff & out$log2FoldChange <= -lfc_cutoff ~ "Down",
    TRUE ~ "Not significant"
  )
  first <- intersect(c("contrast", "gene", "direction"), names(out))
  remaining <- setdiff(names(out), first)
  out[, c(first, remaining), drop = FALSE]
}

top_deg_genes <- function(deg, n = 50, contrast = NULL) {
  out <- deg
  if ("contrast" %in% names(out) && !is.null(contrast)) {
    out <- out[out$contrast == contrast, , drop = FALSE]
  }
  out |>
    dplyr::filter(!is.na(.data$padj)) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange))) |>
    dplyr::slice_head(n = n) |>
    dplyr::pull(.data$gene)
}

parse_gene_list <- function(text) {
  genes <- unlist(strsplit(text, "[,;\\n\\r\\t ]+"))
  unique(trimws(genes[nzchar(trimws(genes))]))
}

expression_long <- function(vsd, metadata, genes, sample_col = "sample") {
  genes <- intersect(genes, rownames(vsd))
  validate(need(length(genes) > 0, "No requested genes were found in the expression matrix."))

  expr <- as.data.frame(t(vsd[genes, , drop = FALSE]))
  expr[[sample_col]] <- rownames(expr)
  merged <- dplyr::left_join(expr, metadata, by = sample_col)
  tidyr::pivot_longer(
    merged,
    cols = dplyr::all_of(genes),
    names_to = "gene",
    values_to = "expression"
  )
}

leading_edge_matches <- function(gsea, gene) {
  if (!"leadingEdge" %in% names(gsea)) return(gsea[0, , drop = FALSE])
  pattern <- paste0("(^|;)", gene, "(;|$)")
  gsea[grepl(pattern, gsea$leadingEdge), , drop = FALSE]
}

add_pathway_overlap_stats <- function(gsea) {
  if (!"leadingEdge" %in% names(gsea)) {
    gsea$leading_edge_genes <- NA_integer_
    return(gsea)
  }
  gsea$leading_edge_genes <- vapply(
    strsplit(ifelse(is.na(gsea$leadingEdge), "", gsea$leadingEdge), ";", fixed = TRUE),
    function(x) sum(nzchar(x)),
    integer(1)
  )
  if ("size" %in% names(gsea)) {
    gsea$leading_edge_fraction <- gsea$leading_edge_genes / gsea$size
  }
  gsea
}

safe_neg_log10 <- function(x) {
  x <- ifelse(is.na(x) | x <= 0, NA_real_, x)
  -log10(x)
}
