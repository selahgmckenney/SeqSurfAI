required_packages <- c(
  "shiny", "ggplot2", "dplyr", "tidyr", "DT", "pheatmap", "bslib", "plotly",
  "enrichR", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install required packages before running this app: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(bslib)

# enrichR 3.4 initializes these options only when the package is attached.
# The app calls the namespace directly to avoid checking every organism endpoint at startup.
options(
  enrichR.base.address = "https://maayanlab.cloud/Enrichr/",
  enrichR.live = TRUE,
  enrichR.quiet = FALSE,
  enrichR.sites = c(
    "Enrichr", "FlyEnrichr", "WormEnrichr",
    "YeastEnrichr", "FishEnrichr", "OxEnrichr"
  ),
  enrichR.sites.base.address = "https://maayanlab.cloud/"
)

source("R/helper_functions.R")
source("R/plot_pca.R")
source("R/plot_volcano.R")
source("R/plot_heatmap.R")
source("R/enrichment_functions.R")
source("R/ai_assistant_functions.R")
source("R/geo_import_functions.R")
source("R/upload_functions.R")

datasets <- get("available_datasets", mode = "function")("data")
active_dataset_choices <- function(dataset_ids) {
  c("No active dataset - select saved results or import/create new results" = "", dataset_choice_list(dataset_ids, "data"))
}
dataset_choices <- active_dataset_choices(datasets)
selected_dataset <- ""

workflow_story_panel <- function(active = "start") {
  steps <- list(
    start = list(label = "Start", title = "Assess", body = "Enter GEO, score reanalysis value, and decide the route."),
    train = list(label = "Prepare Evidence", title = "Study Memory", body = "Collect papers, claims, methods, uncertainty, and citations."),
    validate = list(label = "Validate", title = "Reanalyze", body = "Run the app workflow and score paper agreement."),
    discovery = list(label = "Discovery", title = "Next Analyses", body = "Rank follow-up analyses and export a report.")
  )
  order <- names(steps)
  active_idx <- match(active, order)
  tags$div(
    class = "story-panel",
    tags$div(
      class = "story-panel-copy",
      tags$span("SeqSurf guided workflow"),
      tags$strong("From public study to reproducible reanalysis"),
      tags$p("Assess whether the dataset can be analyzed, build paper-grounded study memory, validate results, then decide what to try next.")
    ),
    tags$div(
      class = "story-steps",
      lapply(seq_along(steps), function(i) {
        step <- steps[[i]]
        status <- if (i < active_idx) "done" else if (i == active_idx) "active" else "waiting"
        tags$div(
          class = paste("story-step", status),
          tags$span(step$label),
          tags$strong(step$title),
          tags$p(step$body)
        )
      })
    )
  )
}

ui <- page_navbar(
  id = "main_nav",
  title = tags$div(
    class = "app-brand",
    tags$span("SeqSurf AI"),
    tags$small("Scientific AI for transcriptomic reanalysis")
  ),
  theme = bs_theme(version = 5, bootswatch = "cosmo"),
  header = tags$head(tags$style(HTML("
    body { background: #f6f8fa; }
    .navbar { border-bottom: 1px solid #dfe3e8; box-shadow: none; padding: 0; background: #fff; }
    .navbar .container-fluid { align-items: stretch; gap: 0; flex-direction: column; flex-wrap: nowrap; padding: 0; }
    .navbar-header { width: 100%; padding: 0.78rem 1.35rem 0.64rem; border-bottom: 1px solid #edf0f2; }
    .navbar-brand { display: block; width: 100%; margin: 0; padding: 0; border-right: 0; }
    .navbar-collapse { width: 100%; overflow-x: auto; scrollbar-width: thin; padding: 0 1.2rem; }
    .app-brand { display: flex; flex-direction: column; justify-content: center; line-height: 1.05; letter-spacing: 0; min-width: 0; }
    .app-brand span { color: #12263f; font-size: 1.58rem; font-weight: 850; }
    .app-brand small { color: #59636e; font-size: 0.78rem; font-weight: 720; letter-spacing: 0; text-transform: none; margin-top: 0.18rem; white-space: normal; }
    .navbar-nav { align-items: center; gap: 1.32rem; flex-wrap: nowrap; line-height: 1.1; min-width: max-content; }
    .navbar-nav .nav-link { position: relative; color: #5f666e; font-size: 1.08rem; font-weight: 540; padding: 0.94rem 0 0.9rem; border: 0; border-radius: 0; margin: 0; white-space: nowrap; }
    .navbar-nav .nav-link:hover { color: #111827; background: transparent; }
    .navbar-nav .nav-link::after { content: ''; position: absolute; left: 0; right: 0; bottom: 0.48rem; height: 3px; background: transparent; }
    .navbar-nav .nav-link.active { color: #111827 !important; background: transparent; box-shadow: none; font-weight: 760; }
    .navbar-nav .nav-link.active::after { background: #111827; }
    @media (max-width: 1200px) {
      .navbar-nav { gap: 1rem; }
      .navbar-nav .nav-link { font-size: 1rem; }
    }
    .nav-underline .nav-link.active, .nav-underline .show > .nav-link { border-bottom-color: transparent; }
    .study-heading { max-width: 900px; padding: 0.5rem 0 1rem; }
    .study-heading h2 { margin: 0.2rem 0 0.5rem; font-size: 2rem; }
    .study-heading p { color: #4b5563; font-size: 1.05rem; margin: 0; }
    .study-kicker { color: #0b6b5f; font-size: 0.75rem; font-weight: 700; letter-spacing: 0.08em; }
    .study-stats { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 1px; background: #d9dee3; border: 1px solid #d9dee3; margin-bottom: 1.5rem; }
    .study-stat { display: flex; flex-direction: column; gap: 0.2rem; background: #fff; padding: 1rem; min-width: 0; }
    .study-stat strong { color: #16324f; font-size: 1.25rem; }
    .study-stat span { color: #59636e; font-size: 0.86rem; }
    .study-section-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 2rem; padding: 0.5rem 0 1.5rem; }
    .study-section, .study-section-grid section { border-top: 3px solid #0b6b5f; padding-top: 0.8rem; margin-bottom: 1.5rem; }
    .study-section h3, .study-section-grid h3 { font-size: 1.15rem; margin-bottom: 0.7rem; }
    .study-definitions { display: grid; grid-template-columns: minmax(150px, 0.28fr) 1fr; column-gap: 1.25rem; margin: 0; }
    .study-definitions dt, .study-definitions dd { border-top: 1px solid #e1e5e9; margin: 0; padding: 0.65rem 0; }
    .study-definitions dt { color: #16324f; font-weight: 600; }
    .study-note { border-left: 4px solid #d6862c; background: #fff8ed; padding: 0.8rem 1rem; margin-bottom: 1.5rem; }
    .assistant-workspace { padding: 1rem; background: #f6f8fa; }
    .assistant-workspace .card { border: 1px solid #d8dee4; box-shadow: 0 1px 2px rgba(16,24,40,0.04); }
    .story-panel { display: grid; grid-template-columns: minmax(240px, 0.35fr) 1fr; gap: 1rem; align-items: stretch; border: 1px solid #d8e3ed; border-left: 5px solid #0b7c86; background: #fff; border-radius: 10px; padding: 1rem; margin-bottom: 1rem; }
    .story-panel-copy span { display: block; color: #0b7c86; font-size: 0.76rem; font-weight: 820; text-transform: uppercase; letter-spacing: 0.04em; margin-bottom: 0.25rem; }
    .story-panel-copy strong { display: block; color: #12263f; font-size: 1.08rem; line-height: 1.2; margin-bottom: 0.35rem; }
    .story-panel-copy p { color: #59636e; margin: 0; line-height: 1.35; }
    .story-steps { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 0.65rem; }
    .story-step { border: 1px solid #e1e6eb; background: #fbfcfd; border-radius: 8px; padding: 0.75rem; min-width: 0; }
    .story-step span { display: block; color: #59636e; font-size: 0.72rem; font-weight: 800; text-transform: uppercase; letter-spacing: 0.04em; margin-bottom: 0.18rem; }
    .story-step strong { display: block; color: #12263f; font-size: 0.92rem; margin-bottom: 0.25rem; }
    .story-step p { color: #59636e; font-size: 0.78rem; line-height: 1.25; margin: 0; }
    .story-step.done { border-color: #2e7d32; background: #f1f8f2; }
    .story-step.active { border-color: #1769aa; background: #eef6fc; box-shadow: 0 0 0 3px rgba(23,105,170,0.12); }
    .assistant-stepper { display: grid; grid-template-columns: repeat(5, minmax(120px, 1fr)); gap: 0.75rem; margin-bottom: 1rem; }
    .assistant-step { background: #fff; border: 1px solid #d8dee4; border-radius: 8px; padding: 0.75rem; min-height: 88px; }
    .assistant-step strong { display: block; color: #16324f; font-size: 0.95rem; margin-bottom: 0.25rem; }
    .assistant-step span { display: block; color: #59636e; font-size: 0.82rem; line-height: 1.25; }
    .assistant-step.done { border-color: #2e7d32; background: #f1f8f2; }
    .assistant-step.active { border-color: #1769aa; box-shadow: 0 0 0 3px rgba(23,105,170,0.16), 0 0 18px rgba(23,105,170,0.18); }
    .assistant-step.waiting { opacity: 0.72; }
    .assistant-action-stack .btn { width: 100%; margin-bottom: 0.55rem; text-align: left; white-space: normal; }
    .assistant-action-stack .next-action { background: #1769aa; border-color: #1769aa; color: #fff; box-shadow: 0 0 0 3px rgba(23,105,170,0.18), 0 0 20px rgba(23,105,170,0.35); }
    .import-workspace { padding: 1rem; background: #f6f8fa; }
    .import-workspace .card { border: 1px solid #d8dee4; box-shadow: 0 1px 2px rgba(16,24,40,0.04); }
    .import-stepper { display: grid; grid-template-columns: repeat(5, minmax(120px, 1fr)); gap: 0.75rem; margin-bottom: 1rem; }
    .import-step { background: #fff; border: 1px solid #d8dee4; border-radius: 8px; padding: 0.75rem; min-height: 88px; }
    .import-step strong { display: block; color: #16324f; font-size: 0.95rem; margin-bottom: 0.25rem; }
    .import-step span { display: block; color: #59636e; font-size: 0.82rem; line-height: 1.25; }
    .import-step.done { border-color: #2e7d32; background: #f1f8f2; }
    .import-step.active { border-color: #1769aa; box-shadow: 0 0 0 3px rgba(23,105,170,0.16), 0 0 18px rgba(23,105,170,0.18); }
    .import-step.waiting { opacity: 0.72; }
    .import-action-stack .btn { width: 100%; margin-bottom: 0.55rem; text-align: left; white-space: normal; }
    .import-action-stack .next-action { background: #1769aa; border-color: #1769aa; color: #fff; box-shadow: 0 0 0 3px rgba(23,105,170,0.18), 0 0 20px rgba(23,105,170,0.35); }
    .assistant-hint { border-left: 4px solid #1769aa; background: #eef6fc; padding: 0.85rem 1rem; margin-bottom: 1rem; border-radius: 6px; color: #16324f; }
    .assistant-empty { color: #6b7280; padding: 0.8rem 0; }
    .seq-avatar-panel { display: flex; align-items: center; gap: 0.9rem; border: 1px solid #cdd9e5; background: #ffffff; border-radius: 8px; padding: 0.9rem 1rem; margin-bottom: 1rem; }
    #ai_status_panel.recalculating { opacity: 1 !important; }
    #ai_status_panel.recalculating * { opacity: 1 !important; }
    .seq-avatar-panel.recalculating { opacity: 1 !important; }
    .seq-avatar-panel.recalculating .seq-avatar, .seq-avatar-panel.recalculating .seq-avatar-copy { opacity: 1 !important; }
    .seq-avatar { width: 48px; height: 48px; border-radius: 50%; display: grid; place-items: center; color: #fff; background: linear-gradient(135deg, #1769aa, #0b6b5f); box-shadow: 0 0 0 4px rgba(23,105,170,0.12); flex: 0 0 auto; }
    .seq-avatar i { font-size: 1.15rem; }
    .seq-avatar-copy { display: flex; flex-direction: column; min-width: 0; position: relative; }
    .seq-avatar-copy strong { color: #12263f; font-size: 0.98rem; }
    .seq-avatar-copy span { color: #1769aa; font-weight: 750; font-size: 0.86rem; }
    .seq-avatar-copy p { color: #4b5563; margin: 0.15rem 0 0; font-size: 0.9rem; line-height: 1.25; }
    .ai-guide-card { border-left: 5px solid #1769aa !important; background: #ffffff; }
    .ai-guide-card .card-header { background: #f0f7fc; color: #12263f; font-weight: 760; }
    .ai-guide-body { color: #26313d; font-size: 0.98rem; line-height: 1.45; }
    .ai-guide-actions { display: flex; flex-wrap: wrap; justify-content: center; gap: 0.55rem; margin-top: 0.75rem; }
    .ai-guide-actions .btn { white-space: normal; }
    .ai-narrative { color: #26313d; font-size: 0.98rem; line-height: 1.45; }
    .ai-narrative h1, .ai-narrative h2, .ai-narrative h3 { display: none; }
    .ai-narrative h4, .ai-narrative h5 { color: #12263f; font-size: 1rem; font-weight: 780; margin: 1rem 0 0.3rem; }
    .ai-narrative h4:first-child, .ai-narrative h5:first-child { margin-top: 0; }
    .ai-narrative p { margin: 0 0 0.72rem; max-width: 1100px; }
    .ai-narrative ul, .ai-narrative ol { margin: 0.2rem 0 0.8rem 1.25rem; padding: 0; }
    .ai-narrative li { margin-bottom: 0.28rem; }
    .ai-narrative strong { color: #12263f; font-weight: 780; }
    .ai-loading { display: flex; align-items: center; gap: 0.8rem; color: #26313d; background: #fbfcfd; border: 1px solid #e1e6eb; border-radius: 8px; padding: 0.9rem 1rem; max-width: 760px; margin: 0 auto; }
    .ai-loading-pulse { width: 0.82rem; height: 0.82rem; border-radius: 50%; background: #1769aa; box-shadow: 0 0 0 rgba(23,105,170,0.55); animation: aiPulse 1.45s infinite; flex: 0 0 auto; }
    @keyframes aiPulse {
      0% { box-shadow: 0 0 0 0 rgba(23,105,170,0.45); }
      70% { box-shadow: 0 0 0 10px rgba(23,105,170,0); }
      100% { box-shadow: 0 0 0 0 rgba(23,105,170,0); }
    }
    .ai-recommendation-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 0.8rem; }
    .ai-recommendation-section { border: 1px solid #e1e6eb; border-radius: 8px; background: #fff; padding: 0.95rem 1rem; min-width: 0; }
    .ai-recommendation-section:first-child { grid-column: 1 / -1; }
    .ai-recommendation-section h4 { color: #12263f; font-size: 0.98rem; font-weight: 820; margin: 0 0 0.42rem; }
    .ai-recommendation-section .ai-section-body { color: #26313d; line-height: 1.45; }
    .ai-recommendation-section .ai-section-body p:last-child,
    .ai-recommendation-section .ai-section-body ul:last-child,
    .ai-recommendation-section .ai-section-body ol:last-child { margin-bottom: 0; }
    .ranked-discovery-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 0.85rem; }
    .ranked-discovery-card { border: 1px solid #dfe7ee; border-radius: 10px; background: #fff; padding: 0.95rem; min-width: 0; border-top: 4px solid #1769aa; }
    .ranked-discovery-card .rank { display: inline-flex; align-items: center; justify-content: center; width: 1.8rem; height: 1.8rem; border-radius: 999px; background: #eef6fc; color: #1769aa; font-weight: 850; margin-bottom: 0.55rem; }
    .ranked-discovery-card h4 { color: #12263f; font-size: 1rem; line-height: 1.22; margin: 0 0 0.45rem; }
    .ranked-discovery-card p { color: #3b4652; font-size: 0.88rem; line-height: 1.38; margin: 0 0 0.55rem; }
    .ranked-discovery-card dl { margin: 0.55rem 0 0; display: grid; grid-template-columns: 92px 1fr; gap: 0.32rem 0.55rem; }
    .ranked-discovery-card dt { color: #59636e; font-size: 0.72rem; font-weight: 820; text-transform: uppercase; letter-spacing: 0.03em; }
    .ranked-discovery-card dd { color: #26313d; font-size: 0.82rem; margin: 0; overflow-wrap: anywhere; }
    .decision-prompts { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 0.55rem; margin-bottom: 0.85rem; }
    .decision-prompt { border: 1px solid #e1e6eb; background: #fbfcfd; border-radius: 8px; padding: 0.68rem 0.75rem; color: #26313d; min-width: 0; }
    .decision-prompt span { display: block; color: #59636e; font-size: 0.7rem; font-weight: 780; text-transform: uppercase; letter-spacing: 0.03em; margin-bottom: 0.2rem; }
    .decision-prompt strong { display: block; color: #12263f; font-size: 0.86rem; line-height: 1.25; }
    .route-summary { display: grid; grid-template-columns: minmax(110px, 0.28fr) 1fr; gap: 0.45rem 1rem; align-items: center; color: #26313d; }
    .route-summary-label { color: #59636e; font-size: 0.78rem; font-weight: 760; text-transform: uppercase; letter-spacing: 0.03em; }
    .route-summary-value { font-size: 0.95rem; line-height: 1.35; }
    .route-badge { display: inline-flex; align-items: center; border-radius: 999px; padding: 0.22rem 0.62rem; font-weight: 780; font-size: 0.86rem; }
    .route-badge.yes { background: #f1f8f2; color: #1b5e20; }
    .route-badge.partial, .route-badge.maybe { background: #fff8ed; color: #8a4b00; }
    .route-badge.no { background: #fbeff0; color: #8a1f2d; }
    .study-snapshot { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 0.8rem; }
    .study-snapshot-header { grid-column: 1 / -1; border-bottom: 1px solid #e5e9ee; padding-bottom: 0.7rem; margin-bottom: 0.1rem; }
    .study-snapshot-header strong { display: block; color: #12263f; font-size: 1.05rem; line-height: 1.25; }
    .study-snapshot-header span { display: block; color: #59636e; font-size: 0.86rem; margin-top: 0.25rem; }
    .study-snapshot-item { border: 1px solid #e1e6eb; background: #fff; border-radius: 8px; padding: 0.85rem; min-width: 0; }
    .study-snapshot-item span { display: block; color: #59636e; font-size: 0.74rem; font-weight: 760; text-transform: uppercase; letter-spacing: 0.03em; margin-bottom: 0.3rem; }
    .study-snapshot-item strong { display: block; color: #26313d; font-size: 0.96rem; line-height: 1.3; overflow-wrap: anywhere; }
    .study-snapshot-note { grid-column: 1 / -1; color: #4b5563; line-height: 1.4; margin: 0.25rem 0 0; }
    .global-sidebar-note { color: #59636e; font-size: 0.82rem; line-height: 1.25; margin-top: 0.35rem; }
    .global-sidebar-note strong { color: #26313d; }
    #download_figures_pdf { width: 100%; white-space: normal; line-height: 1.2; padding: 0.55rem 0.7rem; font-size: 0.92rem; }
    #dataset_id + .selectize-control .selectize-input { min-height: 2.35rem; align-items: center; }
    .assistant-tabset { background: #fff; border: 1px solid #d8dee4; border-radius: 8px; padding: 0.75rem; }
    .assistant-tabset .tab-content { padding-top: 1rem; }
    .import-tabset { background: #fff; border: 1px solid #d8dee4; border-radius: 8px; padding: 0.75rem; }
    .import-tabset .tab-content { padding-top: 1rem; }
    @media (max-width: 800px) {
      .study-stats, .study-section-grid { grid-template-columns: 1fr 1fr; }
      .study-snapshot { grid-template-columns: 1fr 1fr; }
      .ai-recommendation-grid { grid-template-columns: 1fr; }
      .ranked-discovery-grid { grid-template-columns: 1fr; }
      .story-panel { grid-template-columns: 1fr; }
      .story-steps { grid-template-columns: 1fr 1fr; }
      .decision-prompts { grid-template-columns: 1fr 1fr; }
      .assistant-stepper, .import-stepper { grid-template-columns: 1fr; }
      .study-definitions { grid-template-columns: 1fr; }
      .study-definitions dt { padding-bottom: 0.15rem; }
      .study-definitions dd { border-top: 0; padding-top: 0; }
    }
    @media (max-width: 520px) {
      .study-stats, .study-section-grid { grid-template-columns: 1fr; }
      .study-snapshot { grid-template-columns: 1fr; }
      .decision-prompts { grid-template-columns: 1fr; }
      .story-steps { grid-template-columns: 1fr; }
    }
  "))),
  sidebar = sidebar(
    passwordInput("openai_api_key", "OpenAI API key", value = "", placeholder = "Paste once for this session"),
    actionButton("save_openai_api_key", "Save key for session", icon = icon("key")),
    textOutput("openai_api_key_status"),
    textInput("openai_model", "LLM model", value = "gpt-4o-mini"),
    selectInput("dataset_id", "Active dataset", choices = dataset_choices, selected = selected_dataset),
    tags$div(class = "global-sidebar-note", tags$strong("Validate first."), " This selector switches active results after GEO reanalysis or Start-page upload."),
    downloadButton("download_figures_pdf", "Figures PDF"),
    tags$div(class = "global-sidebar-note", "Saved key/model are shared across this Shiny session.")
  ),

  nav_panel(
    "Start",
    layout_sidebar(
      sidebar = sidebar(
        radioButtons(
          "start_input_mode",
          "What do you want to analyze?",
          choices = c("GEO accession" = "geo", "Upload dataset" = "upload"),
          selected = "geo"
        ),
        conditionalPanel(
          "input.start_input_mode == 'geo'",
          textInput("start_geo_accession", "GEO accession", value = "", placeholder = "e.g. GSE19804"),
          actionButton("start_load_demo", "Load local demo", icon = icon("bolt"), class = "btn-outline-primary"),
          textAreaInput(
            "start_own_dataset_notes",
            "Relevant paper/method notes",
            value = "",
            rows = 5,
            placeholder = "Optional: paste abstracts, methods, biomarkers, pathway notes, or supplementary method text."
          ),
          actionButton("start_get_study_info", "Assess GEO dataset", icon = icon("circle-info"))
        ),
        conditionalPanel(
          "input.start_input_mode == 'upload'",
          fileInput("upload_expr", "Expression matrix CSV/TSV", accept = c(".csv", ".tsv", ".txt")),
          fileInput("upload_meta", "Metadata CSV/TSV", accept = c(".csv", ".tsv", ".txt")),
          textInput("upload_dataset_id", "Dataset ID", value = "uploaded_demo"),
          textInput("upload_display_name", "Display name", value = "Uploaded demo dataset"),
          uiOutput("upload_mapping_controls"),
          actionButton("upload_create_dataset", "Create app dataset", icon = icon("file-import"))
        ),
        helpText("Start with GEO for paper-aware reanalysis, or upload an expression matrix and metadata for direct exploration.")
      ),
      tags$div(
        class = "assistant-workspace",
        workflow_story_panel("start"),
        uiOutput("ai_status_panel"),
        card(
          class = "ai-guide-card",
          card_header(info_title("SeqSurfer Recommendation", "AI-guided interpretation of whether this dataset is worth reanalyzing and what to do next.")),
          tags$div(
            class = "ai-guide-body",
            tags$div(
              class = "decision-prompts",
              tags$div(class = "decision-prompt", tags$span("Value"), tags$strong("Why reanalyze this dataset?")),
              tags$div(class = "decision-prompt", tags$span("Evidence"), tags$strong("What has already been shown?")),
              tags$div(class = "decision-prompt", tags$span("Hypotheses"), tags$strong("What should we test next?")),
              tags$div(class = "decision-prompt", tags$span("Risk"), tags$strong("What must be checked first?"))
            ),
            uiOutput("start_reanalysis_narrative"),
            tags$div(
              class = "ai-guide-actions",
              actionButton("start_generate_narrative", "Generate AI recommendation", icon = icon("wand-magic-sparkles")),
              uiOutput("start_next_step_buttons")
            )
          )
        ),
        card(card_header(info_title("Data Quality Signals", "RIN score detection, sample count, expression matrix availability, and PCA-based outlier screening to inform whether this dataset is suitable for reanalysis.")), uiOutput("data_quality_signals_panel")),
        card(
          card_header(info_title("Can SeqSurf Analyze This In-App?", "First-pass route classifier: processed-matrix in-app analysis, count-matrix path, raw-data handoff, or manual curation.")),
          uiOutput("geo_analyzability_route_summary")
        ),
        card(
          card_header(info_title("Dataset Snapshot", "Compact description of the imported study: organism, assay type, sample count, metadata, expression files, and paper context.")),
          uiOutput("start_study_snapshot")
        ),
        tags$div(
          class = "assistant-hint",
          tags$strong("Status: "),
          uiOutput("start_study_status", inline = TRUE)
        ),
        tags$div(
          class = "assistant-tabset",
          navset_tab(
            nav_panel(
              "Evidence",
              layout_columns(
                card(card_header(info_title("Reanalysis Assessment", "Scores for data availability, original-method reproducibility, reanalysis usefulness, novelty opportunity, technical risk, and proceed recommendation.")), DTOutput("start_reanalysis_assessment_table")),
                card(card_header(info_title("Recommended Reanalysis Targets", "Suggested paper-style reproduction checks and modern follow-up analyses derived from GEO/PubMed context and optional paper notes.")), DTOutput("start_reanalysis_recommendations_table"))
              ),
              card(card_header(info_title("Local Learning Prior", "How previous local dataset feedback changes the recommendation for this study route, assay, and analysis pattern.")), DTOutput("start_learning_prior_table")),
              card(card_header(info_title("What Has Been Done Before", "Compact prior-work report from GEO, PubMed, extracted methods/results, dataset reuse clues, and limitations.")), DTOutput("prior_work_report_table"))
            ),
            nav_panel(
              "Analyzability Details",
              layout_columns(
                card(card_header(info_title("Analyzability Checks", "Signals used to classify processed matrices, count matrices, raw files, supplementary links, method detail, and sample-group readiness.")), DTOutput("geo_analyzability_checks_table"))
              ),
              card(card_header(info_title("Raw And Supplementary File Triage", "Practical triage for series matrices, count matrices, raw sequencing files, supplementary archives, and pipeline handoff.")), DTOutput("raw_supplementary_triage_table"))
            ),
            nav_panel(
              "Plan",
              card(card_header(info_title("Exact-Method Reproduction Planner", "Feasibility assessment distinguishing exact reproduction, paper-style approximation, and cases where public files are insufficient.")), DTOutput("exact_reproduction_table")),
              layout_columns(
                card(card_header(info_title("Count-Matrix Automation Plan", "How SeqSurf handles uploaded count matrices and when external DESeq2/edgeR reproduction is needed.")), verbatimTextOutput("count_matrix_pipeline_plan")),
                card(card_header(info_title("FASTQ/SRA External Handoff", "External workflow specification for raw sequencing studies that should not run inside interactive Shiny.")), verbatimTextOutput("raw_sra_handoff_plan"))
              ),
              card(
                card_header(info_title("No-Code AI Handoff", "Plain-language instructions for users who need an AI coding assistant to run external dataset processing and return files to SeqSurf.")),
                uiOutput("no_code_ai_handoff_panel")
              ),
              card(
                card_header(info_title("Runnable Workflow Folder", "Create a local project folder with scripts, preflight checks, prompts, and a double-click starter for Positron/Codex/Claude handoff.")),
                tags$div(
                  class = "ai-guide-body",
                  actionButton("create_external_workspace", "Create workflow folder", icon = icon("folder-plus")),
                  uiOutput("external_workspace_status")
                )
              ),
              card(
                card_header(info_title("Prompt For AI Assistant", "Copy or download this prompt for Claude Code, Positron Assistant, Codex, or another coding assistant that can run R and terminal commands.")),
                verbatimTextOutput("external_ai_analysis_prompt")
              ),
              card(
                card_header(info_title("Next Actions", "Where to go after study lookup.")),
                tags$ol(
                  tags$li("Prepare Evidence: collect study memory from GEO, linked papers, related papers, and any added notes."),
                  tags$li("Use Validate to reproduce the paper-style analysis and score agreement with paper claims."),
                  tags$li("Use Discovery to identify better follow-up analyses, interpret surprising results, and inspect PCA, DEG, volcano, heatmap, GSEA, enrichment, and gene-search views.")
                )
              ),
              layout_columns(
                downloadButton("download_count_matrix_plan", "Download count-matrix plan"),
                downloadButton("download_raw_sra_plan", "Download FASTQ/SRA handoff"),
                downloadButton("download_raw_sra_bundle", "Download raw/SRA workflow bundle"),
                downloadButton("download_no_code_ai_guide", "Download no-code AI guide"),
                downloadButton("download_external_ai_prompt", "Download AI assistant prompt"),
                downloadButton("download_demo_walkthrough", "Download GSE19804 demo walkthrough")
              ),
              card(card_header(info_title("Selected Dataset Study Context", "Study-level background, cohort description, and analysis context for the dataset selected in the global sidebar.")), uiOutput("study_information"))
            ),
            nav_panel(
              "Upload Preview",
              card(card_header(info_title("Upload Status", "Progress and validation messages for user-uploaded expression and metadata files.")), verbatimTextOutput("upload_status")),
              layout_columns(
                card(card_header(info_title("Expression Preview", "First rows and columns from the uploaded expression matrix.")), DTOutput("upload_expr_preview")),
                card(card_header(info_title("Metadata Preview", "First rows and columns from the uploaded metadata table.")), DTOutput("upload_meta_preview"))
              )
            ),
            nav_panel("Count Matrix Import",
              card(
                card_header(info_title("Import Count Matrix", "Drop a raw count matrix (genes x samples) and a metadata table to create an app-ready dataset contract without re-downloading from GEO.")),
                fileInput("count_matrix_file", "Count matrix (CSV, genes as rows)", accept = ".csv"),
                fileInput("count_metadata_file", "Sample metadata (CSV)", accept = ".csv"),
                textInput("count_accession_label", "Dataset label (e.g. GSE12345 or my_study)", placeholder = "Used as the dataset ID"),
                helpText("The count matrix should have a gene symbol or ID column as the first column, then one column per sample. The metadata should have one row per sample with a column matching the count matrix column names."),
                actionButton("import_count_matrix", "Import and create dataset", icon = icon("upload"), class = "btn-primary"),
                uiOutput("count_matrix_import_status")
              )
            )
          )
        )
      )
    )
  ),

  nav_panel(
    "Prepare Evidence",
    layout_sidebar(
      sidebar = sidebar(
        textInput("assistant_geo", "GEO accession", value = "", placeholder = "Carried from Start, or enter e.g. GSE19804"),
        textAreaInput(
          "assistant_extra_papers",
          "Additional papers to add to the study memory",
          value = "",
          rows = 8,
          placeholder = "The app automatically searches PubMed for papers mentioning this GEO accession, then you can paste extra abstracts, methods, key claims, biomarkers, pathway notes, or supplementary method text here."
        ),
        uiOutput("assistant_action_buttons"),
        uiOutput("assistant_kb_picker"),
        downloadButton("download_paper_to_code", "Export paper-to-code brief"),
        helpText("Prepare Evidence first. The same study memory is used by Validate and Discovery.")
      ),
      tags$div(
        class = "assistant-workspace",
        workflow_story_panel("train"),
        uiOutput("assistant_workflow_banner"),
        tags$div(
          class = "assistant-tabset",
          navset_tab(
            id = "train_results_tab",
            nav_panel(
              "Overview",
              card(card_header(info_title("AI Readiness", "Optional dependencies and API key status for GEO/PubMed retrieval, PDF extraction, LLM extraction/chat, survival, and immune deconvolution.")), DTOutput("assistant_setup_table")),
              card(card_header(info_title("Study Memory 2.0 Virtual Lab", "Five specialist agent reviews: literature curator, methods extractor, claims extractor, reproducibility judge, and discovery scientist.")), DTOutput("virtual_lab_summary_table")),
              card(card_header(info_title("Study Knowledge Base", "Paper-grounded extracted study objective, cohort details, methods, findings, biomarkers, pathways, and limitations for the selected GEO accession.")), DTOutput("assistant_kb_table"))
            ),
            nav_panel(
              "Evidence",
              card(card_header(info_title("Knowledge Sources", "Original publication, supplementary files, PMC full text, and indexed local text sources discovered or downloaded for this study.")), DTOutput("assistant_documents_table")),
              card(card_header(info_title("PaperQA Corpus", "Prepared text corpus and manifest for optional citation-grounded PaperQA indexing/chat outside the core Shiny app.")), DTOutput("paperqa_corpus_table")),
              card(card_header(info_title("Related Papers Mentioning GEO", "PubMed papers automatically discovered by searching for the GEO accession. These abstracts are included in the study knowledge base when available.")), DTOutput("assistant_related_papers_table")),
              card(card_header(info_title("Paper Sections", "Section-aware excerpts from downloaded or pasted paper context, including methods, results, discussion, limitations, and supplementary/data-availability text.")), DTOutput("assistant_sections_table")),
              card(card_header(info_title("Claim Candidates", "Rule-detected sentences from paper sections that look like methods/results claims and can be promoted into structured claim extraction.")), DTOutput("assistant_claim_candidates_table")),
              card(card_header(info_title("Structured Publication Claims", "LLM- or rule-extracted paper claims with entity, comparison, direction, significance, evidence sentence, and confidence.")), DTOutput("assistant_claims_table"))
            ),
            nav_panel(
              "Virtual Lab",
              card(card_header(info_title("Agent Reviews", "SeqSurf Study Memory 2.0 agent summaries and evidence-grounded uncertainty notes.")), DTOutput("virtual_lab_summary_table_full")),
              card(card_header(info_title("Agent Recommendations", "Role-specific recommendations from the virtual lab, used to guide Validate and Discovery.")), DTOutput("virtual_lab_recommendation_table")),
              card(card_header(info_title("Citation Evidence", "Claim sentences, agent-cited evidence, and local PaperQA corpus sources used to justify the study memory.")), DTOutput("citation_evidence_table"))
            ),
            nav_panel(
              "Report Brief",
              card(card_header(info_title("Paper-To-Code Brief", "A coding-agent handoff that turns the trained paper knowledge into reproducible analysis-script requirements for this app's dataset contract.")), verbatimTextOutput("paper_to_code_brief"))
            )
          )
        )
      )
    )
  ),

  nav_panel(
    "Validate",
    layout_sidebar(
      sidebar = sidebar(
        textInput("geo_import_accession", "GEO accession", value = "", placeholder = "Carried from Prepare Evidence"),
        uiOutput("geo_import_action_buttons"),
        uiOutput("geo_design_editor"),
        sliderInput("publication_padj", "Publication-check adjusted p-value cutoff", min = 0, max = 0.25, value = 0.1, step = 0.005),
        numericInput("publication_lfc", "Publication-check absolute log2FC cutoff", value = 0.5, min = 0, step = 0.25),
        helpText("Run reanalysis, then validate paper agreement from the same tab.")
      ),
      tags$div(
        class = "import-workspace",
        workflow_story_panel("validate"),
        uiOutput("geo_import_workflow_banner"),
        uiOutput("reproducibility_score_card"),
        tags$div(
          class = "import-tabset",
          navset_tab(
            id = "geo_import_tabset",
            nav_panel(
              "Status",
              card(
                card_header(info_title("Validate Workflow Status", "Current status for GEO fetching, paper training, design proposal, reanalysis, dataset creation, and report generation.")),
                verbatimTextOutput("geo_import_status")
              )
            ),
            nav_panel(
              "Design",
              card(
                card_header(info_title("Proposed Analysis Design", "AI- or rule-proposed sample grouping, contrast, covariates, platform type, and rationale derived from GEO metadata and paper methods.")),
                DTOutput("geo_design_table")
              )
            ),
            nav_panel(
              "Metadata",
              card(
                card_header(info_title("Phenotype Metadata Preview", "GEO sample-level phenotype columns used to choose sample groups, primary contrast, covariates, and sample identifiers.")),
                DTOutput("geo_pheno_preview")
              )
            ),
            nav_panel(
              "Count Files",
              card(
                card_header(info_title("Supplementary Count-Matrix Candidates", "GEO supplementary files classified as likely count matrices, possible rectangular tables, raw files, or non-count archives.")),
                DTOutput("geo_count_candidate_table"),
                downloadButton("download_count_pipeline_script", "Download generated count pipeline script")
              )
            ),
            nav_panel(
              "Report",
              card(
                card_header(info_title("Generated Reproducibility Report", "HTML report comparing the imported dataset's results against paper-extracted claims, plus gap analysis and suggested next analyses.")),
                uiOutput("geo_report_link"),
                downloadButton("download_geo_report", "Download reproducibility report")
              )
            ),
            nav_panel(
              "Paper Agreement",
              uiOutput("dataset_match_warning"),
              card(card_header(info_title("Dataset Verification", "Checks whether the active app dataset is traceable to the trained GEO record and flags sample-count or provenance differences for review.")), DTOutput("assistant_verification_table")),
              card(card_header(info_title("Compare To Original Publication", "Reproducibility report comparing reported biomarkers and pathways from the trained paper knowledge base against the active DEG and GSEA results.")), uiOutput("publication_comparison_table"))
            ),
            nav_panel(
              "Score",
              layout_columns(
                card(card_header(info_title("Post-Analysis Reproducibility Score", "Weighted score summarizing supported, weak, missing, and direction-mismatched publication claims after reanalysis.")), DTOutput("post_reproducibility_score_table")),
                card(card_header(info_title("Relevant Or Different Results", "Triage of reproduced paper findings, findings that diverge from the paper, weak/missing claims, and high-signal app results not found in extracted paper claims.")), DTOutput("relevant_different_results_table"))
              )
            )
          )
        )
      )
    )
  ),

  nav_panel(
    "Discovery",
    tags$div(
      class = "assistant-workspace",
      workflow_story_panel("discovery"),
      card(
        card_header(info_title("Ranked Next Analyses", "Deterministic prioritized recommendations using the active dataset, DEG results, GSEA results, metadata, and paper-agreement status.")),
        uiOutput("ranked_discovery_cards")
      ),
      card(
        class = "ai-guide-card",
        card_header(info_title("SeqSurfer's Next-Step Recommendation", "AI-guided suggestions for better reanalysis strategies and follow-up questions.")),
        tags$div(
          class = "ai-guide-body",
          tags$p("Start here after validation. SeqSurfer uses the trained study memory, active dataset, DEG results, GSEA results, and paper-agreement score to recommend feasible next analyses."),
          textAreaInput("assistant_question", NULL, value = "Based only on the available metadata and app results, what are the three most useful next analyses, which app tabs should I open, and what would each analysis prove or disprove?", rows = 4),
          tags$div(
            class = "ai-guide-actions",
            actionButton("ask_assistant", "Ask Discovery", icon = icon("comment-dots"))
          ),
          uiOutput("assistant_chat_response")
        )
      ),
      card(card_header(info_title("Better Reanalysis Ideas", "Suggested modern or higher-value follow-up analyses based on the active dataset, trained study memory, DEG table, and GSEA results.")), DTOutput("assistant_suggestions_table")),
      layout_columns(
        card(card_header(info_title("Suggested New Comparisons", "Metadata-derived group comparisons that are not clearly represented in existing DEG contrasts and may be worth testing next.")), DTOutput("comparison_suggestions_table")),
        card(card_header(info_title("Gap Analysis", "Analyses available in the app plus unexplored methods that were not clearly detected in the trained publication knowledge.")), DTOutput("gap_analysis_table"))
      ),
      card(card_header(info_title("Standardized Reanalysis Checklist", "Protocol-style checkpoints for making public sequencing reanalysis comparable across datasets.")), DTOutput("standardized_pipeline_table")),
      layout_columns(
        card(card_header(info_title("Common Sequencing Pitfalls", "Route-, metadata-, method-, and learning-aware risks that should be checked before drawing conclusions.")), DTOutput("sequencing_pitfall_table")),
        card(card_header(info_title("Epithelial/EMT Signature Screen", "A fixed first signature layer for standardized interpretation across datasets, including epithelial and mesenchymal coverage.")), DTOutput("signature_screen_table"))
      ),
      card(
        card_header(info_title("Teach SeqSurf", "Save whether this dataset/workflow was useful so the app can build a local learning memory across datasets.")),
        tags$div(
          class = "ai-guide-body",
          layout_columns(
            selectInput(
              "learning_outcome",
              "Outcome",
              choices = c(
                "Useful reanalysis candidate",
                "Reproduced paper well",
                "Interesting disagreement",
                "Needed manual curation",
                "Low-value dataset",
                "Workflow failed"
              ),
              selected = "Useful reanalysis candidate"
            ),
            sliderInput("learning_rating", "Usefulness rating", min = 1, max = 5, value = 4, step = 1)
          ),
          textAreaInput("learning_notes", "What should SeqSurf remember?", value = "", rows = 3, placeholder = "Example: count-matrix route worked, but metadata groups needed manual cleanup; MYCN pathway disagreement was scientifically interesting."),
          actionButton("save_learning_event", "Save learning event", icon = icon("graduation-cap")),
          uiOutput("learning_event_status")
        )
      ),
      layout_columns(
        card(card_header(info_title("Learning Summary", "Local aggregate memory from user feedback and dataset outcomes.")), DTOutput("learning_summary_table")),
        card(card_header(info_title("Learned Lessons", "Route-level patterns SeqSurf can use to guide future recommendations.")), DTOutput("learning_lessons_table"))
      ),
      card(
        card_header(info_title("Shared Learning With Consent", "Export anonymized aggregate learning or import a shared packet from another SeqSurf install. No GEO accession, dataset ID, or free-text notes are shared.")),
        tags$div(
          class = "ai-guide-body",
          tags$p("Shared learning is opt-in. Exported packets contain route, assay signal, sample-size bin, method tags, outcome label, usefulness/reproducibility scores, and claim-count summaries."),
          checkboxInput("shared_learning_consent", "I understand this exports anonymized aggregate learning only, not raw data or paper text.", value = FALSE),
          tags$div(
            class = "ai-guide-actions",
            downloadButton("download_shared_learning_packet", "Export shared learning packet"),
            fileInput("shared_learning_import", "Import shared learning packet", accept = c(".zip", ".csv")),
            actionButton("import_shared_learning", "Import packet", icon = icon("file-import"))
          ),
          uiOutput("shared_learning_status")
        )
      ),
      card(
        card_header(info_title("Training Data Export", "Export structured learning events — analysis contracts, recommendations, warnings, reproducibility scores, and signature plans — for review or use in retrieval-augmented workflows.")),
        tags$div(
          class = "ai-guide-body",
          checkboxInput("training_include_shared", "Include imported shared learning events", value = TRUE),
          checkboxInput("training_include_private", "Include private study context excerpts for local-only model development", value = FALSE),
          DTOutput("training_examples_summary_table"),
          downloadButton("download_training_examples", "Export training examples")
        )
      ),
      layout_columns(
        card(
          card_header(info_title("Survival Screen", "Automatic survival check when the selected dataset contains usable time/event metadata.")),
          actionButton("run_survival_screen", "Run survival screen", icon = icon("chart-line")),
          DTOutput("survival_screen_table")
        ),
        card(
          card_header(info_title("Immune Deconvolution", "Automatic immune-deconvolution check when expression format and optional packages support it.")),
          actionButton("run_immune_screen", "Run immune deconvolution", icon = icon("users")),
          DTOutput("immune_screen_table")
        )
      )
    )
  ),

  nav_panel(
    "Compare Studies",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("compare_dataset_a", "Dataset A", choices = dataset_choices),
        selectInput("compare_dataset_b", "Dataset B", choices = dataset_choices),
        actionButton("run_comparison", "Compare datasets", icon = icon("code-compare"), class = "btn-primary"),
        helpText("Select two datasets that have completed reanalysis to compare DEG overlap, pathway agreement, and claimed biomarkers.")
      ),
      card(
        card_header(info_title("Comparison Snapshot", "High-level counts for significant DEGs, overlap, pathway direction agreement, and claimed biomarkers.")),
        uiOutput("compare_summary_cards")
      ),
      layout_columns(
        card(card_header(info_title("DEG Overlap Direction", "Shared significant genes plotted by log2 fold-change in each dataset. Points in matching quadrants move in the same direction.")), plotOutput("compare_deg_overlap_plot", height = "360px")),
        card(card_header(info_title("Pathway Agreement", "Shared Hallmark pathways plotted by normalized enrichment score in each dataset. Concordant pathways fall in matching-sign quadrants.")), plotOutput("compare_pathway_plot", height = "360px"))
      ),
      card(card_header(info_title("Claimed Biomarker Presence", "Clean claimed biomarkers from either study knowledge base, showing whether each marker appears in each DEG result.")), plotOutput("compare_biomarker_plot", height = "300px"))
    )
  ),

  nav_panel(
    "PCA Explorer",
    layout_sidebar(
      sidebar = sidebar(
        uiOutput("pca_controls")
      ),
      card(card_header(info_title("PCA", "Principal component analysis summarizes the largest expression differences among samples. Nearby samples have more similar global expression profiles; axis percentages show variance explained.")), plotOutput("pca_plot", height = "620px"))
    )
  ),

  nav_panel(
    "Top Loading Genes",
    layout_sidebar(
      sidebar = sidebar(
        uiOutput("loading_controls"),
        downloadButton("download_loadings", "Download table")
      ),
      card(card_header(info_title("Top PCA Loading Genes", "PCA loadings show which genes contribute most strongly to a principal component. Positive and negative loadings drive samples in opposite directions along that axis.")), plotOutput("loadings_plot", height = "620px")),
      card(card_header(info_title("Loading Table", "Downloadable gene-level loading values for the selected principal component and direction.")), DTOutput("loadings_table"))
    )
  ),

  nav_panel(
    "DEG Explorer",
    layout_sidebar(
      sidebar = sidebar(
        uiOutput("deg_controls"),
        downloadButton("download_deg", "Download filtered DEGs")
      ),
      layout_columns(
        value_box("Up", textOutput("deg_up_count"), showcase = NULL),
        value_box("Down", textOutput("deg_down_count"), showcase = NULL),
        value_box("Tested genes", textOutput("deg_total_count"), showcase = NULL)
      ),
      card(card_header(info_title("Differential Expression Table", "Gene-level effect sizes and statistical evidence for the selected contrast. Direction is assigned using the selected adjusted p-value and absolute log2 fold-change cutoffs.")), DTOutput("deg_table"))
    )
  ),

  nav_panel(
    "Volcano Plot",
    layout_sidebar(
      sidebar = sidebar(uiOutput("volcano_controls")),
      card(card_header(info_title("Volcano Plot", "Each point is a gene. The x-axis shows log2 fold change and the y-axis shows statistical significance. Hover over a point for gene-level details.")), uiOutput("volcano_ui"))
    )
  ),

  nav_panel(
    "Heatmap",
    layout_sidebar(
      sidebar = sidebar(uiOutput("heatmap_controls")),
      card(card_header(info_title("Expression Heatmap", "Rows are genes and columns are samples. Values are scaled by gene to emphasize relative high and low expression patterns across the cohort.")), plotOutput("heatmap_plot", height = "720px"))
    )
  ),

  nav_panel(
    "GSEA Explorer",
    layout_sidebar(
      sidebar = sidebar(
        uiOutput("gsea_controls"),
        sliderInput("pathway_padj", "GSEA adjusted p-value cutoff", min = 0, max = 0.25, value = 0.1, step = 0.005),
        downloadButton("download_gsea", "Download pathway table"),
        helpText("Positive NES: enriched toward the positive group. Negative NES: enriched toward the negative group. Hover over bars for pathway size and leading-edge overlap.")
      ),
      card(card_header(info_title("GSEA Normalized Enrichment Scores", "Precomputed preranked GSEA across Hallmark, Reactome, GO Biological Process, or WikiPathways. NES gives enrichment direction and strength; leading-edge genes drive the signal.")), plotly::plotlyOutput("pathway_plot", height = "720px")),
      card(card_header(info_title("GSEA Table", "Pathway-level enrichment statistics, collection and MSigDB version, pathway size, and leading-edge genes for the selected contrast.")), DTOutput("gsea_table"))
    )
  ),

  nav_panel(
    "Enrichr",
    layout_sidebar(
      sidebar = sidebar(
        uiOutput("enrichr_controls"),
        actionButton("run_enrichr", "Run Enrichr", icon = icon("play")),
        downloadButton("download_enrichr", "Download results")
      ),
      card(card_header(info_title("Enrichr Results", "Over-representation analysis of the selected DEG list against an Enrichr library. This analysis sends gene symbols to the public Enrichr service.")), plotOutput("enrichr_plot", height = "620px")),
      card(card_header(info_title("Enrichr Table", "Terms are ranked by adjusted p-value. Overlap reports how many submitted genes occur in each annotated gene set.")), DTOutput("enrichr_table"))
    )
  ),

  nav_panel(
    "clusterProfiler",
    layout_sidebar(
      sidebar = sidebar(
        uiOutput("clusterprofiler_controls"),
        actionButton("run_clusterprofiler", "Run GO enrichment", icon = icon("play")),
        downloadButton("download_clusterprofiler", "Download results")
      ),
      card(card_header(info_title("GO Enrichment", "clusterProfiler tests whether selected DEGs are over-represented in Gene Ontology terms using all tested genes as the background universe.")), plotOutput("clusterprofiler_plot", height = "680px")),
      card(card_header(info_title("clusterProfiler Table", "Gene Ontology terms, enrichment ratios, adjusted p-values, and contributing genes.")), DTOutput("clusterprofiler_table"))
    )
  ),

  nav_panel(
    "Gene Search",
    layout_sidebar(
      sidebar = sidebar(
        textInput("gene_query", "Gene", value = "OXTR"),
        uiOutput("gene_group_control")
      ),
      card(card_header(info_title("DEG Stats", "Differential-expression statistics for the searched gene across available contrasts.")), DTOutput("gene_deg_table")),
      card(card_header(info_title("Expression by Group", "Normalized expression of the searched gene, summarized and displayed for each selected metadata group.")), plotOutput("gene_expression_plot", height = "440px")),
      card(card_header(info_title("Leading-edge Pathways", "GSEA pathways whose leading-edge gene subset contains the searched gene.")), DTOutput("gene_pathway_table"))
    )
  )
)

server <- function(input, output, session) {
  dataset <- reactive({
    validate(need(!is.null(input$dataset_id) && nzchar(input$dataset_id), "No active dataset selected. Import/create a dataset in Validate or choose saved results from the Active dataset menu."))
    load_rnaseq_dataset(input$dataset_id, "data")
  })

  session_openai_api_key <- reactiveVal(Sys.getenv("OPENAI_API_KEY"))

  active_openai_api_key <- reactive({
    first_nonempty(session_openai_api_key(), input$openai_api_key, "")
  })

  observeEvent(input$save_openai_api_key, {
    key <- trimws(first_nonempty(input$openai_api_key, ""))
    if (!nzchar(key)) {
      showNotification("Paste an OpenAI API key before saving it for this session.", type = "warning")
      return()
    }
    session_openai_api_key(key)
    updateTextInput(session, "openai_api_key", value = "")
    showNotification("OpenAI API key saved for this Shiny session.", type = "message")
  }, ignoreInit = TRUE)

  output$openai_api_key_status <- renderText({
    if (nzchar(session_openai_api_key())) {
      "API key: saved for this session"
    } else if (nzchar(trimws(first_nonempty(input$openai_api_key, "")))) {
      "API key: pasted, click Save key for session"
    } else {
      "API key: not saved"
    }
  })

  dataset_matches_accession <- function(dataset_id, accession, data_root = "data") {
    accession <- normalize_geo_accession(accession)
    if (is.null(dataset_id) || !nzchar(dataset_id) || !nzchar(accession)) return(FALSE)
    if (identical(dataset_id, geo_dataset_id(accession))) return(TRUE)
    info_path <- file.path(data_root, dataset_id, "dataset_info.rds")
    if (!file.exists(info_path)) return(FALSE)
    info <- readRDS(info_path)
    haystack <- tolower(paste(
      dataset_id,
      first_nonempty(info$dataset_id, ""),
      first_nonempty(info$display_name, ""),
      first_nonempty(info$notes, ""),
      collapse = " "
    ))
    grepl(tolower(accession), haystack, fixed = TRUE)
  }

  find_dataset_for_accession <- function(accession, data_root = "data") {
    ids <- available_datasets(data_root)
    matches <- ids[vapply(ids, dataset_matches_accession, logical(1), accession = accession, data_root = data_root)]
    if (length(matches) == 0) return(NULL)
    matches[[1]]
  }

  select_dataset_for_accession <- function(accession, quiet = FALSE) {
    match <- find_dataset_for_accession(accession, "data")
    if (!is.null(match)) {
      updateSelectInput(session, "dataset_id", choices = dataset_choices, selected = match)
      choice_label <- names(dataset_choices)[match(match, unname(dataset_choices))]
      choice_label <- first_nonempty(choice_label, match)
      if (!quiet) showNotification(paste("Selected matching dataset:", choice_label), type = "message", duration = 5)
      return(match)
    }
    if (!quiet) {
      showNotification(
        paste("No matching app dataset found for", normalize_geo_accession(accession), "yet. Use Validate first, or select a matching dataset manually."),
        type = "warning",
        duration = 8
      )
    }
    NULL
  }

  assistant_kbs <- reactiveVal(available_study_kbs("data"))
  kb_version <- reactiveVal(0)
  start_study_status <- reactiveVal("Enter a GEO accession or custom dataset description, then click Get study info.")
  start_study_info <- reactiveVal(data.frame(Field = character(), Value = character()))
  start_study_kb <- reactiveVal(NULL)
  start_study_narrative <- reactiveVal("Enter a GEO accession and click Assess GEO dataset. SeqSurfer will score whether the dataset is worth reanalyzing, explain what has already been done, and recommend the next step.")
  ai_status <- reactiveVal(list(status = "Ready", detail = "Enter a GEO accession or upload a dataset to begin."))
  geo_expression_set_status <- reactiveVal(NULL)
  geo_import_status <- reactiveVal("Ready. Enter a GEO accession to begin.")
  geo_design <- reactiveVal(NULL)
  geo_eset <- reactiveVal(NULL)
  geo_report <- reactiveVal(NULL)
  geo_imported_dataset <- reactiveVal(NULL)
  geo_count_candidates <- reactiveVal(data.frame(File = character(), URL = character(), Score = numeric(), Candidate_type = character(), Status = character(), check.names = FALSE))
  upload_status <- reactiveVal("Upload expression and metadata files to begin.")
  external_workspace_info <- reactiveVal(NULL)
  learning_event_status <- reactiveVal("Save feedback after validation to help SeqSurf learn which routes, scores, and findings were actually useful.")
  shared_learning_status <- reactiveVal("No shared learning packet imported this session.")

  uploaded_expr_table <- reactive({
    req(input$upload_expr)
    read_uploaded_table(input$upload_expr$datapath)
  })

  uploaded_meta_table <- reactive({
    req(input$upload_meta)
    read_uploaded_table(input$upload_meta$datapath)
  })

  output$assistant_kb_picker <- renderUI({
    choices <- assistant_kbs()
    selectInput(
      "assistant_kb",
      "Saved study knowledge base",
      choices = c("Choose after entering a GEO accession" = "", stats::setNames(choices, choices)),
      selected = ""
    )
  })

  assistant_kb_choices <- function() {
    choices <- assistant_kbs()
    c("Choose after entering a GEO accession" = "", stats::setNames(choices, choices))
  }

  select_assistant_kb <- function(accession) {
    updateSelectInput(
      session,
      "assistant_kb",
      choices = assistant_kb_choices(),
      selected = normalize_geo_accession(accession)
    )
  }

  output$ai_status_panel <- renderUI({
    state <- ai_status()
    ai_status_badge(state$status, state$detail)
  })

  output$upload_mapping_controls <- renderUI({
    req(input$upload_expr, input$upload_meta)
    expr <- uploaded_expr_table()
    meta <- uploaded_meta_table()
    meta_cols <- names(meta)
    group_cols <- meta_cols[vapply(meta, function(x) {
      n <- length(unique(stats::na.omit(as.character(x))))
      n >= 2 && n <= max(12, floor(nrow(meta) * 0.75))
    }, logical(1))]
    if (length(group_cols) == 0) group_cols <- meta_cols
    tagList(
      selectInput("upload_gene_col", "Gene column", choices = names(expr), selected = names(expr)[[1]]),
      selectInput("upload_sample_col", "Metadata sample ID column", choices = meta_cols, selected = meta_cols[[1]]),
      selectInput("upload_group_col", "Group column", choices = group_cols, selected = group_cols[[1]]),
      uiOutput("upload_group_levels"),
      selectInput("upload_patient_col", "Patient column", choices = c("Use sample IDs" = "", meta_cols), selected = ""),
      selectInput("upload_display_sample_col", "Display sample column", choices = c("Use sample IDs" = "", meta_cols), selected = ""),
      textInput("upload_contrast_label", "Contrast label", value = "")
    )
  })

  output$upload_group_levels <- renderUI({
    req(input$upload_meta, input$upload_group_col)
    meta <- uploaded_meta_table()
    levels <- names(sort(table(as.character(meta[[input$upload_group_col]]), useNA = "no"), decreasing = TRUE))
    validate(need(length(levels) >= 2, "Selected group column needs at least two levels."))
    tagList(
      selectInput("upload_negative_group", "Negative/reference group", choices = levels, selected = levels[[1]]),
      selectInput("upload_positive_group", "Positive/comparison group", choices = levels, selected = levels[[2]])
    )
  })

  output$upload_expr_preview <- renderDT({
    req(input$upload_expr)
    datatable(utils::head(uploaded_expr_table(), 20), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })

  output$upload_meta_preview <- renderDT({
    req(input$upload_meta)
    datatable(utils::head(uploaded_meta_table(), 20), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })

  output$upload_status <- renderText(upload_status())

  sync_geo_accession <- function(accession, source = "") {
    accession <- normalize_geo_accession(accession)
    if (!is_geo_accession(accession)) return(invisible(FALSE))
    if (!identical(source, "start") && !identical(normalize_geo_accession(input$start_geo_accession), accession)) {
      updateTextInput(session, "start_geo_accession", value = accession)
    }
    if (!identical(source, "train") && !identical(normalize_geo_accession(input$assistant_geo), accession)) {
      updateTextInput(session, "assistant_geo", value = accession)
    }
    if (!identical(source, "validate") && !identical(normalize_geo_accession(input$geo_import_accession), accession)) {
      updateTextInput(session, "geo_import_accession", value = accession)
    }
    invisible(TRUE)
  }

  observeEvent(input$start_geo_accession, {
    sync_geo_accession(input$start_geo_accession, "start")
  }, ignoreInit = TRUE)

  observeEvent(input$assistant_geo, {
    sync_geo_accession(input$assistant_geo, "train")
  }, ignoreInit = TRUE)

  observeEvent(input$geo_import_accession, {
    sync_geo_accession(input$geo_import_accession, "validate")
  }, ignoreInit = TRUE)

  output$start_next_step_buttons <- renderUI({
    kb <- start_study_kb()
    if (is.null(kb) || !is_geo_accession(kb$accession)) return(NULL)
    tagList(
      actionButton("start_go_train", "Continue to Prepare Evidence", icon = icon("brain"))
    )
  })

  observeEvent(input$start_go_train, {
    kb <- start_study_kb()
    req(kb)
    sync_geo_accession(kb$accession)
    select_assistant_kb(kb$accession)
    bslib::nav_select("main_nav", selected = "Prepare Evidence", session = session)
  }, ignoreInit = TRUE)

  observeEvent(input$start_go_validate, {
    kb <- start_study_kb()
    req(kb)
    sync_geo_accession(kb$accession)
    select_assistant_kb(kb$accession)
    bslib::nav_select("main_nav", selected = "Validate", session = session)
  }, ignoreInit = TRUE)

  observeEvent(input$start_go_discovery, {
    kb <- start_study_kb()
    req(kb)
    sync_geo_accession(kb$accession)
    select_assistant_kb(kb$accession)
    bslib::nav_select("main_nav", selected = "Discovery", session = session)
  }, ignoreInit = TRUE)

  observeEvent(input$upload_create_dataset, {
    req(input$upload_expr, input$upload_meta, input$upload_gene_col, input$upload_sample_col, input$upload_group_col, input$upload_positive_group, input$upload_negative_group)
    ai_status(list(status = "Creating uploaded dataset", detail = "Building PCA, DEG table, and app dataset contract from uploaded files."))
    upload_status("Creating app-ready dataset from uploaded files...")
    settings <- list(
      dataset_id = input$upload_dataset_id,
      display_name = input$upload_display_name,
      gene_col = input$upload_gene_col,
      sample_col = input$upload_sample_col,
      group_col = input$upload_group_col,
      positive_group = input$upload_positive_group,
      negative_group = input$upload_negative_group,
      patient_id_col = input$upload_patient_col,
      display_sample_col = input$upload_display_sample_col,
      contrast_label = input$upload_contrast_label
    )
    result <- tryCatch(
      write_uploaded_app_dataset(uploaded_expr_table(), uploaded_meta_table(), settings, "data"),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      upload_status(conditionMessage(result))
      ai_status(list(status = "Upload needs review", detail = conditionMessage(result)))
      showNotification(conditionMessage(result), type = "error", duration = 12)
      return()
    }
    datasets <<- get("available_datasets", mode = "function")("data")
    dataset_choices <<- active_dataset_choices(datasets)
    updateSelectInput(session, "dataset_id", choices = dataset_choices, selected = result$dataset_id)
    upload_status(paste("Created", result$dataset_id, "with", result$rows, "samples and", result$genes, "genes. It is now selected for exploration."))
    ai_status(list(status = "Uploaded dataset ready", detail = paste(result$dataset_id, "is available in the explorer tabs.")))
    showNotification("Uploaded dataset created and selected.", type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$import_count_matrix, {
    req(input$count_matrix_file, input$count_metadata_file, nzchar(input$count_accession_label))

    output$count_matrix_import_status <- renderUI(tags$p("Importing...", class = "global-sidebar-note"))

    tryCatch({
      counts_raw <- read.csv(input$count_matrix_file$datapath, row.names = 1, check.names = FALSE)
      meta_raw <- read.csv(input$count_metadata_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)

      label <- normalize_geo_accession(gsub("[^A-Za-z0-9_]", "_", input$count_accession_label))
      dataset_id <- tolower(paste0("cm_", label))
      dataset_dir <- file.path("data", dataset_id)
      dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

      # Detect gene column
      expr <- as.matrix(counts_raw)
      mode(expr) <- "numeric"
      expr <- expr[rowSums(is.na(expr)) < ncol(expr), , drop = FALSE]

      # Basic PCA
      var_genes <- head(order(apply(expr, 1, var, na.rm = TRUE), decreasing = TRUE), min(500, nrow(expr)))
      pca <- prcomp(t(expr[var_genes, , drop = FALSE]), scale. = TRUE)
      pca_df <- as.data.frame(pca$x[, 1:min(10, ncol(pca$x)), drop = FALSE])
      pca_df$sample_id <- rownames(pca_df)
      pca_variance <- summary(pca)$importance[2, 1:min(10, ncol(pca$x))]
      pca_loadings <- as.data.frame(pca$rotation[, 1:min(10, ncol(pca$x)), drop = FALSE])
      pca_loadings$gene <- rownames(pca_loadings)

      # Dataset info
      dataset_info <- list(
        dataset_id = dataset_id,
        label = input$count_accession_label,
        source = "count_matrix_import",
        expression_data_type = "raw_counts",
        sample_id_col = names(meta_raw)[1],
        group_col = if (ncol(meta_raw) > 1) names(meta_raw)[2] else names(meta_raw)[1],
        n_genes = nrow(expr),
        n_samples = ncol(expr),
        created = Sys.time()
      )

      saveRDS(as.data.frame(meta_raw), file.path(dataset_dir, "metadata.rds"))
      saveRDS(expr, file.path(dataset_dir, "vsd_matrix.rds"))
      saveRDS(pca_df, file.path(dataset_dir, "pca_df.rds"))
      saveRDS(pca_variance, file.path(dataset_dir, "pca_variance.rds"))
      saveRDS(pca_loadings, file.path(dataset_dir, "pca_loadings.rds"))
      saveRDS(dataset_info, file.path(dataset_dir, "dataset_info.rds"))
      write.csv(data.frame(gene = character(), log2FC = numeric(), pvalue = numeric(), padj = numeric()), file.path(dataset_dir, "deg_results.csv"), row.names = FALSE)
      write.csv(data.frame(pathway = character(), ES = numeric(), NES = numeric(), pval = numeric()), file.path(dataset_dir, "gsea_hallmark.csv"), row.names = FALSE)
      write.csv(data.frame(pathway = character(), score = numeric()), file.path(dataset_dir, "pathway_scores.csv"), row.names = FALSE)

      datasets <<- get("available_datasets", mode = "function")("data")
      dataset_choices <<- active_dataset_choices(datasets)
      updateSelectInput(session, "dataset_id", choices = dataset_choices, selected = dataset_id)

      output$count_matrix_import_status <- renderUI(
        tags$p(paste0("Dataset '", dataset_id, "' created with ", nrow(expr), " genes and ", ncol(expr), " samples. Select it in the global sidebar to analyze."), class = "global-sidebar-note")
      )
    }, error = function(e) {
      output$count_matrix_import_status <- renderUI(tags$p(paste("Import failed:", conditionMessage(e)), style = "color:red;"))
    })
  }, ignoreInit = TRUE)

  observeEvent(input$start_load_demo, {
    demo_script <- file.path("scripts", "create_compare_demo_datasets.R")
    if (file.exists(demo_script)) {
      tryCatch(sys.source(demo_script, envir = new.env(parent = globalenv())), error = function(e) {
        showNotification(paste("Demo fixture prep skipped:", conditionMessage(e)), type = "warning", duration = 8)
      })
    }

    datasets <<- get("available_datasets", mode = "function")("data")
    dataset_choices <<- active_dataset_choices(datasets)
    updateSelectInput(session, "dataset_id", choices = dataset_choices, selected = "geo_gse19804")
    updateSelectInput(session, "compare_dataset_a", choices = dataset_choices, selected = "geo_gse16476")
    updateSelectInput(session, "compare_dataset_b", choices = dataset_choices, selected = "geo_gse60450")
    updateTextInput(session, "start_geo_accession", value = "GSE19804")

    kb <- read_study_kb("GSE19804", "data")
    if (!is.null(kb)) {
      start_study_kb(kb)
      start_study_info(study_kb_summary_table(kb))
      geo_expression_set_status("available from local prepared demo dataset")
      geo_eset(NULL)
      geo_count_candidates(geo_supplementary_count_candidates(kb, "GSE19804"))
      sync_geo_accession("GSE19804", "start")
      select_assistant_kb("GSE19804")
      assessment <- reanalysis_assessment_table(kb)
      proceed <- assessment$Value[assessment$Score == "Proceed recommendation"]
      start_study_narrative("Local demo loaded. GSE19804 study memory and prepared app dataset are ready; use Validate, Discovery, Compare Studies, and PCA Explorer without live downloads.")
      start_study_status(paste("Local demo ready. GSE19804 loaded with proceed score", proceed, "/ 100. Compare Studies is preloaded with GSE16476 vs GSE60450."))
      ai_status(list(status = "Local demo ready", detail = "Prepared GSE19804, GSE16476, and GSE60450 demo datasets are available without live download."))
      showNotification("Local SeqSurf demo loaded.", type = "message", duration = 6)
    } else {
      start_study_status("Demo data folders are present, but GSE19804 study memory was not found.")
      showNotification("GSE19804 study memory not found in data/ai_knowledge_base.", type = "warning", duration = 8)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$start_get_study_info, {
    accession <- normalize_geo_accession(input$start_geo_accession)
    notes <- trimws(input$start_own_dataset_notes)
    if (is_geo_accession(accession)) {
      ai_status(list(status = "Reading study context", detail = paste("Fetching GEO/PubMed context for", accession)))
      start_study_status(paste("Retrieving GEO/PubMed study information for", accession, "..."))
      withProgress(message = paste("Getting study info for", accession), value = 0.3, {
        result <- tryCatch({
          kb <- get_or_build_study_knowledge_base(accession, notes, "data")
          eset <- tryCatch(fetch_geo_expression_set(accession, "data"), error = function(e) NULL)
          list(kb = kb, eset = eset)
        }, error = function(e) e)
      })
      if (inherits(result, "error")) {
        start_study_status(conditionMessage(result))
        showNotification(conditionMessage(result), type = "error", duration = 12)
        return()
      }
      kb <- result$kb
      sync_geo_accession(accession, "start")
      if (!is.null(result$eset)) {
        geo_eset(result$eset)
        geo_expression_set_status("available")
      } else {
        geo_eset(NULL)
        geo_expression_set_status("unavailable")
      }
      assistant_kbs(available_study_kbs("data"))
      kb_version(kb_version() + 1)
      start_study_kb(kb)
      geo_count_candidates(geo_supplementary_count_candidates(kb, accession))
      ai_status(list(status = "Study scored", detail = paste("Scorecard, prior-work report, and reproduction planner are ready for", accession)))
      start_study_narrative("Scorecard ready. Click Generate AI narrative for a written interpretation.")
      start_study_info(study_kb_summary_table(kb))
      assessment <- reanalysis_assessment_table(kb)
      proceed <- assessment$Value[assessment$Score == "Proceed recommendation"]
      start_study_status(paste("Study info loaded for", accession, ". Proceed recommendation score:", proceed, "/ 100. Send this study to Prepare Evidence or Validate."))
      return()
    }
    if (nzchar(notes)) {
      start_study_kb(NULL)
      geo_expression_set_status(NULL)
      ai_status(list(status = "Custom notes loaded", detail = "Paste/upload expression data or use GEO import for automated study context."))
      start_study_narrative("Custom notes loaded. GEO-backed AI narrative currently requires a trained study knowledge base.")
      start_study_info(data.frame(
        Field = c("Dataset type", "User-provided study notes"),
        Value = c("Custom/private dataset", notes),
        check.names = FALSE
      ))
      start_study_status("Custom dataset notes loaded. Validate automation currently works best for public GEO accessions; use these notes in Prepare Evidence, Discovery, or the paper-to-code brief.")
      return()
    }
    start_study_status("Enter a valid GEO accession such as GSE16476, or paste custom dataset notes.")
  }, ignoreInit = TRUE)

  observeEvent(input$start_generate_narrative, {
    kb <- start_study_kb()
    if (is.null(kb)) {
      showNotification("Get GEO study info first, then generate the narrative.", type = "warning")
      return()
    }
    if (!openai_available(active_openai_api_key())) {
      showNotification("Paste an OpenAI API key and click Save key for session before generating the AI recommendation.", type = "warning", duration = 10)
      start_study_narrative("Scorecard ready. Add a saved OpenAI API key to generate the written SeqSurfer recommendation.")
      return()
    }
    ai_status(list(status = "Writing narrative", detail = "Interpreting scorecard, prior work, and recommended reanalysis targets."))
    start_study_narrative("Generating AI narrative...")
    withProgress(message = "Writing AI reanalysis narrative", value = 0.5, {
      narrative <- tryCatch(
        generate_reanalysis_narrative(kb, active_openai_api_key(), input$openai_model),
        error = function(e) e
      )
    })
    if (inherits(narrative, "error")) {
      start_study_narrative(conditionMessage(narrative))
      showNotification(conditionMessage(narrative), type = "error", duration = 12)
      return()
    }
    ai_status(list(status = "Narrative ready", detail = "Review the written interpretation before deciding whether to import or compare results."))
    start_study_narrative(narrative)
  }, ignoreInit = TRUE)

  output$start_study_status <- renderUI({
    tags$span(start_study_status())
  })

  output$start_study_snapshot <- renderUI({
    kb <- start_study_kb()
    route <- geo_analyzability_route(kb, geo_expression_set_status())
    info <- start_study_info()

    if (is.null(kb)) {
      custom_notes <- if (nrow(info) > 0 && "Field" %in% names(info) && "Value" %in% names(info)) {
        first_nonempty(info$Value[info$Field == "User-provided study notes"], "")
      } else {
        ""
      }
      title <- if (nzchar(custom_notes)) "Custom/private dataset notes" else "No dataset assessed yet"
      note <- if (nzchar(custom_notes)) {
        shorten_text(custom_notes, 260)
      } else {
        "Enter a GEO accession or custom dataset notes, then assess the dataset."
      }
      return(tags$div(
        class = "study-snapshot",
        tags$div(
          class = "study-snapshot-header",
          tags$strong(title),
          tags$span(note)
        ),
        tags$div(class = "study-snapshot-item", tags$span("Organism"), tags$strong("Pending")),
        tags$div(class = "study-snapshot-item", tags$span("Assay type"), tags$strong("Pending")),
        tags$div(class = "study-snapshot-item", tags$span("N samples"), tags$strong("Pending")),
        tags$div(class = "study-snapshot-item", tags$span("Metadata"), tags$strong("Pending")),
        tags$div(class = "study-snapshot-item", tags$span("Expression files"), tags$strong("Pending")),
        tags$div(class = "study-snapshot-item", tags$span("Paper context"), tags$strong(if (nzchar(custom_notes)) "User notes provided" else "Pending"))
      ))
    }

    evidence_text <- paste(
      first_nonempty(kb$geo$title, ""),
      first_nonempty(kb$geo$summary, ""),
      first_nonempty(kb$geo$overall_design, ""),
      first_nonempty(kb$pubmed$title, ""),
      first_nonempty(kb$pubmed$abstract, ""),
      first_nonempty(kb$extracted$statistical_methods, ""),
      paste(kb$geo$platform_accessions, collapse = " "),
      sep = " "
    )
    evidence_low <- tolower(evidence_text)
    assay_type <- if (grepl("rna-seq|rnaseq|rna sequencing|transcriptome sequencing|hiseq|novaseq|nextseq|count matrix|read counts", evidence_low)) {
      "RNA-seq"
    } else if (grepl("microarray|affymetrix|agilent|beadchip|illumina array|expression array|probe", evidence_low)) {
      "Microarray / expression array"
    } else {
      "Not detected yet"
    }
    metadata_status <- if (nzchar(first_nonempty(kb$extracted$sample_groups, ""))) {
      "Sample groups detected"
    } else if (nzchar(first_nonempty(kb$geo$overall_design, ""))) {
      "Study design text available"
    } else {
      "Needs metadata review"
    }
    expression_status <- if (nrow(route) > 0) {
      paste(first_nonempty(route$Route, "Unknown route"), "-", first_nonempty(route$Can_run_in_app, "No"))
    } else {
      "Not checked"
    }
    related_count <- if (!is.null(kb$related_publications) && nrow(kb$related_publications) > 0 && "PubMed_ID" %in% names(kb$related_publications)) {
      sum(nzchar(kb$related_publications$PubMed_ID))
    } else {
      0
    }
    paper_context <- paste(
      if (nzchar(first_nonempty(kb$geo$pubmed_id, ""))) "Original paper linked" else "Original paper not linked",
      paste0(related_count, " citing/related papers found"),
      sep = "; "
    )
    platform_text <- if (length(kb$geo$platform_accessions) > 0) paste(kb$geo$platform_accessions, collapse = "; ") else "No platform parsed"

    tags$div(
      class = "study-snapshot",
      tags$div(
        class = "study-snapshot-header",
        tags$strong(first_nonempty(kb$geo$title, kb$accession)),
        tags$span(paste(first_nonempty(kb$accession, "GEO study"), "|", platform_text))
      ),
      tags$div(class = "study-snapshot-item", tags$span("Organism"), tags$strong(first_nonempty(kb$geo$organism, "Not reported"))),
      tags$div(class = "study-snapshot-item", tags$span("Assay type"), tags$strong(assay_type)),
      tags$div(class = "study-snapshot-item", tags$span("N samples"), tags$strong(first_nonempty(as.character(kb$geo$sample_count), "Unknown"))),
      tags$div(class = "study-snapshot-item", tags$span("Metadata"), tags$strong(metadata_status)),
      tags$div(class = "study-snapshot-item", tags$span("Expression files"), tags$strong(expression_status)),
      tags$div(class = "study-snapshot-item", tags$span("Paper context"), tags$strong(paper_context)),
      tags$p(class = "study-snapshot-note", shorten_text(first_nonempty(kb$extracted$cohort_characteristics, kb$geo$summary), 300))
    )
  })

  output$start_study_info_table <- renderDT({
    datatable(start_study_info(), options = list(pageLength = 12, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$data_quality_signals_panel <- renderUI({
    kb <- start_study_kb()

    rows <- list()

    # RIN score from GEO metadata
    if (!is.null(kb)) {
      rin_text <- tryCatch({
        summary_text <- paste(kb$geo$summary, kb$geo$overall_design, collapse = " ")
        if (grepl("RIN|RNA integrity", summary_text, ignore.case = TRUE)) {
          "RIN score mentioned in GEO metadata"
        } else {
          "No RIN score detected in GEO metadata"
        }
      }, error = function(e) "Unable to check")
      rows <- c(rows, list(data.frame(Signal = "RNA integrity", Finding = rin_text, Status = ifelse(grepl("mentioned", rin_text), "Present", "Review"), stringsAsFactors = FALSE)))

      n_samples <- if (!is.null(kb$geo$sample_count) && !is.na(kb$geo$sample_count)) kb$geo$sample_count else NA
      rows <- c(rows, list(data.frame(Signal = "Sample count", Finding = if (!is.na(n_samples)) paste(n_samples, "samples") else "Unknown", Status = if (!is.na(n_samples) && n_samples >= 6) "Sufficient" else "Review", stringsAsFactors = FALSE)))
    } else {
      rows <- c(rows, list(data.frame(Signal = "Study memory", Finding = "Load local demo or assess a GEO accession", Status = "Pending", stringsAsFactors = FALSE)))
      rows <- c(rows, list(data.frame(Signal = "Sample count", Finding = "Pending study lookup", Status = "Pending", stringsAsFactors = FALSE)))
    }

    eset_status <- geo_expression_set_status()
    active_dat <- tryCatch({
      if (!is.null(input$dataset_id) && nzchar(input$dataset_id)) load_rnaseq_dataset(input$dataset_id, "data") else NULL
    }, error = function(e) NULL)

    matrix_finding <- if (!is.null(active_dat)) {
      paste(nrow(active_dat$vsd), "genes x", ncol(active_dat$vsd), "samples")
    } else if (!is.null(eset_status) && nzchar(eset_status)) {
      eset_status
    } else {
      "Not loaded"
    }
    rows <- c(rows, list(data.frame(Signal = "Expression matrix", Finding = matrix_finding, Status = if (!is.null(active_dat) || (!is.null(eset_status) && grepl("available", eset_status, ignore.case = TRUE))) "Available" else "Pending", stringsAsFactors = FALSE)))

    pca_finding <- "Load expression matrix to run PCA check"
    pca_status <- "Pending"
    if (!is.null(active_dat) && !is.null(active_dat$pca) && nrow(active_dat$pca) > 2) {
      pcs <- intersect(c("PC1", "PC2"), names(active_dat$pca))
      if (length(pcs) >= 2) {
        scores <- active_dat$pca[, pcs, drop = FALSE]
        distances <- sqrt(rowSums(scale(scores)^2))
        outliers <- rownames(scores)[distances > 3]
        if (length(outliers) == 0) {
          pca_finding <- "No obvious PC1/PC2 outliers"
          pca_status <- "Pass"
        } else {
          pca_finding <- paste(length(outliers), "potential PCA outlier(s)")
          pca_status <- "Review"
        }
      } else {
        pca_finding <- "PCA columns unavailable"
        pca_status <- "Review"
      }
    } else {
      eset <- tryCatch(geo_eset(), error = function(e) NULL)
      if (!is.null(eset) && requireNamespace("Biobase", quietly = TRUE)) {
      expr <- tryCatch(Biobase::exprs(eset), error = function(e) NULL)
      if (!is.null(expr) && nrow(expr) > 1 && ncol(expr) > 2) {
        pca_result <- tryCatch({
          top_var <- head(order(apply(expr, 1, var, na.rm = TRUE), decreasing = TRUE), 500)
          pca <- prcomp(t(expr[top_var, , drop = FALSE]), scale. = TRUE)
          scores <- pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE]
          distances <- sqrt(rowSums(scale(scores)^2))
          outliers <- names(distances)[distances > 3]
          if (length(outliers) == 0) {
            list(finding = "No obvious outliers detected on PC1/PC2 (z > 3)", status = "Pass")
          } else {
            list(finding = paste("Potential outliers:", paste(head(outliers, 5), collapse = ", ")), status = "Review")
          }
        }, error = function(e) list(finding = "PCA check failed", status = "Unknown"))
          pca_finding <- pca_result$finding
          pca_status <- pca_result$status
      } else {
          pca_finding <- "Load expression matrix to run PCA check"
          pca_status <- "Pending"
      }
      }
    }
    rows <- c(rows, list(data.frame(Signal = "PCA outlier screen", Finding = pca_finding, Status = pca_status, stringsAsFactors = FALSE)))

    result <- dplyr::bind_rows(rows)
    status_color <- function(status) {
      if (status %in% c("Pass", "Available", "Sufficient", "Present")) return("#168a50")
      if (status %in% c("Review", "Unknown")) return("#c47f00")
      "#667085"
    }
    tags$div(
      style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 0.65rem;",
      lapply(seq_len(nrow(result)), function(i) {
        color <- status_color(result$Status[[i]])
        tags$div(
          style = paste0("border: 1px solid #d6dee8; border-left: 5px solid ", color, "; border-radius: 6px; padding: 0.65rem 0.75rem; background: #fbfcfe; min-width: 0;"),
          tags$div(style = "font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0; font-weight: 800; color: #667085;", result$Signal[[i]]),
          tags$div(style = "font-weight: 800; color: #17233a; margin-top: 0.2rem;", result$Status[[i]]),
          tags$div(style = "font-size: 0.86rem; color: #59636e; line-height: 1.25; overflow-wrap: anywhere;", result$Finding[[i]])
        )
      })
    )
  })

  output$prior_work_report_table <- renderDT({
    datatable(prior_work_report_table(start_study_kb()), options = list(pageLength = 6, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$exact_reproduction_table <- renderDT({
    datatable(exact_reproduction_feasibility_table(start_study_kb()), options = list(pageLength = 6, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$count_matrix_pipeline_plan <- renderText({
    count_matrix_pipeline_plan(start_study_kb())
  })

  output$raw_sra_handoff_plan <- renderText({
    raw_sra_handoff_plan(start_study_kb())
  })

  output$no_code_ai_handoff_panel <- renderUI({
    kb <- start_study_kb()
    accession <- if (is.null(kb)) "this GEO study" else normalize_geo_accession(kb$accession)
    tags$div(
      class = "ai-guide-body",
      tags$p("For datasets that cannot run fully inside SeqSurf, download a guided workflow and let an AI coding assistant help run it outside Shiny."),
      tags$ol(
        tags$li(tags$strong("Install tools: "), "R from CRAN, Positron or another R/Python IDE, and an AI coding assistant such as Claude Code, Codex, or Positron Assistant."),
        tags$li(tags$strong("Download from SeqSurf: "), "raw/SRA workflow bundle, no-code AI guide, and AI assistant prompt."),
        tags$li(tags$strong("Ask the AI to run the workflow: "), "it should resolve/download data, generate counts, write metadata, and document assumptions."),
        tags$li(tags$strong("Return to SeqSurf: "), "upload outputs/seqsurf_gene_counts.csv and outputs/seqsurf_metadata.csv, then use Validate and Discovery.")
      ),
      tags$p(class = "global-sidebar-note", paste("Current study:", accession))
    )
  })

  output$external_ai_analysis_prompt <- renderText({
    external_ai_analysis_prompt(start_study_kb())
  })

  observeEvent(input$create_external_workspace, {
    kb <- start_study_kb()
    if (is.null(kb)) {
      showNotification("Assess a GEO dataset first, then create the workflow folder.", type = "warning", duration = 8)
      return()
    }
    result <- tryCatch(write_external_workflow_workspace(kb, "data"), error = function(e) e)
    if (inherits(result, "error")) {
      external_workspace_info(NULL)
      showNotification(conditionMessage(result), type = "error", duration = 12)
      return()
    }
    external_workspace_info(result)
    showNotification("Workflow folder created. Open it in Positron, Codex, Claude Code, or double-click START_HERE.command.", type = "message", duration = 8)
  }, ignoreInit = TRUE)

  output$external_workspace_status <- renderUI({
    info <- external_workspace_info()
    if (is.null(info)) {
      return(tags$p(class = "global-sidebar-note", "No workflow folder created yet. Assess a GEO dataset, then create a local folder with scripts, prompts, and preflight checks."))
    }
    tags$div(
      tags$p(tags$strong("Folder ready: "), info$directory),
      tags$p(tags$strong("Double-click starter: "), basename(info$starter)),
      tags$p(tags$strong("Prompt file: "), basename(info$prompt)),
      tags$pre(info$command)
    )
  })

  output$download_count_matrix_plan <- downloadHandler(
    filename = function() {
      accession <- if (is.null(start_study_kb())) "study" else normalize_geo_accession(start_study_kb()$accession)
      sanitize_filename(paste0(accession, "_count_matrix_plan.md"))
    },
    content = function(file) {
      writeLines(count_matrix_pipeline_plan(start_study_kb()), file)
    }
  )

  output$download_raw_sra_plan <- downloadHandler(
    filename = function() {
      accession <- if (is.null(start_study_kb())) "study" else normalize_geo_accession(start_study_kb()$accession)
      sanitize_filename(paste0(accession, "_fastq_sra_handoff.md"))
    },
    content = function(file) {
      writeLines(raw_sra_handoff_plan(start_study_kb()), file)
    }
  )

  output$download_raw_sra_bundle <- downloadHandler(
    filename = function() {
      accession <- if (is.null(start_study_kb())) "study" else normalize_geo_accession(start_study_kb()$accession)
      sanitize_filename(paste0(accession, "_raw_sra_workflow_bundle.zip"))
    },
    content = function(file) {
      write_raw_sra_workflow_bundle(start_study_kb(), file)
    },
    contentType = "application/zip"
  )

  output$download_no_code_ai_guide <- downloadHandler(
    filename = function() {
      accession <- if (is.null(start_study_kb())) "study" else normalize_geo_accession(start_study_kb()$accession)
      sanitize_filename(paste0(accession, "_no_code_ai_handoff_guide.md"))
    },
    content = function(file) {
      writeLines(no_code_ai_handoff_guide(start_study_kb()), file)
    }
  )

  output$download_external_ai_prompt <- downloadHandler(
    filename = function() {
      accession <- if (is.null(start_study_kb())) "study" else normalize_geo_accession(start_study_kb()$accession)
      sanitize_filename(paste0(accession, "_prompt_for_ai_assistant.md"))
    },
    content = function(file) {
      writeLines(external_ai_analysis_prompt(start_study_kb()), file)
    }
  )

  output$download_demo_walkthrough <- downloadHandler(
    filename = function() "SeqSurf_GSE19804_demo_walkthrough.md",
    content = function(file) {
      source_path <- file.path("docs", "demo_walkthrough_GSE19804.md")
      if (file.exists(source_path)) {
        file.copy(source_path, file, overwrite = TRUE)
      } else {
        writeLines("Demo walkthrough file is missing from docs/demo_walkthrough_GSE19804.md.", file)
      }
    }
  )

  output$raw_supplementary_triage_table <- renderDT({
    datatable(raw_supplementary_triage_table(start_study_kb()), options = list(pageLength = 5, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$geo_analyzability_route_table <- renderDT({
    datatable(geo_analyzability_route(start_study_kb(), geo_expression_set_status()), options = list(pageLength = 3, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$geo_analyzability_route_summary <- renderUI({
    route <- geo_analyzability_route(start_study_kb(), geo_expression_set_status())
    can_run <- tolower(first_nonempty(route$Can_run_in_app, "No"))
    badge_class <- if (can_run %in% c("yes", "partial", "maybe")) can_run else "no"
    tags$div(
      class = "route-summary",
      tags$div(class = "route-summary-label", "Answer"),
      tags$div(class = "route-summary-value", tags$span(class = paste("route-badge", badge_class), first_nonempty(route$Can_run_in_app, "No"))),
      tags$div(class = "route-summary-label", "Route"),
      tags$div(class = "route-summary-value", first_nonempty(route$Route, "Not ready")),
      tags$div(class = "route-summary-label", "Confidence"),
      tags$div(class = "route-summary-value", paste0(first_nonempty(route$Confidence, 0), "%")),
      tags$div(class = "route-summary-label", "Next step"),
      tags$div(class = "route-summary-value", first_nonempty(route$Recommended_next_step, "Enter a GEO accession and assess the dataset."))
    )
  })

  output$geo_analyzability_checks_table <- renderDT({
    datatable(geo_analyzability_checks(start_study_kb()), options = list(pageLength = 6, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$start_reanalysis_assessment_table <- renderDT({
    datatable(
      reanalysis_assessment_table(start_study_kb()),
      options = list(pageLength = 8, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$start_reanalysis_recommendations_table <- renderDT({
    datatable(
      adaptive_reanalysis_recommendation_table(start_study_kb(), "data"),
      options = list(pageLength = 8, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$start_learning_prior_table <- renderDT({
    datatable(
      learning_prior_for_study(start_study_kb(), "data"),
      options = list(pageLength = 4, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  render_markdown_fragment <- function(text) {
    if (requireNamespace("commonmark", quietly = TRUE)) {
      HTML(commonmark::markdown_html(text))
    } else {
      tags$p(text)
    }
  }

  render_reanalysis_narrative <- function(text) {
    text <- first_nonempty(text, "")
    if (!nzchar(text)) {
      return(tags$div(class = "ai-narrative", tags$p("No recommendation generated yet.")))
    }
    if (grepl("^Generating AI narrative", text)) {
      return(tags$div(
        class = "ai-loading",
        tags$span(class = "ai-loading-pulse"),
        tags$div(
          tags$strong("Generating recommendation"),
          tags$p("Reading the scorecard, study context, and reanalysis route.", style = "margin: 0.12rem 0 0; color: #59636e;")
        )
      ))
    }

    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    heading_idx <- grep("^#{2,6}\\s+", lines)
    if (length(heading_idx) == 0) {
      paragraphs <- strsplit(text, "\n\n+", perl = TRUE)[[1]]
      return(tags$div(class = "ai-narrative", lapply(paragraphs, function(paragraph) tags$p(paragraph))))
    }

    sections <- lapply(seq_along(heading_idx), function(i) {
      start <- heading_idx[[i]]
      end <- if (i < length(heading_idx)) heading_idx[[i + 1]] - 1 else length(lines)
      title <- trimws(sub("^#{2,6}\\s+", "", lines[[start]]))
      body <- if (end > start) paste(lines[(start + 1):end], collapse = "\n") else ""
      list(title = title, body = trimws(body))
    })

    tags$div(
      class = "ai-narrative ai-recommendation-grid",
      lapply(sections, function(section) {
        tags$section(
          class = "ai-recommendation-section",
          tags$h4(section$title),
          tags$div(class = "ai-section-body", render_markdown_fragment(section$body))
        )
      })
    )
  }

  render_discovery_narrative <- function(text) {
    text <- first_nonempty(text, "")
    if (!nzchar(text)) {
      return(tags$div(
        class = "ai-narrative ai-recommendation-grid",
        tags$section(
          class = "ai-recommendation-section",
          tags$h4("What Discovery does"),
          tags$p("After Validate, Discovery recommends higher-value follow-up analyses, new comparisons, risks, and concrete tabs to open next.")
        ),
        tags$section(
          class = "ai-recommendation-section",
          tags$h4("Best first prompt"),
          tags$p("Ask what modern analyses are worth running next, what results are surprising, or what comparison would best strengthen the reanalysis.")
        )
      ))
    }
    if (grepl("^Generating", text, ignore.case = TRUE)) {
      return(tags$div(
        class = "ai-loading",
        tags$span(class = "ai-loading-pulse"),
        tags$div(
          tags$strong("Thinking through next analyses"),
          tags$p("Using study memory, DEG results, GSEA results, metadata, and paper agreement.", style = "margin: 0.12rem 0 0; color: #59636e;")
        )
      ))
    }
    if (grepl("^#{2,6}\\s+", text, perl = TRUE)) {
      return(render_reanalysis_narrative(text))
    }

    lines <- trimws(strsplit(text, "\n", fixed = TRUE)[[1]])
    lines <- lines[nzchar(lines)]
    numbered <- grep("^\\d+\\.\\s+", lines, value = TRUE)
    if (length(numbered) >= 2) {
      cards <- lapply(numbered, function(line) {
        clean <- sub("^\\d+\\.\\s+", "", line)
        title <- sub("^\\*\\*([^*]+)\\*\\*:?\\s*.*$", "\\1", clean)
        body <- if (!identical(title, clean)) sub("^\\*\\*[^*]+\\*\\*:?\\s*", "", clean) else clean
        tags$section(
          class = "ai-recommendation-section",
          tags$h4(title),
          tags$div(class = "ai-section-body", render_markdown_fragment(body))
        )
      })
      return(tags$div(class = "ai-narrative ai-recommendation-grid", cards))
    }

    paragraphs <- strsplit(text, "\n\n+", perl = TRUE)[[1]]
    tags$div(
      class = "ai-narrative ai-recommendation-grid",
      lapply(seq_along(paragraphs), function(i) {
        tags$section(
          class = "ai-recommendation-section",
          tags$h4(if (i == 1) "Discovery recommendation" else paste("Detail", i)),
          tags$div(class = "ai-section-body", render_markdown_fragment(paragraphs[[i]]))
        )
      })
    )
  }

  output$start_reanalysis_narrative <- renderUI({
    render_reanalysis_narrative(start_study_narrative())
  })

  geo_import_workflow_state <- reactive({
    if (is.null(geo_eset())) return("fetch")
    if (is.null(geo_design())) return("design")
    if (is.null(geo_imported_dataset())) return("run")
    if (is.null(geo_report())) return("report")
    "done"
  })

  import_button_class <- function(target, state) {
    if (identical(target, state)) "btn-primary next-action" else "btn-dark"
  }

  output$geo_import_action_buttons <- renderUI({
    state <- geo_import_workflow_state()
    tags$div(
      class = "import-action-stack",
      actionButton(
        "geo_fetch_train",
        "1. Fetch GEO",
        icon = icon("cloud-arrow-down"),
        class = import_button_class("fetch", state)
      ),
      actionButton(
        "geo_propose_design",
        "2. Propose design",
        icon = icon("diagram-project"),
        class = import_button_class("design", state)
      ),
      actionButton(
        "geo_run_import",
        "3. Run reanalysis",
        icon = icon("play"),
        class = import_button_class("run", state)
      ),
      actionButton(
        "geo_scan_counts",
        "Find count matrices",
        icon = icon("table"),
        class = "btn-dark"
      ),
      actionButton(
        "geo_run_count_import",
        "Run count-matrix import",
        icon = icon("calculator"),
        class = "btn-dark"
      ),
      actionButton(
        "compare_publication",
        "4. Validate paper agreement",
        icon = icon("scale-balanced"),
        class = import_button_class("report", state)
      )
    )
  })

  output$geo_design_editor <- renderUI({
    if (is.null(geo_design())) {
      return(tags$div(class = "assistant-empty", "Design fields appear after step 2."))
    }
    tagList(
      tags$hr(),
      tags$strong("Review analysis design"),
      selectInput(
        "geo_expression_type",
        "Expression data type",
        choices = c("raw_counts", "normalized_expression", "microarray_or_processed_expression"),
        selected = "normalized_expression"
      ),
      textInput("geo_group_col", "Group column", value = ""),
      textInput("geo_positive_group", "Positive group", value = ""),
      textInput("geo_negative_group", "Negative group", value = ""),
      textInput("geo_covariates", "Covariates", value = ""),
      textInput("geo_sample_id_col", "Sample ID column", value = "geo_accession"),
      textInput("geo_patient_id_col", "Patient ID column", value = "geo_accession"),
      textInput("geo_display_sample_col", "Display sample column", value = "title"),
      textInput("geo_timepoint_col", "Timepoint/design column", value = ""),
      textInput("geo_contrast_label", "Contrast label", value = "")
    )
  })

  output$geo_import_workflow_banner <- renderUI({
    state <- geo_import_workflow_state()
    step_class <- function(step) {
      if (identical(step, state)) return("import-step active")
      order <- c("fetch", "design", "run", "report", "done")
      if (match(step, order) < match(state, order)) "import-step done" else "import-step waiting"
    }
    hint <- switch(
      state,
      fetch = "Start here: fetch the GEO expression/phenotype data and train the study memory.",
      design = "Next: ask the AI to propose the primary comparison from paper methods and GEO metadata.",
      run = "Next: review/edit the proposed design in the sidebar, then run reanalysis.",
      report = "Next: validate paper agreement against the active reanalysis results.",
      done = "Ready: review the reproducibility score, then move to Discovery for next-step ideas.",
      "Follow the highlighted next step."
    )
    tagList(
      tags$div(class = "assistant-hint", tags$strong("Next step: "), hint),
      tags$div(
        class = "import-stepper",
        tags$div(class = step_class("fetch"), tags$strong("1. Fetch"), tags$span("Download GEO expression and phenotype.")),
        tags$div(class = step_class("design"), tags$strong("2. Design"), tags$span("Propose groups, contrast, covariates.")),
        tags$div(class = step_class("run"), tags$strong("3. Reanalyze"), tags$span("Run DEG, PCA, GSEA, pathway scores.")),
        tags$div(class = step_class("report"), tags$strong("4. Paper Agreement"), tags$span("Compare paper claims to app results.")),
        tags$div(class = step_class("done"), tags$strong("5. Score"), tags$span("Review reproducibility and differences."))
      )
    )
  })

  observeEvent(input$train_assistant, {
    accession <- normalize_geo_accession(input$assistant_geo)
    withProgress(message = paste("Training assistant for", accession), value = 0.2, {
      kb <- tryCatch(
        get_or_build_study_knowledge_base(accession, input$assistant_extra_papers, "data"),
        error = function(e) e
      )
    })
    if (inherits(kb, "error")) {
      showNotification(conditionMessage(kb), type = "error", duration = 10)
      return()
    }
    assistant_kbs(available_study_kbs("data"))
    kb_version(kb_version() + 1)
    select_assistant_kb(accession)
    showNotification(paste("Assistant trained for", accession), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$geo_fetch_train, {
    accession <- normalize_geo_accession(input$geo_import_accession)
    ai_status(list(status = "Fetching GEO", detail = paste("Downloading GEO metadata, expression matrix if available, and study context for", accession)))
    geo_import_status(paste("Fetching GEO expression, supplementary candidates, and study memory for", accession, "..."))
    withProgress(message = paste("Fetching", accession), value = 0.2, {
      result <- tryCatch({
        kb <- get_or_build_study_knowledge_base(accession, input$assistant_extra_papers, "data")
        eset <- tryCatch(fetch_geo_expression_set(accession, "data"), error = function(e) e)
        candidates <- geo_supplementary_count_candidates(kb, accession)
        list(eset = eset, kb = kb, candidates = candidates)
      }, error = function(e) e)
    })
    if (inherits(result, "error")) {
      geo_expression_set_status("unavailable")
      geo_import_status(conditionMessage(result))
      showNotification(conditionMessage(result), type = "error", duration = 12)
      return()
    }
    if (inherits(result$eset, "error")) {
      geo_eset(NULL)
      geo_expression_set_status("unavailable")
    } else {
      geo_eset(result$eset)
      geo_expression_set_status("available")
    }
    geo_count_candidates(result$candidates)
    count_ready <- nrow(result$candidates) > 0 && any(result$candidates$Status %in% c("Likely count matrix", "Possible table"))
    ai_status(list(status = "GEO fetched", detail = paste(accession, "study context is ready.", if (count_ready) "Supplementary count candidates were found." else "No count candidates detected.")))
    assistant_kbs(available_study_kbs("data"))
    kb_version(kb_version() + 1)
    select_assistant_kb(accession)
    updateTabsetPanel(session, "geo_import_tabset", selected = if (count_ready) "Count Files" else "Metadata")
    route_note <- if (identical(geo_expression_set_status(), "available")) {
      "Processed expression matrix available. Now propose/review the design."
    } else if (count_ready) {
      "No processed series matrix was imported, but supplementary count candidates were found. Review Count Files, propose design, then run count-matrix import."
    } else {
      "No processed series matrix or count-matrix candidate was found. Use the FASTQ/SRA handoff or add files manually."
    }
    geo_import_status(paste("Fetched study memory for", accession, ".", route_note))
  }, ignoreInit = TRUE)

  observeEvent(input$geo_propose_design, {
    accession <- normalize_geo_accession(input$geo_import_accession)
    ai_status(list(status = "Designing analysis", detail = paste("Choosing contrast, groups, and covariates for", accession)))
    geo_import_status(paste("Proposing analysis design for", accession, "..."))
    withProgress(message = "Reading methods and proposing design", value = 0.5, {
      design <- tryCatch({
        proposed <- tryCatch(
          propose_geo_analysis_design(accession, read_study_kb(accession, "data"), active_openai_api_key(), input$openai_model, "data"),
          error = function(e) NULL
        )
        if (!is.null(proposed)) {
          proposed
        } else {
          metadata <- fetch_geo_sample_metadata(accession, "data")
          rule_based_geo_metadata_design(accession, metadata)
        }
      }, error = function(e) e)
    })
    if (inherits(design, "error")) {
      geo_import_status(conditionMessage(design))
      showNotification(conditionMessage(design), type = "error", duration = 12)
      return()
    }
    geo_design(design)
    session$onFlushed(function() {
      updateSelectInput(session, "geo_expression_type", selected = design$expression_data_type)
      updateTextInput(session, "geo_group_col", value = design$group_col)
      updateTextInput(session, "geo_positive_group", value = design$positive_group)
      updateTextInput(session, "geo_negative_group", value = design$negative_group)
      updateTextInput(session, "geo_covariates", value = paste(design$covariates, collapse = ", "))
      updateTextInput(session, "geo_sample_id_col", value = design$sample_id_col)
      updateTextInput(session, "geo_patient_id_col", value = design$patient_id_col)
      updateTextInput(session, "geo_display_sample_col", value = design$display_sample_col)
      updateTextInput(session, "geo_timepoint_col", value = design$default_timepoint_col)
      updateTextInput(session, "geo_contrast_label", value = design$contrast_label)
    }, once = TRUE)
    updateTabsetPanel(session, "geo_import_tabset", selected = "Design")
    ai_status(list(status = "Design ready", detail = "Review the proposed contrast and edit fields before running reanalysis."))
    geo_import_status("Design proposed. Review/edit fields, then run reanalysis.")
  }, ignoreInit = TRUE)

  observeEvent(input$geo_run_import, {
    accession <- normalize_geo_accession(input$geo_import_accession)
    ai_status(list(status = "Running reanalysis", detail = paste("Building PCA, DEG, GSEA, pathway scores, and report for", accession)))
    design <- parse_design_from_inputs(
      accession,
      input$geo_expression_type,
      input$geo_expression_type,
      input$geo_group_col,
      input$geo_positive_group,
      input$geo_negative_group,
      input$geo_covariates,
      input$geo_sample_id_col,
      input$geo_patient_id_col,
      input$geo_display_sample_col,
      input$geo_timepoint_col,
      input$geo_contrast_label
    )
    geo_import_status(paste("Running reanalysis for", accession, "..."))
    withProgress(message = "Running GEO reanalysis and writing app dataset", value = 0.4, {
      result <- tryCatch(
        write_geo_app_dataset(accession, design, "data", read_study_kb(accession, "data")),
        error = function(e) e
      )
    })
    if (inherits(result, "error")) {
      geo_import_status(conditionMessage(result))
      showNotification(conditionMessage(result), type = "error", duration = 15)
      return()
    }
    geo_report(result$report_path)
    geo_imported_dataset(result$dataset_id)
    datasets <<- get("available_datasets", mode = "function")("data")
    dataset_choices <<- active_dataset_choices(datasets)
    updateSelectInput(session, "dataset_id", choices = dataset_choices, selected = result$dataset_id)
    select_assistant_kb(accession)
    updateTabsetPanel(session, "geo_import_tabset", selected = "Report")
    ai_status(list(status = "Reanalysis ready", detail = paste(result$dataset_id, "is imported. Click Validate paper agreement, then review the Score tab.")))
    geo_import_status(paste("Imported", result$dataset_id, "with", result$rows, "samples and", result$genes, "genes. Report:", result$report_path))
  }, ignoreInit = TRUE)

  output$geo_import_status <- renderText(geo_import_status())

  output$geo_design_table <- renderDT({
    datatable(design_to_table(geo_design()), options = list(pageLength = 20, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$geo_pheno_preview <- renderDT({
    eset <- geo_eset()
    if (is.null(eset)) {
      accession <- normalize_geo_accession(input$geo_import_accession)
      metadata <- tryCatch(fetch_geo_sample_metadata(accession, "data"), error = function(e) NULL)
      if (is.null(metadata)) {
        return(datatable(data.frame(Status = "Fetch GEO first."), options = list(dom = "t"), rownames = FALSE))
      }
      return(datatable(utils::head(round_for_display(metadata), 25), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE))
    }
    datatable(utils::head(round_for_display(Biobase::pData(eset)), 25), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$geo_count_candidate_table <- renderDT({
    candidates <- geo_count_candidates()
    if (is.null(candidates) || nrow(candidates) == 0) {
      candidates <- data.frame(File = "No supplementary count candidates scanned yet.", URL = "", Score = 0, Candidate_type = "", Status = "Not scanned", check.names = FALSE)
    }
    datatable(candidates, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  observeEvent(input$geo_scan_counts, {
    accession <- normalize_geo_accession(input$geo_import_accession)
    kb <- read_study_kb(accession, "data")
    if (is.null(kb)) {
      withProgress(message = paste("Building study context for", accession), value = 0.3, {
        kb <- tryCatch(get_or_build_study_knowledge_base(accession, input$assistant_extra_papers, "data"), error = function(e) e)
      })
      if (inherits(kb, "error")) {
        showNotification(conditionMessage(kb), type = "error", duration = 12)
        return()
      }
      assistant_kbs(available_study_kbs("data"))
      kb_version(kb_version() + 1)
    }
    candidates <- geo_supplementary_count_candidates(kb, accession)
    geo_count_candidates(candidates)
    updateTabsetPanel(session, "geo_import_tabset", selected = "Count Files")
    ready <- sum(candidates$Status %in% c("Likely count matrix", "Possible table"))
    geo_import_status(paste("Scanned", nrow(candidates), "supplementary files for", accession, "and found", ready, "count/table candidates."))
  }, ignoreInit = TRUE)

  observeEvent(input$geo_run_count_import, {
    accession <- normalize_geo_accession(input$geo_import_accession)
    design <- parse_design_from_inputs(
      accession,
      input$geo_expression_type,
      "raw_counts",
      input$geo_group_col,
      input$geo_positive_group,
      input$geo_negative_group,
      input$geo_covariates,
      input$geo_sample_id_col,
      input$geo_patient_id_col,
      input$geo_display_sample_col,
      input$geo_timepoint_col,
      input$geo_contrast_label
    )
    ai_status(list(status = "Importing count matrix", detail = paste("Scanning GEO supplements and building an app-ready count dataset for", accession)))
    geo_import_status(paste("Importing supplementary count matrix for", accession, "..."))
    withProgress(message = "Importing GEO supplementary count matrix", value = 0.4, {
      result <- tryCatch(
        write_geo_supplementary_count_dataset(accession, design, "data", read_study_kb(accession, "data")),
        error = function(e) e
      )
    })
    if (inherits(result, "error")) {
      geo_import_status(conditionMessage(result))
      showNotification(conditionMessage(result), type = "error", duration = 15)
      return()
    }
    geo_report(result$report_path)
    geo_imported_dataset(result$dataset_id)
    datasets <<- get("available_datasets", mode = "function")("data")
    dataset_choices <<- active_dataset_choices(datasets)
    updateSelectInput(session, "dataset_id", choices = dataset_choices, selected = result$dataset_id)
    select_assistant_kb(accession)
    updateTabsetPanel(session, "geo_import_tabset", selected = "Report")
    ai_status(list(status = "Count import ready", detail = paste(result$dataset_id, "was imported from", result$candidate_file, "using", result$count_method)))
    geo_import_status(paste("Imported", result$dataset_id, "from supplementary count candidate", result$candidate_file, "with", result$rows, "samples and", result$genes, "genes using", result$count_method, ". Script:", result$pipeline_script, "Report:", result$report_path))
  }, ignoreInit = TRUE)

  output$download_count_pipeline_script <- downloadHandler(
    filename = function() {
      dataset_id <- first_nonempty(geo_imported_dataset(), input$dataset_id, "seqsurf_count_matrix")
      sanitize_filename(paste0(dataset_id, "_count_matrix_pipeline.R"))
    },
    content = function(file) {
      dataset_id <- first_nonempty(geo_imported_dataset(), input$dataset_id, "")
      script_path <- file.path("data", dataset_id, "count_matrix_pipeline.R")
      if (nzchar(dataset_id) && file.exists(script_path)) {
        file.copy(script_path, file, overwrite = TRUE)
      } else {
        writeLines("# No generated count-matrix pipeline script yet. Run count-matrix import first.", file)
      }
    },
    contentType = "text/x-r-source"
  )

  output$geo_report_link <- renderUI({
    path <- geo_report()
    if (is.null(path) || !file.exists(path)) return(tags$p("No generated report yet."))
    tags$div(
      tags$p(tags$strong("Report ready."), " Use the download button below to open the generated HTML report."),
      tags$p(class = "global-sidebar-note", basename(path))
    )
  })

  output$download_geo_report <- downloadHandler(
    filename = function() {
      path <- geo_report()
      if (is.null(path) || !file.exists(path)) return("seqsurf_reproducibility_report.html")
      basename(path)
    },
    content = function(file) {
      path <- geo_report()
      if (is.null(path) || !file.exists(path)) {
        writeLines("<html><body><p>No generated report yet.</p></body></html>", file)
      } else {
        file.copy(path, file, overwrite = TRUE)
      }
    },
    contentType = "text/html"
  )

  selected_kb <- reactive({
    kb_version()
    if (is.null(input$assistant_kb) || !nzchar(input$assistant_kb)) return(NULL)
    if (!nzchar(input$assistant_kb)) return(NULL)
    read_study_kb(input$assistant_kb, "data")
  })

  assistant_workflow_state <- reactive({
    kb <- selected_kb()
    if (is.null(kb)) return("train")
    docs_ready <- !is.null(kb$documents) && nrow(kb$documents) > 0 &&
      any(grepl("indexed", kb$documents$Status, ignore.case = TRUE))
    claims_ready <- !is.null(kb$llm) && !is.null(kb$claims) && nrow(normalize_claim_table(kb$claims)) > 0
    if (!docs_ready) return("documents")
    if (!claims_ready) return("claims")
    "ready"
  })

  workflow_button_class <- function(target, state) {
    if (identical(target, state)) "btn-primary next-action" else "btn-dark"
  }

  output$assistant_action_buttons <- renderUI({
    state <- assistant_workflow_state()
    tags$div(
      class = "assistant-action-stack",
      actionButton(
        "train_assistant",
        "1. Train from original paper",
        icon = icon("brain"),
        class = workflow_button_class("train", state)
      ),
      actionButton(
        "fetch_study_documents",
        "2. Fetch evidence text",
        icon = icon("file-arrow-down"),
        class = workflow_button_class("documents", state)
      ),
      actionButton(
        "llm_extract_claims",
        "3. Build Study Memory 2.0",
        icon = icon("wand-magic-sparkles"),
        class = workflow_button_class("claims", state)
      )
    )
  })

  output$assistant_workflow_banner <- renderUI({
    state <- assistant_workflow_state()
    step_class <- function(step) {
      if (identical(step, state)) return("assistant-step active")
      order <- c("train", "documents", "claims", "ready")
      if (match(step, order) < match(state, order)) "assistant-step done" else "assistant-step waiting"
    }
    hint <- switch(
      state,
      train = "Start here: train study memory on GEO metadata, the original linked paper, related PubMed papers mentioning the accession, and any papers you add.",
      documents = "Next: fetch paper evidence with Jina Reader first, then direct PDF/text fallback. This skips large GEO data archives.",
      claims = "Next: run the Study Memory 2.0 virtual lab to extract claims, methods, evidence, reproducibility risks, and discovery ideas.",
      ready = "Study memory is ready. Move to Validate for reanalysis and paper-agreement scoring.",
      "Follow the highlighted next step."
    )
    tagList(
      tags$div(class = "assistant-hint", tags$strong("Next step: "), hint),
      tags$div(
        class = "assistant-stepper",
        tags$div(class = step_class("train"), tags$strong("1. Train"), tags$span("GEO, linked paper, related papers, and added notes.")),
        tags$div(class = step_class("documents"), tags$strong("2. Evidence"), tags$span("Jina Reader text plus PDF/text fallback.")),
        tags$div(class = step_class("claims"), tags$strong("3. Memory"), tags$span("Virtual lab extracts claims and study logic.")),
        tags$div(class = step_class("ready"), tags$strong("4. Ready"), tags$span("Next tab: Validate."))
      ),
      if (identical(state, "ready")) {
        tags$div(
          class = "ai-guide-actions",
          actionButton("train_go_validate", "Go to Validate", icon = icon("arrow-right"), class = "btn-primary")
        )
      }
    )
  })

  observeEvent(input$train_go_validate, {
    kb <- selected_kb()
    if (!is.null(kb) && !is.null(kb$accession)) {
      updateTextInput(session, "geo_import_accession", value = normalize_geo_accession(kb$accession))
    }
    bslib::nav_select("main_nav", selected = "Validate", session = session)
  }, ignoreInit = TRUE)

  output$assistant_setup_table <- renderDT({
    datatable(
      assistant_setup_status(active_openai_api_key()),
      options = list(pageLength = 10, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$assistant_kb_table <- renderDT({
    datatable(
      study_kb_summary_table(selected_kb()),
      options = list(pageLength = 11, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$virtual_lab_summary_table <- renderDT({
    datatable(
      virtual_lab_summary_table(selected_kb()),
      options = list(pageLength = 5, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$virtual_lab_summary_table_full <- renderDT({
    datatable(
      virtual_lab_summary_table(selected_kb()),
      options = list(pageLength = 5, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$virtual_lab_recommendation_table <- renderDT({
    datatable(
      virtual_lab_recommendation_table(selected_kb()),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$citation_evidence_table <- renderDT({
    datatable(
      citation_evidence_table(selected_kb()),
      options = list(pageLength = 12, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$assistant_documents_table <- renderDT({
    kb <- selected_kb()
    docs <- if (is.null(kb)) {
      data.frame(Source = "No study knowledge base selected.", URL = "", Local_path = "", Text_path = "", Ingestion = "", Status = "Needs training")
    } else if (!is.null(kb$documents) && nrow(kb$documents) > 0) {
      kb$documents
    } else {
      discover_document_sources(kb)
    }
    datatable(docs, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$paperqa_corpus_table <- renderDT({
    kb <- selected_kb()
    corpus <- if (is.null(kb)) {
      data.frame(Source = "No study knowledge base selected.", URL = "", Text_file = "", Status = "Needs training", check.names = FALSE)
    } else if (!is.null(kb$paperqa) && !is.null(kb$paperqa$manifest_path) && file.exists(kb$paperqa$manifest_path)) {
      read.csv(kb$paperqa$manifest_path, check.names = FALSE)
    } else {
      data.frame(Source = "No PaperQA corpus prepared yet", URL = "", Text_file = "", Status = "Fetch evidence text first", check.names = FALSE)
    }
    datatable(corpus, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  output$assistant_related_papers_table <- renderDT({
    kb <- selected_kb()
    pubs <- if (is.null(kb) || is.null(kb$related_publications)) {
      empty_related_publications("Train a GEO study first.")
    } else {
      kb$related_publications
    }
    datatable(pubs, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  output$assistant_sections_table <- renderDT({
    datatable(
      summarize_document_sections(selected_kb()),
      options = list(pageLength = 8, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$assistant_claim_candidates_table <- renderDT({
    kb <- selected_kb()
    candidates <- if (is.null(kb)) {
      data.frame(Claim_candidate = "Prepare Evidence or select a study knowledge base first.", Evidence_type = "Not ready", Status = "Not ready", check.names = FALSE)
    } else if (!is.null(kb$claim_candidates) && nrow(kb$claim_candidates) > 0) {
      kb$claim_candidates
    } else {
      section_claim_candidate_table(paste(summarize_document_sections(kb)$Excerpt, collapse = " "))
    }
    datatable(candidates, options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  observeEvent(input$fetch_study_documents, {
    withProgress(message = "Fetching study documents", value = 0.4, {
      kb <- tryCatch(fetch_and_ingest_documents(selected_kb(), "data"), error = function(e) e)
    })
    if (inherits(kb, "error")) {
      showNotification(conditionMessage(kb), type = "error", duration = 10)
      return()
    }
    kb_version(kb_version() + 1)
    showNotification("Study evidence fetched with Jina/direct fallback and PaperQA corpus prepared.", type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$llm_extract_claims, {
    withProgress(message = "Running SeqSurf Study Memory 2.0 virtual lab", value = 0.4, {
      kb <- tryCatch(
        run_study_memory_virtual_lab(selected_kb(), active_openai_api_key(), input$openai_model, "data"),
        error = function(e) e
      )
    })
    if (inherits(kb, "error")) {
      showNotification(conditionMessage(kb), type = "error", duration = 12)
      return()
    }
    kb_version(kb_version() + 1)
    showNotification("Study Memory 2.0 virtual lab complete.", type = "message")
  }, ignoreInit = TRUE)

  output$assistant_claims_table <- renderDT({
    kb <- selected_kb()
    claims <- if (is.null(kb)) empty_claim_table() else normalize_claim_table(kb$claims)
    datatable(claims, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$assistant_verification_table <- renderDT({
    datatable(
      verify_dataset_against_kb(dataset(), selected_kb()),
      options = list(pageLength = 10, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$dataset_match_warning <- renderUI({
    kb <- selected_kb()
    if (is.null(kb)) {
      return(tags$div(class = "assistant-hint", tags$strong("Paper/dataset check: "), "Prepare Evidence or select a study knowledge base first."))
    }
    dat <- dataset()
    matches <- dataset_matches_accession(dat$id, kb$accession, "data")
    if (matches) {
      return(tags$div(
        class = "assistant-hint",
        style = "border-left-color:#2e7d32;background:#f1f8f2;",
        tags$strong("Paper/dataset match: "),
        kb$accession,
        " is being compared against ",
        dat$info$display_name,
        "."
      ))
    }
    tags$div(
      class = "assistant-hint",
      style = "border-left-color:#b26a00;background:#fff8ed;",
      tags$strong("Paper/dataset mismatch: "),
      "You trained the assistant on ",
      kb$accession,
      " but the selected dataset is ",
      dat$info$display_name,
      ". Import/select the matching GEO dataset before trusting the verification."
    )
  })

  publication_comparison <- eventReactive(input$compare_publication, {
    compare_findings_against_publication(
      dataset(),
      selected_kb(),
      input$publication_padj,
      input$publication_lfc
    )
  }, ignoreInit = FALSE)

  observeEvent(input$compare_publication, {
    updateTabsetPanel(session, "geo_import_tabset", selected = "Score")
    showNotification("Comparison updated. Review the Score tab for reproducibility and differences.", type = "message", duration = 6)
  }, ignoreInit = TRUE)

  output$reproducibility_score_card <- renderUI({
    comparison <- tryCatch(publication_comparison(), error = function(e) NULL)
    if (is.null(comparison) || nrow(comparison) == 0) return(NULL)

    status_col <- intersect(c("Status", "Match", "Outcome", "Result"), names(comparison))[1]
    if (is.na(status_col)) return(NULL)

    statuses <- comparison[[status_col]]
    confirmed <- sum(grepl("confirm|support|agree|match", statuses, ignore.case = TRUE))
    total <- nrow(comparison)
    pct <- round(100 * confirmed / max(total, 1))

    color <- if (pct >= 60) "#1a7f37" else if (pct >= 30) "#b08800" else "#cf222e"

    tags$div(
      style = "margin-bottom: 1rem;",
      card(
        style = paste0("border-top: 4px solid ", color, ";"),
        card_body(
          tags$div(
            style = "display: flex; align-items: center; gap: 1.5rem;",
            tags$div(
              tags$span(style = paste0("font-size: 2.5rem; font-weight: 900; color: ", color, ";"), paste0(pct, "%")),
              tags$div(style = "font-size: 1rem; font-weight: 700; color: #12263f;", "Reproducibility Score"),
              tags$div(style = "font-size: 0.85rem; color: #59636e;", paste0(confirmed, " of ", total, " paper claims confirmed in standardized reanalysis"))
            ),
            tags$div(
              style = "flex: 1; font-size: 0.82rem; color: #3b4652; line-height: 1.5;",
              tags$strong("What this means: "), "Claims are matched by direction of effect and statistical significance. Mismatches may reflect preprocessing differences, cohort subsets, or threshold choices — not necessarily errors."
            )
          )
        )
      )
    )
  })

  html_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }

  simple_html_table <- function(df, max_rows = 30) {
    if (is.null(df) || nrow(df) == 0) {
      return(tags$p(class = "assistant-empty", "No rows to display. Run LLM extract claims first, then compare to paper."))
    }
    shown <- utils::head(df, max_rows)
    tags$table(
      class = "table table-sm table-striped table-bordered",
      tags$thead(tags$tr(lapply(names(shown), function(name) tags$th(html_escape(name))))),
      tags$tbody(lapply(seq_len(nrow(shown)), function(i) {
        tags$tr(lapply(shown[i, , drop = TRUE], function(value) tags$td(html_escape(value))))
      }))
    )
  }

  output$publication_comparison_table <- renderUI({
    dat <- publication_comparison()
    status_counts <- if ("Interpretation" %in% names(dat)) {
      as.data.frame(table(dat$Interpretation), responseName = "n")
    } else {
      data.frame(Interpretation = "Rows", n = nrow(dat))
    }
    tagList(
      tags$div(
        class = "assistant-hint",
        tags$strong("Comparison generated: "),
        nrow(dat),
        " paper-derived claims checked against ",
        dataset()$info$display_name,
        "."
      ),
      simple_html_table(status_counts, max_rows = 12),
      simple_html_table(dat, max_rows = 40)
    )
  })

  output$post_reproducibility_score_table <- renderDT({
    datatable(
      score_reproducibility_from_comparison(publication_comparison()),
      options = list(pageLength = 6, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$relevant_different_results_table <- renderDT({
    datatable(
      summarize_relevant_different_results(
        dataset(),
        selected_kb(),
        publication_comparison(),
        input$publication_padj,
        input$publication_lfc
      ),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  comparison_results <- eventReactive(input$run_comparison, {
    req(nzchar(input$compare_dataset_a), nzchar(input$compare_dataset_b))
    req(input$compare_dataset_a != input$compare_dataset_b)

    load_deg <- function(dataset_id) {
      path <- file.path("data", dataset_id, "deg_results.csv")
      if (!file.exists(path)) return(NULL)
      tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
    }
    load_pathway <- function(dataset_id) {
      path <- file.path("data", dataset_id, "gsea_hallmark.csv")
      if (!file.exists(path)) return(NULL)
      tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
    }

    deg_a <- load_deg(input$compare_dataset_a)
    deg_b <- load_deg(input$compare_dataset_b)
    pw_a  <- load_pathway(input$compare_dataset_a)
    pw_b  <- load_pathway(input$compare_dataset_b)

    list(deg_a = deg_a, deg_b = deg_b, pw_a = pw_a, pw_b = pw_b,
         id_a = input$compare_dataset_a, id_b = input$compare_dataset_b)
  })

  comparison_tables <- reactive({
    res <- comparison_results()
    deg_a <- res$deg_a; deg_b <- res$deg_b

    empty <- list(
      summary = data.frame(Metric = character(), Value = numeric(), Detail = character()),
      deg_overlap = data.frame(),
      pathway = data.frame(),
      biomarkers = data.frame(),
      id_a = res$id_a,
      id_b = res$id_b
    )

    if (is.null(deg_a) || is.null(deg_b)) return(empty)

    gene_col_a <- intersect(c("gene", "Gene", "gene_symbol", "symbol"), names(deg_a))[1]
    gene_col_b <- intersect(c("gene", "Gene", "gene_symbol", "symbol"), names(deg_b))[1]
    padj_col_a <- intersect(c("padj", "adj.P.Val", "FDR", "qvalue"), names(deg_a))[1]
    padj_col_b <- intersect(c("padj", "adj.P.Val", "FDR", "qvalue"), names(deg_b))[1]
    lfc_col_a <- intersect(c("log2FoldChange", "logFC", "log2FC"), names(deg_a))[1]
    lfc_col_b <- intersect(c("log2FoldChange", "logFC", "log2FC"), names(deg_b))[1]

    sig_a <- character()
    sig_b <- character()
    deg_overlap <- data.frame()
    if (!is.na(gene_col_a) && !is.na(gene_col_b)) {
      sig_a <- if (!is.na(padj_col_a)) deg_a[!is.na(deg_a[[padj_col_a]]) & deg_a[[padj_col_a]] < 0.05, gene_col_a] else deg_a[[gene_col_a]]
      sig_b <- if (!is.na(padj_col_b)) deg_b[!is.na(deg_b[[padj_col_b]]) & deg_b[[padj_col_b]] < 0.05, gene_col_b] else deg_b[[gene_col_b]]
      overlap <- intersect(sig_a, sig_b)
      if (length(overlap) > 0 && !is.na(lfc_col_a) && !is.na(lfc_col_b)) {
        merge_a <- deg_a[deg_a[[gene_col_a]] %in% overlap, c(gene_col_a, lfc_col_a, padj_col_a), drop = FALSE]
        merge_b <- deg_b[deg_b[[gene_col_b]] %in% overlap, c(gene_col_b, lfc_col_b, padj_col_b), drop = FALSE]
        names(merge_a) <- c("Gene", "log2FC_A", "padj_A")
        names(merge_b) <- c("Gene", "log2FC_B", "padj_B")
        deg_overlap <- merge(merge_a, merge_b, by = "Gene")
        deg_overlap$Direction <- ifelse(sign(deg_overlap$log2FC_A) == sign(deg_overlap$log2FC_B), "Same direction", "Opposite direction")
        deg_overlap$Max_abs_log2FC <- pmax(abs(deg_overlap$log2FC_A), abs(deg_overlap$log2FC_B), na.rm = TRUE)
      }
    }

    pathway <- data.frame()
    pw_a <- res$pw_a; pw_b <- res$pw_b
    if (!is.null(pw_a) && !is.null(pw_b)) {
      pw_col_a <- intersect(c("pathway", "Pathway", "Description", "Term"), names(pw_a))[1]
      pw_col_b <- intersect(c("pathway", "Pathway", "Description", "Term"), names(pw_b))[1]
      nes_col_a <- intersect(c("NES", "nes", "ES"), names(pw_a))[1]
      nes_col_b <- intersect(c("NES", "nes", "ES"), names(pw_b))[1]
      if (!is.na(pw_col_a) && !is.na(pw_col_b) && !is.na(nes_col_a) && !is.na(nes_col_b)) {
        shared <- intersect(pw_a[[pw_col_a]], pw_b[[pw_col_b]])
        if (length(shared) > 0) {
          sub_a <- pw_a[pw_a[[pw_col_a]] %in% shared, c(pw_col_a, nes_col_a), drop = FALSE]
          sub_b <- pw_b[pw_b[[pw_col_b]] %in% shared, c(pw_col_b, nes_col_b), drop = FALSE]
          names(sub_a) <- c("Pathway", "NES_A")
          names(sub_b) <- c("Pathway", "NES_B")
          pathway <- merge(sub_a, sub_b, by = "Pathway")
          pathway$Direction <- ifelse(sign(pathway$NES_A) == sign(pathway$NES_B), "Concordant", "Discordant")
          pathway$Label <- gsub("^HALLMARK_", "", pathway$Pathway)
        }
      }
    }

    biomarkers <- data.frame()
    kb_a <- tryCatch(read_study_kb(gsub("^geo_", "", res$id_a), "data"), error = function(e) NULL)
    kb_b <- tryCatch(read_study_kb(gsub("^geo_", "", res$id_b), "data"), error = function(e) NULL)
    claims_a <- if (!is.null(kb_a)) claim_terms(kb_a, "biomarker") else character()
    claims_b <- if (!is.null(kb_b)) claim_terms(kb_b, "biomarker") else character()
    all_claimed <- union(claims_a, claims_b)
    if (length(all_claimed) > 0) {
      biomarkers <- data.frame(
        Biomarker = rep(all_claimed, each = 4),
        Layer = rep(c("Claimed in A", "Claimed in B", "DEG in A", "DEG in B"), times = length(all_claimed)),
        Present = c(rbind(
          all_claimed %in% claims_a,
          all_claimed %in% claims_b,
          all_claimed %in% sig_a,
          all_claimed %in% sig_b
        )),
        stringsAsFactors = FALSE
      )
    }

    summary <- data.frame(
      Metric = c("Significant DEGs A", "Significant DEGs B", "Shared DEGs", "Concordant pathways", "Claimed markers seen"),
      Value = c(length(sig_a), length(sig_b), nrow(deg_overlap), sum(pathway$Direction == "Concordant", na.rm = TRUE), length(unique(biomarkers$Biomarker[biomarkers$Layer %in% c("DEG in A", "DEG in B") & biomarkers$Present]))),
      Detail = c(res$id_a, res$id_b, "padj < 0.05 in both", "shared Hallmark NES sign", "claimed and detected in either DEG list"),
      stringsAsFactors = FALSE
    )

    list(summary = summary, deg_overlap = deg_overlap, pathway = pathway, biomarkers = biomarkers, id_a = res$id_a, id_b = res$id_b)
  })

  output$compare_summary_cards <- renderUI({
    tables <- comparison_tables()
    if (nrow(tables$summary) == 0) {
      return(tags$p("Click Compare datasets after selecting two completed reanalysis datasets."))
    }
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 0.65rem; align-items: stretch; max-width: 100%;",
      lapply(seq_len(nrow(tables$summary)), function(i) {
        tags$div(
          style = "flex: 1 1 180px; border: 1px solid #d6dee8; border-left: 5px solid #0f8c8f; border-radius: 6px; padding: 0.6rem 0.7rem; background: #f8fbfd; min-width: 160px;",
          tags$div(style = "font-size: 1.3rem; font-weight: 800; color: #17233a; line-height: 1;", tables$summary$Value[[i]]),
          tags$div(style = "font-weight: 800; font-size: 0.92rem; line-height: 1.15; margin-top: 0.3rem; overflow-wrap: anywhere;", tables$summary$Metric[[i]]),
          tags$div(style = "font-size: 0.78rem; color: #667085; line-height: 1.2; margin-top: 0.2rem; overflow-wrap: anywhere;", tables$summary$Detail[[i]])
        )
      })
    )
  })

  output$compare_deg_overlap_plot <- renderPlot({
    dat <- comparison_tables()$deg_overlap
    validate(need(nrow(dat) > 0, "No significant DEG overlap detected for this pair."))
    label_dat <- dat[order(dat$Max_abs_log2FC, decreasing = TRUE), , drop = FALSE]
    label_dat <- utils::head(label_dat, 10)
    ggplot(dat, aes(x = log2FC_A, y = log2FC_B, color = Direction, size = -log10(pmax(padj_A, padj_B)))) +
      geom_hline(yintercept = 0, color = "#aab2bd", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "#aab2bd", linewidth = 0.4) +
      geom_point(alpha = 0.82) +
      geom_text(data = label_dat, aes(label = Gene), size = 3, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
      scale_color_manual(values = c("Same direction" = "#168a50", "Opposite direction" = "#c2410c")) +
      labs(x = "Dataset A log2FC", y = "Dataset B log2FC", color = NULL, size = "Significance") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  output$compare_pathway_plot <- renderPlot({
    dat <- comparison_tables()$pathway
    validate(need(nrow(dat) > 0, "No shared pathway results detected for this pair."))
    label_dat <- dat[order(abs(dat$NES_A - dat$NES_B), decreasing = TRUE), , drop = FALSE]
    label_dat <- utils::head(label_dat, 8)
    ggplot(dat, aes(x = NES_A, y = NES_B, color = Direction)) +
      geom_hline(yintercept = 0, color = "#aab2bd", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "#aab2bd", linewidth = 0.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#8792a2") +
      geom_point(size = 4, alpha = 0.86) +
      geom_text(data = label_dat, aes(label = Label), size = 3, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
      scale_color_manual(values = c("Concordant" = "#2563eb", "Discordant" = "#c2410c")) +
      labs(x = "Dataset A NES", y = "Dataset B NES", color = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  output$compare_biomarker_plot <- renderPlot({
    dat <- comparison_tables()$biomarkers
    validate(need(nrow(dat) > 0, "No claimed biomarkers found in the selected study memories."))
    dat$Layer <- factor(dat$Layer, levels = c("Claimed in A", "Claimed in B", "DEG in A", "DEG in B"))
    dat$Present_label <- ifelse(dat$Present, "Present", "Absent")
    ggplot(dat, aes(x = Layer, y = Biomarker, fill = Present_label)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = ifelse(Present, "yes", "")), color = "white", fontface = "bold", size = 3.5) +
      scale_fill_manual(values = c("Present" = "#0f8c8f", "Absent" = "#d9e2ec")) +
      labs(x = NULL, y = NULL, fill = NULL) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", panel.grid = element_blank())
  })

  output$gap_analysis_table <- renderDT({
    datatable(
      gap_analysis_report(dataset(), selected_kb()),
      options = list(pageLength = 12, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$assistant_suggestions_table <- renderDT({
    datatable(
      suggest_ai_analyses(dataset(), selected_kb()),
      options = list(pageLength = 10, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$comparison_suggestions_table <- renderDT({
    datatable(
      suggest_unexplored_comparisons(dataset(), selected_kb()),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$standardized_pipeline_table <- renderDT({
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    datatable(
      standardized_pipeline_checklist(active_dat, selected_kb()),
      options = list(pageLength = 8, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$sequencing_pitfall_table <- renderDT({
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    datatable(
      sequencing_pitfall_table(active_dat, selected_kb(), "data"),
      options = list(pageLength = 8, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$signature_screen_table <- renderDT({
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    datatable(
      signature_screen_table(active_dat, selected_kb()),
      options = list(pageLength = 5, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  observeEvent(input$save_learning_event, {
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    kb <- selected_kb()
    comparison <- tryCatch(publication_comparison(), error = function(e) NULL)
    accession <- first_nonempty(
      if (!is.null(kb)) kb$accession else "",
      input$geo_import_accession,
      input$assistant_geo,
      input$start_geo_accession,
      ""
    )
    row <- tryCatch(
      record_learning_event(
        accession = accession,
        dataset_id = first_nonempty(input$dataset_id, ""),
        event_type = "user_feedback",
        outcome = input$learning_outcome,
        usefulness_rating = input$learning_rating,
        notes = input$learning_notes,
        kb = kb,
        dat = active_dat,
        comparison = comparison,
        data_root = "data"
      ),
      error = function(e) e
    )
    if (inherits(row, "error")) {
      learning_event_status(conditionMessage(row))
      showNotification(conditionMessage(row), type = "error", duration = 10)
      return()
    }
    learning_event_status(paste("Saved learning event for", first_nonempty(row$accession, row$dataset_id, "this dataset"), "at", row$timestamp))
    showNotification("SeqSurf learning event saved locally.", type = "message")
  }, ignoreInit = TRUE)

  output$learning_event_status <- renderUI({
    tags$p(class = "global-sidebar-note", learning_event_status())
  })

  output$learning_summary_table <- renderDT({
    input$save_learning_event
    input$import_shared_learning
    datatable(
      learning_summary_table("data"),
      options = list(pageLength = 6, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$learning_lessons_table <- renderDT({
    input$save_learning_event
    input$import_shared_learning
    datatable(
      learning_lessons_table("data"),
      options = list(pageLength = 8, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$download_shared_learning_packet <- downloadHandler(
    filename = function() {
      paste0("seqsurf_shared_learning_packet_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(file) {
      if (!isTRUE(input$shared_learning_consent)) {
        stop("Check the consent box before exporting a shared learning packet.", call. = FALSE)
      }
      write_shared_learning_export(file, "data")
    },
    contentType = "application/zip"
  )

  observeEvent(input$import_shared_learning, {
    if (is.null(input$shared_learning_import)) {
      shared_learning_status("Choose a shared learning packet first.")
      showNotification("Choose a shared learning packet first.", type = "warning")
      return()
    }
    result <- tryCatch(
      import_shared_learning_packet(input$shared_learning_import$datapath, "data"),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      shared_learning_status(conditionMessage(result))
      showNotification(conditionMessage(result), type = "error", duration = 10)
      return()
    }
    shared_learning_status(paste0(result$Status[[1]], " Imported rows: ", result$Imported_rows[[1]], "; shared pool rows: ", result$Total_shared_rows[[1]], "."))
    showNotification("Shared learning packet imported.", type = "message")
  }, ignoreInit = TRUE)

  output$shared_learning_status <- renderUI({
    tags$p(class = "global-sidebar-note", shared_learning_status())
  })

  output$training_examples_summary_table <- renderDT({
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    datatable(
      training_examples_summary_table(
        selected_kb(),
        active_dat,
        "data",
        include_shared = isTRUE(input$training_include_shared),
        include_private = isTRUE(input$training_include_private)
      ),
      options = list(pageLength = 5, scrollX = TRUE, dom = "t"),
      rownames = FALSE
    )
  })

  output$download_training_examples <- downloadHandler(
    filename = function() {
      paste0("seqsurf_training_examples_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(file) {
      active_dat <- tryCatch(dataset(), error = function(e) NULL)
      write_training_data_export(
        file,
        selected_kb(),
        active_dat,
        "data",
        include_shared = isTRUE(input$training_include_shared),
        include_private = isTRUE(input$training_include_private)
      )
    },
    contentType = "application/zip"
  )

  output$ranked_discovery_cards <- renderUI({
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    comparison <- tryCatch(publication_comparison(), error = function(e) NULL)
    plan <- adaptive_ranked_discovery_plan(active_dat, selected_kb(), comparison, "data")
    tags$div(
      class = "ranked-discovery-grid",
      lapply(seq_len(nrow(plan)), function(i) {
        row <- plan[i, , drop = FALSE]
        tags$section(
          class = "ranked-discovery-card",
          tags$span(class = "rank", row$Rank),
          tags$h4(row$Analysis),
          tags$p(row$Why),
          tags$dl(
            tags$dt("Open"),
            tags$dd(row$Where_to_open),
            tags$dt("Evidence"),
            tags$dd(row$Evidence),
            tags$dt("Trust check"),
            tags$dd(row$Trust_check)
          )
        )
      })
    )
  })

  assistant_answer <- eventReactive(input$ask_assistant, {
    tryCatch(
      chat_with_study_assistant(dataset(), selected_kb(), input$assistant_question, active_openai_api_key(), input$openai_model),
      error = function(e) conditionMessage(e)
    )
  }, ignoreInit = TRUE)

  output$assistant_chat_response <- renderUI({
    if (input$ask_assistant == 0) {
      return(render_discovery_narrative(""))
    } else {
      render_discovery_narrative(assistant_answer())
    }
  })

  survival_screen <- eventReactive(input$run_survival_screen, {
    run_survival_summary(dataset())
  }, ignoreInit = TRUE)

  output$survival_screen_table <- renderDT({
    out <- if (is.null(input$run_survival_screen) || input$run_survival_screen == 0) {
      data.frame(Status = "Click Run survival screen when ready.")
    } else {
      survival_screen()
    }
    datatable(out, options = list(pageLength = 5, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  immune_screen <- eventReactive(input$run_immune_screen, {
    run_immune_deconvolution_summary(dataset())
  }, ignoreInit = TRUE)

  output$immune_screen_table <- renderDT({
    out <- if (is.null(input$run_immune_screen) || input$run_immune_screen == 0) {
      data.frame(Status = "Click Run immune deconvolution when ready.")
    } else {
      immune_screen()
    }
    datatable(out, options = list(pageLength = 5, scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  output$paper_to_code_brief <- renderText({
    active_dat <- tryCatch(dataset(), error = function(e) NULL)
    paper_to_code_brief(active_dat, selected_kb())
  })

  output$download_paper_to_code <- downloadHandler(
    filename = function() {
      accession <- if (is.null(input$assistant_kb) || !nzchar(input$assistant_kb)) first_nonempty(input$dataset_id, "study") else input$assistant_kb
      sanitize_filename(paste0(accession, "_paper_to_code_brief.md"))
    },
    content = function(file) {
      active_dat <- tryCatch(dataset(), error = function(e) NULL)
      writeLines(paper_to_code_brief(active_dat, selected_kb()), file)
    }
  )

  pca_figure <- reactive({
    req(input$pca_x, input$pca_y, input$pca_color)
    plot_pca(
      dataset()$pca,
      input$pca_x,
      input$pca_y,
      input$pca_color,
      input$pca_shape,
      input$pca_label,
      dataset()$pca_variance
    )
  })

  output$dataset_summary <- renderTable({
    metadata_summary(dataset()$metadata, dataset()$info)
  })

  output$study_information <- renderUI({
    dat <- dataset()
    dataset_study_content(dat$id, dat$metadata, dat$info)
  })

  output$response_summary <- renderTable({
    dat <- dataset()
    col <- dat$info$default_group_col
    validate(need(col %in% names(dat$metadata), "No default response column found."))
    as.data.frame(table(dat$metadata[[col]], useNA = "ifany"), responseName = "n") |>
      dplyr::rename(Group = Var1)
  })

  output$sex_summary <- renderTable({
    dat <- dataset()
    sex_col <- if ("inferred_sex" %in% names(dat$metadata)) {
      "inferred_sex"
    } else if ("sex" %in% names(dat$metadata)) {
      "sex"
    } else {
      NULL
    }
    validate(need(!is.null(sex_col), "No sex column found."))
    as.data.frame(table(dat$metadata[[sex_col]], useNA = "ifany"), responseName = "n") |>
      dplyr::rename(Sex = Var1)
  })

  output$metadata_table <- renderDT({
    datatable(round_for_display(dataset()$metadata), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$pca_controls <- renderUI({
    dat <- dataset()
    pca_df <- dat$pca[, !duplicated(names(dat$pca)), drop = FALSE]
    pcs <- pc_columns(pca_df)
    color_choices <- unique(c(
      dat$info$default_group_col,
      "inferred_sex",
      "rin",
      "epithelial_score",
      "mitochondrial_score",
      names(pca_df)
    ))
    color_choices <- color_choices[color_choices %in% names(pca_df)]
    if (length(color_choices) == 0) color_choices <- names(pca_df)
    shape_choices <- c("None" = "", categorical_columns(pca_df))
    label_choices <- c("None" = "", "Display sample" = dat$info$display_sample_col, "Patient" = dat$info$patient_id_col, names(pca_df))
    label_choices <- label_choices[label_choices == "" | label_choices %in% names(pca_df)]
    selected_color <- if (dat$info$default_group_col %in% color_choices) dat$info$default_group_col else color_choices[[1]]
    selected_shape <- if (dat$info$default_timepoint_col %in% unname(shape_choices)) dat$info$default_timepoint_col else ""

    tagList(
      selectInput("pca_x", "PC x-axis", choices = pcs, selected = pcs[[1]]),
      selectInput("pca_y", "PC y-axis", choices = pcs, selected = pcs[[2]]),
      selectInput("pca_color", "Color by", choices = color_choices, selected = selected_color),
      selectInput("pca_shape", "Shape by", choices = shape_choices, selected = selected_shape),
      selectInput("pca_label", "Label by", choices = label_choices, selected = "")
    )
  })

  output$pca_plot <- renderPlot({
    tryCatch(
      pca_figure(),
      error = function(e) {
        plot.new()
        text(0.5, 0.55, "PCA could not render for this dataset.", cex = 1.2)
        text(0.5, 0.45, conditionMessage(e), cex = 0.9)
      }
    )
  })

  output$loading_controls <- renderUI({
    dat <- dataset()
    tagList(
      selectInput("loading_pc", "PC number", choices = unique(dat$loadings$PC), selected = "PC1"),
      selectInput("loading_n", "Number of genes", choices = c(10, 20, 50), selected = 20),
      selectInput("loading_direction", "Loading direction", choices = c("absolute", "positive", "negative"), selected = "absolute")
    )
  })

  selected_loadings <- reactive({
    req(input$loading_pc, input$loading_n, input$loading_direction)
    top_loading_genes(dataset()$loadings, input$loading_pc, as.integer(input$loading_n), input$loading_direction)
  })

  loadings_figure <- reactive({
    plot_top_loadings(selected_loadings())
  })

  output$loadings_plot <- renderPlot({
    loadings_figure()
  })

  output$loadings_table <- renderDT({
    datatable(round_for_display(selected_loadings()), options = list(pageLength = 20, scrollX = TRUE))
  })

  output$download_loadings <- downloadHandler(
    filename = function() sanitize_filename(paste0(
      input$dataset_id, "_", input$loading_pc, "_top", input$loading_n,
      "_", input$loading_direction, "_loadings.csv"
    )),
    content = function(file) write.csv(selected_loadings(), file, row.names = FALSE)
  )

  output$deg_controls <- renderUI({
    dat <- dataset()
    contrasts <- if ("contrast" %in% names(dat$deg)) unique(dat$deg$contrast) else dat$info$default_contrast
    tagList(
      selectInput("deg_contrast", "Contrast", choices = contrasts, selected = dat$info$default_contrast),
      sliderInput("deg_padj", "Adjusted p-value cutoff", min = 0, max = 0.25, value = 0.05, step = 0.005),
      numericInput("deg_lfc", "Absolute log2FC cutoff", value = 1, min = 0, step = 0.25)
    )
  })

  filtered_deg <- reactive({
    req(input$deg_contrast, input$deg_padj, input$deg_lfc)
    filter_deg(dataset()$deg, input$deg_contrast, input$deg_padj, input$deg_lfc)
  })

  output$deg_up_count <- renderText(sum(filtered_deg()$direction == "Up", na.rm = TRUE))
  output$deg_down_count <- renderText(sum(filtered_deg()$direction == "Down", na.rm = TRUE))
  output$deg_total_count <- renderText(nrow(filtered_deg()))

  output$deg_table <- renderDT({
    datatable(round_for_display(filtered_deg()), options = list(pageLength = 15, scrollX = TRUE), filter = "top")
  })

  output$download_deg <- downloadHandler(
    filename = function() sanitize_filename(paste0(
      input$dataset_id, "_", input$deg_contrast,
      "_padj", cutoff_label(input$deg_padj),
      "_abslog2FC", cutoff_label(input$deg_lfc),
      "_DEG.csv"
    )),
    content = function(file) write.csv(filtered_deg(), file, row.names = FALSE)
  )

  output$volcano_controls <- renderUI({
    dat <- dataset()
    contrasts <- if ("contrast" %in% names(dat$deg)) unique(dat$deg$contrast) else dat$info$default_contrast
    tagList(
      selectInput("volcano_contrast", "Contrast", choices = contrasts, selected = dat$info$default_contrast),
      sliderInput("volcano_padj", "Adjusted p-value cutoff", min = 0, max = 0.25, value = 0.05, step = 0.005),
      numericInput("volcano_lfc", "Absolute log2FC cutoff", value = 1, min = 0, step = 0.25),
      helpText(if (requireNamespace("plotly", quietly = TRUE)) "Hover over points to see gene stats." else "Install plotly for hover/click interactivity.")
    )
  })

  output$volcano_ui <- renderUI({
    if (requireNamespace("plotly", quietly = TRUE)) {
      plotly::plotlyOutput("volcano_plot", height = "700px")
    } else {
      plotOutput("volcano_plot", height = "700px")
    }
  })

  volcano_data <- reactive({
    req(input$volcano_contrast)
    dat <- dataset()$deg
    if ("contrast" %in% names(dat)) {
      dat <- dat[dat$contrast == input$volcano_contrast, , drop = FALSE]
    }
    dat
  })

  if (requireNamespace("plotly", quietly = TRUE)) {
    output$volcano_plot <- plotly::renderPlotly({
      plot_volcano(volcano_data(), input$volcano_padj, input$volcano_lfc)
    })
  } else {
    output$volcano_plot <- renderPlot({
      plot_volcano(volcano_data(), input$volcano_padj, input$volcano_lfc)
    })
  }

  volcano_figure <- reactive({
    req(input$volcano_contrast, input$volcano_padj, input$volcano_lfc)
    plot_volcano_ggplot(volcano_data(), input$volcano_padj, input$volcano_lfc)
  })

  output$heatmap_controls <- renderUI({
    dat <- dataset()
    contrasts <- if ("contrast" %in% names(dat$deg)) unique(dat$deg$contrast) else dat$info$default_contrast
    annotation_choices <- intersect(
      c(dat$info$default_group_col, "inferred_sex", dat$info$default_timepoint_col, "rin"),
      names(dat$metadata)
    )
    tagList(
      radioButtons("heatmap_mode", "Genes", choices = c("Top DEGs" = "top", "Manual gene list" = "manual"), selected = "top"),
      selectInput("heatmap_contrast", "DEG contrast", choices = contrasts, selected = dat$info$default_contrast),
      selectInput("heatmap_top_n", "Top DEGs", choices = c(25, 50, 100), selected = 50),
      textAreaInput("heatmap_genes", "Manual genes", value = "OXTR\nAQP6\nSLC44A4", rows = 5),
      checkboxGroupInput("heatmap_annotations", "Sample annotations", choices = annotation_choices, selected = annotation_choices)
    )
  })

  heatmap_genes <- reactive({
    req(input$heatmap_mode)
    if (input$heatmap_mode == "manual") {
      parse_gene_list(input$heatmap_genes)
    } else {
      req(input$heatmap_contrast)
      top_deg_genes(dataset()$deg, as.integer(input$heatmap_top_n), input$heatmap_contrast)
    }
  })

  output$heatmap_plot <- renderPlot({
    dat <- dataset()
    plot_expression_heatmap(
      dat$vsd,
      dat$metadata,
      heatmap_genes(),
      input$heatmap_annotations,
      sample_col = dat$info$sample_id_col
    )
  })

  heatmap_figure <- reactive({
    dat <- dataset()
    plot_expression_heatmap(
      dat$vsd,
      dat$metadata,
      heatmap_genes(),
      input$heatmap_annotations,
      sample_col = dat$info$sample_id_col,
      silent = TRUE
    )
  })

  output$gsea_controls <- renderUI({
    dat <- dataset()
    contrasts <- if ("contrast" %in% names(dat$gsea)) unique(dat$gsea$contrast) else dat$info$default_contrast
    collections <- if ("collection" %in% names(dat$gsea)) unique(dat$gsea$collection) else "Hallmark"
    tagList(
      selectInput("pathway_contrast", "Contrast", choices = contrasts, selected = dat$info$default_contrast),
      selectInput("gsea_collection", "Gene-set collection", choices = collections, selected = "Hallmark"),
      selectInput("gsea_direction", "Enrichment direction", choices = c("Both", "Positive", "Negative"), selected = "Both"),
      selectInput("gsea_top_n", "Pathways in plot", choices = c(10, 20, 50), selected = 20)
    )
  })

  pathway_data <- reactive({
    req(input$pathway_contrast, input$gsea_collection)
    dat <- dataset()$gsea
    if ("contrast" %in% names(dat)) {
      dat <- dat[dat$contrast == input$pathway_contrast, , drop = FALSE]
    }
    if ("collection" %in% names(dat)) {
      dat <- dat[dat$collection == input$gsea_collection, , drop = FALSE]
    }
    dat
  })

  pathway_figure <- reactive({
    req(input$pathway_contrast, input$pathway_padj, input$gsea_direction, input$gsea_top_n)
    plot_pathway_nes_ggplot(
      pathway_data(), input$pathway_padj, input$gsea_direction, as.integer(input$gsea_top_n)
    )
  })

  filtered_gsea <- reactive({
    req(input$pathway_padj, input$gsea_direction)
    out <- add_pathway_overlap_stats(pathway_data())
    out <- out[!is.na(out$padj) & out$padj <= input$pathway_padj, , drop = FALSE]
    if (input$gsea_direction == "Positive") {
      out <- out[out$NES > 0, , drop = FALSE]
    } else if (input$gsea_direction == "Negative") {
      out <- out[out$NES < 0, , drop = FALSE]
    }
    out[order(out$padj, -abs(out$NES)), , drop = FALSE]
  })

  output$pathway_plot <- plotly::renderPlotly({
    plot_pathway_nes(
      pathway_data(), input$pathway_padj, input$gsea_direction, as.integer(input$gsea_top_n)
    )
  })

  output$gsea_table <- renderDT({
    datatable(
      round_for_display(filtered_gsea()),
      options = list(pageLength = 15, scrollX = TRUE),
      filter = "top"
    )
  })

  output$download_gsea <- downloadHandler(
    filename = function() {
      sanitize_filename(paste0(
        input$dataset_id, "_", input$pathway_contrast,
        "_", input$gsea_collection,
        "_", input$gsea_direction,
        "_GSEA_padj", cutoff_label(input$pathway_padj), ".csv"
      ))
    },
    content = function(file) {
      write.csv(filtered_gsea(), file, row.names = FALSE)
    }
  )

  output$gene_group_control <- renderUI({
    dat <- dataset()
    group_choices <- intersect(
      c(dat$info$default_group_col, "inferred_sex", dat$info$default_timepoint_col, categorical_columns(dat$metadata)),
      names(dat$metadata)
    )
    selectInput("gene_group", "Expression grouped by", choices = unique(group_choices), selected = dat$info$default_group_col)
  })

  gene_query <- reactive({
    req(input$gene_query)
    toupper(trimws(input$gene_query))
  })

  gene_expression_figure <- reactive({
    req(input$gene_group)
    dat <- dataset()
    expr <- expression_long(dat$vsd, dat$metadata, gene_query(), sample_col = dat$info$sample_id_col)
    plot_gene_expression(expr, input$gene_group)
  })

  output$gene_deg_table <- renderDT({
    dat <- dataset()
    out <- dat$deg[toupper(dat$deg$gene) == gene_query(), , drop = FALSE]
    datatable(round_for_display(out), options = list(dom = "t", scrollX = TRUE))
  })

  output$gene_expression_plot <- renderPlot({
    gene_expression_figure()
  })

  output$gene_pathway_table <- renderDT({
    dat <- dataset()
    matches <- add_pathway_overlap_stats(leading_edge_matches(dat$gsea, gene_query()))
    datatable(round_for_display(matches), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$enrichr_controls <- renderUI({
    dat <- dataset()
    contrasts <- if ("contrast" %in% names(dat$deg)) unique(dat$deg$contrast) else dat$info$default_contrast
    tagList(
      selectInput("enrichr_contrast", "Contrast", choices = contrasts, selected = dat$info$default_contrast),
      selectInput("enrichr_direction", "Gene direction", choices = c("Up", "Down", "Both"), selected = "Up"),
      sliderInput("enrichr_padj", "DEG adjusted p-value cutoff", min = 0, max = 0.25, value = 0.05, step = 0.005),
      numericInput("enrichr_lfc", "Absolute log2FC cutoff", value = 1, min = 0, step = 0.25),
      selectInput(
        "enrichr_database",
        "Enrichr library",
        choices = c(
          "GO Biological Process 2025" = "GO_Biological_Process_2025",
          "GO Molecular Function 2025" = "GO_Molecular_Function_2025",
          "GO Cellular Component 2025" = "GO_Cellular_Component_2025",
          "Reactome Pathways 2024" = "Reactome_Pathways_2024",
          "MSigDB Hallmark 2020" = "MSigDB_Hallmark_2020",
          "Drug perturbations from GEO" = "Drug_Perturbations_from_GEO_2014"
        )
      ),
      numericInput("enrichr_max_genes", "Maximum submitted genes", value = 500, min = 10, max = 2000, step = 50)
    )
  })

  enrichr_signature <- reactive({
    req(
      input$enrichr_contrast, input$enrichr_direction, input$enrichr_padj,
      input$enrichr_lfc, input$enrichr_max_genes
    )
    select_signature_genes(
      dataset()$deg,
      input$enrichr_contrast,
      input$enrichr_padj,
      input$enrichr_lfc,
      input$enrichr_direction,
      input$enrichr_max_genes
    )
  })

  enrichr_results <- eventReactive(input$run_enrichr, {
    req(input$enrichr_database)
    genes <- enrichr_signature()$gene
    validate(need(length(genes) >= 5, "At least five genes pass the selected thresholds."))
    withProgress(message = "Querying Enrichr", value = 0.5, {
      result <- tryCatch(
        enrichR::enrichr(genes, input$enrichr_database)[[input$enrichr_database]],
        error = function(e) e
      )
    })
    if (inherits(result, "error")) {
      validate(need(
        FALSE,
        paste(
          "Enrichr request failed. Check the internet connection or try again later:",
          conditionMessage(result)
        )
      ))
    }
    result
  }, ignoreInit = TRUE)

  output$enrichr_plot <- renderPlot({
    out <- enrichr_results()
    out <- out[!is.na(out$Adjusted.P.value) & out$Adjusted.P.value <= 0.1, , drop = FALSE]
    enrichment_barplot(
      out,
      "Term",
      "Adjusted.P.value",
      score_col = "Combined.Score",
      title = input$enrichr_database
    )
  })

  output$enrichr_table <- renderDT({
    datatable(
      round_for_display(enrichr_results()),
      options = list(pageLength = 15, scrollX = TRUE),
      filter = "top"
    )
  })

  output$download_enrichr <- downloadHandler(
    filename = function() sanitize_filename(paste0(
      input$dataset_id, "_", input$enrichr_contrast, "_", input$enrichr_direction,
      "_", input$enrichr_database, "_Enrichr.csv"
    )),
    content = function(file) write.csv(enrichr_results(), file, row.names = FALSE)
  )

  output$clusterprofiler_controls <- renderUI({
    dat <- dataset()
    contrasts <- if ("contrast" %in% names(dat$deg)) unique(dat$deg$contrast) else dat$info$default_contrast
    tagList(
      selectInput("cp_contrast", "Contrast", choices = contrasts, selected = dat$info$default_contrast),
      selectInput("cp_direction", "Gene direction", choices = c("Up", "Down", "Both"), selected = "Up"),
      sliderInput("cp_deg_padj", "DEG adjusted p-value cutoff", min = 0, max = 0.25, value = 0.05, step = 0.005),
      numericInput("cp_lfc", "Absolute log2FC cutoff", value = 1, min = 0, step = 0.25),
      selectInput(
        "cp_ontology",
        "Gene Ontology domain",
        choices = c(
          "Biological Process" = "BP",
          "Molecular Function" = "MF",
          "Cellular Component" = "CC"
        )
      ),
      sliderInput("cp_result_padj", "GO adjusted p-value cutoff", min = 0, max = 0.25, value = 0.05, step = 0.005),
      numericInput("cp_max_genes", "Maximum selected genes", value = 1000, min = 10, max = 5000, step = 100)
    )
  })

  clusterprofiler_signature <- reactive({
    req(input$cp_contrast, input$cp_direction, input$cp_deg_padj, input$cp_lfc, input$cp_max_genes)
    select_signature_genes(
      dataset()$deg,
      input$cp_contrast,
      input$cp_deg_padj,
      input$cp_lfc,
      input$cp_direction,
      input$cp_max_genes
    )
  })

  clusterprofiler_results <- eventReactive(input$run_clusterprofiler, {
    req(input$cp_ontology, input$cp_result_padj)
    withProgress(message = "Running clusterProfiler", value = 0.5, {
      run_clusterprofiler_go(
        clusterprofiler_signature(),
        unique(dataset()$deg$gene),
        input$cp_ontology,
        input$cp_result_padj
      )
    })
  }, ignoreInit = TRUE)

  output$clusterprofiler_plot <- renderPlot({
    enrichment_barplot(
      clusterprofiler_results(),
      "Description",
      "p.adjust",
      title = paste("Gene Ontology", input$cp_ontology)
    )
  })

  output$clusterprofiler_table <- renderDT({
    datatable(
      round_for_display(clusterprofiler_results()),
      options = list(pageLength = 15, scrollX = TRUE),
      filter = "top"
    )
  })

  output$download_clusterprofiler <- downloadHandler(
    filename = function() sanitize_filename(paste0(
      input$dataset_id, "_", input$cp_contrast, "_", input$cp_direction,
      "_GO_", input$cp_ontology, "_padj", cutoff_label(input$cp_result_padj), ".csv"
    )),
    content = function(file) write.csv(clusterprofiler_results(), file, row.names = FALSE)
  )

  output$download_figures_pdf <- downloadHandler(
    filename = function() {
      sanitize_filename(paste0(
        input$dataset_id, "_current_figures_",
        format(Sys.Date(), "%Y-%m-%d"), ".pdf"
      ))
    },
    content = function(file) {
      grDevices::pdf(file, width = 11, height = 8.5, onefile = TRUE)
      on.exit(grDevices::dev.off(), add = TRUE)

      graphics::plot.new()
      graphics::text(0.5, 0.62, dataset()$info$display_name, cex = 1.8, font = 2)
      graphics::text(0.5, 0.52, "RNA-seq Explorer: current figure selections", cex = 1.3)
      graphics::text(0.5, 0.44, paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")), cex = 1)

      print(pca_figure())
      print(loadings_figure())
      print(volcano_figure())
      grid::grid.newpage()
      grid::grid.draw(heatmap_figure()$gtable)
      print(pathway_figure())
      print(gene_expression_figure())
    }
  )
}

shinyApp(ui, server)
