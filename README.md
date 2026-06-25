# SeqSurf AI

[![R](https://img.shields.io/badge/R-Shiny-276DC3)](https://shiny.posit.co/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Release](https://img.shields.io/github/v/release/selahgmckenney/SeqSurfAI?include_prereleases)](https://github.com/selahgmckenney/SeqSurfAI/releases)
[![Status](https://img.shields.io/badge/status-research%20prototype-blue)](#current-scope)

**SeqSurf AI is a Shiny application for paper-grounded transcriptomic reanalysis.**
It helps researchers start from a GEO accession, understand what has already
been done with a public dataset, classify whether reanalysis is feasible,
validate app results against publication claims, and identify scientifically
useful next analyses.

Public transcriptomic datasets are reused constantly, but conclusions are often
difficult to compare because preprocessing choices, phenotype definitions,
contrasts, and validation standards vary across studies. SeqSurf AI is designed
as a reproducibility layer between public repositories, papers, and exploratory
bioinformatics workflows.

![SeqSurf AI interface](docs/screenshots/00_current_app_view.png)

## Why This Exists

| Problem | SeqSurf AI response |
|---|---|
| Public GEO studies are reused inconsistently. | Classifies each dataset into an analyzable route before reanalysis. |
| AI tools can summarize papers but do not enforce analysis contracts. | Connects GEO metadata, paper claims, methods, app results, and reproducibility scores. |
| Exact reproduction is often impossible from public files alone. | Separates exact reproduction, paper-style approximation, count-matrix import, and raw/SRA handoff. |
| Reanalysis value is hard to judge up front. | Scores data availability, reproducibility potential, usefulness, novelty, and technical risk. |
| Cross-study comparison is usually manual. | Provides visual comparison of DEG overlap, pathway agreement, and claimed biomarkers across prepared studies. |

## Core Capabilities

SeqSurf AI currently supports:

- **GEO triage:** enter a GEO accession and classify the study as in-app
  processed-matrix reanalysis, count-matrix workflow, raw/SRA handoff, or manual
  curation.
- **Study memory:** build local study context from GEO, PubMed, related papers,
  optional notes, and extracted claims.
- **Paper-aware scoring:** estimate whether reanalysis is useful before spending
  time on deeper analysis.
- **Processed-matrix validation:** generate or load an app-ready dataset
  contract with metadata, expression matrix, PCA, DEG results, GSEA/pathway
  summaries, and report files.
- **Publication claim comparison:** compare paper-reported biomarkers and
  pathways with app DEG/GSEA results.
- **Reproducibility score:** summarize supported, weak, missing, and
  direction-mismatched claims.
- **Discovery recommendations:** rank next analyses using metadata, DEG/GSEA
  signals, paper agreement, and learned local feedback.
- **Compare Studies:** visually compare prepared studies by DEG overlap,
  pathway direction agreement, and shared claimed biomarkers.
- **External workflow handoff:** export raw/SRA or count-matrix plans for
  Positron, Codex, Claude Code, or another coding assistant.
- **Learning export:** save local feedback and export structured examples for
  future SeqSurf-specific model development.

## App Workflow

1. **Start**
   - Enter a GEO accession or load the local demo.
   - Review data quality signals, analyzability route, prior work, and
     reanalysis scores.

2. **Prepare Evidence**
   - Build study memory from GEO/PubMed and optional paper or method notes.
   - Extract structured claims, methods, and citation evidence.

3. **Validate**
   - Fetch or load a prepared dataset.
   - Propose a contrast, run/reuse PCA, DEG, GSEA, pathway scores, and report
     generation.
   - Compare results to publication claims and calculate a reproducibility score.

4. **Discovery**
   - Identify relevant, different, weak, missing, and high-signal results.
   - Surface sequencing pitfalls and standardized signature checks.
   - Save feedback so future recommendations can adapt locally.

5. **Compare Studies**
   - Compare two completed app-contract datasets.
   - View DEG overlap direction, pathway NES agreement, and claimed biomarker
     presence as figures rather than tables.

## Quick Start

Clone the repository:

```bash
git clone https://github.com/selahgmckenney/SeqSurfAI.git
cd SeqSurfAI
```

Install core R packages:

```r
install.packages(c(
  "shiny", "bslib", "ggplot2", "dplyr", "tidyr", "DT",
  "plotly", "pheatmap", "httr2", "jsonlite"
))
```

Run the app:

```r
shiny::runApp(host = "127.0.0.1", port = 4248)
```

Optional Bioconductor/GEO features use packages such as `GEOquery`, `Biobase`,
`limma`, `edgeR`, `fgsea`, `msigdbr`, `AnnotationDbi`, and `org.Hs.eg.db`.

## Fast Local Demo

For a presentation, use the built-in demo path:

1. Open the app.
2. Go to **Start**.
3. Click **Load local demo**.
4. Show **Data Quality Signals**, analyzability route, and the study snapshot.
5. Move to **Validate**, **Discovery**, **PCA Explorer**, and **Compare Studies**.

The local demo is designed to avoid slow live downloads. It prepares:

- `GSE19804` as the primary processed-matrix demo.
- `GSE16476` and `GSE60450` as lightweight comparison fixtures for the
  **Compare Studies** tab.

The comparison fixtures can also be recreated directly:

```r
source("scripts/create_compare_demo_datasets.R")
```

## Dataset Contract

Imported or uploaded datasets live under `data/<dataset_id>/` and follow a
standard app contract:

```text
metadata.rds
vsd_matrix.rds
pca_df.rds
pca_variance.rds
pca_loadings.rds
deg_results.csv
gsea_hallmark.csv
pathway_scores.csv
dataset_info.rds
```

This contract allows the same PCA, DEG, volcano, heatmap, GSEA, gene search,
paper-agreement, discovery, and cross-study comparison views to work across
different studies.

## AI And Privacy Notes

OpenAI-backed extraction and chat are optional. If an API key is provided, it is
held in Shiny session memory and is not written to project files. Local learning
logs, generated workflow folders, downloaded papers, and app-generated GEO data
are ignored by git by default.

Shared learning export is opt-in and intended to omit GEO accession, dataset ID,
raw data, paper text, and free-text notes.

## Current Scope

SeqSurf AI is a research prototype. It is strongest for GEO studies with usable
processed series-matrix expression data. Raw FASTQ/SRA processing is not run
inside Shiny; instead, SeqSurf AI exports workflow plans and handoff prompts so
external tools can generate import-ready files.

Current limitations:

- Exact paper reproduction depends on public file completeness and method detail.
- Full-text/PDF extraction is section-aware but not perfect for every article.
- Upload-mode DEG uses a simple prototype baseline and should be replaced with a
  study-appropriate model for publication-grade analysis.
- Shared/global learning is file-based; a consent-aware backend is future work.
- The SeqSurf-specific ML model is planned, but this repository currently exports
  structured training examples rather than fine-tuning a model in-app.

## Documentation

- [Demo walkthrough: GSE19804](docs/demo_walkthrough_GSE19804.md)
- [Demo dataset matrix](docs/demo_dataset_matrix.md)
- [Presentation demo script](docs/presentation_demo_script.md)
- [Healthcare AI overview](docs/SeqSurf_AI.qmd)

## Citation

If you use SeqSurf AI, cite the repository release used for your analysis.

```text
McKenney S. SeqSurf AI: paper-grounded transcriptomic reanalysis for GEO
dataset triage, claim validation, reproducibility scoring, and adaptive
learning. GitHub: https://github.com/selahgmckenney/SeqSurfAI
```

For manuscript text:

```text
Source code for SeqSurf AI is available at
https://github.com/selahgmckenney/SeqSurfAI.
```

