# SeqSurf AI

SeqSurf AI is a Shiny prototype for AI-assisted transcriptomic reanalysis.
It starts from a GEO accession or custom study notes, builds study context from
GEO/PubMed, PubMed papers that mention the GEO accession, and optional paper
text, proposes a reproducible analysis design, writes an app-ready dataset, and
then lets the user explore and verify results.

This app is intentionally separate from `RNAseq_Shiny_App/`, which remains the
research dashboard for existing osteosarcoma, neuroblastoma, and ANBL datasets.

## What SeqSurf Does

SeqSurf is designed to answer four questions before and after reanalysis:

1. Can this public dataset be analyzed in the app?
2. What has already been done with this dataset?
3. How well does a public reanalysis reproduce or approximate the paper?
4. What should be analyzed next?

For studies that cannot run fully inside Shiny, SeqSurf exports a reproducible
external workflow bundle and a no-code AI prompt so a user can ask Claude Code,
Codex, Positron Assistant, or another coding assistant to run the raw/count
workflow and return import-ready files.

## Product Workflow

1. Start
   - Enter a GEO accession or paste custom/private dataset notes.
   - Retrieve GEO/PubMed study context.
   - Review what has already been done, detected methods/results, prior reuse
     clues, limitations, exact-reproduction feasibility, and raw/supplementary
     file triage.
   - Classify whether the study can run in-app, needs a count-matrix pipeline,
     needs a raw sequencing pipeline, or requires manual curation.
   - Send the study to Prepare Evidence or Validate.

2. Prepare Evidence
   - Build study memory from GEO/PubMed, related papers mentioning the GEO
     accession, and optional user-provided text.
   - Fetch paper-like PDFs/HTML when available.
   - Extract structured publication claims from the study context.
   - Export a paper-to-code brief.

3. Validate
   - Fetch a GEO ExpressionSet when a processed series matrix is available.
   - Propose a primary contrast and analysis design.
   - Run PCA, DEG, GSEA, and pathway-score generation.
   - Write the app dataset contract under `data/geo_<accession>/`.
   - Verify whether the imported dataset matches the trained GEO accession.
   - Compare reproduced DEG/GSEA results against paper claims.
   - Score post-analysis reproducibility from supported, weak, missing, and
     direction-mismatched claims.

4. Discovery
   - Rank next analyses using the active dataset, DEG/GSEA results, metadata,
     trained study memory, and paper-agreement status.
   - Triage reproduced findings, divergent findings, weak/missing claims, and
     high-signal new DEG/pathway results.
   - Report gaps, suggested comparisons, and follow-up analysis ideas.
   - Explore dataset overview, PCA, top loading genes, DEG, volcano, heatmap,
     GSEA, Enrichr, clusterProfiler, and gene search.

The paper-to-code brief is designed for Claude Code, Codex, or another coding
agent. It includes section-aware paper excerpts, the analyzability route,
exact-method feasibility, raw/supplementary triage, required app-contract files,
and acceptance criteria for generated analysis scripts.

6. No-code external handoff
   - Download the raw/SRA workflow bundle.
   - Download the no-code AI guide and AI assistant prompt.
   - Run the external workflow with help from an AI coding assistant.
   - Upload `outputs/seqsurf_gene_counts.csv` and
     `outputs/seqsurf_metadata.csv` back into SeqSurf.

5. Upload mode on Start
   - Upload an expression matrix and metadata table from the Start page.
   - Map gene, sample, group, patient, and display columns.
   - Create an app-ready dataset contract for PCA, DEG, volcano, heatmap, and
     gene-level exploration.
   - Uploaded DEG currently uses a simple Welch t-test baseline suitable for
     demos/prototypes; publishable analyses should replace this with a
     study-appropriate model.

## Current Scope

The app currently works best for GEO studies that expose a usable processed
series-matrix expression table through `GEOquery::getGEO(..., GSEMatrix = TRUE)`.
Raw FASTQ/count-matrix reconstruction is handled as a triage and pipeline-handoff
plan rather than an in-Shiny FASTQ workflow. The upload workflow supports
app-ready expression matrices plus metadata.

SeqSurf AI includes an analyzability classifier so the app can say whether a GEO
study is currently suitable for in-app processed-matrix reanalysis, partially
supported through a count-matrix workflow, or outside the Shiny runtime because
it requires raw FASTQ/SRA/BAM-style processing.

The first version of the scoring layer separates:

- reanalysis usefulness score
- data availability score
- reproducibility score
- novelty/opportunity score
- technical risk score

These scores are deterministic heuristics from GEO/PubMed metadata, extracted
paper methods/findings, optional user-provided paper notes, sample count,
supplementary-file availability, reported claims, and platform/method language.

The Start page also includes an optional LLM narrative layer. When an OpenAI API
key is available, the app can write a concise scientist-facing interpretation of
what has already been done with the dataset, whether reanalysis is worth doing,
what to reproduce first, what modern follow-up analysis is most promising, and
what manual checks are needed before proceeding.

## Run

From this folder:

```r
shiny::runApp()
```

Or navigate to the cloned directory first:

```r
setwd("path/to/SeqSurf_AI")
shiny::runApp()
```

Recommended local run:

```r
shiny::runApp(host = "127.0.0.1", port = 4248)
```

## Optional Dependencies

Core app packages are checked at startup. GEO/AI features also benefit from:

- `GEOquery`
- `Biobase`
- `rentrez`
- `pdftools`
- `httr2`
- `jsonlite`
- `limma`
- `edgeR`
- `matrixStats`
- `fgsea`
- `msigdbr`
- `survival`
- `immunedeconv`
- `commonmark`

OpenAI-backed extraction/chat requires an API key pasted into the app or set as
`OPENAI_API_KEY`. In the app, paste the key once in the global sidebar and click
`Save key for session`; SeqSurf AI keeps it in Shiny session memory for Start,
Prepare Evidence, Validate, and Discovery actions. The key is not written to
project files.

## Installation Notes

For a clean demo machine, install R first from CRAN:

- https://cran.r-project.org/

Then install an IDE. Positron is a good choice for users who want an AI-assisted
R/Python workflow:

- https://positron.posit.co/download.html

For no-code external workflow help, use an AI coding assistant that can inspect a
folder and run terminal/R commands. SeqSurf exports prompts for Claude Code,
Codex, Positron Assistant, or similar tools.

## Dataset Contract

Imported or uploaded datasets should live under `data/<dataset_id>/` and contain:

- `metadata.rds`
- `vsd_matrix.rds`
- `pca_df.rds`
- `pca_variance.rds`
- `pca_loadings.rds`
- `deg_results.csv`
- `gsea_hallmark.csv`
- `pathway_scores.csv`
- `dataset_info.rds`

The expression matrix should have genes as rows and sample IDs as columns.
`metadata.rds` should include the sample ID column named in `dataset_info.rds`.

## Relationship To The Original App

`RNAseq_Shiny_App/` keeps the existing precomputed research datasets and
dataset-specific study pages. `SeqSurf_AI/` starts clean, with no bundled
personal datasets, so it can become the general GEO/upload reanalysis product.

## Demo Notes

For a demo or abstract, present automated GEO import as strongest for processed
series-matrix datasets. For raw/count-only studies, use the raw/supplementary
triage and paper-to-code brief as the honest handoff into a reproducible external
pipeline, then load the resulting app contract back into SeqSurf AI.

Recommended demo materials:

- `docs/presentation_demo_script.md`
- `docs/demo_walkthrough_GSE19804.md`
- `docs/demo_dataset_matrix.md`

Recommended live demo accession:

- `GSE19804`

Route examples:

- Processed matrix: `GSE19804`
- Count-matrix path: `GSE60450`
- Raw/SRA handoff: `GSE16476`

## Current Limitations

- Full-text/PDF extraction is section-aware but not perfect for every PDF layout.
- Exact paper reproduction is only possible when public files and methods are
  sufficiently complete.
- Raw FASTQ/SRA processing is intentionally exported as an external workflow
  bundle instead of running inside Shiny.
- PaperQA support currently prepares a corpus/manifest; full citation-grounded
  chat remains a future integration target.
