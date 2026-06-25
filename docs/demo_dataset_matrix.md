# SeqSurf AI Demo Dataset Matrix

Use this matrix to test the three major product paths before a presentation or paper demo.

| Route | GEO | Expected Classification | Demo Purpose | What To Show |
|---|---:|---|---|---|
| Processed matrix | `GSE19804` | In-app GEO processed-matrix reanalysis | Main live demo | Start assessment, Train memory, Validate reanalysis, Score, Discovery, Report |
| Count matrix | `GSE60450` | Count-matrix pipeline then app import | Shows semi-automated count support | Count candidate detection, count-matrix plan, generated count pipeline script |
| Raw/SRA handoff | `GSE16476` | Raw sequencing pipeline handoff | Shows honest external workflow route | Analyzability “No”, raw/SRA workflow bundle, no-code AI prompt |

## Acceptance Checks

### `GSE19804`

- Start page loads title, organism, assay, sample count, and paper context.
- Analyzability says the dataset can run in-app.
- Validate can produce an app dataset.
- Score and report tabs produce a reproducibility report.
- Discovery shows ranked next analyses.

### `GSE60450`

- Start page identifies count-matrix or partial route.
- Validate -> Count Files finds likely count/table candidates.
- Run count-matrix import either succeeds or clearly explains why candidate parsing failed.
- Generated script is available from Validate -> Count Files.

### `GSE16476`

- Start page identifies raw/SRA or external handoff route.
- No-code AI handoff panel is available in Start -> Plan.
- Raw/SRA workflow bundle downloads successfully.
- Prompt for AI assistant explains how to generate `outputs/seqsurf_gene_counts.csv` and `outputs/seqsurf_metadata.csv`.

## Demo Reliability Rule

For the live presentation, use `GSE19804` as the main path. Use `GSE60450` and `GSE16476` as route-classification examples unless you have already pre-run the external/count workflow.
