# SeqSurf AI Presentation Demo Script

Use `GSE19804` for the main live demo because it follows the cleanest path: GEO accession, processed matrix, in-app reanalysis, paper-agreement scoring, Discovery, and report export.

## Slide Story

1. Problem: public transcriptomic datasets are valuable but hard to reanalyze reproducibly.
2. SeqSurf input: a GEO accession plus optional paper notes.
3. SeqSurf decision: in-app analysis, count-matrix workflow, raw/SRA handoff, or manual curation.
4. Study Memory 2.0: paper-aware memory with claims, methods, limitations, and citations.
5. Validate: reanalyze and score agreement with the original paper.
6. Discovery: rank next analyses and identify what is worth doing next.
7. Export: reproducibility report plus no-code AI handoff for datasets requiring external workflows.

## Live Click Path

### 1. Start

Screenshot target: `docs/screenshots/01_start_assessment.png`

- Open SeqSurf AI.
- Go to **Start**.
- Enter `GSE19804`.
- Click **Assess GEO dataset**.
- Narration: “SeqSurf first asks whether this dataset can be analyzed in the app, should use a count-matrix workflow, or needs raw-data handoff.”
- Show:
  - SeqSurfer Recommendation
  - Can SeqSurf Analyze This In-App?
  - Dataset Snapshot

### 2. Train

Screenshot target: `docs/screenshots/02_train_memory.png`

- Click **Continue to Train** or open **Train**.
- Click **1. Train from original paper** if needed.
- Click **2. Fetch evidence text**.
- Click **3. Build Study Memory 2.0**.
- Open **Virtual Lab**.
- Narration: “This is the intelligence layer. It separates literature curation, methods extraction, claim extraction, reproducibility judgment, and discovery planning.”

### 3. Validate

Screenshot target: `docs/screenshots/03_validate_workflow.png`

- Click **Go to Validate**.
- Click **1. Fetch GEO**.
- Review **Design**.
- Click **2. Propose design** if needed.
- Click **3. Run reanalysis**.
- Click **4. Validate paper agreement**.
- Open **Score**.
- Narration: “The app does not just make plots. It asks whether the public reanalysis supports or diverges from paper claims.”

### 4. Discovery

Screenshot target: `docs/screenshots/04_discovery_ranked.png`

- Open **Discovery**.
- Show **Ranked Next Analyses**.
- Ask: “Based only on the available metadata and app results, what are the three most useful next analyses?”
- Narration: “Discovery turns the validation result into a prioritized next-step plan.”

### 5. Export

Screenshot target: `docs/screenshots/05_report_export.png`

- Return to **Validate -> Report**.
- Click **Download reproducibility report**.
- Narration: “The report is the artifact for a meeting, abstract, paper supplement, or GitHub demo.”

## Backup Datasets

- Processed matrix demo: `GSE19804`
- Count-matrix route demo: `GSE60450`
- Raw/SRA handoff demo: `GSE16476`

If network retrieval fails during a live meeting, use saved `geo_gse19804` results and show the exported report plus Discovery cards.

## Key Lines To Say

- “SeqSurf is deliberately honest about analyzability. It does not pretend arbitrary raw FASTQ processing belongs inside Shiny.”
- “The app first reproduces or approximates the paper, then recommends more modern analyses.”
- “For non-coders, SeqSurf exports an AI-assisted workflow bundle and prompt so Claude/Positron/Codex can help run the external workflow.”
- “The final output is not only a plot. It is a reproducibility score, agreement table, discovery plan, and report.”
