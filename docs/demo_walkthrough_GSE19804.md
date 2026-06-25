# SeqSurf AI Demo Walkthrough: GSE19804

Use this dataset for a live demo because SeqSurf can score it, import processed expression data, run a paper-style reanalysis, compare claims, and export a reproducibility report inside the app.

## Dataset

- GEO accession: `GSE19804`
- Title: Genome-wide screening of transcriptional modulation in non-smoking female lung cancer in Taiwan
- Organism: Homo sapiens
- Assay: Microarray / expression array
- Samples: 120
- Demo route: In-app GEO processed-matrix reanalysis

## Story To Tell

SeqSurf starts from a GEO accession, asks whether reanalysis is feasible, builds a study memory from paper evidence, reanalyzes the public expression matrix, compares results with published claims, then proposes modern follow-up analyses.

## Click Script

1. Open SeqSurf AI and go to **Start**.
2. Enter `GSE19804`.
3. Click **Assess GEO dataset**.
4. Point out the **Can SeqSurf Analyze This In-App?** card.
   - Expected answer: `Yes`
   - Expected route: `In-app GEO processed-matrix reanalysis`
5. Review the **Dataset Snapshot**.
   - Emphasize organism, assay type, sample count, metadata, expression files, and paper context.
6. Go to **Train**.
7. Confirm `GSE19804` carried into the GEO field.
8. Click **1. Train from original paper** or **Fetch GEO / evidence text** if the local state is blank.
9. Click **Fetch paper PDFs** if available.
10. Click **Build Study Memory 2.0**.
11. Open the **Virtual Lab** tab.
12. Show:
   - Literature Curator
   - Methods Extractor
   - Claims Extractor
   - Reproducibility Judge
   - Discovery Scientist
   - Citation Evidence
13. Go to **Validate**.
14. Click **1. Fetch GEO**.
15. Review **Design** and confirm the group/contrast.
16. Click **2. Propose design** if needed.
17. Click **3. Run reanalysis**.
18. Open **Paper Agreement** and show claim agreement/disagreement.
19. Open **Score** and show the post-analysis reproducibility score.
20. Open **Report** and click **Download reproducibility report**.
21. Go to **Discovery**.
22. Show modern follow-up ideas, suggested new comparisons, and gap analysis.

## Main Demo Points

- SeqSurf does not pretend every GEO dataset is immediately analyzable.
- It classifies the route first: in-app processed matrix, count-matrix pipeline, or raw FASTQ/SRA handoff.
- The study memory is not just a chat summary; it has claims, source text, related papers, agent reviews, and citation evidence.
- Validation is separate from discovery: first reproduce or approximate the paper, then ask what should be done next.
- The exported report is the artifact for a meeting, abstract, paper supplement, or GitHub demo.

## Backup Line If A Network Step Fails

GEO/PubMed retrieval depends on network access. If live retrieval fails, select the saved `geo_gse19804` dataset if it exists, then show the prepared report, score, and explorer tabs.
