# wk_260615 — TMAO Pathway Analysis: Weekly Working Directory

**Week of:** 2026-06-15  
**Handed off from:** `TMAO_pathway_analysis/`  
**Purpose:** Re-run and extend the TMAO pathway analysis pipeline from scratch. All input data is pre-prepared; just run the scripts in order.

---

## Background (read this first)

This project investigates how the gut microbiome metabolite **TMAO (trimethylamine N-oxide)** influences plasma metabolomics. The experiment has two perturbation periods applied sequentially to the same subjects:

- **Perturbation X ("Early period"):** Timepoints T1 → T2 → T3, where T1 is the baseline.
- **Perturbation Y ("Late period"):** Timepoints T4 → T5 → T6, where T4 is the baseline.

The analysis runs in two stages:

1. **Upstream (`src/data_analysis_r/`)** — Starts from raw MS-DIAL feature tables, runs statistics, PLS-DA, and classifies both subjects and metabolites. Generates the CSV files used downstream.
2. **Downstream** — Two independent pipelines that consume the upstream outputs:
   - **`src/pathway_analysis_r/`** — Pathway ORA, KEGG analysis, Cytoscape export, bile acid × TMAO correlation.
   - **`src/python/`** — Ontology enrichment and trajectory visualization.

---

## Directory Structure

```
wk_260615/
├── src/
│   ├── data_analysis_r/          # ← UPSTREAM: raw data → stats → PLS-DA → classification
│   │   ├── load_data.R                     # Load dummy data into R environment
│   │   ├── load_data_production.R          # Load production data into R environment
│   │   ├── boxplot_pairwise_dummy.R        # Paired Wilcoxon tests + boxplots (dummy)
│   │   ├── boxplot_pairwise_production.R   # Paired Wilcoxon tests + boxplots (production)
│   │   ├── enrich_pairwise_results.R       # Add HMDB/MS-FINDER annotation to stats output
│   │   ├── classify_subjects_by_response.R # Group subjects by TMAO response pattern
│   │   ├── pca_per_timepoint.R             # PCA score plots per timepoint
│   │   ├── plot_timeseries.R               # Time-series plots (dummy data)
│   │   ├── plot_timeseries_production.R    # Time-series plots (production data)
│   │   ├── plot_selected_features.R        # Boxplots for specific Alignment_IDs
│   │   ├── pls_alignment414.R              # PLS with TMAO (ID 414) as response
│   │   ├── pls_alignment414_per_timepoint.R# Same PLS, each subject×TP as observation
│   │   ├── plsda_per_timepoint.R           # PLS-DA score plots per timepoint
│   │   ├── plsda_loadings_vip.R            # VIP scores + loading plots from PLS-DA
│   │   ├── plsda_mechanism_candidates.R    # Classify metabolites into mechanism classes
│   │   └── classify_metabolites_metaboanalystr.R # Alternative classification for MetaboAnalystR
│   │
│   ├── pathway_analysis_r/       # ← DOWNSTREAM: pathway ORA + KEGG + Cytoscape + H1-A
│   │   ├── 01_preprocess.R
│   │   ├── 02_enrichment_ora.R
│   │   ├── 03_pathway_analysis.R
│   │   ├── 04_cytoscape_export.R
│   │   ├── 05_H1A_scatter.R
│   │   ├── 06_group_composite_scatter.R
│   │   └── 07_bothup_report.R
│   │
│   └── python/                   # ← DOWNSTREAM: ontology enrichment + trajectories
│       ├── 01_ontology_enrichment.py
│       ├── 02_make_excel.py
│       └── 04_trajectory_analysis.py
│
├── data/
│   ├── production/               # Production data (excluded from Git)
│   │   ├── raw/                  #   Raw MS-DIAL output (df_conf.csv etc.)
│   │   └── processed/            #   Processed/intermediate files
│   └── ref/                      # Reference files (small, manually curated)
│       └── tmao_subject_groups.csv
│
├── output/
│   ├── figures/
│   ├── scatter_xy/
│   └── reports/
└── docs/
```

---

## Data Handling Rules

| Folder | Git Tracked | Purpose |
|---|---|---|
| `data/dummy/` | ✅ Yes | Dummy data for development and testing |
| `data/production/` | ❌ No | Actual production data (excluded via `.gitignore`) |
| `data/ref/` | ✅ Yes | Small curated reference files (safe to commit) |

---

## Input Data: What to Copy

### For the upstream pipeline (`data_analysis_r/`)

The upstream scripts expect the raw MS-DIAL output and processed files to be in `data/production/`. Place the following in the indicated subdirectories:

| File | Where to put it | Used by |
|---|---|---|
| `df_conf.csv` (MS-DIAL raw output) | `data/production/raw/` | `enrich_pairwise_results.R` |
| `samplesheet.csv` | `data/production/processed/` | `load_data_production.R` |
| Feature matrix and metadata CSVs | `data/production/processed/` | `load_data_production.R` |

### For the downstream pipelines

These files are outputs of the upstream pipeline (or copies of them from `TMAO_pathway_analysis/`). If you are skipping the upstream re-run and starting from pre-computed results, copy them now:

```bash
# From inside wk_260615/, run once:
cp ../TMAO_pathway_analysis/data/plsda_mechanism_candidates.csv   data/production/processed/
cp ../TMAO_pathway_analysis/data/plsda_metaboanalyst_ora.csv      data/production/processed/
cp ../TMAO_pathway_analysis/data/plsda_metaboanalyst_msea.csv     data/production/processed/
cp ../TMAO_pathway_analysis/data/log2FC_values_production.csv     data/production/processed/
cp ../TMAO_pathway_analysis/data/hmdb_results_260515.csv          data/production/processed/
cp ../TMAO_pathway_analysis/data/pairwise_test_results_log2_production.csv data/production/processed/
```

### Reference file: `tmao_subject_groups.csv`

This file (`data/ref/tmao_subject_groups.csv`) contains curated TMAO response group assignments per subject. Columns: `Subject`, `X` (mean log2FC in X period), `Y` (mean log2FC in Y period), `n_detected`, `tmao_group`.

It is a refined version of the output of `classify_subjects_by_response.R` (which produces the simpler 4-group `subject_response_groups_414.csv`). The `tmao_group` values used in scripts 05–07 are:

| Group | Meaning |
|---|---|
| `Both_up` | TMAO rises in both X and Y periods — high producer |
| `Y_only` | TMAO rises only in Y period — carnitine-responsive |
| `X_down_Y_up` | TMAO falls in X, rises in Y — primed suppression |
| `X_up_Y_down` | TMAO rises in X, falls in Y — X-period responder |
| `X_only` | TMAO rises only in X period |
| `Non_producer` | No TMAO response in either period |
| `Sparse` / `Both_down` / `Weak` | Ambiguous; merged into "Weak" by scripts 05–07 |

The file is already in place at `data/ref/tmao_subject_groups.csv`. **No copy needed.**

---

## ⚠️ Path Fixes Required

### Python scripts

The Python scripts have absolute paths hardcoded to an old session (`/sessions/charming-youthful-pasteur/...`). Fix them:

```bash
# Run from inside wk_260615/
REPO="$(pwd)"
OLD="/sessions/charming-youthful-pasteur/mnt/TMAO_pathway_analysis"

for f in src/python/01_ontology_enrichment.py src/python/02_make_excel.py src/python/04_trajectory_analysis.py; do
  sed -i '' "s|${OLD}/data|${REPO}/data/production/processed|g" "$f"
  sed -i '' "s|${OLD}/output|${REPO}/output|g" "$f"
done
```

### R scripts 05–07 (`pathway_analysis_r/`)

Scripts 05–07 reference `output/tmao_subject_groups.csv` and `data/log2FC_values_production.csv`, but the files are now at `data/ref/tmao_subject_groups.csv` and `data/production/processed/log2FC_values_production.csv`. Update the path variables near the top of each script:

| Script | Variable | Old value | New value |
|---|---|---|---|
| 05, 06, 07 | `GROUPS_CSV` | `"output/tmao_subject_groups.csv"` | `"data/ref/tmao_subject_groups.csv"` |
| 05, 06, 07 | `LFC_CSV` | `"data/log2FC_values_production.csv"` | `"data/production/processed/log2FC_values_production.csv"` |

---

## Part 1: Upstream Analysis — `src/data_analysis_r/`

**Run all scripts from the project root** (`wk_260615/`), either via `Rscript` or by sourcing in RStudio after setting working directory to the project root.

The scripts must be run in the order below. Steps marked **[prereq]** must be run (or their outputs must exist) before the following scripts will work.

---

### Step 0 — Load data into R environment `[prereq for everything]`

#### `load_data_production.R`

Loads the three core data objects into the R environment: `samplesheet`, `feat_meta` (feature metadata), and `feat_mat` (feature intensity matrix). This script does not save any files — it must be sourced at the start of each RStudio session before running other `data_analysis_r/` scripts.

**Input:** `data/production/processed/samplesheet.csv` (+ feature matrix and metadata CSVs in the same folder)  
**Output:** R objects `samplesheet`, `feat_meta`, `feat_mat` in memory  
**Run:** `source("src/data_analysis_r/load_data_production.R")`

> Also available: `load_data.R` — same script but reads from `data/dummy/` for development with dummy data.

---

### Step 1 — Pairwise statistics

#### `boxplot_pairwise_production.R`

Runs paired Wilcoxon signed-rank tests for all annotated metabolites across all 6 timepoint comparisons (T1vT2, T1vT3, T2vT3 for the Early period; T4vT5, T4vT6, T5vT6 for the Late period). Applies BH FDR correction across all tests. Classifies each metabolite as `both`, `early_only`, `late_only`, or `ns`. Saves a stats results CSV and generates boxplot PDFs.

Key parameters at the top of the script: `TRANSFORMS` (default `"log2FC"`), `ALPHA` (default `0.05`), `INCLUDE_MSFINDER`.

**Input:** R objects from `load_data_production.R`  
**Output:**
- `data/production/processed/pairwise_test_results_log2FC.csv` — per-feature stats with classification
- `output/figures/boxplots_*.pdf` — boxplot PDFs

**Run:** `source("src/data_analysis_r/boxplot_pairwise_production.R")`

---

### Step 2 — Annotation enrichment

#### `enrich_pairwise_results.R`

Enriches the pairwise stats CSV with HMDB annotation (InChIKey, MS-FINDER structure, etc.) from the MS-DIAL raw output. Produces the `_enriched.csv` variant used by `pathway_analysis_r/01_preprocess.R`.

**Input:**
- `data/production/processed/pairwise_test_results_{TRANSFORM}.csv` (from step 1)
- `data/production/raw/df_conf.csv` (MS-DIAL raw output)
- `data/production/processed/feature_metadata_msfinder_2090.csv` (optional; adds MS-FINDER annotations for unknowns)

**Output:** `data/production/processed/pairwise_test_results_{TRANSFORM}_enriched.csv`

**Run:** `source("src/data_analysis_r/enrich_pairwise_results.R")`

> **Note on filename in `pathway_analysis_r/01_preprocess.R`:** That script references `pairwise_test_results_none_enriched.csv`. Either rename the file or update `STATISTICS_CSV` in `01_preprocess.R` to match the actual filename.

---

### Step 3 — Subject classification `[prereq for PLS scripts]`

#### `classify_subjects_by_response.R`

Classifies each subject into one of four response groups based on whether TMAO (Alignment ID 414 by default) rises in the Early period, Late period, both, or neither. This classification is used as the response variable (`Y`) in all PLS-DA analyses.

**Input:** R objects from `load_data_production.R`  
**Output:** `data/production/processed/subject_response_groups_414.csv`
- Columns: `Subject`, `group` (`both` / `X_only` / `Y_only` / `other`), log2FC values per timepoint, raw intensities

**Run:** `source("src/data_analysis_r/classify_subjects_by_response.R")`

> **Relationship to `tmao_subject_groups.csv`:** The `data/ref/tmao_subject_groups.csv` file is a more refined version of this output, with 7 group labels instead of 4 (`Both_up`, `Y_only`, `X_up_Y_down`, etc.). It was generated from a previous analysis run and is provided as a reference file. Scripts 05–07 use this refined version.

---

### Step 4 — PLS analyses `[prereq: steps 0 + 3]`

The following scripts all require both `load_data_production.R` and `classify_subjects_by_response.R` to have been run first.

#### `pls_alignment414.R`

PLS regression with TMAO (ID 414) log2FC as the continuous response variable. Fits two independent models: one for the Early period (per-subject average of T1–T3) and one for the Late period (T4–T6). Color-codes score plots by the 4 response groups.

**Output:** `output/figures/pls_early_414.pdf`, `output/figures/pls_late_414.pdf`

#### `pls_alignment414_per_timepoint.R`

Same as above but treats each subject × timepoint as an independent observation (n = subjects × 3 per period), giving higher statistical power at the cost of within-subject correlation.

**Output:** `output/figures/pls_early_414_pertp.pdf`, `output/figures/pls_late_414_pertp.pdf`

#### `plsda_per_timepoint.R`

PLS-DA with the 4 response groups as the categorical Y. Each subject × timepoint is an observation. Produces score plots (Comp1 vs Comp2) with 95% confidence ellipses per group.

**Output:** `output/figures/plsda_early_pertp.pdf`, `output/figures/plsda_late_pertp.pdf`

#### `plsda_loadings_vip.R`

Extracts VIP scores and Comp1 loadings from the PLS-DA model. Generates VIP bar plots (top N features), a VIP_early vs VIP_late scatter (identifies period-specific vs shared drivers), and a loading scatter (shows direction of effect). **This script produces the key CSV needed by `plsda_mechanism_candidates.R`.**

**Output:**
- `data/production/processed/plsda_vip_loadings.csv` — all features with VIP_early, VIP_late, loading_early, loading_late
- `output/figures/plsda_vip_early_top.pdf/.png`
- `output/figures/plsda_vip_late_top.pdf/.png`
- `output/figures/plsda_vip_scatter.pdf/.png`
- `output/figures/plsda_loading_scatter.pdf/.png`

---

### Step 5 — Metabolite classification `[prereq: step 4]`

#### `plsda_mechanism_candidates.R`

Classifies metabolites into mechanism classes using VIP scores and Comp1 loading quadrants. **This is where `plsda_mechanism_candidates.csv` comes from** — the master input for the Python pipeline and `pathway_analysis_r/` scripts 01–04.

Classification logic:

| `mechanism_class` | Condition |
|---|---|
| `Early_specific` | VIP_early ≥ threshold, VIP_late < threshold |
| `Late_specific` | VIP_late ≥ threshold, VIP_early < threshold |
| `Reversed_Q2` | Both VIP pass + Early loading (−), Late loading (+) |
| `Reversed_Q4` | Both VIP pass + Early loading (+), Late loading (−) |
| `Shared_Q1` / `Shared_Q3` | Both VIP pass, same loading direction |

Annotation confidence determines VIP threshold: well-annotated metabolites use a more lenient threshold; unknowns use a stricter one.

**Input:**
- `data/production/processed/plsda_vip_loadings.csv` (from `plsda_loadings_vip.R`)
- `data/production/processed/feature_metadata_msfinder_2090.csv`
- `data/production/processed/pairwise_test_results_{TRANSFORM}.csv` (optional, for p-value annotation)

**Output:**
- `plsda_mechanism_candidates.csv` — all candidates with mechanism_class, VIP, loading, InChIKey, Ontology
- `plsda_metaboanalyst_ora.csv` — stripped to ID + mechanism_class for MetaboAnalystR ORA
- `plsda_metaboanalyst_msea.csv` — ID + ranking score for MSEA

#### `classify_metabolites_metaboanalystr.R`

Alternative metabolite classification using a slightly different 5-class scheme (`early_specific`, `late_specific`, `shared_concordant`, `shared_discordant`, `ns`). Can be used as a cross-check against `plsda_mechanism_candidates.R`.

---

### Visualization scripts (run anytime after steps 1–3)

#### `plot_timeseries_production.R`

Time-series plots for all annotated metabolites with 1–5 annotation counts. Each feature gets a small panel showing individual-subject lines across T1–T6, with the early/late split marked. Generates a multi-page PDF.

Key settings: `TRANSFORM` (default `"log2FC"`), `EXCLUDE_LOW_SCORE`, `INCLUDE_MSFINDER`.

**Output:** `output/figures/timeseries_MSFINDER_log2FC_production.pdf`

#### `plot_selected_features.R`

Boxplots for a hand-specified list of `Alignment_IDs` (edit `TARGET_IDS` in the script). Reads p-values from a pre-computed pairwise stats CSV. Useful for zooming in on a handful of features of interest.

**Output:** `output/figures/selected_features_{transform}.pdf`

#### `pca_per_timepoint.R`

PCA score plots treating each subject × timepoint as an observation. Color = 4 response groups, shape = timepoint. Two plots: Early period and Late period.

**Output:** `output/figures/pca_early_pertp.pdf`, `output/figures/pca_late_pertp.pdf`

---

## Part 2: Downstream Analysis — `src/pathway_analysis_r/`

> **Working directory note:** Scripts 01–04 are designed to run from `data/production/processed/` (relative paths to bare filenames). Scripts 05–07 run from the **project root** (`wk_260615/`). See the run commands below.

---

### Scripts 01–04: Pathway ORA and Cytoscape export

These four scripts form a sequential pipeline. Each feeds the next.

```bash
cd wk_260615/data/production/processed
```

#### `01_preprocess.R`

Joins the HMDB annotation table with pairwise statistical results, filters by annotation match type and HMDB ID availability, deduplicates on HMDB ID (keeping highest-scoring annotation), computes per-perturbation log2FC and significance flags, and writes MetaboAnalystR-ready input files.

**Input (relative to `data/production/processed/`):**
- `hmdb_results_260515.csv` — HMDB annotation (InChIKey, primary pathway, match type, S/N)
- `pairwise_test_results_none_enriched.csv` ← rename the enriched file to match, or change `STATISTICS_CSV` in the script

**Key parameters:**

| Variable | Default | What it controls |
|---|---|---|
| `MATCH_TYPES` | `c("prefix")` | Annotation confidence filter. Use `c("exact")` for strictest. |
| `PADJ_THRESHOLD` | `0.05` | Significance cutoff for ORA inputs |
| `MIN_LOG2FC` | `0` | Min \|log2FC\| for the QEA ranked list |

**Output (in working directory):**
- `filtered_features.csv` — Annotated feature table with `log2FC_X`, `log2FC_Y`, `padj_X`, `padj_Y`. **Key input for scripts 02–04.**
- `input_ora_X.csv`, `input_ora_Y.csv` — Significant HMDB IDs for ORA
- `input_pathway_X.csv`, `input_pathway_Y.csv` — Ranked HMDB ID + log2FC for QEA

**Run:** `Rscript ../../../src/pathway_analysis_r/01_preprocess.R`

---

#### `02_enrichment_ora.R`

Runs MetaboAnalystR ORA on the significant metabolite lists from script 01, separately for perturbation X and Y. Outputs bar charts, dot plots, and a cross-comparison table.

> **macOS Apple Silicon note:** Cairo must be installed and is loaded first. The script auto-patches a MetaboAnalystR `default.dpi` crash. See `TMAO_pathway_analysis/docs/Cairo_issue.md` if needed.

**Input:** `input_ora_X.csv`, `input_ora_Y.csv`, `filtered_features.csv`

**Key parameters:** `MSET_LIBRARY` (default `"smpdb_pathway"`), `MIN_SET_SIZE` (default `2`), `USE_CUSTOM_BG` (default `FALSE`)

**Output:**
- `pertX/ora_results_pertX.csv`, `pertX/ora_bar_pertX.png`, `pertX/ora_dot_pertX.png`
- `pertY/` — same set for perturbation Y
- `ora_comparison_XY.csv` — X-only / Y-only / Both / Neither labels
- `ora_dotplot_comparison.png`

**Run:** `Rscript ../../../src/pathway_analysis_r/02_enrichment_ora.R`

---

#### `03_pathway_analysis.R`

Independent KEGG pathway ORA using `KEGGREST` (no MetaboAnalystR). Maps HMDB → KEGG compound IDs, fetches all human pathway memberships, runs hypergeometric tests. Results are cached in `.kegg_cache/` after the first run.

> **First run takes ~10–15 min** (KEGG API calls). Subsequent runs are fast.

**Input:** `filtered_features.csv`

**Output:**
- `pertX/pathway_results_pertX.csv`, `pertX/pathway_bubble_pertX.png`
- `pertY/` — same for perturbation Y
- `pathway_comparison_XY.csv`
- `pathway_scatter_XvsY.png`

**Run:** `Rscript ../../../src/pathway_analysis_r/03_pathway_analysis.R`

---

#### `04_cytoscape_export.R`

Builds Cytoscape-ready node and edge tables from pathway results. Metabolites and pathways become nodes; "member of" edges connect them.

> **Path note:** Edit `PROCESSED_DIR`, `PATHWAY_DIR`, and `OUT_DIR` at the top to match actual file locations before running.

**Input:** `filtered_features.csv`, `pathway_comparison_XY.csv`, optionally MetaboAnalystR `.rda` files

**Key parameters:** `PATHWAY_FILTER` (default `"sig_either"`), `MAX_PATHWAYS` (default `50`)

**Output (in `output/cytoscape/`):**
- `nodes_metabolites.csv`, `nodes_pathways.csv`, `nodes_combined.csv`
- `edges.csv`
- `network.sif` — import directly into Cytoscape

**Cytoscape import:** File → Import → Network from File → `network.sif`, then import `nodes_combined.csv` and `edges.csv`.

**Run:** `Rscript ../../../src/pathway_analysis_r/04_cytoscape_export.R`

---

### Scripts 05–07: TMAO Group × Bile Acid Analysis

> **Run from the project root** (`wk_260615/`). These scripts do not depend on scripts 01–04.

> **⚠️ Path fix needed before running:** Open each script and update two variables:
> - `GROUPS_CSV` → `"data/ref/tmao_subject_groups.csv"`
> - `LFC_CSV` → `"data/production/processed/log2FC_values_production.csv"`

---

#### `05_H1A_scatter.R`

Tests **Hypothesis 1-A**: secondary bile acid levels are positively correlated with TMAO production. Plots 16 pre-selected CDCA-conjugate and secondary bile acids against TMAO log2FC (IDs 414 + 415, averaged), broken down by TMAO response group. Produces box plots, a heatmap, scatter plots, and a correlation summary CSV.

**Input:** `log2FC_values_production.csv`, `tmao_subject_groups.csv`

**Output:**
- `output/figures/H1A_boxplot.png` — 4×4 facet: Y-period mean log2FC per group per bile acid
- `output/figures/H1A_heatmap.png` — Subject × bile acid heatmap, rows sorted by TMAO group
- `output/figures/H1A_scatter_all.png` — 4×4 scatter: each bile acid vs TMAO, Y-period mean, Pearson r shown
- `output/figures/H1A_scatter_tp_split.png` — Same scatter with T4vsT5 (●) and T4vsT6 (▲) as separate shapes
- `output/figures/H1A_scatter_tp_split_Xperiod.png` — Same split-TP scatter for X period
- `output/H1A_correlation_table.csv` — Per-bile-acid Pearson r for T4vsT5, T4vsT6, and pooled Y period

**Run:** `cd wk_260615 && Rscript src/pathway_analysis_r/05_H1A_scatter.R`

---

#### `06_group_composite_scatter.R`

Reduces the 16 bile acids to a single **composite z-score** per subject per timepoint pair (z-standardize all 16 × all subjects × both TP pairs globally, then average per subject per TP). Plots composite vs TMAO log2FC — one two-panel figure per TMAO group. Focal group highlighted in color against grey background of all subjects.

**Input:** `log2FC_values_production.csv`, `tmao_subject_groups.csv`

**Output:**
- `output/scatter_xy/{group}_scatter.png` — 7 PNGs (one per group), Y period left panel / X period right panel
- `output/scatter_xy/r_summary.csv` — Pearson r and p-value per group per period

**Run:** `cd wk_260615 && Rscript src/pathway_analysis_r/06_group_composite_scatter.R`

---

#### `07_bothup_report.R`

Generates a **17-page A4 PDF** for the `Both_up` group (n = 5 high TMAO producers).

- **Page 1:** TMAO time-series for all 5 subjects (IDs 414 and 415 as separate lines, not averaged), T1–T6.
- **Pages 2–17:** 4 period comparisons (T1vsT2, T1vsT3, T4vsT5, T4vsT6) × 4 bile acid chunks (4 acids per chunk) = 16 pages. Each page is a 4×2 grid (4 bile acids × 2 TMAO IDs), points labeled by subject, Pearson r annotated.

**Input:** `log2FC_values_production.csv`, `tmao_subject_groups.csv`

**Output:** `output/Both_up_detailed_report.pdf` (17 pages, A4 portrait, `cairo_pdf`)

**Run:** `cd wk_260615 && Rscript src/pathway_analysis_r/07_bothup_report.R`

---

## Part 3: Python Pipeline — `src/python/`

**Run from the project root** (`wk_260615/`). Independent of the R pipelines.

---

#### `01_ontology_enrichment.py`

Fisher's exact test on chemical Ontology classes, separately for Early_specific, Late_specific, and Reversed metabolites vs the full PLS-DA background. Also does a direct Early vs Late comparison. Appends KEGG pathway annotations via a hand-curated map.

**Input:** `plsda_mechanism_candidates.csv`, `plsda_metaboanalyst_ora.csv`, `plsda_metaboanalyst_msea.csv`

**Output:**
- `output/01_ontology_ORA_all_groups.csv`
- `output/02_early_vs_late_ontology_comparison.csv`
- `output/03_Early_specific_metabolites.csv`
- `output/04_Late_specific_metabolites.csv`
- `output/05_Reversed_metabolites.csv`

**Run:** `python src/python/01_ontology_enrichment.py`

---

#### `02_make_excel.py`

Same enrichment logic as `01_ontology_enrichment.py` but outputs a polished 8-sheet Excel workbook (Summary, Early vs Late comparison, ORA per group, metabolite lists).

**Output:** `output/TMAO_pathway_results.xlsx`

**Run:** `python src/python/02_make_excel.py`

---

#### `04_trajectory_analysis.py`

Spaghetti plots (individual-subject log2FC across T1–T6) for 15 pre-selected carnitines and bile acids. Grey lines = individual subjects, colored line = mean ± SD. Color = mechanism class.

**Input:** `log2FC_values_production.csv`

**Output:**
- `output/fig1_carnitine_trajectories.png` — 2×4: 8 acylcarnitines + free carnitine
- `output/fig2_bileacid_trajectories.png` — 2×4: 7 bile acids
- `output/fig3_representative_comparison.png` — 1×4: one representative per mechanism class
- `output/fig4_mean_response_summary.png` — horizontal bar: mean within-period response ± SEM

**Run:** `python src/python/04_trajectory_analysis.py`

---

## Full Run Order (summary)

```bash
cd wk_260615

# ── 0. Setup ──────────────────────────────────────────────────────────────────
# Copy pre-computed downstream inputs (if skipping upstream re-run):
cp ../TMAO_pathway_analysis/data/plsda_mechanism_candidates.csv   data/production/processed/
cp ../TMAO_pathway_analysis/data/plsda_metaboanalyst_ora.csv      data/production/processed/
cp ../TMAO_pathway_analysis/data/plsda_metaboanalyst_msea.csv     data/production/processed/
cp ../TMAO_pathway_analysis/data/log2FC_values_production.csv     data/production/processed/
cp ../TMAO_pathway_analysis/data/hmdb_results_260515.csv          data/production/processed/
cp ../TMAO_pathway_analysis/data/pairwise_test_results_log2_production.csv data/production/processed/
# (tmao_subject_groups.csv is already in data/ref/ — no copy needed)

# Fix Python paths (see "Path Fixes" section above)
# Fix paths in pathway_analysis_r/05–07 (GROUPS_CSV and LFC_CSV)

# ── 1. Upstream R pipeline (from project root, sequential) ────────────────────
# Open RStudio, set WD to project root, then source in order:
#   source("src/data_analysis_r/load_data_production.R")      # load data
#   source("src/data_analysis_r/boxplot_pairwise_production.R") # stats
#   source("src/data_analysis_r/enrich_pairwise_results.R")    # annotate
#   source("src/data_analysis_r/classify_subjects_by_response.R") # group subjects
#   source("src/data_analysis_r/plsda_loadings_vip.R")         # VIP + loadings
#   source("src/data_analysis_r/plsda_mechanism_candidates.R") # classify metabolites

# ── 2. Downstream R 01-04 (from data/production/processed/, sequential) ───────
cd data/production/processed
Rscript ../../../src/pathway_analysis_r/01_preprocess.R
Rscript ../../../src/pathway_analysis_r/02_enrichment_ora.R
Rscript ../../../src/pathway_analysis_r/03_pathway_analysis.R   # slow on first run
Rscript ../../../src/pathway_analysis_r/04_cytoscape_export.R

# ── 3. Downstream R 05-07 (from project root, any order) ──────────────────────
cd ../../..   # back to wk_260615/
Rscript src/pathway_analysis_r/05_H1A_scatter.R
Rscript src/pathway_analysis_r/06_group_composite_scatter.R
Rscript src/pathway_analysis_r/07_bothup_report.R

# ── 4. Python pipeline (from project root, any order) ─────────────────────────
python src/python/01_ontology_enrichment.py
python src/python/02_make_excel.py
python src/python/04_trajectory_analysis.py
```

---

## Expected Outputs at a Glance

| Script | Key output to look at first |
|---|---|
| `data_analysis_r/plsda_loading_scatter` | `output/figures/plsda_loading_scatter.png` — quadrant plot of Early vs Late loadings |
| `data_analysis_r/plsda_mechanism_candidates` | `plsda_mechanism_candidates.csv` — check group sizes and annotation quality |
| `pathway_analysis_r/02` | `ora_comparison_XY.csv` — which SMPDB pathway sets are X-only / Y-only / Both |
| `pathway_analysis_r/03` | `pathway_scatter_XvsY.png` — KEGG pathway enrichment overview |
| `pathway_analysis_r/04` | Load `output/cytoscape/network.sif` into Cytoscape, then node/edge tables |
| `pathway_analysis_r/05` | `output/figures/H1A_scatter_all.png` — do any bile acids correlate with TMAO? |
| `pathway_analysis_r/06` | `output/scatter_xy/r_summary.csv` — which group drives the correlation? |
| `pathway_analysis_r/07` | `output/Both_up_detailed_report.pdf` — deep dive into the 5 high-producer subjects |
| `python/01` | `output/02_early_vs_late_ontology_comparison.csv` — chemical class differences |
| `python/02` | `output/TMAO_pathway_results.xlsx` → **Summary** sheet |
| `python/04` | `output/fig1_carnitine_trajectories.png`, `fig2_bileacid_trajectories.png` |

---

## Dependencies

**Python:** `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `openpyxl`

```bash
pip install pandas numpy scipy statsmodels matplotlib openpyxl --break-system-packages
```

**R:** `tidyverse`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `ggrepel`, `patchwork`, `pheatmap`, `scales`, `Cairo`, `MetaboAnalystR`, `KEGGREST`, `tibble`

```r
install.packages(c("tidyverse", "ggrepel", "patchwork", "pheatmap", "scales", "Cairo", "tibble"))
BiocManager::install("KEGGREST")
# MetaboAnalystR: see TMAO_pathway_analysis/docs/MetaboAnalystR_README.md
```

---

## Contact / Handoff Notes

- Handed off by: sano (ducky@keio.jp), 2026-06-12
- Source project: `TMAO_pathway_analysis/`
- For PLS-DA classification logic: `TMAO_pathway_analysis/docs/downstream_metaboanalystr.md` and `analysis_roadmap.md`
- For Cairo install issues on macOS Apple Silicon: `TMAO_pathway_analysis/docs/Cairo_issue.md`
