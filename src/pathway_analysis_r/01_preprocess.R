# ==============================================================================
# 01_preprocess.R
# Purpose : Load, join, filter, and prepare MetaboAnalystR-ready input tables
# Inputs  : data/production/annotation.csv   — HMDB ID / metadata
#           data/production/statistics.csv   — log2FC / padj per comparison
# Outputs : data/production/processed/filtered_features.csv
#           data/production/processed/input_ora_X.csv
#           data/production/processed/input_ora_Y.csv
#           data/production/processed/input_pathway_X.csv
#           data/production/processed/input_pathway_Y.csv
#
# Experimental design
#   Perturbation X : T1 (ctrl) vs T2, T1 vs T3
#   Perturbation Y : T4 (ctrl) vs T5, T4 vs T6
#   T2vT3 / T5vT6  : within-perturbation replicate consistency (reference only)
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)

# ------------------------------------------------------------------------------
# Configuration — edit paths and thresholds here
# ------------------------------------------------------------------------------
ANNOTATION_CSV <- "hmdb_results_260515.csv"
STATISTICS_CSV <- "pairwise_test_results_none_enriched.csv"
OUT_DIR        <- "."

# match_type filter: "exact" only is recommended for high-confidence results.
# Set to c("exact", "prefix") to include near-exact (same formula, ambiguous stereo).
MATCH_TYPES    <- c("prefix")

# Adjusted p-value threshold for ORA input (significant metabolite list)
PADJ_THRESHOLD <- 0.05

# Minimum mean log2FC absolute value to include in the ranked list for pathway QEA
# Set to 0 to include all annotated features
MIN_LOG2FC     <- 0

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
anno  <- read_csv(ANNOTATION_CSV, show_col_types = FALSE)
stats <- read_csv(STATISTICS_CSV, show_col_types = FALSE)

cat("Annotation rows :", nrow(anno),  "\n")
cat("Statistics rows :", nrow(stats), "\n")

# ------------------------------------------------------------------------------
# 2. Join on Alignment_ID
# ------------------------------------------------------------------------------
# Select only the columns needed from annotation to avoid duplication
anno_slim <- anno %>%
  select(
    Alignment_ID,
    hmdb_id,
    match_type,
    hmdb_inchikey,
    INCHIKEY_final,
    classification,
    primary_pathway,
    pathway_count,
    total_score  = `Total score`,
    sn_average   = `S/N average`,
    Matched_peaks_percentage = `Matched peaks percentage`
  )

merged <- stats %>%
  inner_join(anno_slim, by = "Alignment_ID")

cat("After join        :", nrow(merged), "features\n")

# ------------------------------------------------------------------------------
# 3. Filter by match_type and HMDB ID availability
# ------------------------------------------------------------------------------
filtered <- merged %>%
  filter(
    match_type %in% MATCH_TYPES,
    !is.na(hmdb_id),
    hmdb_id != ""
  )

cat("After match filter:", nrow(filtered), "features (match_type =",
    paste(MATCH_TYPES, collapse = "/"), ")\n")

# Flag duplicate HMDB IDs (multiple features mapping to same metabolite)
dup_ids <- filtered %>%
  group_by(hmdb_id) %>%
  filter(n() > 1) %>%
  nrow()
cat("Rows with duplicated HMDB ID:", dup_ids, "\n")

# When the same HMDB ID appears multiple times, keep the row with the highest
# Total score (most confident annotation).
filtered <- filtered %>%
  group_by(hmdb_id) %>%
  slice_max(total_score, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("After deduplication:", nrow(filtered), "unique features\n")

# ------------------------------------------------------------------------------
# 4. Compute summary statistics per perturbation
# ------------------------------------------------------------------------------
filtered <- filtered %>%
  mutate(
    # Mean log2FC across the two ctrl-vs-perturbation comparisons
    log2FC_X     = rowMeans(cbind(T1vT2, T1vT3), na.rm = TRUE),
    log2FC_Y     = rowMeans(cbind(T4vT5, T4vT6), na.rm = TRUE),

    # Within-perturbation replicate delta (|T2 vs T3|, |T5 vs T6|)
    # Small values indicate high replicate consistency
    rep_delta_X  = abs(T2vT3),
    rep_delta_Y  = abs(T5vT6),

    # Most significant padj across the two ctrl comparisons
    padj_X       = pmin(padj_T1vT2, padj_T1vT3, na.rm = TRUE),
    padj_Y       = pmin(padj_T4vT5, padj_T4vT6, na.rm = TRUE),

    # Significance flags
    sig_X        = padj_X < PADJ_THRESHOLD,
    sig_Y        = padj_Y < PADJ_THRESHOLD
  )

cat("\nSignificant in X (padj <", PADJ_THRESHOLD, "):",
    sum(filtered$sig_X, na.rm = TRUE), "features\n")
cat("Significant in Y (padj <", PADJ_THRESHOLD, "):",
    sum(filtered$sig_Y, na.rm = TRUE), "features\n")

# ------------------------------------------------------------------------------
# 5. Save full filtered feature table
# ------------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

write_csv(filtered, file.path(OUT_DIR, "filtered_features.csv"))
cat("\nSaved:", file.path(OUT_DIR, "filtered_features.csv"), "\n")

# ------------------------------------------------------------------------------
# 6. Prepare MetaboAnalystR input tables
#
# Format for ORA (msetora / pathora):
#   A plain text file with one HMDB ID per line (significant metabolites only).
#
# Format for Pathway QEA (pathqea):
#   Two-column CSV: HMDB ID, log2FC score (all annotated features, ranked).
#   MetaboAnalystR expects the first column to be compound IDs and the second
#   to be numeric scores.
# ------------------------------------------------------------------------------

# --- ORA input: significant HMDB IDs ---
ora_X <- filtered %>%
  filter(sig_X) %>%
  select(hmdb_id, Metabolite_name, log2FC_X, padj_X) %>%
  arrange(padj_X)

ora_Y <- filtered %>%
  filter(sig_Y) %>%
  select(hmdb_id, Metabolite_name, log2FC_Y, padj_Y) %>%
  arrange(padj_Y)

write_csv(ora_X, file.path(OUT_DIR, "input_ora_X.csv"))
write_csv(ora_Y, file.path(OUT_DIR, "input_ora_Y.csv"))
cat("Saved ORA input X:", nrow(ora_X), "metabolites\n")
cat("Saved ORA input Y:", nrow(ora_Y), "metabolites\n")

# --- Pathway QEA input: ranked full list with log2FC score ---
# Filter by minimum absolute log2FC if set
pathway_X <- filtered %>%
  filter(abs(log2FC_X) >= MIN_LOG2FC | is.na(log2FC_X)) %>%
  select(hmdb_id, log2FC_X) %>%
  arrange(desc(log2FC_X)) %>%
  rename(Score = log2FC_X)

pathway_Y <- filtered %>%
  filter(abs(log2FC_Y) >= MIN_LOG2FC | is.na(log2FC_Y)) %>%
  select(hmdb_id, log2FC_Y) %>%
  arrange(desc(log2FC_Y)) %>%
  rename(Score = log2FC_Y)

write_csv(pathway_X, file.path(OUT_DIR, "input_pathway_X.csv"))
write_csv(pathway_Y, file.path(OUT_DIR, "input_pathway_Y.csv"))
cat("Saved pathway input X:", nrow(pathway_X), "metabolites\n")
cat("Saved pathway input Y:", nrow(pathway_Y), "metabolites\n")

# ------------------------------------------------------------------------------
# 7. Quick summary table (console)
# ------------------------------------------------------------------------------
cat("\n--- Summary ---\n")
summary_tbl <- filtered %>%
  summarise(
    n_total      = n(),
    n_sig_X      = sum(sig_X, na.rm = TRUE),
    n_sig_Y      = sum(sig_Y, na.rm = TRUE),
    n_sig_both   = sum(sig_X & sig_Y, na.rm = TRUE),
    n_sig_X_only = sum(sig_X & !sig_Y, na.rm = TRUE),
    n_sig_Y_only = sum(!sig_X & sig_Y, na.rm = TRUE)
  )
print(summary_tbl)
