# =============================================================================
# classify_subjects_by_response_4263.R
# Classify subjects into 4 response groups for Alignment ID 4263
#
# Identical logic to classify_subjects_by_response.R (TARGET_ID = 414),
# but with TARGET_ID = 4263.
#
# Classification (per subject):
#   early_up : log2FC_T2 > 0  OR  log2FC_T3 > 0   (T2 vs T1, or T3 vs T1)
#   late_up  : log2FC_T5 > 0  OR  log2FC_T6 > 0   (T5 vs T4, or T6 vs T4)
#
#   "both"   : early_up AND  late_up
#   "X_only" : early_up AND !late_up
#   "Y_only" : !early_up AND  late_up
#   "other"  : !early_up AND !late_up
#
# Output:
#   data/production/processed/subject_response_groups_4263.csv
#
# Prerequisites:
#   Run load_data_production.R first (samplesheet / feat_meta / feat_mat).
# =============================================================================

TARGET_ID  <- 4263
OUT_DIR    <- "data/production/processed"

# -----------------------------------------------------------------------------
# 0. Check required objects
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
}

aid_char <- as.character(TARGET_ID)
if (!aid_char %in% rownames(feat_mat)) {
  stop(sprintf("Alignment_ID %s not found in feat_mat.", aid_char))
}

# -----------------------------------------------------------------------------
# 1. Package
# -----------------------------------------------------------------------------
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# -----------------------------------------------------------------------------
# 2. Average intensities over repeats → one value per Subject × Timepoint
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

raw_vals <- sapply(group_cols, function(cols) {
  mean(as.numeric(feat_mat[aid_char, cols]), na.rm = TRUE)
})

subjects   <- unique(bio$Subject)
timepoints <- 1:6

raw_mat <- matrix(
  NA_real_,
  nrow     = length(subjects),
  ncol     = length(timepoints),
  dimnames = list(subjects, paste0("raw_T", timepoints))
)
for (subj in subjects) {
  for (tp in timepoints) {
    key <- paste0(subj, "-", tp)
    if (key %in% names(raw_vals)) raw_mat[subj, paste0("raw_T", tp)] <- raw_vals[[key]]
  }
}

# -----------------------------------------------------------------------------
# 3. log2 transform
# -----------------------------------------------------------------------------
all_vals <- as.vector(raw_mat)
nz       <- all_vals[all_vals > 0 & !is.na(all_vals)]
offset   <- if (length(nz) > 0) min(nz) / 2 else 1
log2_mat <- log2(raw_mat + offset)

# -----------------------------------------------------------------------------
# 4. Period-specific log2FC
# -----------------------------------------------------------------------------
fc_mat <- matrix(
  NA_real_,
  nrow     = length(subjects),
  ncol     = length(timepoints),
  dimnames = list(subjects, paste0("log2FC_T", timepoints))
)
for (subj in subjects) {
  ref_early <- log2_mat[subj, "raw_T1"]
  for (tp in 1:3) {
    fc_mat[subj, paste0("log2FC_T", tp)] <- log2_mat[subj, paste0("raw_T", tp)] - ref_early
  }
  ref_late <- log2_mat[subj, "raw_T4"]
  for (tp in 4:6) {
    fc_mat[subj, paste0("log2FC_T", tp)] <- log2_mat[subj, paste0("raw_T", tp)] - ref_late
  }
}

early_up <- fc_mat[, "log2FC_T2"] > 0 | fc_mat[, "log2FC_T3"] > 0
late_up  <- fc_mat[, "log2FC_T5"] > 0 | fc_mat[, "log2FC_T6"] > 0
early_up[is.na(early_up)] <- FALSE
late_up [is.na(late_up)]  <- FALSE

# -----------------------------------------------------------------------------
# 5. Classify
# -----------------------------------------------------------------------------
group <- ifelse(
   early_up &  late_up, "both",
  ifelse( early_up & !late_up, "X_only",
  ifelse(!early_up &  late_up, "Y_only",
                               "other"))
)

# -----------------------------------------------------------------------------
# 6. Assemble and save
# -----------------------------------------------------------------------------
result_df <- data.frame(
  Subject  = subjects,
  group    = group,
  early_up = early_up,
  late_up  = late_up,
  stringsAsFactors = FALSE
)
fc_df  <- as.data.frame(round(fc_mat, 4));  fc_df$Subject  <- rownames(fc_df)
raw_df <- as.data.frame(round(raw_mat, 2)); raw_df$Subject <- rownames(raw_df)
result_df <- merge(result_df, fc_df,  by = "Subject")
result_df <- merge(result_df, raw_df, by = "Subject")
result_df <- result_df[order(result_df$group, result_df$Subject), ]
rownames(result_df) <- NULL

meta_row <- feat_meta[feat_meta$Alignment_ID == TARGET_ID, ]
name_str <- if (nrow(meta_row) > 0) meta_row$Metabolite_name[1] else "Unknown"

cat(sprintf("\n=== Alignment_ID %d : %s ===\n", TARGET_ID, name_str))
cat(sprintf("  Offset (log2 transform) : %.4g\n", offset))
cat(sprintf("  Subjects analysed       : %d\n",   length(subjects)))
cat(sprintf("\n  Group summary:\n"))
print(table(result_df$group))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(OUT_DIR, sprintf("subject_response_groups_%d.csv", TARGET_ID))
write.csv(result_df, out_csv, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", out_csv))
