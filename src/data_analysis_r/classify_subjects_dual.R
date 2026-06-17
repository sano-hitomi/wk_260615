# =============================================================================
# classify_subjects_dual.R
# Combine 414 (TMAO) and 4263 response classifications into 3 groups
#
# Group assignment (per subject):
#   resp_414  := group_414  != "other"   (any upregulation in TMAO)
#   resp_4263 := group_4263 != "other"   (any upregulation in 4263)
#
#   "dual"      :  resp_414 &  resp_4263
#   "TMAO_only" :  resp_414 & !resp_4263
#   "tgt_only"  : !resp_414 &  resp_4263
#   excluded    : !resp_414 & !resp_4263  (non-responders)
#   excluded    : outliers (Mahalanobis distance > threshold)
#
# Outlier detection:
#   Features used: log2FC_T2, log2FC_T3, log2FC_T5, log2FC_T6 from BOTH
#   metabolites (8-dim vector per subject). Mahalanobis distance computed
#   with robust covariance (MASS::cov.rob); subjects beyond chi2(df=8, p=0.01)
#   are flagged.
#
# Outputs:
#   data/production/processed/subject_dual_groups.csv
#     Subject, dual_group, resp_414, resp_4263,
#     group_414, group_4263,
#     outlier_flag, mahal_dist,
#     log2FC columns for both metabolites,
#     raw intensity columns for both metabolites
#
# Prerequisites:
#   Run load_data_production.R first.
#   Run classify_subjects_by_response.R   (generates subject_response_groups_414.csv)
#   Run classify_subjects_by_response_4263.R (generates subject_response_groups_4263.csv)
# =============================================================================

ID_A   <- 414
ID_B   <- 4263
OUT_DIR <- "data/production/processed"
OUT_CSV <- file.path(OUT_DIR, "subject_dual_groups.csv")

OUTLIER_P <- 0.01   # chi2 tail probability threshold for Mahalanobis

# -----------------------------------------------------------------------------
# 0. Prerequisites
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0)
  stop("Missing: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")

csv_a <- file.path(OUT_DIR, sprintf("subject_response_groups_%d.csv", ID_A))
csv_b <- file.path(OUT_DIR, sprintf("subject_response_groups_%d.csv", ID_B))
if (!file.exists(csv_a))
  stop("Not found: ", csv_a, "\nRun classify_subjects_by_response.R first.")
if (!file.exists(csv_b))
  stop("Not found: ", csv_b, "\nRun classify_subjects_by_response_4263.R first.")

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
for (p in c("dplyr", "MASS")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Load both classifications
# -----------------------------------------------------------------------------
df_a <- read.csv(csv_a, stringsAsFactors = FALSE)
df_b <- read.csv(csv_b, stringsAsFactors = FALSE)

df_a$Subject <- as.character(df_a$Subject)
df_b$Subject <- as.character(df_b$Subject)

# Suffix columns to avoid collision
rename_fc <- function(df, suffix) {
  fc_cols  <- grep("^log2FC_T", colnames(df), value = TRUE)
  raw_cols <- grep("^raw_T",    colnames(df), value = TRUE)
  colnames(df)[colnames(df) %in% fc_cols]  <- paste0(fc_cols,  "_", suffix)
  colnames(df)[colnames(df) %in% raw_cols] <- paste0(raw_cols, "_", suffix)
  colnames(df)[colnames(df) == "group"]    <- paste0("group_", suffix)
  colnames(df)[colnames(df) == "early_up"] <- paste0("early_up_", suffix)
  colnames(df)[colnames(df) == "late_up"]  <- paste0("late_up_", suffix)
  df
}

df_a <- rename_fc(df_a, ID_A)
df_b <- rename_fc(df_b, ID_B)

# Inner join on Subject
merged <- merge(df_a, df_b, by = "Subject")
cat(sprintf("Subjects in both analyses: %d\n", nrow(merged)))

# -----------------------------------------------------------------------------
# 3. Responder flags
# -----------------------------------------------------------------------------
merged$resp_414  <- merged[[paste0("group_", ID_A)]]  != "other"
merged$resp_4263 <- merged[[paste0("group_", ID_B)]] != "other"

# 3-group classification (pre-outlier removal)
merged$dual_group_raw <- with(merged, ifelse(
   resp_414 &  resp_4263, "dual",
  ifelse( resp_414 & !resp_4263, "TMAO_only",
  ifelse(!resp_414 &  resp_4263, "tgt_only",
                                 "non_responder"))))

cat("\nGroup counts (before outlier exclusion):\n")
print(table(merged$dual_group_raw))

# -----------------------------------------------------------------------------
# 4. Outlier detection via Mahalanobis distance
#    Input: 8 log2FC values per subject (T2,T3,T5,T6 × 2 metabolites)
# -----------------------------------------------------------------------------
fc_cols_a <- paste0("log2FC_T", c(2,3,5,6), "_", ID_A)
fc_cols_b <- paste0("log2FC_T", c(2,3,5,6), "_", ID_B)
fc_use    <- c(fc_cols_a, fc_cols_b)

fc_mat <- as.matrix(merged[, fc_use])
rownames(fc_mat) <- merged$Subject

# Fill NA with column means (rare; should not happen in practice)
for (j in seq_len(ncol(fc_mat))) {
  na_idx <- is.na(fc_mat[, j])
  if (any(na_idx)) fc_mat[na_idx, j] <- mean(fc_mat[, j], na.rm = TRUE)
}

# Robust covariance → Mahalanobis
rob  <- tryCatch(
  MASS::cov.rob(fc_mat),
  error = function(e) {
    warning("cov.rob failed; falling back to standard cov. e: ", conditionMessage(e))
    list(center = colMeans(fc_mat), cov = cov(fc_mat))
  }
)
mahal <- mahalanobis(fc_mat, center = rob$center, cov = rob$cov)

# Chi-squared threshold (df = 8 variables)
chi2_thresh <- qchisq(1 - OUTLIER_P, df = ncol(fc_mat))
cat(sprintf("\nMahalanobis chi2 threshold (df=%d, p=%.3f): %.2f\n",
            ncol(fc_mat), OUTLIER_P, chi2_thresh))

merged$mahal_dist   <- round(mahal, 3)
merged$outlier_flag <- mahal > chi2_thresh

cat(sprintf("Outliers detected: %d\n", sum(merged$outlier_flag)))
if (any(merged$outlier_flag)) {
  cat("  Subjects flagged:\n")
  out_subj <- merged[merged$outlier_flag, c("Subject", "dual_group_raw", "mahal_dist")]
  print(out_subj)
}

# -----------------------------------------------------------------------------
# 5. Final 3-group assignment
# -----------------------------------------------------------------------------
merged$dual_group <- ifelse(
  merged$outlier_flag | merged$dual_group_raw == "non_responder",
  NA_character_,
  merged$dual_group_raw
)

cat("\nFinal group counts (NA = excluded):\n")
print(table(merged$dual_group, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(merged, OUT_CSV, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 7. Summary plot (boxplot of FC per group)
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  library(tidyr)

  plot_df <- merged[!is.na(merged$dual_group), ]

  # Pivot long: T2/T3/T5/T6 FC for each metabolite
  long_a <- plot_df[, c("Subject", "dual_group",
                         paste0("log2FC_T", c(2,3,5,6), "_", ID_A))]
  long_a <- tidyr::pivot_longer(long_a, cols = starts_with("log2FC"),
                                 names_to = "timepoint", values_to = "log2FC")
  long_a$metabolite <- paste0("ID ", ID_A, " (TMAO)")
  long_a$timepoint  <- sub(paste0("_", ID_A), "", long_a$timepoint)

  long_b <- plot_df[, c("Subject", "dual_group",
                         paste0("log2FC_T", c(2,3,5,6), "_", ID_B))]
  long_b <- tidyr::pivot_longer(long_b, cols = starts_with("log2FC"),
                                 names_to = "timepoint", values_to = "log2FC")
  long_b$metabolite <- paste0("ID ", ID_B)
  long_b$timepoint  <- sub(paste0("_", ID_B), "", long_b$timepoint)

  long_all <- rbind(long_a, long_b)
  long_all$dual_group <- factor(long_all$dual_group,
                                 levels = c("dual", "TMAO_only", "tgt_only"))

  p_box <- ggplot(long_all, aes(x = dual_group, y = log2FC, fill = dual_group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60", linetype = "dashed") +
    geom_boxplot(outlier.size = 1.5, alpha = 0.8, width = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
    facet_grid(metabolite ~ timepoint, scales = "free_y") +
    scale_fill_manual(values = c("dual"      = "#E41A1C",
                                  "TMAO_only" = "#377EB8",
                                  "tgt_only"  = "#4DAF4A")) +
    labs(
      title    = "log2FC by dual-response group",
      subtitle = sprintf("Outliers excluded (Mahalanobis p < %.2f); non-responders excluded",
                         OUTLIER_P),
      x        = NULL, y = "log2FC (period-specific reference)",
      fill     = "Group"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 10, face = "bold"),
      legend.position  = "right"
    )

  out_dir_fig <- "output/figures"
  dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
  base_name <- file.path(out_dir_fig, "subject_dual_groups_fc")
  ggsave(paste0(base_name, ".pdf"), plot = p_box, width = 10, height = 5)
  ggsave(paste0(base_name, ".png"), plot = p_box, width = 10, height = 5, dpi = 150)
  cat(sprintf("Plot saved: %s\n", base_name))
}

# -----------------------------------------------------------------------------
# 8. Print subject list per group
# -----------------------------------------------------------------------------
cat("\n=== Subjects per group ===\n")
for (g in c("dual", "TMAO_only", "tgt_only")) {
  subj_g <- merged$Subject[!is.na(merged$dual_group) & merged$dual_group == g]
  cat(sprintf("  [%s] n=%d : %s\n", g, length(subj_g), paste(sort(subj_g), collapse = ", ")))
}
exc <- merged$Subject[is.na(merged$dual_group)]
cat(sprintf("  [excluded] n=%d : %s\n", length(exc), paste(sort(exc), collapse = ", ")))
