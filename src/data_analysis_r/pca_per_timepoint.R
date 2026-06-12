# =============================================================================
# pca_per_timepoint.R
# PCA score plots — each subject × timepoint as an independent observation
#
# Design:
#   X : all metabolite features, log2-transformed, averaged over repeats
#       Rows = subject-timepoint pairs (e.g. "001-T1", "001-T2", ...)
#   Two independent analyses:
#     Early (Perturbation X) : T1, T2, T3  → n_subjects × 3 observations
#     Late  (Perturbation Y) : T4, T5, T6  → n_subjects × 3 observations
#
# Coloring / shape:
#   Color : 4-group classification from GROUPS_CSV
#           (output of classify_subjects_by_response.R)
#   Shape : timepoint within the period
#
# Output:
#   output/figures/pca_early_pertp.pdf
#   output/figures/pca_late_pertp.pdf
#
# Prerequisites:
#   1. load_data_production.R  (samplesheet, feat_meta, feat_mat)
#   2. classify_subjects_by_response.R  (generates GROUPS_CSV)
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
TARGET_ID  <- 414
GROUPS_CSV <- sprintf("data/production/processed/subject_response_groups_%d.csv",
                      TARGET_ID)
OUT_DIR    <- "output/figures"

GROUP_COLORS <- c(
  both   = "#E41A1C",
  X_only = "#377EB8",
  Y_only = "#4DAF4A",
  other  = "#984EA3"
)

TP_SHAPES <- c("1" = 16, "2" = 17, "3" = 15,
               "4" = 16, "5" = 17, "6" = 15)

# -----------------------------------------------------------------------------
# 0. Check prerequisites
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
}
if (!file.exists(GROUPS_CSV)) {
  stop("Groups CSV not found: ", GROUPS_CSV,
       "\nRun classify_subjects_by_response.R first.")
}

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# -----------------------------------------------------------------------------
# 2. Load group classifications
# -----------------------------------------------------------------------------
groups_df <- read.csv(GROUPS_CSV, stringsAsFactors = FALSE)
groups_df$Subject <- as.character(groups_df$Subject)

# -----------------------------------------------------------------------------
# 3. Build per-subject-timepoint log2 matrix (averaged over repeats)
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

avg_list <- lapply(group_cols, function(cols) {
  rowMeans(feat_mat[, cols, drop = FALSE], na.rm = TRUE)
})
avg_mat_full <- do.call(cbind, avg_list)   # features × subject-timepoints

# log2 transform with per-feature offset (min non-zero / 2)
log2_mat_full <- t(apply(avg_mat_full, 1, function(x) {
  nz     <- x[x > 0 & !is.na(x)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2(x + offset)
}))
# log2_mat_full: features × subject-timepoints → transpose to obs × features
log2_mat_full <- t(log2_mat_full)
colnames(log2_mat_full) <- rownames(feat_mat)

row_meta <- data.frame(
  SubjTP    = names(group_cols),
  Subject   = sub("-[0-9]+$", "", names(group_cols)),
  Timepoint = as.character(as.integer(sub(".*-", "", names(group_cols)))),
  stringsAsFactors = FALSE
)
rownames(log2_mat_full) <- row_meta$SubjTP

# -----------------------------------------------------------------------------
# 4. Subset to period, impute, drop constant columns
# -----------------------------------------------------------------------------
make_X <- function(tp_range) {
  idx  <- which(as.integer(row_meta$Timepoint) %in% tp_range)
  mat  <- log2_mat_full[idx, , drop = FALSE]
  meta <- row_meta[idx, ]

  # Drop constant columns
  keep <- apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)
  mat  <- mat[, keep, drop = FALSE]

  # Impute NA with column mean
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) mat[na_idx, j] <- mean(mat[, j], na.rm = TRUE)
  }

  list(X = mat, meta = meta)
}

early_data <- make_X(1:3)
late_data  <- make_X(4:6)

cat(sprintf("Early X: %d obs × %d features\n", nrow(early_data$X), ncol(early_data$X)))
cat(sprintf("Late  X: %d obs × %d features\n", nrow(late_data$X),  ncol(late_data$X)))

# -----------------------------------------------------------------------------
# 5. PCA (center + scale)
# -----------------------------------------------------------------------------
early_pca <- prcomp(early_data$X, center = TRUE, scale. = TRUE)
late_pca  <- prcomp(late_data$X,  center = TRUE, scale. = TRUE)

pct_var <- function(pca_obj) {
  v <- pca_obj$sdev^2
  round(v / sum(v) * 100, 1)
}
early_var <- pct_var(early_pca)
late_var  <- pct_var(late_pca)

# -----------------------------------------------------------------------------
# 6. Build score data frames
# -----------------------------------------------------------------------------
make_score_df <- function(pca_obj, meta, groups_df) {
  scores           <- as.data.frame(pca_obj$x[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$SubjTP    <- meta$SubjTP
  scores$Subject   <- meta$Subject
  scores$Timepoint <- meta$Timepoint

  scores <- merge(scores, groups_df[, c("Subject", "group")],
                  by = "Subject", all.x = TRUE)
  scores$group[is.na(scores$group)] <- "other"
  scores$group <- factor(scores$group,
                         levels = c("both", "X_only", "Y_only", "other"))
  scores
}

early_scores <- make_score_df(early_pca, early_data$meta, groups_df)
late_scores  <- make_score_df(late_pca,  late_data$meta,  groups_df)

# -----------------------------------------------------------------------------
# 7. Plot function
# -----------------------------------------------------------------------------
make_pca_plot <- function(scores, var_expl, period_label, tp_range) {
  x_lab <- sprintf("PC1 (%.1f%%)", var_expl[1])
  y_lab <- sprintf("PC2 (%.1f%%)", var_expl[2])

  tp_levels <- as.character(tp_range)
  tp_shapes <- TP_SHAPES[tp_levels]

  ggplot(scores, aes(x = PC1, y = PC2,
                     color = group, shape = Timepoint)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(size = 2.8, alpha = 0.80) +
    scale_color_manual(
      values = GROUP_COLORS,
      name   = sprintf("Response\n(ID %d)", TARGET_ID),
      labels = c(both   = "Both (X & Y)",
                 X_only = "X only (early↑)",
                 Y_only = "Y only (late↑)",
                 other  = "Other")
    ) +
    scale_shape_manual(
      values = tp_shapes,
      name   = "Timepoint"
    ) +
    labs(
      title    = sprintf("PCA score plot — %s", period_label),
      subtitle = "All features, log2-transformed; each point = one subject × timepoint",
      x = x_lab,
      y = y_lab
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.title       = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )
}

p_early <- make_pca_plot(early_scores, early_var, "Early (T1–T3)", 1:3)
p_late  <- make_pca_plot(late_scores,  late_var,  "Late  (T4–T6)", 4:6)

# -----------------------------------------------------------------------------
# 8. Save PDFs
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

early_pdf <- file.path(OUT_DIR, "pca_early_pertp.pdf")
late_pdf  <- file.path(OUT_DIR, "pca_late_pertp.pdf")

early_png <- file.path(OUT_DIR, "pca_early_pertp.png")
late_png  <- file.path(OUT_DIR, "pca_late_pertp.png")

ggsave(early_pdf, plot = p_early, width = 7.5, height = 5.5)
ggsave(late_pdf,  plot = p_late,  width = 7.5, height = 5.5)
ggsave(early_png, plot = p_early, width = 7.5, height = 5.5, dpi = 150)
ggsave(late_png,  plot = p_late,  width = 7.5, height = 5.5, dpi = 150)

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n  %s\n",
            early_pdf, late_pdf, early_png, late_png))

# -----------------------------------------------------------------------------
# 9. Variance explained summary
# -----------------------------------------------------------------------------
cat("\nEarly PCA — variance explained:\n")
for (i in 1:min(5, length(early_var)))
  cat(sprintf("  PC%d: %.1f%%\n", i, early_var[i]))
cat("\nLate PCA — variance explained:\n")
for (i in 1:min(5, length(late_var)))
  cat(sprintf("  PC%d: %.1f%%\n", i, late_var[i]))
