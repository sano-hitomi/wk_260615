# =============================================================================
# plsda_per_timepoint.R
# PLS-DA with 4-group response classification as Y
# Each subject × timepoint treated as an independent observation
#
# Design:
#   Predictors : all metabolite features, log2-transformed, averaged over repeats
#                Rows = subject-timepoint pairs
#   Response   : 4-group label (both / X_only / Y_only / other)
#                loaded from GROUPS_CSV (output of classify_subjects_by_response.R)
#   Two independent analyses:
#     Early (Perturbation X) : T1, T2, T3
#     Late  (Perturbation Y) : T4, T5, T6
#
# Plot:
#   Score plot (comp1 vs comp2)
#     Color : 4 groups
#     Shape : timepoint within the period
#   Ellipses: 95% confidence ellipse per group (if ≥3 observations)
#
# Output:
#   output/figures/plsda_early_pertp.pdf
#   output/figures/plsda_late_pertp.pdf
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
NCOMP      <- 2
OUT_DIR    <- "output/figures"

GROUP_COLORS <- c(
  both   = "#E41A1C",
  X_only = "#377EB8",
  Y_only = "#4DAF4A",
  other  = "#984EA3"
)
GROUP_LABELS <- c(
  both   = "Both (X & Y)",
  X_only = "X only (early↑)",
  Y_only = "Y only (late↑)",
  other  = "Other"
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
pkgs <- c("ggplot2", "mixOmics")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "mixOmics") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("mixOmics", ask = FALSE)
    } else {
      install.packages(p)
    }
  }
  library(p, character.only = TRUE)
}

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

log2_mat_full <- t(apply(avg_mat_full, 1, function(x) {
  nz     <- x[x > 0 & !is.na(x)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2(x + offset)
}))
log2_mat_full <- t(log2_mat_full)          # obs × features
colnames(log2_mat_full) <- rownames(feat_mat)

row_meta <- data.frame(
  SubjTP    = names(group_cols),
  Subject   = sub("-[0-9]+$", "", names(group_cols)),
  Timepoint = as.character(as.integer(sub(".*-", "", names(group_cols)))),
  stringsAsFactors = FALSE
)
rownames(log2_mat_full) <- row_meta$SubjTP

# -----------------------------------------------------------------------------
# 4. Subset to period, attach group labels, drop constant / impute
# -----------------------------------------------------------------------------
make_data <- function(tp_range) {
  idx  <- which(as.integer(row_meta$Timepoint) %in% tp_range)
  mat  <- log2_mat_full[idx, , drop = FALSE]
  meta <- row_meta[idx, ]

  # Attach group label (per subject → applied to all its timepoints)
  meta <- merge(meta, groups_df[, c("Subject", "group")],
                by = "Subject", all.x = TRUE)
  meta$group[is.na(meta$group)] <- "other"
  # Restore original row order
  meta <- meta[match(row_meta$SubjTP[idx], meta$SubjTP), ]

  # Drop constant columns
  keep <- apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)
  mat  <- mat[, keep, drop = FALSE]

  # Impute NA with column mean
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) mat[na_idx, j] <- mean(mat[, j], na.rm = TRUE)
  }

  Y <- factor(meta$group, levels = c("both", "X_only", "Y_only", "other"))

  list(X = mat, Y = Y, meta = meta)
}

early_data <- make_data(1:3)
late_data  <- make_data(4:6)

cat(sprintf("Early — obs: %d, features: %d\n", nrow(early_data$X), ncol(early_data$X)))
cat(sprintf("Late  — obs: %d, features: %d\n", nrow(late_data$X),  ncol(late_data$X)))
cat("\nEarly group counts:\n"); print(table(early_data$Y))
cat("\nLate  group counts:\n"); print(table(late_data$Y))

# -----------------------------------------------------------------------------
# 5. PLS-DA
# -----------------------------------------------------------------------------
ncomp_use <- min(NCOMP, nlevels(early_data$Y) - 1,
                 nrow(early_data$X) - 1)
ncomp_use <- max(ncomp_use, 2)   # need at least 2 for a 2D score plot

early_plsda <- plsda(X = early_data$X, Y = early_data$Y,
                     ncomp = ncomp_use, scale = TRUE)
late_plsda  <- plsda(X = late_data$X,  Y = late_data$Y,
                     ncomp = ncomp_use, scale = TRUE)

# Variance explained by each component
# mixOmics >= 6.x uses $prop_expl_var$X; older versions use $explained_variance$X
var_exp <- function(obj) {
  ve <- obj$prop_expl_var$X
  if (is.null(ve)) ve <- obj$explained_variance$X
  if (is.null(ve)) return(rep(NA_real_, obj$ncomp))
  round(ve * 100, 1)
}
early_var <- var_exp(early_plsda)
late_var  <- var_exp(late_plsda)

# -----------------------------------------------------------------------------
# 6. Build score data frames
# -----------------------------------------------------------------------------
make_score_df <- function(obj, data) {
  scores           <- as.data.frame(obj$variates$X[, 1:2])
  colnames(scores) <- c("Comp1", "Comp2")
  scores$SubjTP    <- data$meta$SubjTP
  scores$Subject   <- data$meta$Subject
  scores$Timepoint <- data$meta$Timepoint
  scores$group     <- factor(data$meta$group,
                             levels = c("both", "X_only", "Y_only", "other"))
  scores
}

early_scores <- make_score_df(early_plsda, early_data)
late_scores  <- make_score_df(late_plsda,  late_data)

# -----------------------------------------------------------------------------
# 7. Plot function
# -----------------------------------------------------------------------------
make_plsda_plot <- function(scores, var_expl, period_label, tp_range) {
  x_lab     <- sprintf("Component 1 (%s%%)", var_expl[1])
  y_lab     <- sprintf("Component 2 (%s%%)", var_expl[2])
  tp_levels <- as.character(tp_range)
  tp_shapes <- TP_SHAPES[tp_levels]

  ggplot(scores, aes(x = Comp1, y = Comp2,
                     color = group, shape = Timepoint)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    stat_ellipse(aes(group = group, color = group),
                 type = "norm", level = 0.95,
                 linewidth = 0.7, linetype = "dashed",
                 show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(
      values = GROUP_COLORS,
      name   = "Group",
      labels = GROUP_LABELS
    ) +
    scale_shape_manual(
      values = tp_shapes,
      name   = "Timepoint"
    ) +
    labs(
      title    = sprintf("PLS-DA score plot — %s", period_label),
      subtitle = sprintf(
        "Response: 4-group label (ID %d); Predictors: all features, log2-transformed",
        TARGET_ID),
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

p_early <- make_plsda_plot(early_scores, early_var, "Early (T1–T3)", 1:3)
p_late  <- make_plsda_plot(late_scores,  late_var,  "Late  (T4–T6)", 4:6)

# -----------------------------------------------------------------------------
# 8. Save PDFs
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

early_pdf <- file.path(OUT_DIR, "plsda_early_pertp.pdf")
late_pdf  <- file.path(OUT_DIR, "plsda_late_pertp.pdf")
early_png <- file.path(OUT_DIR, "plsda_early_pertp.png")
late_png  <- file.path(OUT_DIR, "plsda_late_pertp.png")

ggsave(early_pdf, plot = p_early, width = 7.5, height = 5.5)
ggsave(late_pdf,  plot = p_late,  width = 7.5, height = 5.5)
ggsave(early_png, plot = p_early, width = 7.5, height = 5.5, dpi = 150)
ggsave(late_png,  plot = p_late,  width = 7.5, height = 5.5, dpi = 150)

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n  %s\n",
            early_pdf, late_pdf, early_png, late_png))

# -----------------------------------------------------------------------------
# 9. Variance explained summary
# -----------------------------------------------------------------------------
cat("\nEarly PLS-DA — variance explained (Predictors):\n")
for (i in seq_along(early_var))
  cat(sprintf("  Comp%d: %.1f%%\n", i, early_var[i]))
cat("\nLate PLS-DA — variance explained (Predictors):\n")
for (i in seq_along(late_var))
  cat(sprintf("  Comp%d: %.1f%%\n", i, late_var[i]))
