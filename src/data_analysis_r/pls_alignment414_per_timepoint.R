# =============================================================================
# pls_alignment414_per_timepoint.R
# PLS analysis with Alignment ID 414 as continuous response (Y)
# Each subject × timepoint is treated as an independent observation
#
# Design:
#   X : all metabolite features except TARGET_ID, log2-transformed
#       Rows = subject-timepoint pairs (e.g. "001-T1", "001-T2", ...)
#   Y : log2-transformed intensity of TARGET_ID at each subject-timepoint
#   Two independent analyses:
#     Early (Perturbation X) : T1, T2, T3  → n_subjects × 3 observations
#     Late  (Perturbation Y) : T4, T5, T6  → n_subjects × 3 observations
#
# Coloring:
#   4-group classification from GROUPS_CSV (per subject;
#   the same color is applied to all timepoints of a subject)
#
# Output:
#   output/figures/pls_early_414_pertp.pdf
#   output/figures/pls_late_414_pertp.pdf
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

# Point shape per timepoint
TP_SHAPES <- c("1" = 16, "2" = 17, "3" = 15,   # early: circle, triangle, square
               "4" = 16, "5" = 17, "6" = 15)   # late : same set

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
# 3. Build per-subject-timepoint intensity matrix (averaged over repeats)
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]
subjects   <- sort(unique(bio$Subject))
all_tp     <- 1:6

# Average over repeats for each Subject × Timepoint
group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

avg_list <- lapply(group_cols, function(cols) {
  rowMeans(feat_mat[, cols, drop = FALSE], na.rm = TRUE)
})
avg_mat_full <- do.call(cbind, avg_list)   # features × subject-timepoints

# log2 transform: per-feature offset
log2_mat_full <- apply(avg_mat_full, 1, function(x) {
  nz     <- x[x > 0 & !is.na(x)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2(x + offset)
})
# log2_mat_full: subject-timepoints × features (after apply transpose)
colnames(log2_mat_full) <- rownames(feat_mat)

# Row metadata: Subject and Timepoint for each row
row_meta <- data.frame(
  SubjTP    = names(group_cols),
  Subject   = sub("-[0-9]+$", "", names(group_cols)),
  Timepoint = as.integer(sub(".*-", "", names(group_cols))),
  stringsAsFactors = FALSE
)
rownames(log2_mat_full) <- row_meta$SubjTP

# -----------------------------------------------------------------------------
# 4. Subset to early / late and split X / Y
# -----------------------------------------------------------------------------
aid_char <- as.character(TARGET_ID)

make_XY <- function(tp_range) {
  idx <- which(row_meta$Timepoint %in% tp_range)
  mat <- log2_mat_full[idx, , drop = FALSE]
  meta_sub <- row_meta[idx, ]

  Y <- mat[, aid_char]
  X <- mat[, colnames(mat) != aid_char, drop = FALSE]

  # Drop constant columns
  keep <- apply(X, 2, function(x) var(x, na.rm = TRUE) > 0)
  X    <- X[, keep, drop = FALSE]

  # Impute remaining NA with column mean
  for (j in seq_len(ncol(X))) {
    na_idx <- is.na(X[, j])
    if (any(na_idx)) X[na_idx, j] <- mean(X[, j], na.rm = TRUE)
  }
  Y[is.na(Y)] <- mean(Y, na.rm = TRUE)

  list(X = X, Y = Y, meta = meta_sub)
}

early_data <- make_XY(1:3)
late_data  <- make_XY(4:6)

cat(sprintf("Early X: %d obs × %d features\n", nrow(early_data$X), ncol(early_data$X)))
cat(sprintf("Late  X: %d obs × %d features\n", nrow(late_data$X),  ncol(late_data$X)))

# -----------------------------------------------------------------------------
# 5. PLS (continuous Y)
# -----------------------------------------------------------------------------
run_pls <- function(data, ncomp) {
  ncomp_use <- min(ncomp, nrow(data$X) - 1)
  pls(X = data$X, Y = as.matrix(data$Y),
      ncomp = ncomp_use, scale = TRUE, mode = "regression")
}

early_pls <- run_pls(early_data, NCOMP)
late_pls  <- run_pls(late_data,  NCOMP)

var_exp <- function(pls_obj) {
  ve <- pls_obj$prop_expl_var$X
  if (is.null(ve)) ve <- pls_obj$explained_variance$X
  if (is.null(ve)) return(rep(NA_real_, pls_obj$ncomp))
  round(ve * 100, 1)
}
early_var <- var_exp(early_pls)
late_var  <- var_exp(late_pls)

# -----------------------------------------------------------------------------
# 6. Build score data frames
# -----------------------------------------------------------------------------
make_score_df <- function(pls_obj, meta, groups_df) {
  scores           <- as.data.frame(pls_obj$variates$X)
  colnames(scores) <- paste0("Comp", seq_len(ncol(scores)))
  scores$SubjTP    <- meta$SubjTP
  scores$Subject   <- meta$Subject
  scores$Timepoint <- as.character(meta$Timepoint)

  scores <- merge(scores, groups_df[, c("Subject", "group")],
                  by = "Subject", all.x = TRUE)
  scores$group[is.na(scores$group)] <- "other"
  scores$group <- factor(scores$group,
                         levels = c("both", "X_only", "Y_only", "other"))
  scores
}

early_scores <- make_score_df(early_pls, early_data$meta, groups_df)
late_scores  <- make_score_df(late_pls,  late_data$meta,  groups_df)

# -----------------------------------------------------------------------------
# 7. Plot function
# -----------------------------------------------------------------------------
make_pls_plot <- function(scores, var_expl, period_label, target_id, tp_range) {
  x_lab <- sprintf("Component 1 (%s%%)", var_expl[1])
  y_lab <- sprintf("Component 2 (%s%%)", var_expl[2])

  tp_levels  <- as.character(tp_range)
  tp_shapes  <- TP_SHAPES[tp_levels]

  ggplot(scores, aes(x = Comp1, y = Comp2,
                     color = group, shape = Timepoint)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(size = 2.8, alpha = 0.80) +
    scale_color_manual(
      values = GROUP_COLORS,
      name   = sprintf("Response\n(ID %d)", target_id),
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
      title    = sprintf("PLS score plot — %s  (Y = Alignment ID %d)",
                         period_label, target_id),
      subtitle = sprintf(
        "X: all features except ID %d, log2-transformed; each point = one subject × timepoint",
        target_id),
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

p_early <- make_pls_plot(early_scores, early_var,
                         "Early (T1–T3)", TARGET_ID, 1:3)
p_late  <- make_pls_plot(late_scores,  late_var,
                         "Late  (T4–T6)", TARGET_ID, 4:6)

# -----------------------------------------------------------------------------
# 8. Save PDFs
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

early_pdf <- file.path(OUT_DIR, sprintf("pls_early_%d_pertp.pdf", TARGET_ID))
late_pdf  <- file.path(OUT_DIR, sprintf("pls_late_%d_pertp.pdf",  TARGET_ID))

early_png <- file.path(OUT_DIR, sprintf("pls_early_%d_pertp.png", TARGET_ID))
late_png  <- file.path(OUT_DIR, sprintf("pls_late_%d_pertp.png",  TARGET_ID))

ggsave(early_pdf, plot = p_early, width = 7.5, height = 5.5)
ggsave(late_pdf,  plot = p_late,  width = 7.5, height = 5.5)
ggsave(early_png, plot = p_early, width = 7.5, height = 5.5, dpi = 150)
ggsave(late_png,  plot = p_late,  width = 7.5, height = 5.5, dpi = 150)

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n  %s\n",
            early_pdf, late_pdf, early_png, late_png))

# -----------------------------------------------------------------------------
# 9. Variance explained summary
# -----------------------------------------------------------------------------
cat(sprintf("\nEarly PLS — variance explained (X):\n"))
for (i in seq_along(early_var))
  cat(sprintf("  Comp%d: %.1f%%\n", i, early_var[i]))
cat(sprintf("\nLate PLS — variance explained (X):\n"))
for (i in seq_along(late_var))
  cat(sprintf("  Comp%d: %.1f%%\n", i, late_var[i]))
