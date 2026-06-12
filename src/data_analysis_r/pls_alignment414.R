# =============================================================================
# pls_alignment414.R
# PLS analysis with Alignment ID 414 as continuous response (Y)
#
# Design:
#   X : all metabolite features except TARGET_ID, log2-transformed,
#       per-subject average within each period
#   Y : log2-transformed intensity of TARGET_ID,
#       per-subject average within each period
#   Two independent analyses:
#     Early (Perturbation X) : average of T1, T2, T3 per subject
#     Late  (Perturbation Y) : average of T4, T5, T6 per subject
#
# Coloring:
#   4-group classification loaded from GROUPS_CSV
#   (output of classify_subjects_by_response.R)
#
# Output:
#   output/figures/pls_early_414.pdf
#   output/figures/pls_late_414.pdf
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
NCOMP      <- 2          # number of PLS components to compute
OUT_DIR    <- "output/figures"

GROUP_COLORS <- c(
  both   = "#E41A1C",   # red
  X_only = "#377EB8",   # blue
  Y_only = "#4DAF4A",   # green
  other  = "#984EA3"    # purple
)

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
# 3. Build per-subject period-averaged intensity matrix (log2)
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]
subjects <- sort(unique(bio$Subject))

# Helper: average over repeats for given subject × timepoints
period_avg <- function(subj_vec, tp_vec) {
  # Returns matrix: rows = subjects, cols = feature rows in feat_mat
  # averaged over the given timepoints
  mat <- sapply(subj_vec, function(subj) {
    cols <- bio$label[bio$Subject == subj & bio$Timepoint %in% tp_vec]
    cols <- intersect(cols, colnames(feat_mat))
    if (length(cols) == 0) return(rep(NA_real_, nrow(feat_mat)))
    rowMeans(feat_mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  # mat: features × subjects → transpose to subjects × features
  t(mat)
}

# Raw averaged matrices
early_raw <- period_avg(subjects, 1:3)   # subjects × features
late_raw  <- period_avg(subjects, 4:6)

rownames(early_raw) <- subjects
rownames(late_raw)  <- subjects
colnames(early_raw) <- rownames(feat_mat)
colnames(late_raw)  <- rownames(feat_mat)

# log2 transform: per-feature offset = min non-zero / 2 across all data
log2_transform <- function(mat) {
  # Apply per-column (feature) offset
  result <- mat
  for (j in seq_len(ncol(mat))) {
    x  <- mat[, j]
    nz <- x[x > 0 & !is.na(x)]
    offset <- if (length(nz) > 0) min(nz) / 2 else 1
    result[, j] <- log2(x + offset)
  }
  result
}

early_log2 <- log2_transform(early_raw)
late_log2  <- log2_transform(late_raw)

# -----------------------------------------------------------------------------
# 4. Split X and Y
# -----------------------------------------------------------------------------
aid_char <- as.character(TARGET_ID)

# Y: log2 intensity of target feature
early_Y <- early_log2[, aid_char]
late_Y  <- late_log2[, aid_char]

# X: all other features; drop constant columns
early_X <- early_log2[, colnames(early_log2) != aid_char, drop = FALSE]
late_X  <- late_log2[,  colnames(late_log2)  != aid_char, drop = FALSE]

drop_constant <- function(mat) {
  keep <- apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)
  mat[, keep, drop = FALSE]
}
early_X <- drop_constant(early_X)
late_X  <- drop_constant(late_X)

# Replace any remaining NA with column mean
impute_colmean <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) mat[na_idx, j] <- mean(mat[, j], na.rm = TRUE)
  }
  mat
}
early_X <- impute_colmean(early_X)
late_X  <- impute_colmean(late_X)
early_Y[is.na(early_Y)] <- mean(early_Y, na.rm = TRUE)
late_Y[is.na(late_Y)]   <- mean(late_Y,  na.rm = TRUE)

cat(sprintf("Early X: %d subjects × %d features\n", nrow(early_X), ncol(early_X)))
cat(sprintf("Late  X: %d subjects × %d features\n", nrow(late_X),  ncol(late_X)))

# -----------------------------------------------------------------------------
# 5. PLS (continuous Y)
# -----------------------------------------------------------------------------
ncomp_use <- min(NCOMP, nrow(early_X) - 1)

early_pls <- pls(X = early_X, Y = as.matrix(early_Y),
                 ncomp = ncomp_use, scale = TRUE, mode = "regression")
late_pls  <- pls(X = late_X,  Y = as.matrix(late_Y),
                 ncomp = ncomp_use, scale = TRUE, mode = "regression")

# Variance explained by each component
var_exp <- function(pls_obj) {
  # sum of squared scores / total variance in X
  expl <- pls_obj$explained_variance$X
  round(expl * 100, 1)
}
early_var <- var_exp(early_pls)
late_var  <- var_exp(late_pls)

# -----------------------------------------------------------------------------
# 6. Build score data frames for plotting
# -----------------------------------------------------------------------------
make_score_df <- function(pls_obj, subjects, groups_df, var_expl) {
  scores <- as.data.frame(pls_obj$variates$X)
  colnames(scores) <- paste0("Comp", seq_len(ncol(scores)))
  scores$Subject <- subjects
  scores <- merge(scores, groups_df[, c("Subject", "group")],
                  by = "Subject", all.x = TRUE)
  scores$group[is.na(scores$group)] <- "other"
  scores$group <- factor(scores$group,
                         levels = c("both", "X_only", "Y_only", "other"))
  scores
}

early_scores <- make_score_df(early_pls, subjects, groups_df, early_var)
late_scores  <- make_score_df(late_pls,  subjects, groups_df, late_var)

# -----------------------------------------------------------------------------
# 7. Plot function
# -----------------------------------------------------------------------------
make_pls_plot <- function(scores, var_expl, period_label, target_id) {
  x_lab <- sprintf("Component 1 (%s%%)", var_expl[1])
  y_lab <- sprintf("Component 2 (%s%%)", var_expl[2])

  ggplot(scores, aes(x = Comp1, y = Comp2, color = group, label = Subject)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(size = 3, alpha = 0.85) +
    geom_text(size = 2.5, vjust = -0.8, hjust = 0.5, show.legend = FALSE) +
    scale_color_manual(
      values = GROUP_COLORS,
      name   = sprintf("Response\n(ID %d)", target_id),
      labels = c(both   = "Both (X & Y)",
                 X_only = "X only (early↑)",
                 Y_only = "Y only (late↑)",
                 other  = "Other")
    ) +
    labs(
      title    = sprintf("PLS score plot — %s period  (Y = Alignment ID %d)",
                         period_label, target_id),
      subtitle = sprintf("X: all features except ID %d, log2-transformed, per-subject average",
                         target_id),
      x = x_lab,
      y = y_lab
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(size = 10, face = "bold"),
      plot.subtitle     = element_text(size = 8, color = "grey40")
    )
}

p_early <- make_pls_plot(early_scores, early_var, "Early (T1–T3)", TARGET_ID)
p_late  <- make_pls_plot(late_scores,  late_var,  "Late  (T4–T6)", TARGET_ID)

# -----------------------------------------------------------------------------
# 8. Save PDFs
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

early_pdf <- file.path(OUT_DIR, sprintf("pls_early_%d.pdf", TARGET_ID))
late_pdf  <- file.path(OUT_DIR, sprintf("pls_late_%d.pdf",  TARGET_ID))

ggsave(early_pdf, plot = p_early, width = 7, height = 5.5)
ggsave(late_pdf,  plot = p_late,  width = 7, height = 5.5)

cat(sprintf("\nSaved:\n  %s\n  %s\n", early_pdf, late_pdf))

# -----------------------------------------------------------------------------
# 9. Print variance explained summary
# -----------------------------------------------------------------------------
cat(sprintf("\nEarly PLS — variance explained by X:\n"))
for (i in seq_along(early_var))
  cat(sprintf("  Comp%d: %.1f%%\n", i, early_var[i]))
cat(sprintf("\nLate PLS — variance explained by X:\n"))
for (i in seq_along(late_var))
  cat(sprintf("  Comp%d: %.1f%%\n", i, late_var[i]))
