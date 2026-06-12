# =============================================================================
# plot_selected_features.R
# Plot boxplots for user-specified Alignment_IDs
#
# Reads p-values from a pre-computed pairwise_test_results CSV,
# and draws the same boxplots as boxplot_pairwise_production.R.
#
# Prerequisites:
#   Run load_data_production.R first (samplesheet / feat_meta / feat_mat).
#
# Usage:
#   1. Edit the Settings section below (TARGET_IDS, TRANSFORM, etc.)
#   2. source("src/r/plot_selected_features.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Settings — edit here
# -----------------------------------------------------------------------------
TARGET_IDS <- c(723, 724, 414, 415)          # Alignment_IDs to plot (integer or character)

TRANSFORM  <- "log2FC"               # "none" | "log2" | "log2FC"

# Path to the pairwise test CSV for the chosen transform
PAIRWISE_CSV <- sprintf(
  "data/production/processed/pairwise_test_results_%s_production.csv",
  TRANSFORM
)

OUT_PDF <- sprintf(
  "output/figures/selected_features_%s.pdf",
  TRANSFORM
)

# -----------------------------------------------------------------------------
# 0. Check required objects
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
}

if (!file.exists(PAIRWISE_CSV)) {
  stop("CSV not found: ", PAIRWISE_CSV,
       "\nRun boxplot_pairwise_production.R first to generate it.")
}

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork", "ggsignif")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# -----------------------------------------------------------------------------
# 2. Load pairwise test results
# -----------------------------------------------------------------------------
pairwise_res <- read.csv(PAIRWISE_CSV, stringsAsFactors = FALSE)
pairwise_res$Alignment_ID <- as.character(pairwise_res$Alignment_ID)

target_ids <- as.character(TARGET_IDS)

# Validate
missing_ids <- setdiff(target_ids, pairwise_res$Alignment_ID)
if (length(missing_ids) > 0) {
  warning("These Alignment_IDs are not in the CSV and will be skipped: ",
          paste(missing_ids, collapse = ", "))
  target_ids <- intersect(target_ids, pairwise_res$Alignment_ID)
}
if (length(target_ids) == 0) stop("No valid Alignment_IDs found.")

cat(sprintf("Plotting %d feature(s): %s\n", length(target_ids),
            paste(target_ids, collapse = ", ")))

# -----------------------------------------------------------------------------
# 3. Build averaged intensity matrix (same logic as boxplot_pairwise_production.R)
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

sub_mat  <- feat_mat[target_ids, , drop = FALSE]
avg_mat0 <- sapply(group_cols, function(cols) {
  rowMeans(sub_mat[, cols, drop = FALSE], na.rm = TRUE)
})
if (length(target_ids) == 1) {
  avg_mat0 <- matrix(avg_mat0, nrow = 1,
                     dimnames = list(target_ids, names(group_cols)))
}

# -----------------------------------------------------------------------------
# 4. Apply transform
# -----------------------------------------------------------------------------
avg_mat <- avg_mat0

if (TRANSFORM == "log2") {
  offsets <- apply(avg_mat, 1, function(x) {
    nz <- x[x > 0 & !is.na(x)]; if (length(nz) == 0) return(1); min(nz) / 2
  })
  avg_mat <- log2(avg_mat + offsets)
  y_label <- "log2 Intensity (avg)"

} else if (TRANSFORM == "log2FC") {
  offsets <- apply(avg_mat, 1, function(x) {
    nz <- x[x > 0 & !is.na(x)]; if (length(nz) == 0) return(1); min(nz) / 2
  })
  avg_mat <- log2(avg_mat + offsets)
  subjects <- unique(sub("-[0-9]+$", "", colnames(avg_mat)))
  for (subj in subjects) {
    for (ref_tp in c(1L, 4L)) {
      ref_col  <- paste0(subj, "-", ref_tp)
      if (!ref_col %in% colnames(avg_mat)) next
      ref_vals  <- avg_mat[, ref_col]
      tp_range  <- if (ref_tp == 1L) 1:3 else 4:6
      subj_cols <- intersect(paste0(subj, "-", tp_range), colnames(avg_mat))
      avg_mat[, subj_cols] <- sweep(avg_mat[, subj_cols, drop = FALSE], 1, ref_vals, "-")
    }
  }
  y_label <- "log2FC (T1–T3: vs T1 / T4–T6: vs T4)"

} else {
  y_label <- "Intensity (avg)"
}

# -----------------------------------------------------------------------------
# 5. Convert to long format
# -----------------------------------------------------------------------------
df_avg              <- as.data.frame(avg_mat)
df_avg$Alignment_ID <- rownames(df_avg)

df_long <- df_avg %>%
  pivot_longer(-Alignment_ID, names_to = "SubjTP", values_to = "intensity") %>%
  mutate(
    Subject   = sub("-[0-9]+$", "", SubjTP),
    Timepoint = as.integer(sub(".*-", "", SubjTP))
  )

# -----------------------------------------------------------------------------
# 6. Plot helpers (identical to boxplot_pairwise_production.R)
# -----------------------------------------------------------------------------
pval_to_stars <- function(p) {
  if (is.na(p))  return("n.s.")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")
}

make_boxplot_panel <- function(df_feat, timepoints, pairs, raw_pvals,
                               panel_title, y_label) {
  df_sub           <- df_feat[df_feat$Timepoint %in% timepoints, ]
  df_sub$Timepoint <- factor(df_sub$Timepoint, levels = timepoints)

  y_max   <- max(df_sub$intensity, na.rm = TRUE)
  y_min   <- min(df_sub$intensity, na.rm = TRUE)
  y_range <- max(y_max - y_min, .Machine$double.eps)
  step    <- y_range * 0.20
  y_top   <- y_max + step * (length(pairs) + 1.5)

  comparisons <- lapply(pairs, function(p) as.character(p))
  annotations <- sapply(raw_pvals, pval_to_stars)
  y_positions <- y_max + step * seq_along(pairs)

  ggplot(df_sub, aes(x = Timepoint, y = intensity)) +
    geom_boxplot(fill = "#E8F4FD", outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.0, alpha = 0.5, color = "#2C7BB6") +
    ggsignif::geom_signif(
      comparisons = comparisons,
      annotations = annotations,
      y_position  = y_positions,
      tip_length  = 0.01,
      textsize    = 3.5,
      vjust       = 0.4
    ) +
    scale_y_continuous(limits = c(y_min - y_range * 0.05, y_top)) +
    labs(title = panel_title, x = "Timepoint", y = y_label) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
}

make_feature_page <- function(aid, df_long, pairwise_res, y_label) {
  meta_row <- feat_meta[feat_meta$Alignment_ID == as.integer(aid), ]
  res_row  <- pairwise_res[pairwise_res$Alignment_ID == aid, ]

  # Display label
  is_msfinder <- isTRUE(meta_row$MSFINDER_annotated)
  if (is_msfinder) {
    display_name <- sprintf("%s [MSFINDER, score=%.2f]",
                            meta_row$MSFINDER_structure,
                            meta_row$MSFINDER_total_score)
  } else {
    display_name <- meta_row$Metabolite_name
  }

  cls <- if (nrow(res_row) > 0) res_row$classification else "N/A"

  title_str <- sprintf(
    "%s  |  ID=%s  Rt=%.2f  m/z=%.4f  [%s]",
    display_name, aid,
    round(as.numeric(meta_row$Rt_min), 2),
    round(as.numeric(meta_row$Mz), 4),
    cls
  )

  df_feat <- df_long[df_long$Alignment_ID == aid, ]

  if (nrow(res_row) > 0) {
    early_pvals <- unlist(res_row[, c("T1vT2", "T1vT3", "T2vT3")])
    late_pvals  <- unlist(res_row[, c("T4vT5", "T4vT6", "T5vT6")])
  } else {
    early_pvals <- rep(NA_real_, 3)
    late_pvals  <- rep(NA_real_, 3)
  }

  p_early <- make_boxplot_panel(
    df_feat, 1:3, list(c(1, 2), c(1, 3), c(2, 3)),
    early_pvals, "Timepoints 1-3 (early)", y_label
  )
  p_late <- make_boxplot_panel(
    df_feat, 4:6, list(c(4, 5), c(4, 6), c(5, 6)),
    late_pvals, "Timepoints 4-6 (late)", y_label
  )

  (p_early | p_late) +
    plot_annotation(
      title = title_str,
      theme = theme(plot.title = element_text(size = 9.5, face = "bold"))
    )
}

# -----------------------------------------------------------------------------
# 7. Generate PDF
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Saving PDF → %s\n", OUT_PDF))
pdf(OUT_PDF, width = 10, height = 5)
for (aid in target_ids) {
  p <- make_feature_page(aid, df_long, pairwise_res, y_label)
  print(p)
}
dev.off()
cat(sprintf("Done: %s\n", OUT_PDF))
