# =============================================================================
# plot_timeseries.R
# Time-series plots for metabolite features with 1–5 annotations
#
# Transformation options (TRANSFORM):
#   "none"    : raw intensity
#   "log2"    : log2(x + offset), offset = min non-zero value / 2 per feature
#   "log2FC"  : period-specific log2FC per subject
#               Early (T1–T3): T1 = 0 reference
#               Late  (T4–T6): T4 = 0 reference (independent intervention)
#
# Output:
#   output/figures/timeseries_identified_{transform}_dummy.pdf
#
# Usage (RStudio):
#   1. Run load_data.R first (or have samplesheet / feat_meta / feat_mat loaded)
#   2. Set TRANSFORM as desired
#   3. source("src/r/plot_timeseries.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
EXCLUDE_LOW_SCORE <- TRUE       # exclude "low score:" annotations
COUNT_MIN         <- 1          # minimum annotation count (inclusive)
COUNT_MAX         <- 5          # maximum annotation count (inclusive)
PLOTS_PER_PAGE    <- 6          # layout: 2 rows x 3 cols
TRANSFORM         <- "log2FC"   # "none" | "log2" | "log2FC"
INCLUDE_MSFINDER  <- TRUE       # TRUE: also plot MS-FINDER annotated Unknown features
OUT_PDF           <- sprintf("output/figures/timeseries_identified_%s_dummy.pdf", TRANSFORM)

# -----------------------------------------------------------------------------
# 0. Check required objects
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun load_data.R first.")
}

# -----------------------------------------------------------------------------
# 1. Select target features
#    (a) MS-DIAL identified features  : COUNT_MIN <= name count <= COUNT_MAX,
#                                       not "Unknown", optionally not "low score:"
#    (b) MS-FINDER annotated features : Unknown features that MS-FINDER annotated
#                                       (MSFINDER_annotated == TRUE), added when
#                                       INCLUDE_MSFINDER is TRUE
# -----------------------------------------------------------------------------
name_counts <- table(feat_meta$Metabolite_name)

target_names <- names(name_counts)[
  name_counts >= COUNT_MIN &
  name_counts <= COUNT_MAX &
  names(name_counts) != "Unknown"
]

if (EXCLUDE_LOW_SCORE) {
  target_names <- target_names[!grepl("^low score:", target_names)]
}

target_features <- feat_meta[feat_meta$Metabolite_name %in% target_names, ]

# Add MS-FINDER annotated Unknown features
if (INCLUDE_MSFINDER && isTRUE("MSFINDER_annotated" %in% colnames(feat_meta))) {
  msfinder_feats <- feat_meta[
    feat_meta$Metabolite_name == "Unknown" &
    !is.na(feat_meta$MSFINDER_annotated) &
    feat_meta$MSFINDER_annotated == TRUE, ]
  target_features <- rbind(target_features, msfinder_feats)
  cat(sprintf("MS-FINDER features added: %d\n", nrow(msfinder_feats)))
}

target_features <- target_features[order(target_features$Metabolite_name,
                                         target_features$Alignment_ID), ]

cat(sprintf("Target MS-DIAL names    : %d\n", length(target_names)))
cat(sprintf("Target features (total) : %d\n", nrow(target_features)))
cat(sprintf("Output PDF              : %s\n\n", OUT_PDF))

# -----------------------------------------------------------------------------
# 2. Prepare averaged matrix: rows = features, cols = Subject-Timepoint
#    Average over technical replicates (Repeat 1–3)
# -----------------------------------------------------------------------------
bio <- samplesheet[samplesheet$type == "biological", ]

group_key   <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols  <- split(bio$label, group_key)

# Keep only columns present in feat_mat
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

avg_mat <- sapply(group_cols, function(cols) {
  rowMeans(feat_mat[, cols, drop = FALSE], na.rm = TRUE)
})

# Ensure column order: Subject × Timepoint
subjects   <- sort(unique(bio$Subject))
timepoints <- sort(unique(bio$Timepoint))
col_order  <- as.vector(outer(subjects, timepoints,
                               FUN = function(s, t) paste(s, t, sep = "-")))
col_order  <- intersect(col_order, colnames(avg_mat))
avg_mat    <- avg_mat[, col_order, drop = FALSE]

# -----------------------------------------------------------------------------
# 2.5 Apply intensity transformation
# -----------------------------------------------------------------------------
if (TRANSFORM == "log2") {
  offsets <- apply(avg_mat, 1, function(x) {
    nz <- x[x > 0 & !is.na(x)]
    if (length(nz) == 0) return(1)
    min(nz) / 2
  })
  avg_mat <- log2(avg_mat + offsets)
  y_label <- "log2 Intensity (avg)"
  cat(sprintf("Transformation : log2(x + offset)\n"))

} else if (TRANSFORM == "log2FC") {
  offsets <- apply(avg_mat, 1, function(x) {
    nz <- x[x > 0 & !is.na(x)]
    if (length(nz) == 0) return(1)
    min(nz) / 2
  })
  avg_mat <- log2(avg_mat + offsets)
  # Period-specific reference:
  #   Early (T1–T3) → subtract T1 per subject
  #   Late  (T4–T6) → subtract T4 per subject (independent intervention)
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
  cat(sprintf("Transformation : log2FC  [T1 ref for T1-T3; T4 ref for T4-T6]\n"))

} else {
  y_label <- "Intensity (avg)"
  cat(sprintf("Transformation : none (raw intensity)\n"))
}

# -----------------------------------------------------------------------------
# 3. Draw time-series plots → PDF
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

n_features <- nrow(target_features)
n_cols_layout <- 3
n_rows_layout <- 2

pdf(OUT_PDF,
    width  = 4 * n_cols_layout,
    height = 3 * n_rows_layout)
par(mfrow = c(n_rows_layout, n_cols_layout),
    mar   = c(4, 4, 3, 1))

for (i in seq_len(n_features)) {

  align_id   <- target_features$Alignment_ID[i]
  metab_name <- target_features$Metabolite_name[i]
  rt         <- round(as.numeric(target_features$Rt_min[i]), 2)
  mz         <- round(as.numeric(target_features$Mz[i]), 4)

  # Display label: use MS-FINDER structure name for annotated Unknowns
  is_msfinder <- isTRUE(target_features$MSFINDER_annotated[i])
  if (is_msfinder) {
    msf_str   <- target_features$MSFINDER_structure[i]
    msf_score <- target_features$MSFINDER_total_score[i]
    display_label <- sprintf("%s [MSFINDER, score=%.2f]", msf_str, msf_score)
  } else {
    display_label <- metab_name
  }

  # Extract averaged values for this feature
  row_idx <- which(rownames(avg_mat) == as.character(align_id))
  if (length(row_idx) == 0) next

  vals <- avg_mat[row_idx, ]  # named vector: "001-1", "001-2", ...

  # Reshape to matrix: rows = subjects, cols = timepoints
  mat_plot <- matrix(NA, nrow = length(subjects), ncol = length(timepoints),
                     dimnames = list(subjects, as.character(timepoints)))
  for (s in subjects) {
    for (t in timepoints) {
      key <- paste(s, t, sep = "-")
      if (key %in% names(vals)) mat_plot[s, as.character(t)] <- vals[key]
    }
  }

  # Y-axis range
  yrange <- range(mat_plot, na.rm = TRUE)
  if (diff(yrange) == 0) yrange <- yrange + c(-1, 1)

  # Title
  title_str <- sprintf("%s\nID=%s  Rt=%.2f  m/z=%.4f", display_label, align_id, rt, mz)

  # Plot empty frame
  plot(NA,
       xlim = c(1, length(timepoints)),
       ylim = yrange,
       xlab = "Timepoint",
       ylab = y_label,
       main = title_str,
       xaxt = "n",
       cex.main = 0.75,
       cex.lab  = 0.85)
  axis(1, at = seq_along(timepoints), labels = timepoints)

  # One line per subject (light grey)
  for (s in subjects) {
    lines(seq_along(timepoints), mat_plot[s, ],
          col = "#AAAAAA", lwd = 0.8, type = "l")
  }

  # Grand mean across subjects (bold black)
  grand_mean <- colMeans(mat_plot, na.rm = TRUE)
  lines(seq_along(timepoints), grand_mean,
        col = "black", lwd = 2.5, type = "b", pch = 16, cex = 0.9)

  # Progress log every 100 features
  if (i %% 100 == 0) cat(sprintf("  %d / %d features plotted...\n", i, n_features))
}

dev.off()
cat(sprintf("\nDone. PDF saved to: %s\n", OUT_PDF))
