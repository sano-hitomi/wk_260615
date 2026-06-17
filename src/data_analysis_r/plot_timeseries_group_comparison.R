# =============================================================================
# plot_timeseries_group_comparison.R
# Time-series comparison of two subject groups for selected alignment IDs
#
# Groups:
#   Group 1 (blue)   : Subjects 14, 18, 37, 9
#   Group 2 (red)    : Subjects  5, 18, 32
#
# Alignments to plot: 414, 4263
#
# Transformation: log2FC (period-specific)
#   Early (T1–T3): reference = T1 per subject
#   Late  (T4–T6): reference = T4 per subject
#
# Output:
#   output/figures/timeseries_group_comparison.pdf
#
# Usage (RStudio):
#   1. Set working directory to project root (TMAO_data_analysis/)
#      Session > Set Working Directory > To Project Directory
#   2. source("src/r/plot_timeseries_group_comparison.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
GROUP1_SUBJECTS  <- c(14, 18, 37, 9)       # Group 1 (blue)
GROUP2_SUBJECTS  <- c(5, 18, 32)            # Group 2 (red)
TARGET_ALIGN_IDS <- c(414, 4263)            # Alignment IDs to plot
TRANSFORM        <- "log2FC"                # "none" | "log2" | "log2FC"
USE_DUMMY        <- TRUE                    # TRUE: dummy data / FALSE: production data

OUT_PDF <- "output/figures/timeseries_group_comparison.pdf"

COLOR1 <- "#2196F3"   # blue  – Group 1
COLOR2 <- "#FF5722"   # red   – Group 2
ALPHA  <- 0.6         # line transparency (approximated via col2rgb)
LWD_SUBJ  <- 1.5     # per-subject line width
LWD_MEAN  <- 3.0     # group mean line width

# -----------------------------------------------------------------------------
# 0. Load data (same as load_data.R, supports dummy / production)
# -----------------------------------------------------------------------------
data_dir <- if (USE_DUMMY) "data/dummy/processed" else "data/production/processed"

ss_path   <- file.path(data_dir, "samplesheet.csv")
meta_path <- file.path(data_dir, "feature_metadata.csv")
mat_path  <- file.path(data_dir, "feature_matrix.csv")

if (!file.exists(ss_path)) {
  stop(
    "Data not found: ", ss_path, "\n",
    "Set the working directory to the project root:\n",
    "  Session > Set Working Directory > To Project Directory"
  )
}

samplesheet <- read.csv(ss_path,  stringsAsFactors = FALSE)
feat_meta   <- read.csv(meta_path, stringsAsFactors = FALSE)
feat_mat    <- read.csv(mat_path,  row.names = 1,
                         check.names = FALSE, stringsAsFactors = FALSE)

# Convert feat_mat values to numeric
feat_mat[] <- lapply(feat_mat, as.numeric)

cat(sprintf("Loaded: %d samples, %d features\n", nrow(samplesheet), nrow(feat_mat)))

# -----------------------------------------------------------------------------
# 1. Restrict to biological samples (exclude rerun, BLK, QC)
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" &
  (is.na(samplesheet$rerun_suffix) | samplesheet$rerun_suffix == ""),
]

bio$Subject    <- as.numeric(bio$Subject)
bio$Timepoint  <- as.numeric(bio$Timepoint)
bio <- bio[!is.na(bio$Subject) & !is.na(bio$Timepoint), ]

# All subjects and timepoints present in data
all_subjects   <- sort(unique(bio$Subject))
all_timepoints <- sort(unique(bio$Timepoint))   # expected: 1–6

# -----------------------------------------------------------------------------
# 2. Helper: build averaged (over repeats) intensity matrix for one alignment
#    Returns named vector: names = "Subject-Timepoint", values = mean intensity
# -----------------------------------------------------------------------------
get_avg_vals <- function(align_id) {
  aid_char <- as.character(align_id)
  if (!aid_char %in% rownames(feat_mat)) {
    warning(sprintf("Alignment_ID %s not found in feat_mat.", aid_char))
    return(NULL)
  }

  row_vals   <- as.numeric(feat_mat[aid_char, ])
  names(row_vals) <- colnames(feat_mat)

  group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
  group_cols <- split(bio$label, group_key)
  group_cols <- lapply(group_cols, function(cols) intersect(cols, names(row_vals)))
  group_cols <- group_cols[sapply(group_cols, length) > 0]

  sapply(group_cols, function(cols) mean(row_vals[cols], na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# 3. Helper: apply log2FC transformation to a raw matrix
#    Input : raw_mat [subjects × timepoints], dimnames set
#    Output: fc_mat  [subjects × timepoints]
# -----------------------------------------------------------------------------
apply_transform <- function(raw_mat, transform = TRANSFORM) {
  nz     <- raw_mat[raw_mat > 0 & !is.na(raw_mat)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2_mat <- log2(raw_mat + offset)

  if (transform == "none") {
    return(list(mat = raw_mat, ylabel = "Intensity (avg)"))
  }

  if (transform == "log2") {
    return(list(mat = log2_mat, ylabel = "log2 Intensity (avg)"))
  }

  # log2FC: period-specific reference
  fc_mat <- log2_mat
  for (subj in rownames(log2_mat)) {
    # Early: T1–T3, reference = T1
    ref_early_col <- which(colnames(log2_mat) == "1")
    if (length(ref_early_col) > 0) {
      ref_early <- log2_mat[subj, ref_early_col]
      for (tp in intersect(c(1, 2, 3), as.numeric(colnames(log2_mat)))) {
        col_name <- as.character(tp)
        fc_mat[subj, col_name] <- log2_mat[subj, col_name] - ref_early
      }
    }
    # Late: T4–T6, reference = T4
    ref_late_col <- which(colnames(log2_mat) == "4")
    if (length(ref_late_col) > 0) {
      ref_late <- log2_mat[subj, ref_late_col]
      for (tp in intersect(c(4, 5, 6), as.numeric(colnames(log2_mat)))) {
        col_name <- as.character(tp)
        fc_mat[subj, col_name] <- log2_mat[subj, col_name] - ref_late
      }
    }
  }
  list(mat = fc_mat, ylabel = "log2FC (T1–T3: vs T1 / T4–T6: vs T4)")
}

# -----------------------------------------------------------------------------
# 4. Helper: semi-transparent colour (base R, no extra packages)
# -----------------------------------------------------------------------------
alpha_col <- function(hex, alpha = 0.5) {
  rgb_vals <- col2rgb(hex) / 255
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha)
}

# -----------------------------------------------------------------------------
# 5. Helper: draw one time-series panel for a given alignment ID
# -----------------------------------------------------------------------------
draw_panel <- function(align_id, avg_vals) {

  # --- build Subject × Timepoint matrix (raw intensities) ---
  subj_tp_mat <- matrix(
    NA_real_,
    nrow     = length(all_subjects),
    ncol     = length(all_timepoints),
    dimnames = list(as.character(all_subjects), as.character(all_timepoints))
  )
  for (subj in all_subjects) {
    for (tp in all_timepoints) {
      key <- paste(subj, tp, sep = "-")
      if (key %in% names(avg_vals))
        subj_tp_mat[as.character(subj), as.character(tp)] <- avg_vals[[key]]
    }
  }

  # --- transform ---
  res    <- apply_transform(subj_tp_mat)
  tmat   <- res$mat
  ylabel <- res$ylabel

  # --- subset to groups ---
  g1_rows <- intersect(as.character(GROUP1_SUBJECTS), rownames(tmat))
  g2_rows <- intersect(as.character(GROUP2_SUBJECTS), rownames(tmat))

  g1_mat <- tmat[g1_rows, , drop = FALSE]
  g2_mat <- tmat[g2_rows, , drop = FALSE]

  # --- y-axis range: union of both groups ---
  yrange <- range(c(g1_mat, g2_mat), na.rm = TRUE)
  if (!is.finite(diff(yrange)) || diff(yrange) == 0)
    yrange <- yrange + c(-1, 1)

  # --- metadata for title ---
  meta_row  <- feat_meta[feat_meta$Alignment_ID == align_id, ]
  name_str  <- if (nrow(meta_row) > 0) meta_row$Metabolite_name[1] else "Unknown"
  rt_str    <- if (nrow(meta_row) > 0) sprintf("%.2f", as.numeric(meta_row$Rt_min[1])) else "NA"
  mz_str    <- if (nrow(meta_row) > 0) sprintf("%.4f", as.numeric(meta_row$Mz[1]))    else "NA"

  title_str <- sprintf(
    "%s\nID=%d  Rt=%s  m/z=%s", name_str, align_id, rt_str, mz_str
  )

  # --- empty plot frame ---
  plot(NA,
       xlim  = c(1, length(all_timepoints)),
       ylim  = yrange,
       xlab  = "Timepoint",
       ylab  = ylabel,
       main  = title_str,
       xaxt  = "n",
       cex.main = 0.85,
       cex.lab  = 0.85)
  axis(1, at = seq_along(all_timepoints), labels = all_timepoints)
  abline(v = 3.5, lty = 2, col = "grey60")   # separator: early | late
  grid(nx = NA, ny = NULL, lty = 1, col = "grey92")

  # --- per-subject lines: Group 1 ---
  c1_alpha <- alpha_col(COLOR1, ALPHA)
  for (i in seq_len(nrow(g1_mat))) {
    lines(seq_along(all_timepoints), g1_mat[i, ],
          col = c1_alpha, lwd = LWD_SUBJ, type = "l")
    points(seq_along(all_timepoints), g1_mat[i, ],
           col = c1_alpha, pch = 16, cex = 0.6)
  }

  # --- per-subject lines: Group 2 ---
  c2_alpha <- alpha_col(COLOR2, ALPHA)
  for (i in seq_len(nrow(g2_mat))) {
    lines(seq_along(all_timepoints), g2_mat[i, ],
          col = c2_alpha, lwd = LWD_SUBJ, type = "l")
    points(seq_along(all_timepoints), g2_mat[i, ],
           col = c2_alpha, pch = 17, cex = 0.6)
  }

  # --- group mean lines ---
  g1_mean <- colMeans(g1_mat, na.rm = TRUE)
  g2_mean <- colMeans(g2_mat, na.rm = TRUE)

  lines(seq_along(all_timepoints), g1_mean,
        col = COLOR1, lwd = LWD_MEAN, type = "b", pch = 16, cex = 1.0)
  lines(seq_along(all_timepoints), g2_mean,
        col = COLOR2, lwd = LWD_MEAN, type = "b", pch = 17, cex = 1.0)

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# 6. Main: build PDF with one panel per alignment ID
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

n_plots  <- length(TARGET_ALIGN_IDS)
n_cols   <- min(n_plots, 2)
n_rows   <- ceiling(n_plots / n_cols)

pdf(OUT_PDF,
    width  = 6.5 * n_cols,
    height = 5.5 * n_rows + 1.0)   # +1 for legend space

layout(
  rbind(seq_len(n_plots), rep(n_plots + 1L, n_cols)),
  heights = c(5.5 * n_rows, 1.0)
)
par(mar = c(4.5, 4.5, 3.5, 1.5))

for (aid in TARGET_ALIGN_IDS) {
  avg_vals <- get_avg_vals(aid)
  if (is.null(avg_vals)) next
  draw_panel(aid, avg_vals)
}

# --- shared legend panel ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = c(
    sprintf("Group 1  (Sub %s)",  paste(sort(GROUP1_SUBJECTS), collapse = ", ")),
    sprintf("Group 2  (Sub %s)",  paste(sort(GROUP2_SUBJECTS), collapse = ", ")),
    "Group mean"
  ),
  col    = c(COLOR1, COLOR2, "black"),
  lwd    = c(LWD_SUBJ, LWD_SUBJ, LWD_MEAN),
  pch    = c(16, 17, NA),
  pt.cex = 0.9,
  lty    = 1,
  horiz  = TRUE,
  cex    = 1.0,
  bty    = "n"
)

dev.off()
cat(sprintf("\nDone. PDF saved to: %s\n", OUT_PDF))
