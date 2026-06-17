# =============================================================================
# plot_timeseries_allsubjects.R
# All subjects (1–40) time-series comparison: Alignment 414 (left) vs 4263 (right)
#
# Layout  : A4 portrait (8.27" × 11.69"), 5 subjects per page × 8 pages
#           Each row = one subject; left panel = Alignment 414, right = Alignment 4263
# Y-axis  : globally unified (same limits across all 80 panels)
# Transform: log2FC (period-specific)
#              Early T1–T3 : reference = T1 per subject
#              Late  T4–T6 : reference = T4 per subject
#
# Output  : output/figures/timeseries_allsubjects_log2FC.pdf
#
# Usage (RStudio):
#   Set working directory to project root, then:
#   source("src/r/plot_timeseries_allsubjects.R")
# =============================================================================

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
TARGET_IDS      <- c(414, 4263)
SUBJECTS        <- 1:40
SUBJECTS_PER_PAGE <- 5
USE_DUMMY       <- TRUE        # FALSE → data/production/processed/
TRANSFORM       <- "log2FC"    # "none" | "log2" | "log2FC"

OUT_PDF <- "output/figures/timeseries_allsubjects_log2FC.pdf"

# A4 portrait dimensions in inches
A4_W <- 8.27
A4_H <- 11.69

# ---------------------------------------------------------------------------
# 0. Load data
# ---------------------------------------------------------------------------
data_dir  <- if (USE_DUMMY) "data/dummy/processed" else "data/production/processed"
ss_path   <- file.path(data_dir, "samplesheet.csv")
meta_path <- file.path(data_dir, "feature_metadata.csv")
mat_path  <- file.path(data_dir, "feature_matrix.csv")

if (!file.exists(ss_path))
  stop("Data not found: ", ss_path,
       "\nSet WD to project root: Session > Set Working Directory > To Project Directory")

samplesheet <- read.csv(ss_path,  stringsAsFactors = FALSE)
feat_meta   <- read.csv(meta_path, stringsAsFactors = FALSE)
feat_mat    <- read.csv(mat_path,  row.names = 1,
                        check.names = FALSE, stringsAsFactors = FALSE)
feat_mat[]  <- lapply(feat_mat, as.numeric)

cat(sprintf("Loaded: %d samples × %d features\n", nrow(samplesheet), nrow(feat_mat)))

# Biological samples, excluding reruns
bio <- samplesheet[
  samplesheet$type == "biological" &
  (samplesheet$rerun_suffix == "" | is.na(samplesheet$rerun_suffix)),
]
bio$Subject   <- as.numeric(bio$Subject)
bio$Timepoint <- as.numeric(bio$Timepoint)
bio <- bio[!is.na(bio$Subject) & !is.na(bio$Timepoint), ]

timepoints <- sort(unique(bio$Timepoint))   # 1–6

# ---------------------------------------------------------------------------
# 1. Build averaged (over repeats) Subject × Timepoint matrix per alignment
#    Returns matrix: rows = subjects, cols = timepoints (column names = tp)
# ---------------------------------------------------------------------------
build_raw_mat <- function(align_id, subjects) {
  aid_char <- as.character(align_id)
  if (!aid_char %in% rownames(feat_mat))
    stop(sprintf("Alignment_ID %d not found in feat_mat.", align_id))

  row_vals        <- as.numeric(feat_mat[aid_char, ])
  names(row_vals) <- colnames(feat_mat)

  mat <- matrix(
    NA_real_,
    nrow     = length(subjects),
    ncol     = length(timepoints),
    dimnames = list(as.character(subjects), as.character(timepoints))
  )
  for (subj in subjects) {
    for (tp in timepoints) {
      cols <- bio$label[bio$Subject == subj & bio$Timepoint == tp]
      cols <- intersect(cols, names(row_vals))
      if (length(cols) > 0)
        mat[as.character(subj), as.character(tp)] <- mean(row_vals[cols], na.rm = TRUE)
    }
  }
  mat
}

# ---------------------------------------------------------------------------
# 2. log2FC transformation (period-specific, per subject)
# ---------------------------------------------------------------------------
apply_log2fc <- function(raw_mat) {
  nz     <- raw_mat[raw_mat > 0 & !is.na(raw_mat)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  lmat   <- log2(raw_mat + offset)

  if (TRANSFORM == "none")  return(raw_mat)
  if (TRANSFORM == "log2")  return(lmat)

  # log2FC
  fc <- lmat
  for (subj in rownames(lmat)) {
    ref1 <- lmat[subj, "1"]
    for (tp in intersect(c(1,2,3), timepoints))
      fc[subj, as.character(tp)] <- lmat[subj, as.character(tp)] - ref1
    if ("4" %in% colnames(lmat)) {
      ref4 <- lmat[subj, "4"]
      for (tp in intersect(c(4,5,6), timepoints))
        fc[subj, as.character(tp)] <- lmat[subj, as.character(tp)] - ref4
    }
  }
  fc
}

# ---------------------------------------------------------------------------
# 3. Build transformed matrices for both alignments
# ---------------------------------------------------------------------------
cat("Computing averaged intensities for Alignment", TARGET_IDS[1], "...\n")
raw414  <- build_raw_mat(TARGET_IDS[1], SUBJECTS)
fc414   <- apply_log2fc(raw414)

cat("Computing averaged intensities for Alignment", TARGET_IDS[2], "...\n")
raw4263 <- build_raw_mat(TARGET_IDS[2], SUBJECTS)
fc4263  <- apply_log2fc(raw4263)

# ---------------------------------------------------------------------------
# 4. Global y-axis limits (same for all 80 panels)
# ---------------------------------------------------------------------------
global_range <- range(c(fc414, fc4263), na.rm = TRUE)
# Add 5% padding
pad <- diff(global_range) * 0.05
ylim_global <- global_range + c(-pad, pad)
if (!is.finite(diff(ylim_global))) ylim_global <- c(-1, 1)

cat(sprintf("Global y-axis (log2FC): [%.3f, %.3f]\n", ylim_global[1], ylim_global[2]))

# ---------------------------------------------------------------------------
# 5. Metadata strings for panel titles
# ---------------------------------------------------------------------------
get_label <- function(align_id) {
  m <- feat_meta[feat_meta$Alignment_ID == align_id, ]
  if (nrow(m) == 0) return(sprintf("ID %d", align_id))
  name <- m$Metabolite_name[1]
  rt   <- sprintf("%.2f", as.numeric(m$Rt_min[1]))
  mz   <- sprintf("%.4f", as.numeric(m$Mz[1]))
  sprintf("Align %d | %s | Rt=%s | m/z=%s", align_id, name, rt, mz)
}
label414  <- get_label(TARGET_IDS[1])
label4263 <- get_label(TARGET_IDS[2])

# ---------------------------------------------------------------------------
# 6. Helper: draw one subject's time-series panel
# ---------------------------------------------------------------------------
y_label <- switch(TRANSFORM,
  none    = "Intensity (avg)",
  log2    = "log2 Intensity",
  log2FC  = "log2FC"
)

draw_subject_panel <- function(fc_mat, subj, title_str) {
  vals <- fc_mat[as.character(subj), ]

  plot(NA,
       xlim  = c(0.8, length(timepoints) + 0.2),
       ylim  = ylim_global,
       xlab  = "Timepoint",
       ylab  = y_label,
       main  = title_str,
       xaxt  = "n",
       cex.main = 0.72,
       cex.lab  = 0.70,
       cex.axis = 0.65)
  axis(1, at = seq_along(timepoints), labels = timepoints, cex.axis = 0.65)

  # Background colour: magenta if any timepoint has log2FC < 0, else blue
  has_negative <- any(vals < 0, na.rm = TRUE)
  bg_col <- if (has_negative) "#FFE0EE" else "#E0EEFF"
  rect(0.8, ylim_global[1], length(timepoints) + 0.2, ylim_global[2],
       col = bg_col, border = NA)

  # Horizontal reference line at 0 (for log2FC)
  if (TRANSFORM == "log2FC") abline(h = 0, lty = 2, col = "grey50", lwd = 0.8)

  # Period separator
  abline(v = 3.5, lty = 3, col = "grey60", lwd = 0.8)

  box()

  # Data line
  lines(seq_along(timepoints), vals,
        col = "#333333", lwd = 1.4, type = "b", pch = 16, cex = 0.55)
}

# ---------------------------------------------------------------------------
# 7. Render PDF
# ---------------------------------------------------------------------------
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

n_pages <- ceiling(length(SUBJECTS) / SUBJECTS_PER_PAGE)

# Page margins (inches): bottom, left, top, right
PAGE_MAR <- c(0.5, 0.3, 0.6, 0.3)
# Panel margins (lines): bottom, left, top, right
PANEL_MAR <- c(3.2, 3.5, 2.2, 0.8)

pdf(OUT_PDF, width = A4_W, height = A4_H, paper = "a4")

for (page in seq_len(n_pages)) {
  # Subjects on this page
  idx_start <- (page - 1) * SUBJECTS_PER_PAGE + 1
  idx_end   <- min(page * SUBJECTS_PER_PAGE, length(SUBJECTS))
  page_subjects <- SUBJECTS[idx_start:idx_end]
  n_rows <- length(page_subjects)

  # Layout: n_rows rows × 2 cols + header strip at top
  # Use layout() for fine control
  layout(
    matrix(1:(n_rows * 2), nrow = n_rows, ncol = 2, byrow = TRUE),
    widths  = c(1, 1),
    heights = rep(1, n_rows)
  )
  par(mar = PANEL_MAR, oma = PAGE_MAR)

  for (subj in page_subjects) {
    # Left panel: Alignment 414
    draw_subject_panel(
      fc414, subj,
      sprintf("Sub %02d  |  %s", subj, label414)
    )
    # Right panel: Alignment 4263
    draw_subject_panel(
      fc4263, subj,
      sprintf("Sub %02d  |  %s", subj, label4263)
    )
  }

  # Page header (outer margin)
  mtext(
    sprintf(
      "Time-series: Alignment %d (left) vs %d (right)   |   Transform: %s   |   Page %d / %d",
      TARGET_IDS[1], TARGET_IDS[2], TRANSFORM, page, n_pages
    ),
    side = 3, outer = TRUE, line = 0.1, cex = 0.65, col = "grey30"
  )
  mtext(
    sprintf("Global y-axis [%.2f, %.2f]   |   Grey = Early (T1–T3)   |   Blue = Late (T4–T6)",
            ylim_global[1], ylim_global[2]),
    side = 1, outer = TRUE, line = 0.0, cex = 0.58, col = "grey40"
  )
}

dev.off()
cat(sprintf("\nDone. %d pages saved to: %s\n", n_pages, OUT_PDF))
