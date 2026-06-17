# =============================================================================
# plot_timeseries_allsubjects_gg.R
# All subjects (1–40): Alignment 414 (left) vs 4263 (right), ggplot2 version
#
# Layout  : A4 portrait, 5 subjects per page × 8 pages
#           facet_grid(subject ~ alignment_id)
# Y-axis  : globally unified across all 80 panels
# Transform: log2FC (Early T1–T3 vs T1; Late T4–T6 vs T4)
# BG colour: light magenta if any timepoint < 0, else light blue
#
# Output  : output/figures/timeseries_allsubjects_gg.pdf
#
# Usage   : source("src/r/plot_timeseries_allsubjects_gg.R")
#           (working directory = project root)
# =============================================================================

# ---------------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------------
pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------
TARGET_IDS        <- c(414L, 4263L)
SUBJECTS          <- 1:40
SUBJECTS_PER_PAGE <- 5
USE_DUMMY         <- TRUE
OUT_PDF           <- "output/figures/timeseries_allsubjects_gg.pdf"
A4_W <- 8.27;  A4_H <- 11.69   # inches, portrait

# ---------------------------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------------------------
data_dir <- if (USE_DUMMY) "data/dummy/processed" else "data/production/processed"
if (!file.exists(file.path(data_dir, "samplesheet.csv")))
  stop("Working directory must be the project root (TMAO_data_analysis/).")

samplesheet <- read.csv(file.path(data_dir, "samplesheet.csv"),
                        stringsAsFactors = FALSE)
feat_mat    <- read.csv(file.path(data_dir, "feature_matrix.csv"),
                        row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
feat_mat[]  <- lapply(feat_mat, as.numeric)

bio <- samplesheet[
  samplesheet$type == "biological" &
  (is.na(samplesheet$rerun_suffix) | samplesheet$rerun_suffix == ""), ]
bio$Subject   <- as.integer(bio$Subject)
bio$Timepoint <- as.integer(bio$Timepoint)
bio <- bio[!is.na(bio$Subject) & !is.na(bio$Timepoint), ]

timepoints <- sort(unique(bio$Timepoint))
cat(sprintf("Loaded: %d biological samples, %d features\n",
            nrow(bio), nrow(feat_mat)))

# ---------------------------------------------------------------------------
# 3. Build long-format data frame with log2FC values
# ---------------------------------------------------------------------------
build_long <- function(align_id) {
  aid <- as.character(align_id)
  row_vals        <- as.numeric(feat_mat[aid, ])
  names(row_vals) <- colnames(feat_mat)

  # Average over technical replicates → one value per Subject × Timepoint
  rows <- lapply(SUBJECTS, function(subj) {
    lapply(timepoints, function(tp) {
      cols <- bio$label[bio$Subject == subj & bio$Timepoint == tp]
      cols <- intersect(cols, names(row_vals))
      val  <- if (length(cols) > 0) mean(row_vals[cols], na.rm = TRUE) else NA_real_
      data.frame(Subject = subj, Timepoint = tp, raw = val)
    })
  })
  df <- do.call(rbind, do.call(c, rows))

  # log2 transform with global offset
  nz     <- df$raw[df$raw > 0 & !is.na(df$raw)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  df$log2val <- log2(df$raw + offset)

  # Period-specific log2FC per subject
  df <- df[order(df$Subject, df$Timepoint), ]
  df$log2FC <- NA_real_
  for (subj in SUBJECTS) {
    idx <- df$Subject == subj
    ref1 <- df$log2val[idx & df$Timepoint == 1]
    ref4 <- df$log2val[idx & df$Timepoint == 4]
    if (length(ref1) == 1 && !is.na(ref1))
      df$log2FC[idx & df$Timepoint %in% 1:3] <-
        df$log2val[idx & df$Timepoint %in% 1:3] - ref1
    if (length(ref4) == 1 && !is.na(ref4))
      df$log2FC[idx & df$Timepoint %in% 4:6] <-
        df$log2val[idx & df$Timepoint %in% 4:6] - ref4
  }

  df$Alignment <- align_id
  df
}

cat("Computing log2FC for both alignments...\n")
long_df <- rbind(build_long(TARGET_IDS[1]), build_long(TARGET_IDS[2]))
long_df$Alignment <- factor(long_df$Alignment, levels = TARGET_IDS)

# ---------------------------------------------------------------------------
# 4. Global y-axis limits
# ---------------------------------------------------------------------------
ylim_global <- range(long_df$log2FC, na.rm = TRUE)
pad <- diff(ylim_global) * 0.05
ylim_global <- ylim_global + c(-pad, pad)
cat(sprintf("Global y-axis: [%.3f, %.3f]\n", ylim_global[1], ylim_global[2]))

# ---------------------------------------------------------------------------
# 5. Background colour per Subject × Alignment
#    light magenta if any timepoint has log2FC < 0, else light blue
# ---------------------------------------------------------------------------
bg_df <- long_df |>
  group_by(Subject, Alignment) |>
  summarise(has_neg = any(log2FC < 0, na.rm = TRUE), .groups = "drop") |>
  mutate(bg_col = ifelse(has_neg, "#FFE0EE", "#E0EEFF"))

# Expand to rect coordinates spanning the full x range
bg_df$xmin <- min(timepoints) - 0.4
bg_df$xmax <- max(timepoints) + 0.4
bg_df$ymin <- ylim_global[1]
bg_df$ymax <- ylim_global[2]

# ---------------------------------------------------------------------------
# 6. Build ggplot for a subset of subjects (one page)
# ---------------------------------------------------------------------------
make_page_plot <- function(page_subjects) {
  ld  <- long_df[long_df$Subject %in% page_subjects, ]
  bld <- bg_df  [bg_df$Subject   %in% page_subjects, ]

  # Subject label: "Sub 01" etc.
  ld$subj_label  <- factor(sprintf("Sub %02d", ld$Subject),
                            levels = sprintf("Sub %02d", page_subjects))
  bld$subj_label <- factor(sprintf("Sub %02d", bld$Subject),
                            levels = sprintf("Sub %02d", page_subjects))

  ggplot(ld, aes(x = Timepoint, y = log2FC)) +

    # Background rectangles
    geom_rect(data = bld,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = bg_col),
              inherit.aes = FALSE, alpha = 1) +
    scale_fill_identity() +

    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 3.5, linetype = "dotted",
               colour = "grey60", linewidth = 0.4) +

    # Time-series line + points
    geom_line(colour = "#333333", linewidth = 0.8) +
    geom_point(colour = "#333333", size = 1.2) +

    # Facets: rows = subject, cols = alignment
    facet_grid(subj_label ~ Alignment, switch = "y") +

    # Axes
    scale_x_continuous(breaks = timepoints) +
    coord_cartesian(ylim = ylim_global, xlim = c(0.6, 6.4)) +

    labs(x = "Timepoint", y = "log2FC") +

    theme_bw(base_size = 8) +
    theme(
      strip.text.x      = element_text(size = 9, face = "bold"),
      strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
      strip.placement   = "outside",
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.spacing     = unit(3, "pt"),
      axis.text         = element_text(size = 7),
      axis.title        = element_text(size = 8),
      plot.margin       = margin(4, 4, 4, 4, "pt"),
      legend.position   = "none"
    )
}

# ---------------------------------------------------------------------------
# 7. Render multi-page PDF
# ---------------------------------------------------------------------------
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

n_pages <- ceiling(length(SUBJECTS) / SUBJECTS_PER_PAGE)

pdf(OUT_PDF, width = A4_W, height = A4_H, paper = "a4")
for (page in seq_len(n_pages)) {
  idx_start     <- (page - 1) * SUBJECTS_PER_PAGE + 1
  idx_end       <- min(page  * SUBJECTS_PER_PAGE, length(SUBJECTS))
  page_subjects <- SUBJECTS[idx_start:idx_end]

  p <- make_page_plot(page_subjects)
  print(p)

  cat(sprintf("  Page %d / %d done\n", page, n_pages))
}
dev.off()

cat(sprintf("\nDone. %d pages → %s\n", n_pages, OUT_PDF))
