# =============================================================================
# boxplot_pairwise_dummy.R
# Pairwise timepoint comparisons (paired Wilcoxon) + box plots (dummy data)
#
# Comparisons:
#   Early group (period = "early"): T1vT2, T1vT3, T2vT3
#   Late  group (period = "late" ): T4vT5, T4vT6, T5vT6
#
# Target metabolites: same filter as plot_timeseries.R
#   (COUNT_MIN <= annotation count <= COUNT_MAX, not "Unknown", not "low score:")
#   When INCLUDE_MSFINDER = TRUE, Unknown features annotated by MS-FINDER are
#   also included (MSFINDER_annotated == TRUE in feat_meta, added by load_data.R)
#
# Statistical method:
#   Paired Wilcoxon signed-rank test, per-feature × per-pair
#   BH (Benjamini-Hochberg) FDR correction across ALL tests
#
# Classification (based on adjusted p-value < ALPHA):
#   "both"        : significant in ANY early pair AND ANY late pair
#   "early_only"  : significant in ANY early pair only
#   "late_only"   : significant in ANY late  pair only
#   "ns"          : not significant in either group
#
# Transformation options (TRANSFORMS):
#   "none"    : raw intensity
#   "log2"    : log2(x + offset), offset = min non-zero value / 2 per feature
#   "log2FC"  : period-specific log2FC per subject
#               Early (T1–T3): T1 = 0 reference
#               Late  (T4–T6): T4 = 0 reference (independent intervention)
#               Note: p-values algebraically identical to "log2"
#   All can be run in one execution by listing them in TRANSFORMS.
#
# Output (per transform):
#   data/dummy/processed/pairwise_test_results_{transform}_dummy.csv
#   output/figures/boxplot_sig_both_{transform}_dummy.pdf
#   output/figures/boxplot_sig_early_only_{transform}_dummy.pdf
#   output/figures/boxplot_sig_late_only_{transform}_dummy.pdf
#
# Usage (RStudio):
#   1. Run load_data.R first (samplesheet / feat_meta / feat_mat must exist)
#   2. source("src/r/boxplot_pairwise_dummy.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
COUNT_MIN         <- 1                  # minimum annotation count (inclusive)
COUNT_MAX         <- 5                  # maximum annotation count (inclusive)
EXCLUDE_LOW_SCORE <- TRUE               # exclude "low score:" annotations
ALPHA             <- 0.05              # significance threshold (BH-adjusted p-value)
TRANSFORMS        <- c("none", "log2", "log2FC") # transformations to run: "none" | "log2" | "log2FC"
INCLUDE_MSFINDER  <- TRUE              # TRUE: also include MS-FINDER annotated Unknown features
OUT_DIR           <- "output/figures"
OUT_PROCESSED     <- "data/dummy/processed"

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
# 1. Package loading
# -----------------------------------------------------------------------------
pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork", "ggsignif")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# -----------------------------------------------------------------------------
# 2. Select target features (same filter as plot_timeseries.R)
#    (a) MS-DIAL identified : COUNT_MIN <= name count <= COUNT_MAX,
#                             not "Unknown", optionally not "low score:"
#    (b) MS-FINDER annotated: Unknown features with MSFINDER_annotated == TRUE,
#                             added when INCLUDE_MSFINDER is TRUE
# -----------------------------------------------------------------------------
name_counts  <- table(feat_meta$Metabolite_name)
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
target_ids <- as.character(target_features$Alignment_ID)

cat(sprintf("Target MS-DIAL names    : %d\n", length(target_names)))
cat(sprintf("Target features (total) : %d\n", length(target_ids)))

# -----------------------------------------------------------------------------
# 3. Build base averaged matrix (raw)
#    Average intensity over Repeats → one value per Subject × Timepoint
#    Biological samples only (reruns excluded)
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

sub_mat  <- feat_mat[target_ids, , drop = FALSE]
avg_mat0 <- sapply(group_cols, function(cols) {         # raw averaged matrix (preserved)
  rowMeans(sub_mat[, cols, drop = FALSE], na.rm = TRUE)
})

# -----------------------------------------------------------------------------
# 4. Pair definitions (fixed across all transforms)
# -----------------------------------------------------------------------------
early_pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
late_pairs  <- list(c(4L, 5L), c(4L, 6L), c(5L, 6L))
all_pairs   <- c(early_pairs, late_pairs)
pair_labels <- c("T1vT2", "T1vT3", "T2vT3", "T4vT5", "T4vT6", "T5vT6")

early_padj_cols <- paste0("padj_", c("T1vT2", "T1vT3", "T2vT3"))
late_padj_cols  <- paste0("padj_", c("T4vT5", "T4vT6", "T5vT6"))

# -----------------------------------------------------------------------------
# 5. Plot helper functions
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

make_metabolite_page <- function(idx, df_long, pval_df, classification, y_label) {
  aid      <- target_ids[idx]
  meta_row <- feat_meta[feat_meta$Alignment_ID == as.integer(aid), ]
  cls      <- classification[idx]

  # Display label: use MS-FINDER structure for annotated Unknowns
  is_msfinder <- isTRUE(meta_row$MSFINDER_annotated)
  if (is_msfinder) {
    display_name <- sprintf("%s [MSFINDER, score=%.2f]",
                            meta_row$MSFINDER_structure,
                            meta_row$MSFINDER_total_score)
  } else {
    display_name <- meta_row$Metabolite_name
  }

  title_str <- sprintf(
    "%s  |  ID=%s  Rt=%.2f  m/z=%.4f  [%s]",
    display_name, aid,
    round(as.numeric(meta_row$Rt_min), 2),
    round(as.numeric(meta_row$Mz), 4),
    cls
  )

  df_feat     <- df_long[df_long$Alignment_ID == aid, ]
  early_pvals <- unlist(pval_df[idx, c("T1vT2", "T1vT3", "T2vT3")])
  late_pvals  <- unlist(pval_df[idx, c("T4vT5", "T4vT6", "T5vT6")])

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

save_pdf_group <- function(indices, filename, group_label,
                           df_long, pval_df, classification, y_label) {
  if (length(indices) == 0) {
    cat(sprintf("  No features in '%s' — skipping.\n", group_label))
    return(invisible(NULL))
  }
  cat(sprintf("  Saving '%s' (%d features) → %s\n",
              group_label, length(indices), filename))
  pdf(filename, width = 10, height = 5)
  for (i in seq_along(indices)) {
    p <- make_metabolite_page(indices[i], df_long, pval_df, classification, y_label)
    print(p)
    if (i %% 50 == 0) cat(sprintf("    %d / %d\n", i, length(indices)))
  }
  dev.off()
  cat(sprintf("  Done: %s\n", filename))
}

# =============================================================================
# Main loop: run pipeline for each transform
# =============================================================================
dir.create(OUT_DIR,       recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PROCESSED, recursive = TRUE, showWarnings = FALSE)

all_summaries <- list()   # collect classification tables for final comparison

for (transform in TRANSFORMS) {

  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  Transform: %s\n", transform))
  cat(sprintf("%s\n", strrep("=", 70)))

  # ---------------------------------------------------------------------------
  # Step A: Apply transformation to avg_mat
  # ---------------------------------------------------------------------------
  avg_mat <- avg_mat0   # start from raw copy each iteration

  if (transform == "log2") {
    offsets <- apply(avg_mat, 1, function(x) {
      nz <- x[x > 0 & !is.na(x)]
      if (length(nz) == 0) return(1)
      min(nz) / 2
    })
    avg_mat <- log2(avg_mat + offsets)
    y_label <- "log2 Intensity (avg)"
    cat("  Transformation : log2(x + offset)  [offset = min non-zero / 2 per feature]\n")

  } else if (transform == "log2FC") {
    # Step 1: log2 transform with per-feature offset
    offsets <- apply(avg_mat, 1, function(x) {
      nz <- x[x > 0 & !is.na(x)]
      if (length(nz) == 0) return(1)
      min(nz) / 2
    })
    avg_mat <- log2(avg_mat + offsets)
    # Step 2: period-specific reference per subject
    #   Early (T1–T3) → subtract T1  (T1 = start of early intervention)
    #   Late  (T4–T6) → subtract T4  (T4 = start of late intervention, independent)
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
    cat("  Transformation : log2FC  [T1 ref for T1-T3; T4 ref for T4-T6]\n")
    cat("  Note: p-values are algebraically identical to 'log2' transform\n")

  } else {
    y_label <- "Intensity (avg)"
    cat("  Transformation : none (raw intensity)\n")
  }

  # ---------------------------------------------------------------------------
  # Step B: Convert to long format
  # ---------------------------------------------------------------------------
  df_avg              <- as.data.frame(avg_mat)
  df_avg$Alignment_ID <- rownames(df_avg)

  df_long <- df_avg %>%
    pivot_longer(-Alignment_ID, names_to = "SubjTP", values_to = "intensity") %>%
    mutate(
      Subject   = sub("-[0-9]+$", "", SubjTP),
      Timepoint = as.integer(sub(".*-", "", SubjTP))
    )

  # ---------------------------------------------------------------------------
  # Step B2: Save log2FC values as CSV (log2FC transform only)
  # ---------------------------------------------------------------------------
  if (transform == "log2FC") {
    logfc_out <- file.path(OUT_PROCESSED, "log2FC_values_dummy.csv")
    logfc_long <- df_long %>%
      left_join(feat_meta[, c("Alignment_ID", "Metabolite_name")],
                by = "Alignment_ID")
    write.csv(logfc_long, logfc_out, row.names = FALSE)
    cat(sprintf("  log2FC values saved: %s\n", logfc_out))
  }

  # ---------------------------------------------------------------------------
  # Step C: Pairwise paired Wilcoxon tests
  # ---------------------------------------------------------------------------
  cat("  Running pairwise Wilcoxon tests ...\n")

  pval_list <- lapply(target_ids, function(aid) {
    df_feat <- df_long[df_long$Alignment_ID == aid, ]
    pvals <- sapply(all_pairs, function(pair) {
      d1 <- df_feat[df_feat$Timepoint == pair[1], ]
      d2 <- df_feat[df_feat$Timepoint == pair[2], ]
      common_subj <- intersect(d1$Subject, d2$Subject)
      if (length(common_subj) < 3) return(NA_real_)
      d1 <- d1[match(common_subj, d1$Subject), ]
      d2 <- d2[match(common_subj, d2$Subject), ]
      tryCatch(
        wilcox.test(d1$intensity, d2$intensity, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    })
    setNames(as.data.frame(t(pvals)), pair_labels)
  })

  pval_df              <- bind_rows(pval_list)
  pval_df$Alignment_ID <- target_ids
  cat(sprintf("  %d tests computed (%d features x %d pairs)\n",
              length(target_ids) * 6, length(target_ids), 6))

  # ---------------------------------------------------------------------------
  # Step D: BH FDR correction
  # ---------------------------------------------------------------------------
  pval_mat <- as.matrix(pval_df[, pair_labels])
  padj_vec <- p.adjust(as.vector(pval_mat), method = "BH")
  padj_mat <- matrix(
    padj_vec,
    nrow     = nrow(pval_mat),
    ncol     = ncol(pval_mat),
    dimnames = list(NULL, paste0("padj_", pair_labels))
  )

  # ---------------------------------------------------------------------------
  # Step E: Classification
  # ---------------------------------------------------------------------------
  sig_early <- apply(padj_mat[, early_padj_cols], 1,
                     function(x) any(x < ALPHA, na.rm = TRUE))
  sig_late  <- apply(padj_mat[, late_padj_cols],  1,
                     function(x) any(x < ALPHA, na.rm = TRUE))

  classification <- ifelse(
    sig_early & sig_late, "both",
    ifelse(sig_early,     "early_only",
    ifelse(sig_late,      "late_only", "ns"))
  )

  tbl <- table(classification)
  cat("\n  Classification summary:\n")
  print(tbl)
  all_summaries[[transform]] <- tbl

  # ---------------------------------------------------------------------------
  # Step F: Save summary CSV
  # ---------------------------------------------------------------------------
  out_csv <- file.path(OUT_PROCESSED,
                       sprintf("pairwise_test_results_%s_dummy.csv", transform))

  idx_in_meta <- match(as.integer(target_ids), feat_meta$Alignment_ID)
  summary_df <- data.frame(
    Alignment_ID       = target_ids,
    Metabolite_name    = feat_meta$Metabolite_name[idx_in_meta],
    MSFINDER_structure = if ("MSFINDER_structure" %in% colnames(feat_meta))
                           feat_meta$MSFINDER_structure[idx_in_meta] else NA,
    MSFINDER_score     = if ("MSFINDER_total_score" %in% colnames(feat_meta))
                           feat_meta$MSFINDER_total_score[idx_in_meta] else NA,
    Rt_min             = feat_meta$Rt_min[idx_in_meta],
    Mz                 = feat_meta$Mz[idx_in_meta],
    classification     = classification,
    stringsAsFactors   = FALSE
  )
  summary_df <- cbind(summary_df,
                      pval_df[, pair_labels],
                      as.data.frame(padj_mat))

  write.csv(summary_df, out_csv, row.names = FALSE)
  cat(sprintf("\n  Results saved: %s\n", out_csv))

  # ---------------------------------------------------------------------------
  # Step G: Save PDFs
  # ---------------------------------------------------------------------------
  idx_both       <- which(classification == "both")
  idx_early_only <- which(classification == "early_only")
  idx_late_only  <- which(classification == "late_only")

  save_pdf_group(
    idx_both,
    file.path(OUT_DIR, sprintf("boxplot_sig_both_%s_dummy.pdf", transform)),
    "both", df_long, pval_df, classification, y_label
  )
  save_pdf_group(
    idx_early_only,
    file.path(OUT_DIR, sprintf("boxplot_sig_early_only_%s_dummy.pdf", transform)),
    "early_only", df_long, pval_df, classification, y_label
  )
  save_pdf_group(
    idx_late_only,
    file.path(OUT_DIR, sprintf("boxplot_sig_late_only_%s_dummy.pdf", transform)),
    "late_only", df_long, pval_df, classification, y_label
  )
}

# =============================================================================
# Final summary across all transforms
# =============================================================================
cat(sprintf("\n%s\n", strrep("=", 70)))
cat("  Final classification summary across transforms\n")
cat(sprintf("%s\n", strrep("=", 70)))
for (transform in TRANSFORMS) {
  cat(sprintf("\n  [%s]\n", transform))
  print(all_summaries[[transform]])
}
cat(sprintf("\n  Outputs:\n"))
cat(sprintf("    CSVs : %s/pairwise_test_results_*_dummy.csv\n", OUT_PROCESSED))
cat(sprintf("    PDFs : %s/boxplot_sig_*_dummy.pdf\n", OUT_DIR))
cat(sprintf("%s\n", strrep("=", 70)))
