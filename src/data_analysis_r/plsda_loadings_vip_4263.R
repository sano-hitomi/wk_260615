# =============================================================================
# plsda_loadings_vip_4263.R
# PLS-DA with Alignment ID 4263 response groups as Y — VIP and loading extraction
#
# Identical logic to plsda_loadings_vip.R, but uses the 4-group classification
# derived from Alignment ID 4263 (subject_response_groups_4263.csv).
#
# Output:
#   data/production/processed/plsda_vip_loadings_4263.csv
#   output/figures/plsda_vip_early_top_4263.pdf / .png
#   output/figures/plsda_vip_late_top_4263.pdf  / .png
#   output/figures/plsda_vip_scatter_4263.pdf   / .png
#   output/figures/plsda_loading_scatter_4263.pdf / .png
#
# Prerequisites:
#   1. load_data_production.R
#   2. classify_subjects_by_response_4263.R  (generates GROUPS_CSV)
# =============================================================================

TARGET_ID  <- 4263
GROUPS_CSV <- sprintf("data/production/processed/subject_response_groups_%d.csv",
                      TARGET_ID)
NCOMP      <- 2
TOP_N      <- 30
OUT_DIR    <- "output/figures"
OUT_CSV    <- sprintf("data/production/processed/plsda_vip_loadings_%d.csv", TARGET_ID)

# -----------------------------------------------------------------------------
# 0. Prerequisites
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0)
  stop("Missing: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
if (!file.exists(GROUPS_CSV))
  stop("Groups CSV not found: ", GROUPS_CSV,
       "\nRun classify_subjects_by_response_4263.R first.")

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
pkgs <- c("ggplot2", "dplyr", "mixOmics")
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
# 3. Build per-subject-timepoint log2 matrix
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
avg_mat_full <- do.call(cbind, avg_list)

log2_mat_full <- t(apply(avg_mat_full, 1, function(x) {
  nz     <- x[x > 0 & !is.na(x)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2(x + offset)
}))
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
# 4. Build X / Y per period
# -----------------------------------------------------------------------------
make_data <- function(tp_range) {
  idx  <- which(as.integer(row_meta$Timepoint) %in% tp_range)
  mat  <- log2_mat_full[idx, , drop = FALSE]
  meta <- row_meta[idx, ]

  meta <- merge(meta, groups_df[, c("Subject", "group")],
                by = "Subject", all.x = TRUE)
  meta$group[is.na(meta$group)] <- "other"
  meta <- meta[match(row_meta$SubjTP[idx], meta$SubjTP), ]

  keep <- apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0)
  mat  <- mat[, keep, drop = FALSE]

  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) mat[na_idx, j] <- mean(mat[, j], na.rm = TRUE)
  }

  Y <- factor(meta$group, levels = c("both", "X_only", "Y_only", "other"))
  list(X = mat, Y = Y, meta = meta)
}

early_data <- make_data(1:3)
late_data  <- make_data(4:6)

# -----------------------------------------------------------------------------
# 5. PLS-DA
# -----------------------------------------------------------------------------
ncomp_use <- min(NCOMP, nlevels(early_data$Y) - 1, nrow(early_data$X) - 1)
ncomp_use <- max(ncomp_use, 2)

early_plsda <- plsda(X = early_data$X, Y = early_data$Y,
                     ncomp = ncomp_use, scale = TRUE)
late_plsda  <- plsda(X = late_data$X,  Y = late_data$Y,
                     ncomp = ncomp_use, scale = TRUE)

cat("PLS-DA complete.\n")

# -----------------------------------------------------------------------------
# 6. Extract VIP and loadings (Comp1)
# -----------------------------------------------------------------------------
early_vip  <- vip(early_plsda)[, 1]
late_vip   <- vip(late_plsda)[, 1]

early_load <- early_plsda$loadings$X[, 1]
late_load  <- late_plsda$loadings$X[, 1]

common_ids <- intersect(names(early_vip), names(late_vip))

# -----------------------------------------------------------------------------
# 7. Sign alignment (anchor on TARGET_ID if present)
# -----------------------------------------------------------------------------
ref_id <- as.character(TARGET_ID)
if (ref_id %in% common_ids) {
  sign_flip_late <- sign(early_load[ref_id]) != sign(late_load[ref_id])
  cat(sprintf("Sign alignment: using ID %s as anchor (early=%.4f, late=%.4f).\n",
              ref_id, early_load[ref_id], late_load[ref_id]))
} else {
  load_cor <- cor(early_load[common_ids], late_load[common_ids], use = "complete.obs")
  sign_flip_late <- load_cor < 0
  cat(sprintf("Sign alignment: TARGET_ID absent; loading correlation r=%.3f.\n", load_cor))
}

if (sign_flip_late) {
  cat("  → Flipping sign of Late Comp1.\n")
  late_load                        <- -late_load
  late_plsda$loadings$X[, 1]      <- -late_plsda$loadings$X[, 1]
  late_plsda$variates$X[, 1]      <- -late_plsda$variates$X[, 1]
} else {
  cat("  → Signs consistent; no flip.\n")
}

# -----------------------------------------------------------------------------
# 8. Name lookup
# -----------------------------------------------------------------------------
id_to_name <- function(ids) {
  idx  <- match(as.integer(ids), feat_meta$Alignment_ID)
  name <- feat_meta$Metabolite_name[idx]
  if ("MSFINDER_annotated" %in% colnames(feat_meta)) {
    is_msf <- !is.na(feat_meta$MSFINDER_annotated[idx]) &
               feat_meta$MSFINDER_annotated[idx]
    name[is_msf] <- feat_meta$MSFINDER_structure[idx[is_msf]]
  }
  name
}

# -----------------------------------------------------------------------------
# 9. Build summary CSV
# -----------------------------------------------------------------------------
summary_df <- data.frame(
  Alignment_ID    = common_ids,
  Metabolite_name = id_to_name(common_ids),
  VIP_early       = round(early_vip[common_ids],  4),
  VIP_late        = round(late_vip[common_ids],   4),
  loading_early   = round(early_load[common_ids], 6),
  loading_late    = round(late_load[common_ids],  6),
  stringsAsFactors = FALSE
)
summary_df <- summary_df[order(-pmax(summary_df$VIP_early, summary_df$VIP_late)), ]

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, OUT_CSV, row.names = FALSE)
cat(sprintf("VIP/loading table saved: %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 10. Helper
# -----------------------------------------------------------------------------
trunc_name <- function(x, n = 35) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n), "…"), x)
}

# -----------------------------------------------------------------------------
# 11. VIP bar plots
# -----------------------------------------------------------------------------
make_vip_bar <- function(df, vip_col, period_label, top_n) {
  df_top    <- df[order(-df[[vip_col]]), ][1:min(top_n, nrow(df)), ]
  raw_name  <- ifelse(is.na(df_top$Metabolite_name) | df_top$Metabolite_name == "",
                      "", df_top$Metabolite_name)
  name_trunc <- trunc_name(raw_name, n = 28)
  full_label  <- ifelse(name_trunc == "",
                        as.character(df_top$Alignment_ID),
                        paste0(df_top$Alignment_ID, " ", name_trunc))
  df_top$label <- trunc_name(full_label, n = 38)
  df_top$label <- make.unique(df_top$label, sep = " #")
  df_top$label <- factor(df_top$label, levels = rev(df_top$label))

  ggplot(df_top, aes(x = .data[[vip_col]], y = label)) +
    geom_col(fill = "#2C7BB6", alpha = 0.85, width = 0.7) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "firebrick", linewidth = 0.5) +
    labs(
      title    = sprintf("Top %d features by VIP (Comp1) — %s [Y = ID %d]",
                         top_n, period_label, TARGET_ID),
      subtitle = "Dashed line: VIP = 1",
      x        = "VIP score (Comp1)",
      y        = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 10, face = "bold"),
      plot.subtitle      = element_text(size = 8, color = "grey40"),
      axis.text.y        = element_text(size = 8)
    )
}

p_vip_early <- make_vip_bar(summary_df, "VIP_early", "Early (T1–T3)", TOP_N)
p_vip_late  <- make_vip_bar(summary_df, "VIP_late",  "Late  (T4–T6)", TOP_N)

# -----------------------------------------------------------------------------
# 12. VIP scatter
# -----------------------------------------------------------------------------
vip_thresh <- 1.0

p_vip_scatter <- ggplot(summary_df, aes(x = VIP_early, y = VIP_late)) +
  geom_hline(yintercept = vip_thresh, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = vip_thresh, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.2, alpha = 0.4, color = "grey70") +
  geom_point(
    data = subset(summary_df, VIP_early > vip_thresh | VIP_late > vip_thresh),
    aes(color = case_when(
      VIP_early > vip_thresh & VIP_late > vip_thresh ~ "Both",
      VIP_early > vip_thresh                         ~ "Early only",
      TRUE                                           ~ "Late only"
    )),
    size = 2.0, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("Both" = "#E41A1C", "Early only" = "#377EB8", "Late only" = "#4DAF4A"),
    name   = "VIP > 1"
  ) +
  labs(
    title    = sprintf("VIP comparison: Early vs Late [Y = ID %d]", TARGET_ID),
    subtitle = "Features above dashed lines have VIP > 1 in that period",
    x        = "VIP — Early (T1–T3)",
    y        = "VIP — Late  (T4–T6)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

# -----------------------------------------------------------------------------
# 13. Loading scatter
# -----------------------------------------------------------------------------
load_df <- summary_df
load_df$vip_flag <- case_when(
  load_df$VIP_early > vip_thresh & load_df$VIP_late > vip_thresh ~ "Both",
  load_df$VIP_early > vip_thresh                                  ~ "Early only",
  load_df$VIP_late  > vip_thresh                                  ~ "Late only",
  TRUE                                                            ~ "Low VIP"
)

p_load_scatter <- ggplot(load_df,
                         aes(x = loading_early, y = loading_late,
                             color = vip_flag)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(
    values = c("Both"       = "#E41A1C",
               "Early only" = "#377EB8",
               "Late only"  = "#4DAF4A",
               "Low VIP"    = "grey70"),
    name = "VIP > 1"
  ) +
  labs(
    title    = sprintf("Comp1 loadings: Early vs Late  [Y = ID %d]", TARGET_ID),
    subtitle = "Q1(++)/Q3(--): same direction | Q2/Q4: reversal between periods",
    x        = "Loading Comp1 — Early (T1–T3)",
    y        = "Loading Comp1 — Late  (T4–T6)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

# -----------------------------------------------------------------------------
# 14. Save plots
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, base_name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(base_name, ".pdf")),
         plot = p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(base_name, ".png")),
         plot = p, width = w, height = h, dpi = 150)
  cat(sprintf("  %s\n", base_name))
}

cat("\nSaving plots:\n")
save_plot(p_vip_early,    sprintf("plsda_vip_early_top_%d",    TARGET_ID), 8,   6)
save_plot(p_vip_late,     sprintf("plsda_vip_late_top_%d",     TARGET_ID), 8,   6)
save_plot(p_vip_scatter,  sprintf("plsda_vip_scatter_%d",      TARGET_ID), 6.5, 5.5)
save_plot(p_load_scatter, sprintf("plsda_loading_scatter_%d",  TARGET_ID), 6.5, 5.5)

# -----------------------------------------------------------------------------
# 15. Summary
# -----------------------------------------------------------------------------
cat(sprintf("\nVIP > 1 summary (Y = ID %d):\n", TARGET_ID))
cat(sprintf("  Early only : %d\n", sum(summary_df$VIP_early > 1 & summary_df$VIP_late <= 1)))
cat(sprintf("  Late  only : %d\n", sum(summary_df$VIP_early <= 1 & summary_df$VIP_late > 1)))
cat(sprintf("  Both       : %d\n", sum(summary_df$VIP_early > 1 & summary_df$VIP_late > 1)))
cat(sprintf("  Neither    : %d\n", sum(summary_df$VIP_early <= 1 & summary_df$VIP_late <= 1)))
