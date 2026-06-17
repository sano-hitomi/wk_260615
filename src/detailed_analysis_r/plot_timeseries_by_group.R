# =============================================================================
# plot_timeseries_by_group.R
# 指定した Alignment_ID の時系列を、TMAOレスポンスグループ色で描画
#
# 前提:
#   load_data_production.R 実行済み
#   classify_subjects_by_response.R 実行済み（GROUPS_CSV が存在する）
#
# 出力:
#   output/figures/timeseries_by_group.pdf
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
TARGET_IDS <- c(414L, 415L, 4263L)
GROUPS_CSV <- "data/production/processed/subject_response_groups_414.csv"
TRANSFORM  <- "log2FC"   # "none" | "log2" | "log2FC"
OUT_PDF    <- "output/figures/timeseries_by_group.pdf"

GROUP_COLORS <- c(
  both   = "#E41A1C",
  X_only = "#377EB8",
  Y_only = "#4DAF4A",
  other  = "#984EA3"
)
GROUP_LABELS <- c(
  both   = "Both (X & Y)",
  X_only = "X only",
  Y_only = "Y only",
  other  = "Other"
)

# -----------------------------------------------------------------------------
# 0. 前提確認
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0)
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
if (!file.exists(GROUPS_CSV))
  stop("Groups CSV not found: ", GROUPS_CSV,
       "\nRun classify_subjects_by_response.R first.")

# -----------------------------------------------------------------------------
# 1. パッケージ
# -----------------------------------------------------------------------------
pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# -----------------------------------------------------------------------------
# 2. グループ読み込み
# -----------------------------------------------------------------------------
groups_df <- read.csv(GROUPS_CSV, stringsAsFactors = FALSE)
groups_df$Subject <- as.character(groups_df$Subject)

# -----------------------------------------------------------------------------
# 3. テクニカルレプリカートの平均 → 被験者×時点行列
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

ids_chr <- as.character(TARGET_IDS)
avg_raw <- sapply(group_cols, function(cols) {
  rowMeans(feat_mat[ids_chr, cols, drop = FALSE], na.rm = TRUE)
})
if (length(ids_chr) == 1)
  avg_raw <- matrix(avg_raw, nrow = 1, dimnames = list(ids_chr, names(group_cols)))

# -----------------------------------------------------------------------------
# 4. 変換
# -----------------------------------------------------------------------------
subjects   <- sort(unique(bio$Subject))
timepoints <- sort(unique(bio$Timepoint))

if (TRANSFORM %in% c("log2", "log2FC")) {
  offsets <- apply(avg_raw, 1, function(x) {
    nz <- x[x > 0 & !is.na(x)]; if (length(nz) == 0) 1 else min(nz) / 2
  })
  avg_mat <- log2(avg_raw + offsets)
} else {
  avg_mat <- avg_raw
}

if (TRANSFORM == "log2FC") {
  for (subj in subjects) {
    for (ref_tp in c(1L, 4L)) {
      ref_col  <- paste0(subj, "-", ref_tp)
      if (!ref_col %in% colnames(avg_mat)) next
      ref_vals <- avg_mat[, ref_col]
      tp_range <- if (ref_tp == 1L) 1:3 else 4:6
      subj_cols <- intersect(paste0(subj, "-", tp_range), colnames(avg_mat))
      avg_mat[, subj_cols] <- sweep(avg_mat[, subj_cols, drop = FALSE], 1, ref_vals, "-")
    }
  }
  y_label <- "log2FC (T1–T3: vs T1 / T4–T6: vs T4)"
} else if (TRANSFORM == "log2") {
  y_label <- "log2 Intensity"
} else {
  y_label <- "Intensity"
}

# -----------------------------------------------------------------------------
# 5. Long format に変換してグループ情報を付与
# -----------------------------------------------------------------------------
df_long <- as.data.frame(avg_mat) %>%
  mutate(Alignment_ID = rownames(avg_mat)) %>%
  pivot_longer(-Alignment_ID, names_to = "SubjTP", values_to = "value") %>%
  mutate(
    Subject   = sub("-[0-9]+$", "", SubjTP),
    Timepoint = as.integer(sub(".*-", "", SubjTP))
  ) %>%
  left_join(groups_df[, c("Subject", "group")], by = "Subject") %>%
  mutate(
    group = ifelse(is.na(group), "other", group),
    group = factor(group, levels = names(GROUP_COLORS))
  )

# -----------------------------------------------------------------------------
# 6. 1 特徴量あたり 1 パネル描画
# -----------------------------------------------------------------------------
make_panel <- function(aid) {
  meta_row <- feat_meta[feat_meta$Alignment_ID == as.integer(aid), ]
  name_str <- if (nrow(meta_row) > 0) meta_row$Metabolite_name[1] else "Unknown"
  rt_str   <- if (nrow(meta_row) > 0) sprintf("Rt=%.2f", meta_row$Rt_min[1]) else ""
  mz_str   <- if (nrow(meta_row) > 0) sprintf("m/z=%.4f", meta_row$Mz[1]) else ""

  df_feat <- df_long[df_long$Alignment_ID == aid, ]

  # グループ別 grand mean
  df_gm <- df_feat %>%
    group_by(group, Timepoint) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  ggplot(df_feat, aes(x = Timepoint, y = value,
                      color = group, group = Subject)) +
    geom_vline(xintercept = 3.5, linetype = "dashed",
               color = "grey70", linewidth = 0.5) +
    # 個別被験者（細線・半透明）
    geom_line(alpha = 0.35, linewidth = 0.6) +
    geom_point(alpha = 0.35, size = 1.2) +
    # グループ mean（太線）
    geom_line(data = df_gm,
              aes(x = Timepoint, y = mean_val, color = group, group = group),
              linewidth = 1.8, inherit.aes = FALSE) +
    geom_point(data = df_gm,
               aes(x = Timepoint, y = mean_val, color = group, group = group),
               size = 3.5, inherit.aes = FALSE) +
    scale_x_continuous(breaks = 1:6) +
    scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                       name = "TMAO response") +
    labs(
      title = sprintf("ID %s  %s  |  %s  %s", aid, name_str, rt_str, mz_str),
      x     = "Timepoint",
      y     = y_label
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 9, face = "bold"),
      legend.position  = "right"
    )
}

panels <- lapply(ids_chr, make_panel)

p_combined <- wrap_plots(panels, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = sprintf("時系列: Alignment ID %s（%s）",
                       paste(TARGET_IDS, collapse = " / "), TRANSFORM),
    subtitle = "細線: 個別被験者  |  太線: グループ mean  |  破線: Early/Late 境界",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40")
    )
  )

# -----------------------------------------------------------------------------
# 7. 保存
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PDF,
       plot   = p_combined,
       width  = 9,
       height = 4.5 * length(TARGET_IDS))
ggsave(sub("\\.pdf$", ".png", OUT_PDF),
       plot   = p_combined,
       width  = 9,
       height = 4.5 * length(TARGET_IDS),
       dpi    = 150)

cat(sprintf("保存完了: %s\n", OUT_PDF))
