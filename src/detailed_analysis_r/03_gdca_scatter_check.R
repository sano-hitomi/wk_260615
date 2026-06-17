# =============================================================================
# 03_gdca_scatter_check.R
#
# ① GDCA × TMAO 期別 log2FC 散布図（X期・Y期）
#    — Y期に着目: GDCA高・TMAO高 の群は誰か
#
# ② ID 4263（GDCA）の時点別ボックスプロット
#    — Early T2 のピークが信頼できる実データかを被験者レベルで確認
#
# 出力:
#   output/figures/gdca_tmao_scatter_XY.pdf/.png    — ①
#   output/figures/gdca_per_timepoint_box.pdf/.png  — ②
#
# 前提:
#   load_data_production.R 実行済み
#   classify_subjects_by_response.R 実行済み（GROUPS_CSV）
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
TMAO_ID    <- 414L
GDCA_ID    <- 4263L
GROUPS_CSV <- "data/production/processed/subject_response_groups_414.csv"
OUT_DIR    <- "output/figures"

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
pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# -----------------------------------------------------------------------------
# 2. 共通: テクニカルレプリカート平均 → log2 変換
# -----------------------------------------------------------------------------
groups_df <- read.csv(GROUPS_CSV, stringsAsFactors = FALSE)
groups_df$Subject <- as.character(groups_df$Subject)

bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]
group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

ids_chr <- as.character(c(TMAO_ID, GDCA_ID))
avg_raw <- sapply(group_cols, function(cols) {
  rowMeans(feat_mat[ids_chr, cols, drop = FALSE], na.rm = TRUE)
})

log2_mat <- t(apply(avg_raw, 1, function(x) {
  nz <- x[x > 0 & !is.na(x)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2(x + offset)
}))

# log2FC（期内ベースライン差分）
subjects <- sort(unique(bio$Subject))
lfc_mat  <- log2_mat
for (subj in subjects) {
  for (ref_tp in c(1L, 4L)) {
    ref_col   <- paste0(subj, "-", ref_tp)
    if (!ref_col %in% colnames(lfc_mat)) next
    ref_vals  <- lfc_mat[, ref_col]
    tp_range  <- if (ref_tp == 1L) 1:3 else 4:6
    subj_cols <- intersect(paste0(subj, "-", tp_range), colnames(lfc_mat))
    lfc_mat[, subj_cols] <- sweep(lfc_mat[, subj_cols, drop = FALSE], 1, ref_vals, "-")
  }
}

# Long format
df_long <- as.data.frame(lfc_mat) %>%
  mutate(Alignment_ID = rownames(lfc_mat)) %>%
  pivot_longer(-Alignment_ID, names_to = "SubjTP", values_to = "lfc") %>%
  mutate(
    Subject   = sub("-[0-9]+$", "", SubjTP),
    Timepoint = as.integer(sub(".*-", "", SubjTP))
  ) %>%
  left_join(groups_df[, c("Subject", "group")], by = "Subject") %>%
  mutate(group = ifelse(is.na(group), "other", group),
         group = factor(group, levels = names(GROUP_COLORS)))

# ─────────────────────────────────────────────────────────────────────────────
# ① 散布図: TMAO log2FC × GDCA log2FC（X期・Y期）
# ─────────────────────────────────────────────────────────────────────────────

# 被験者単位の期別平均 log2FC
df_wide <- df_long %>%
  group_by(Alignment_ID, Subject, group) %>%
  summarise(
    lfc_X = mean(lfc[Timepoint %in% 2:3], na.rm = TRUE),
    lfc_Y = mean(lfc[Timepoint %in% 5:6], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols     = c(Subject, group),
    names_from  = Alignment_ID,
    values_from = c(lfc_X, lfc_Y),
    names_sep   = "_"
  )

tmao_col_X <- paste0("lfc_X_", TMAO_ID)
tmao_col_Y <- paste0("lfc_Y_", TMAO_ID)
gdca_col_X <- paste0("lfc_X_", GDCA_ID)
gdca_col_Y <- paste0("lfc_Y_", GDCA_ID)

make_scatter <- function(df, xcol, ycol, period_label, corr_label_pos = "tl") {
  r  <- cor(df[[xcol]], df[[ycol]], use = "complete.obs")
  pv <- cor.test(df[[xcol]], df[[ycol]])$p.value
  ann <- sprintf("r = %.3f\np = %.3f", r, pv)
  xpos <- if (corr_label_pos == "tl") -Inf else Inf
  hjust <- if (corr_label_pos == "tl") -0.1 else 1.1

  ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]],
                 color = group, label = Subject)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_smooth(method = "lm", se = TRUE, inherit.aes = FALSE,
                aes(x = .data[[xcol]], y = .data[[ycol]]),
                color = "grey50", linewidth = 0.7, linetype = "dashed") +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = 20) +
    annotate("text", x = xpos, y = Inf,
             label = ann, hjust = hjust, vjust = 1.4,
             size = 3.5, color = "grey30") +
    scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                       name = "TMAO response") +
    labs(
      title = period_label,
      x = sprintf("TMAO (ID %d) log2FC", TMAO_ID),
      y = sprintf("GDCA (ID %d) log2FC", GDCA_ID)
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 11, face = "bold"))
}

p_s_X <- make_scatter(df_wide, tmao_col_X, gdca_col_X, "Early (X期, T2–T3 vs T1)")
p_s_Y <- make_scatter(df_wide, tmao_col_Y, gdca_col_Y, "Late  (Y期, T5–T6 vs T4)")

p_scatter <- (p_s_X | p_s_Y) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "① TMAO × GDCA 期別 log2FC 散布図",
    subtitle = "各点 = 被験者（T2-T3 or T5-T6 の平均 log2FC）",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40")
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# ② 4263（GDCA）の時点別ボックスプロット — T2 ピークの信頼性確認
# ─────────────────────────────────────────────────────────────────────────────

df_gdca_lfc <- df_long %>%
  filter(Alignment_ID == as.character(GDCA_ID)) %>%
  mutate(Timepoint = factor(Timepoint, levels = 1:6))

# 外れ値候補（|lfc| > 1.5 × IQR + Q3 を超える点）にラベルを付ける
outlier_thresh <- df_gdca_lfc %>%
  group_by(Timepoint) %>%
  summarise(
    q1  = quantile(lfc, 0.25, na.rm = TRUE),
    q3  = quantile(lfc, 0.75, na.rm = TRUE),
    iqr = IQR(lfc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(upper = q3 + 1.5 * iqr, lower = q1 - 1.5 * iqr)

df_gdca_label <- df_gdca_lfc %>%
  left_join(outlier_thresh[, c("Timepoint", "upper", "lower")],
            by = "Timepoint") %>%
  filter(lfc > upper | lfc < lower)

p_box <- ggplot(df_gdca_lfc, aes(x = Timepoint, y = lfc, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60",
             linetype = "dashed") +
  geom_vline(xintercept = 3.5, linewidth = 0.5, color = "grey60",
             linetype = "dashed") +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.55,
               position = position_dodge(width = 0)) +
  # 全個別点（jitter）
  geom_jitter(aes(color = group), width = 0.2, size = 1.5, alpha = 0.6) +
  # 外れ値候補にラベル
  geom_text_repel(data = df_gdca_label,
                  aes(label = Subject, color = group),
                  size = 2.8, show.legend = FALSE,
                  max.overlaps = 20, segment.size = 0.3) +
  annotate("text", x = 2,   y = -Inf, label = "Early (X期)", vjust = -0.5,
           size = 3.5, color = "grey40") +
  annotate("text", x = 5,   y = -Inf, label = "Late (Y期)",  vjust = -0.5,
           size = 3.5, color = "grey40") +
  scale_fill_manual(values  = GROUP_COLORS, labels = GROUP_LABELS,
                    name = "TMAO response") +
  scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                     name = "TMAO response") +
  labs(
    title    = sprintf("② GDCA（ID %d）時点別 log2FC — T2 ピークの信頼性確認", GDCA_ID),
    subtitle = "外れ値候補（1.5×IQR 基準）に被験者番号を表示",
    x        = "Timepoint",
    y        = "log2FC"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        plot.title       = element_text(size = 11, face = "bold"),
        plot.subtitle    = element_text(size = 9,  color = "grey40"))

# ─────────────────────────────────────────────────────────────────────────────
# 保存
# ─────────────────────────────────────────────────────────────────────────────
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), plot = p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), plot = p, width = w, height = h, dpi = 150)
  cat(sprintf("  %s\n", name))
}

cat("保存完了:\n")
save_plot(p_scatter, "gdca_tmao_scatter_XY",    12, 5.5)
save_plot(p_box,     "gdca_per_timepoint_box",  9,  5.5)
