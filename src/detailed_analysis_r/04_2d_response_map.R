# =============================================================================
# 04_2d_response_map.R
#
# ③ TMAO_Y × GDCA_Y の2次元応答マップ
#
# 目的:
#   TMAOグループ分類（4群, 414 の log2FC に基づく）は1次元の情報しか使っていない。
#   GDCA（4263）の Y期応答を第2軸に加えた2D空間で被験者を分類し、
#   新たな亜群（特に「GDCA↑ TMAO↑」vs「GDCA↑ TMAO↓」）を探索する。
#
# 象限の解釈:
#   Q1 (TMAO↑, GDCA↑): 両代謝物が上昇 → 腸内細菌叢が全般的に活性化
#   Q2 (TMAO↓, GDCA↑): GDCA産生菌は活性だがTMA産生は抑制
#   Q3 (TMAO↓, GDCA↓): 腸内細菌叢全般的に非活性
#   Q4 (TMAO↑, GDCA↓): TMA産生は活性だがGDCA産生菌は少ない
#
# 出力:
#   output/figures/2d_response_map.pdf/.png        — 2D散布図（象限分類）
#   output/figures/2d_group_comparison.pdf/.png    — 元4群 vs 新2D象限の比較
#   data/production/processed/subjects_2d_class.csv  — 被験者ごとの象限分類
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
OUT_CSV    <- "data/production/processed/subjects_2d_class.csv"

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

QUAD_COLORS <- c(
  Q1 = "#C00000",   # TMAO↑ GDCA↑
  Q2 = "#70AD47",   # TMAO↓ GDCA↑
  Q3 = "#4472C4",   # TMAO↓ GDCA↓
  Q4 = "#ED7D31"    # TMAO↑ GDCA↓
)
QUAD_LABELS <- c(
  Q1 = "Q1: TMAO↑ GDCA↑",
  Q2 = "Q2: TMAO↓ GDCA↑",
  Q3 = "Q3: TMAO↓ GDCA↓",
  Q4 = "Q4: TMAO↑ GDCA↓"
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
# 2. データ準備（テクニカルレプリカート平均 → log2FC）
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

# 被験者単位の期別平均 log2FC
df_lfc <- as.data.frame(lfc_mat) %>%
  mutate(Alignment_ID = rownames(lfc_mat)) %>%
  pivot_longer(-Alignment_ID, names_to = "SubjTP", values_to = "lfc") %>%
  mutate(
    Subject   = sub("-[0-9]+$", "", SubjTP),
    Timepoint = as.integer(sub(".*-", "", SubjTP))
  ) %>%
  group_by(Alignment_ID, Subject) %>%
  summarise(
    lfc_X = mean(lfc[Timepoint %in% 2:3], na.rm = TRUE),
    lfc_Y = mean(lfc[Timepoint %in% 5:6], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(id_cols = Subject,
              names_from = Alignment_ID,
              values_from = c(lfc_X, lfc_Y),
              names_sep = "_") %>%
  left_join(groups_df[, c("Subject", "group")], by = "Subject") %>%
  mutate(group = ifelse(is.na(group), "other", group),
         group = factor(group, levels = names(GROUP_COLORS)))

tmao_X <- paste0("lfc_X_", TMAO_ID)
tmao_Y <- paste0("lfc_Y_", TMAO_ID)
gdca_X <- paste0("lfc_X_", GDCA_ID)
gdca_Y <- paste0("lfc_Y_", GDCA_ID)

# -----------------------------------------------------------------------------
# 3. Y期の2次元象限分類
#    閾値: 0（log2FC = 0 が自然なベースライン境界）
# -----------------------------------------------------------------------------
df_lfc <- df_lfc %>%
  mutate(
    quadrant = case_when(
      .data[[tmao_Y]] >= 0 & .data[[gdca_Y]] >= 0 ~ "Q1",
      .data[[tmao_Y]] <  0 & .data[[gdca_Y]] >= 0 ~ "Q2",
      .data[[tmao_Y]] <  0 & .data[[gdca_Y]] <  0 ~ "Q3",
      .data[[tmao_Y]] >= 0 & .data[[gdca_Y]] <  0 ~ "Q4",
      TRUE ~ NA_character_
    ),
    quadrant = factor(quadrant, levels = c("Q1","Q2","Q3","Q4"))
  )

cat("\n=== Y期 2D象限 分類結果 ===\n")
print(table(df_lfc$quadrant, useNA = "ifany"))
cat("\n=== 元4群 × 2D象限 クロス集計 ===\n")
print(table(group = df_lfc$group, quadrant = df_lfc$quadrant, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 4. CSV 保存
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
out_df <- df_lfc %>%
  select(Subject, group, quadrant,
         TMAO_X = all_of(tmao_X), TMAO_Y = all_of(tmao_Y),
         GDCA_X = all_of(gdca_X), GDCA_Y = all_of(gdca_Y))
write.csv(out_df, OUT_CSV, row.names = FALSE)
cat(sprintf("\n2D分類 CSV 保存: %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 5. メインプロット: 2D 応答マップ（Y期）
# -----------------------------------------------------------------------------
# 象限ラベルの配置位置（各象限の端）
quad_ann <- data.frame(
  x     = c( 2.2, -2.2, -2.2,  2.2),
  y     = c( max(df_lfc[[gdca_Y]], na.rm=TRUE) * 0.92,
             max(df_lfc[[gdca_Y]], na.rm=TRUE) * 0.92,
             min(df_lfc[[gdca_Y]], na.rm=TRUE) * 0.92,
             min(df_lfc[[gdca_Y]], na.rm=TRUE) * 0.92),
  label = c("Q1\nTMAO↑ GDCA↑", "Q2\nTMAO↓ GDCA↑",
            "Q3\nTMAO↓ GDCA↓", "Q4\nTMAO↑ GDCA↓"),
  color = c("#C00000","#70AD47","#4472C4","#ED7D31"),
  stringsAsFactors = FALSE
)

p_2d <- ggplot(df_lfc, aes(x = .data[[tmao_Y]], y = .data[[gdca_Y]])) +
  # 象限背景
  annotate("rect", xmin =  0, xmax = Inf, ymin =  0, ymax = Inf,
           fill = "#C00000", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin =  0, ymax = Inf,
           fill = "#70AD47", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#4472C4", alpha = 0.05) +
  annotate("rect", xmin =  0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#ED7D31", alpha = 0.05) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey50") +
  # 象限ラベル
  geom_text(data = quad_ann, aes(x = x, y = y, label = label, color = color),
            size = 3.2, fontface = "bold", inherit.aes = FALSE) +
  scale_color_identity() +
  # 被験者点（元4群の色）
  ggnewscale::new_scale_color() +
  geom_point(aes(color = group), size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = Subject, color = group),
                  size = 2.6, show.legend = FALSE, max.overlaps = 25) +
  scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                     name = "元TMAOグループ") +
  labs(
    title    = "③ 2次元応答マップ（Y期, log2FC）",
    subtitle = sprintf("x軸: TMAO (ID %d) Y期 log2FC  |  y軸: GDCA (ID %d) Y期 log2FC",
                       TMAO_ID, GDCA_ID),
    x = "TMAO log2FC (Y期, T5-T6 vs T4)",
    y = "GDCA log2FC (Y期, T5-T6 vs T4)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 9,  color = "grey40"),
        legend.position  = "right")

# ggnewscale が使えない場合のフォールバック
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
  library(ggnewscale)
  # プロットを再生成（install後に再実行）
  source(sys.frame(1)$ofile)
}

# -----------------------------------------------------------------------------
# 6. サブプロット: 元4群のY期応答分布（ヒートマップ風集計）
# -----------------------------------------------------------------------------
df_cross <- df_lfc %>%
  count(group, quadrant) %>%
  complete(group, quadrant, fill = list(n = 0)) %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_cross <- ggplot(df_cross, aes(x = quadrant, y = group, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08306B",
                      name = "% in group") +
  scale_y_discrete(labels = GROUP_LABELS) +
  scale_x_discrete(labels = QUAD_LABELS) +
  labs(
    title    = "元TMAOグループ × 2D象限 クロス集計",
    subtitle = "各セル: 人数 / グループ内比率",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.x  = element_text(size = 9, angle = 15, hjust = 1),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

# -----------------------------------------------------------------------------
# 7. X期の2D マップ（参考）
# -----------------------------------------------------------------------------
p_2d_X <- ggplot(df_lfc, aes(x = .data[[tmao_X]], y = .data[[gdca_X]],
                               color = group, label = Subject)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                     name = "元TMAOグループ") +
  labs(
    title    = "参考: X期の2D応答（早期介入）",
    x        = "TMAO log2FC (X期, T2-T3 vs T1)",
    y        = "GDCA log2FC (X期, T2-T3 vs T1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

# -----------------------------------------------------------------------------
# 8. 保存
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), plot = p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), plot = p, width = w, height = h, dpi = 150)
  cat(sprintf("  %s\n", name))
}

cat("\n保存完了:\n")
save_plot(p_2d,   "2d_response_map_Y",      8.5, 7)
save_plot(p_cross,"2d_group_comparison",     7,   4)
save_plot(p_2d_X, "2d_response_map_X",      7.5, 6)

# -----------------------------------------------------------------------------
# 9. コンソールサマリー
# -----------------------------------------------------------------------------
cat("\n=== 象限ごとの被験者一覧（Y期）===\n")
for (q in c("Q1","Q2","Q3","Q4")) {
  subjs <- df_lfc$Subject[!is.na(df_lfc$quadrant) & df_lfc$quadrant == q]
  grps  <- df_lfc$group[!is.na(df_lfc$quadrant) & df_lfc$quadrant == q]
  cat(sprintf("\n%s (%s): n=%d\n", q, QUAD_LABELS[[q]], length(subjs)))
  for (i in seq_along(subjs))
    cat(sprintf("  Subject %s  [%s]\n", subjs[i], grps[i]))
}
