# =============================================================================
# 02_tmao_gdca_correlation.R
# TMAO（統合シグナル）× GDCA（Alignment ID 4263）の関連解析
#
# 仮説 H1-A:
#   二次胆汁酸（GDCA = Glycodeoxycholate, ID 4263）は腸内細菌叢活性の
#   サロゲートマーカーであり、TMAOレスポンスが高い群で系統的に高値を示す。
#
# 解析内容:
#   1. TMAO × GDCA の期別 log2FC 散布図（被験者単位、グループで色分け）
#   2. TMAOレスポンスグループ別 GDCA ボックスプロット（X期・Y期）
#   3. 全被験者×全時点での TMAO × GDCA log2 強度相関
#   4. 時系列: TMAO（統合）と GDCA の grand mean を重ねてプロット
#
# 出力:
#   output/figures/tmao_gdca_scatter_lfc.pdf/.png   — 期別 log2FC 散布図
#   output/figures/tmao_gdca_boxplot.pdf/.png        — グループ別ボックスプロット
#   output/figures/tmao_gdca_corr_all.pdf/.png       — 全点 log2 強度相関
#   output/figures/tmao_gdca_timeseries.pdf/.png     — 時系列比較
#
# 前提:
#   1. load_data_production.R 実行済み
#   2. classify_subjects_by_response.R 実行済み（GROUPS_CSV が存在する）
#   3. 01_tmao_dual_peaks.R 実行済み（TMAO_INT_CSV が存在する）
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
GDCA_ID    <- 4263L
TMAO_A_ID  <- 414L
TMAO_B_ID  <- 415L

# classify_subjects_by_response.R の出力
GROUPS_CSV   <- "data/production/processed/subject_response_groups_414.csv"

# 01_tmao_dual_peaks.R の出力（統合 log2 シグナル）
TMAO_INT_CSV <- "data/production/processed/tmao_integrated.csv"

OUT_DIR <- "output/figures"

# グループカラー（analysis_roadmap.md と統一）
GROUP_COLORS <- c(
  both   = "#E41A1C",   # 赤: X期・Y期ともTMAO上昇
  X_only = "#377EB8",   # 青: X期のみ
  Y_only = "#4DAF4A",   # 緑: Y期のみ
  other  = "#984EA3"    # 紫: その他
)
GROUP_LABELS <- c(
  both   = "Both (X & Y)",
  X_only = "X only (early↑)",
  Y_only = "Y only (late↑)",
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

if (!file.exists(TMAO_INT_CSV))
  stop("Integrated TMAO CSV not found: ", TMAO_INT_CSV,
       "\nRun 01_tmao_dual_peaks.R first.")

if (!as.character(GDCA_ID) %in% rownames(feat_mat))
  stop(sprintf("Alignment_ID %d (GDCA) not found in feat_mat.", GDCA_ID))

# -----------------------------------------------------------------------------
# 1. パッケージ
# -----------------------------------------------------------------------------
pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# -----------------------------------------------------------------------------
# 2. データ読み込み
# -----------------------------------------------------------------------------
groups_df  <- read.csv(GROUPS_CSV,   stringsAsFactors = FALSE)
groups_df$Subject <- as.character(groups_df$Subject)

tmao_int   <- read.csv(TMAO_INT_CSV, stringsAsFactors = FALSE)
tmao_int$Subject   <- as.character(tmao_int$Subject)
tmao_int$Timepoint <- as.integer(tmao_int$Timepoint)

# -----------------------------------------------------------------------------
# 3. GDCA の平均行列（テクニカルレプリカートの平均）
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

gdca_raw <- sapply(group_cols, function(cols) {
  mean(as.numeric(feat_mat[as.character(GDCA_ID), cols]), na.rm = TRUE)
})

# log2 変換
nz_gdca    <- gdca_raw[gdca_raw > 0 & !is.na(gdca_raw)]
gdca_offset <- if (length(nz_gdca) > 0) min(nz_gdca) / 2 else 1
log2_gdca  <- log2(gdca_raw + gdca_offset)

df_gdca <- data.frame(
  SubjTP    = names(log2_gdca),
  Subject   = sub("-[0-9]+$", "", names(log2_gdca)),
  Timepoint = as.integer(sub(".*-", "", names(log2_gdca))),
  log2_gdca = log2_gdca,
  stringsAsFactors = FALSE
)

# log2FC（期内ベースライン差分）
df_gdca <- df_gdca %>%
  group_by(Subject) %>%
  mutate(
    ref_T1  = log2_gdca[Timepoint == 1],
    ref_T4  = log2_gdca[Timepoint == 4],
    lfc_gdca = case_when(
      Timepoint %in% 1:3 ~ log2_gdca - ref_T1,
      Timepoint %in% 4:6 ~ log2_gdca - ref_T4,
      TRUE               ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  select(-ref_T1, -ref_T4)

# -----------------------------------------------------------------------------
# 4. 被験者単位の期別 log2FC サマリー（X期平均・Y期平均）
# -----------------------------------------------------------------------------
# TMAO 統合シグナル
tmao_summary <- tmao_int %>%
  group_by(Subject) %>%
  summarise(
    TMAO_X = mean(log2FC[Timepoint %in% 2:3], na.rm = TRUE),  # T2,T3 vs T1
    TMAO_Y = mean(log2FC[Timepoint %in% 5:6], na.rm = TRUE),  # T5,T6 vs T4
    .groups = "drop"
  )

# GDCA
gdca_summary <- df_gdca %>%
  group_by(Subject) %>%
  summarise(
    GDCA_X = mean(lfc_gdca[Timepoint %in% 2:3], na.rm = TRUE),
    GDCA_Y = mean(lfc_gdca[Timepoint %in% 5:6], na.rm = TRUE),
    .groups = "drop"
  )

df_summary <- tmao_summary %>%
  inner_join(gdca_summary, by = "Subject") %>%
  left_join(groups_df[, c("Subject", "group")], by = "Subject") %>%
  mutate(group = ifelse(is.na(group), "other", group),
         group = factor(group, levels = names(GROUP_COLORS)))

cat(sprintf("\n=== 被験者数 ===\n"))
print(table(df_summary$group))

# 期別相関係数
for (period in c("X", "Y")) {
  tmao_col <- paste0("TMAO_", period)
  gdca_col <- paste0("GDCA_", period)
  r <- cor(df_summary[[tmao_col]], df_summary[[gdca_col]], use = "complete.obs")
  p <- cor.test(df_summary[[tmao_col]], df_summary[[gdca_col]])$p.value
  cat(sprintf("TMAO vs GDCA — %s期: r = %.3f, p = %.3f (n = %d)\n",
              period, r, p, sum(!is.na(df_summary[[tmao_col]]) & !is.na(df_summary[[gdca_col]]))))
}

# -----------------------------------------------------------------------------
# 5. 散布図: TMAO log2FC × GDCA log2FC（X期・Y期）
# -----------------------------------------------------------------------------
make_lfc_scatter <- function(df, tmao_col, gdca_col, period_label) {
  r   <- cor(df[[tmao_col]], df[[gdca_col]], use = "complete.obs")
  pv  <- cor.test(df[[tmao_col]], df[[gdca_col]])$p.value
  ann <- sprintf("r = %.3f, p = %.3f", r, pv)

  ggplot(df, aes(x = .data[[tmao_col]], y = .data[[gdca_col]],
                 color = group, label = Subject)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_smooth(method = "lm", se = TRUE, color = "grey50",
                linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE,
                aes(x = .data[[tmao_col]], y = .data[[gdca_col]])) +
    geom_point(size = 3, alpha = 0.85) +
    geom_text_repel(size = 2.8, show.legend = FALSE, max.overlaps = 15) +
    annotate("text", x = -Inf, y = Inf,
             label = ann, hjust = -0.05, vjust = 1.4,
             size = 3.5, color = "grey30") +
    scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                       name = "TMAO response") +
    labs(
      title = sprintf("%s: TMAO log2FC × GDCA log2FC", period_label),
      x     = "TMAO 統合シグナル log2FC",
      y     = "GDCA (ID 4263) log2FC"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = 11, face = "bold"),
          legend.position  = "right")
}

p_scatter_X <- make_lfc_scatter(df_summary, "TMAO_X", "GDCA_X", "Early (X期, T1–T3)")
p_scatter_Y <- make_lfc_scatter(df_summary, "TMAO_Y", "GDCA_Y", "Late  (Y期, T4–T6)")

p_scatter_combined <- (p_scatter_X / p_scatter_Y) +
  plot_annotation(
    title    = "TMAO × GDCA（Glycodeoxycholate）: 期別 log2FC 相関",
    subtitle = "被験者ごとに X期（T2-3 vs T1）・Y期（T5-6 vs T4）の平均 log2FC",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

# -----------------------------------------------------------------------------
# 6. グループ別ボックスプロット: GDCA の期別 log2FC
# -----------------------------------------------------------------------------
df_box <- df_summary %>%
  pivot_longer(cols = c(GDCA_X, GDCA_Y),
               names_to  = "period",
               values_to = "GDCA_lfc") %>%
  mutate(period = recode(period,
    GDCA_X = "Early (X期, T2-3 vs T1)",
    GDCA_Y = "Late  (Y期, T5-6 vs T4)"
  ))

p_boxplot <- ggplot(df_box, aes(x = group, y = GDCA_lfc, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60", linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.55) +
  geom_jitter(aes(color = group), width = 0.15, size = 2.0, alpha = 0.8) +
  facet_wrap(~period, ncol = 2) +
  scale_fill_manual(values  = GROUP_COLORS, labels = GROUP_LABELS, guide = "none") +
  scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS, name = "TMAO response") +
  scale_x_discrete(labels = GROUP_LABELS) +
  labs(
    title    = "TMAOレスポンスグループ別 GDCA（ID 4263）log2FC",
    subtitle = "水平破線: log2FC = 0（ベースライン比変化なし）",
    x        = NULL,
    y        = "GDCA log2FC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(angle = 20, hjust = 1, size = 9),
    legend.position   = "none",
    plot.title        = element_text(size = 11, face = "bold"),
    plot.subtitle     = element_text(size = 9,  color = "grey40"),
    strip.text        = element_text(size = 10, face = "bold")
  )

# -----------------------------------------------------------------------------
# 7. 全被験者×全時点の log2 強度相関（TMAO vs GDCA）
# -----------------------------------------------------------------------------
df_all <- tmao_int %>%
  select(SubjTP, Subject, Timepoint, log2_int) %>%
  inner_join(df_gdca %>% select(SubjTP, log2_gdca), by = "SubjTP") %>%
  left_join(groups_df[, c("Subject", "group")], by = "Subject") %>%
  mutate(group     = ifelse(is.na(group), "other", group),
         group     = factor(group, levels = names(GROUP_COLORS)),
         Timepoint = factor(Timepoint))

r_all  <- cor(df_all$log2_int, df_all$log2_gdca, use = "complete.obs")
pv_all <- cor.test(df_all$log2_int, df_all$log2_gdca)$p.value

p_corr_all <- ggplot(df_all, aes(x = log2_int, y = log2_gdca,
                                  color = group, shape = Timepoint)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50",
              linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE,
              aes(x = log2_int, y = log2_gdca)) +
  geom_point(size = 1.8, alpha = 0.7) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.3f, p = %.2e\n(全被験者×全時点, n = %d)",
                           r_all, pv_all, nrow(df_all)),
           hjust = -0.05, vjust = 1.4, size = 3.5, color = "grey30") +
  scale_color_manual(values = GROUP_COLORS, labels = GROUP_LABELS,
                     name = "TMAO response") +
  scale_shape_manual(values = c("1"=16,"2"=17,"3"=15,"4"=3,"5"=7,"6"=8),
                     name = "Timepoint") +
  labs(
    title    = "TMAO × GDCA log2 強度相関（全被験者×全時点）",
    subtitle = "各点 = 被験者 × 時点（テクニカルレプリカート平均）",
    x        = "TMAO 統合シグナル log2 強度",
    y        = "GDCA (ID 4263) log2 強度"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 11, face = "bold"),
        plot.subtitle    = element_text(size = 9,  color = "grey40"))

# -----------------------------------------------------------------------------
# 8. 時系列比較: TMAO 統合シグナル × GDCA grand mean
# -----------------------------------------------------------------------------
# 0-1 スケーリングで同一軸に重ねる
scale01 <- function(x) (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))

df_ts_tmao <- tmao_int %>%
  group_by(Timepoint) %>%
  summarise(mean_log2 = mean(log2_int, na.rm = TRUE), .groups = "drop") %>%
  mutate(scaled = scale01(mean_log2), signal = "TMAO (統合)")

df_ts_gdca <- df_gdca %>%
  group_by(Timepoint) %>%
  summarise(mean_log2 = mean(log2_gdca, na.rm = TRUE), .groups = "drop") %>%
  mutate(scaled = scale01(mean_log2), signal = "GDCA (ID 4263)")

df_ts <- bind_rows(df_ts_tmao, df_ts_gdca)

TS_COLORS <- c("TMAO (統合)" = "#E41A1C", "GDCA (ID 4263)" = "#377EB8")

p_ts <- ggplot(df_ts, aes(x = Timepoint, y = scaled, color = signal, group = signal)) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = 1:6) +
  scale_color_manual(values = TS_COLORS, name = NULL) +
  annotate("text", x = 2,   y = -0.06, label = "Early (X期)", size = 3.5, color = "grey40") +
  annotate("text", x = 4.8, y = -0.06, label = "Late (Y期)",  size = 3.5, color = "grey40") +
  labs(
    title    = "TMAO × GDCA 時系列（0–1 スケーリング, grand mean）",
    subtitle = "縦破線: Early / Late 期の境界  |  各シグナルを 0–1 正規化して比較",
    x        = "Timepoint",
    y        = "Scaled log2 intensity (0–1)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top",
        plot.title       = element_text(size = 11, face = "bold"),
        plot.subtitle    = element_text(size = 9,  color = "grey40"))

# -----------------------------------------------------------------------------
# 9. 保存
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), plot = p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), plot = p, width = w, height = h, dpi = 150)
  cat(sprintf("  %s.pdf/.png\n", name))
}

cat("\n保存完了:\n")
save_plot(p_scatter_combined, "tmao_gdca_scatter_lfc",  7.5, 11)
save_plot(p_boxplot,          "tmao_gdca_boxplot",       9,   5.5)
save_plot(p_corr_all,         "tmao_gdca_corr_all",      7.5,  6)
save_plot(p_ts,               "tmao_gdca_timeseries",    8,    5)

cat(sprintf("\n=== 解釈の手がかり ===\n"))
cat("  散布図で both 群が右上に集まる → TMAOとGDCAが連動して上昇する群が存在\n")
cat("  ボックスプロットで both 群の GDCA_Y が他群より高い → H1-A 支持\n")
cat("  時系列でTMAO・GDCAの動態が同期 → 腸内細菌叢活性の共通ドライバーを示唆\n")
