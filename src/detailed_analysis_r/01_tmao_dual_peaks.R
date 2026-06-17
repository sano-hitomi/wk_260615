# =============================================================================
# 01_tmao_dual_peaks.R
# TMAO 2峰（Alignment ID 414・415）の相関確認と統合シグナル作成
#
# 目的:
#   MS-DIALでTMAO（HMDB0000925）として検出された2つのピーク（ID 414・415）が
#   同一化合物を捉えているかを確認し、後続解析で使う統合シグナルを決定する。
#
# 出力:
#   output/figures/tmao_dual_peaks_scatter.pdf  — 414 vs 415 相関散布図
#   output/figures/tmao_dual_peaks_timeseries.pdf — 時系列（414・415・統合値）
#   data/production/processed/tmao_integrated.csv — 被験者×時点の統合TMAOシグナル
#
# 統合ロジック:
#   Pearson r ≥ CORR_THRESHOLD → sum（生強度を加算してから log2変換）
#   r < CORR_THRESHOLD         → 414 単独（メッセージで通知）
#
# 前提:
#   load_data_production.R 実行済み（samplesheet / feat_meta / feat_mat）
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
ID_A           <- 414L          # TMAO peak A
ID_B           <- 415L          # TMAO peak B
CORR_THRESHOLD <- 0.80          # Pearson r の閾値
OUT_DIR        <- "output/figures"
OUT_CSV        <- "data/production/processed/tmao_integrated.csv"

# -----------------------------------------------------------------------------
# 0. 前提オブジェクト確認
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
}
for (aid in c(ID_A, ID_B)) {
  if (!as.character(aid) %in% rownames(feat_mat))
    stop(sprintf("Alignment_ID %d not found in feat_mat.", aid))
}

# -----------------------------------------------------------------------------
# 1. パッケージ
# -----------------------------------------------------------------------------
pkgs <- c("ggplot2", "ggrepel", "dplyr", "tidyr", "patchwork")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# -----------------------------------------------------------------------------
# 2. 生物サンプルの平均行列（テクニカルレプリカートの平均）
# -----------------------------------------------------------------------------
bio <- samplesheet[
  samplesheet$type == "biological" & samplesheet$rerun_suffix == "",
]

group_key  <- paste(bio$Subject, bio$Timepoint, sep = "-")
group_cols <- split(bio$label, group_key)
group_cols <- lapply(group_cols, function(cols) intersect(cols, colnames(feat_mat)))
group_cols <- group_cols[sapply(group_cols, length) > 0]

sub_ids <- as.character(c(ID_A, ID_B))
avg_raw <- sapply(group_cols, function(cols) {
  rowMeans(feat_mat[sub_ids, cols, drop = FALSE], na.rm = TRUE)
})
# avg_raw: 2 × n_SubjTP

# -----------------------------------------------------------------------------
# 3. 対数変換
# -----------------------------------------------------------------------------
log2_with_offset <- function(vec) {
  nz <- vec[vec > 0 & !is.na(vec)]
  offset <- if (length(nz) > 0) min(nz) / 2 else 1
  log2(vec + offset)
}

log2_A <- log2_with_offset(avg_raw[1, ])
log2_B <- log2_with_offset(avg_raw[2, ])

# -----------------------------------------------------------------------------
# 4. Pearson 相関の計算
# -----------------------------------------------------------------------------
r_val  <- cor(log2_A, log2_B, use = "complete.obs")
r_sq   <- r_val^2
n_pts  <- sum(!is.na(log2_A) & !is.na(log2_B))
p_val  <- cor.test(log2_A, log2_B, method = "pearson")$p.value

cat(sprintf("\n=== TMAO 2峰 相関 (ID %d vs ID %d) ===\n", ID_A, ID_B))
cat(sprintf("  Pearson r   = %.4f\n", r_val))
cat(sprintf("  R²          = %.4f\n", r_sq))
cat(sprintf("  p-value     = %.2e\n", p_val))
cat(sprintf("  n (SubjTP)  = %d\n\n", n_pts))

# 統合方針を決定
if (r_val >= CORR_THRESHOLD) {
  integrate_method <- "sum"
  cat(sprintf("→ r ≥ %.2f: 統合シグナル = 生強度の和（ID %d + ID %d）を log2変換\n",
              CORR_THRESHOLD, ID_A, ID_B))
} else {
  integrate_method <- "single"
  cat(sprintf("→ r < %.2f: 相関が低いため ID %d のみを使用\n",
              CORR_THRESHOLD, ID_A))
}

# -----------------------------------------------------------------------------
# 5. 統合シグナル作成
# -----------------------------------------------------------------------------
if (integrate_method == "sum") {
  raw_sum   <- avg_raw[1, ] + avg_raw[2, ]   # 生強度の和
  log2_int  <- log2_with_offset(raw_sum)
  int_label <- sprintf("TMAO (ID%d+ID%d, log2 sum)", ID_A, ID_B)
} else {
  log2_int  <- log2_A
  int_label <- sprintf("TMAO (ID%d only, log2)", ID_A)
}

# SubjTP → Subject / Timepoint
subj_tp_names <- names(log2_int)
subject_vec   <- sub("-[0-9]+$", "", subj_tp_names)
timepoint_vec <- as.integer(sub(".*-", "", subj_tp_names))

df_int <- data.frame(
  SubjTP       = subj_tp_names,
  Subject      = subject_vec,
  Timepoint    = timepoint_vec,
  log2_A       = log2_A,
  log2_B       = log2_B,
  log2_int     = log2_int,
  stringsAsFactors = FALSE
)

# log2FC（期内ベースライン差分）
df_int <- df_int %>%
  group_by(Subject) %>%
  mutate(
    ref_T1   = log2_int[Timepoint == 1],
    ref_T4   = log2_int[Timepoint == 4],
    log2FC   = case_when(
      Timepoint %in% 1:3 ~ log2_int - ref_T1,
      Timepoint %in% 4:6 ~ log2_int - ref_T4,
      TRUE               ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  select(-ref_T1, -ref_T4)

# CSV 保存
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(df_int, OUT_CSV, row.names = FALSE)
cat(sprintf("統合シグナル CSV 保存: %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 6. 散布図: ID 414 vs ID 415（log2 強度）
# -----------------------------------------------------------------------------
df_scatter <- data.frame(
  SubjTP    = subj_tp_names,
  Subject   = subject_vec,
  Timepoint = as.factor(timepoint_vec),
  log2_A    = log2_A,
  log2_B    = log2_B,
  stringsAsFactors = FALSE
)

annotation_text <- sprintf("r = %.3f,  R² = %.3f\np = %.2e  (n = %d)", r_val, r_sq, p_val, n_pts)

p_scatter <- ggplot(df_scatter, aes(x = log2_A, y = log2_B, color = Timepoint)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE,
              aes(x = log2_A, y = log2_B)) +
  annotate("text", x = -Inf, y = Inf,
           label = annotation_text,
           hjust = -0.05, vjust = 1.3,
           size = 3.5, color = "grey30") +
  scale_color_brewer(palette = "Set2", name = "Timepoint") +
  labs(
    title    = sprintf("TMAO 2峰の相関: ID %d vs ID %d (log2 強度)", ID_A, ID_B),
    subtitle = sprintf("統合方針: %s", int_label),
    x        = sprintf("log2 Intensity  (ID %d)", ID_A),
    y        = sprintf("log2 Intensity  (ID %d)", ID_B)
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40"))

# -----------------------------------------------------------------------------
# 7. 時系列プロット: 414・415・統合値の grand mean ± 個別被験者
# -----------------------------------------------------------------------------
df_ts <- df_int %>%
  pivot_longer(cols = c(log2_A, log2_B, log2_int),
               names_to = "signal", values_to = "value") %>%
  mutate(signal = recode(signal,
    log2_A   = sprintf("ID %d", ID_A),
    log2_B   = sprintf("ID %d", ID_B),
    log2_int = "統合シグナル"
  ))

# Grand mean
df_gm <- df_ts %>%
  group_by(Timepoint, signal) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

SIGNAL_COLORS <- c(
  setNames(c("#4DAF4A", "#377EB8", "#E41A1C"),
           c(sprintf("ID %d", ID_A), sprintf("ID %d", ID_B), "統合シグナル"))
)

make_ts_panel <- function(df_ts_sub, df_gm_sub, signal_name, tp_range, period_label) {
  df_s  <- df_ts_sub[df_ts_sub$signal == signal_name &
                       df_ts_sub$Timepoint %in% tp_range, ]
  df_gm_s <- df_gm_sub[df_gm_sub$signal == signal_name &
                          df_gm_sub$Timepoint %in% tp_range, ]
  col <- SIGNAL_COLORS[[signal_name]]
  ggplot(df_s, aes(x = Timepoint, y = value)) +
    geom_line(aes(group = Subject), color = col, alpha = 0.25, linewidth = 0.5) +
    geom_line(data = df_gm_s, aes(x = Timepoint, y = mean_val, group = 1),
              color = col, linewidth = 1.8) +
    geom_point(data = df_gm_s, aes(x = Timepoint, y = mean_val),
               color = col, size = 3) +
    scale_x_continuous(breaks = tp_range) +
    labs(title = sprintf("%s — %s", signal_name, period_label),
         x = "Timepoint", y = "log2 Intensity") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 9, face = "bold"))
}

signals <- c(sprintf("ID %d", ID_A), sprintf("ID %d", ID_B), "統合シグナル")
panels <- list()
for (sig in signals) {
  panels[[length(panels) + 1]] <- make_ts_panel(df_ts, df_gm, sig, 1:3, "Early (T1–T3)")
  panels[[length(panels) + 1]] <- make_ts_panel(df_ts, df_gm, sig, 4:6, "Late  (T4–T6)")
}

p_ts <- wrap_plots(panels, ncol = 2, nrow = 3) +
  plot_annotation(
    title    = "TMAO 2峰 + 統合シグナル 時系列（log2 強度）",
    subtitle = "細線: 個別被験者  |  太線: grand mean",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40")
    )
  )

# -----------------------------------------------------------------------------
# 8. 保存
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(OUT_DIR, "tmao_dual_peaks_scatter.pdf"),
       plot = p_scatter, width = 6.5, height = 5.5)
ggsave(file.path(OUT_DIR, "tmao_dual_peaks_scatter.png"),
       plot = p_scatter, width = 6.5, height = 5.5, dpi = 150)

ggsave(file.path(OUT_DIR, "tmao_dual_peaks_timeseries.pdf"),
       plot = p_ts, width = 8, height = 10)
ggsave(file.path(OUT_DIR, "tmao_dual_peaks_timeseries.png"),
       plot = p_ts, width = 8, height = 10, dpi = 150)

cat(sprintf("\n保存完了:\n"))
cat(sprintf("  %s/tmao_dual_peaks_scatter.pdf/.png\n", OUT_DIR))
cat(sprintf("  %s/tmao_dual_peaks_timeseries.pdf/.png\n", OUT_DIR))
cat(sprintf("  %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 9. 後続スクリプトへの引き継ぎ情報
# -----------------------------------------------------------------------------
cat(sprintf("\n=== 後続スクリプトへの引き継ぎ ===\n"))
cat(sprintf("  統合方針        : %s\n", integrate_method))
cat(sprintf("  統合シグナル列  : log2_int  (in tmao_integrated.csv)\n"))
cat(sprintf("  CSV path        : %s\n", OUT_CSV))
cat(sprintf("  → 02_tmao_gdca_correlation.R でこの CSV を読み込んで使用\n"))
