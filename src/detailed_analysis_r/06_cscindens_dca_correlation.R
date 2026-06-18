# ==============================================================================
# 06_cscindens_dca_correlation.R
# Purpose : Clostridium scindens abundance × DCA 強度の相関解析
#
#   腸内細菌 C. scindens は BSH（胆汁酸塩加水分解酵素）および 7α-HSDH を持ち，
#   CA → DCA 変換（7α-脱水酸素化）の主要な実施菌。
#   16S rRNA / メタゲノムデータが利用可能な場合，DCA 強度との相関を確認する。
#
# 想定する 16S rRNA データ形式 (data/ref/microbiome_16S.csv):
#   Subject, Timepoint (or 'baseline'), cscindens_rel_abund (or 相対存在量列)
#   └ 被験者ごとの C. scindens 相対存在量（%）
#
# 実際のファイルが存在しない場合は，期待されるフォーマットと代替戦略を
# コンソールに出力して終了する。
#
# 出力:
#   output/figures/cscindens_dca_scatter.pdf / .png
#   output/reports/cscindens_dca_correlation.csv
#
# 依存: install.packages(c("tidyverse", "patchwork", "ggrepel"))
# ==============================================================================

library(tidyverse)
library(patchwork)

# ── パス設定 ──────────────────────────────────────────────────────────────────
PROJ_ROOT  <- "."
LFC_CSV    <- file.path(PROJ_ROOT, "data/production/processed/log2FC_values_production.csv")
GROUPS_CSV <- file.path(PROJ_ROOT, "data/ref/tmao_subject_groups.csv")
OUT_FIG    <- file.path(PROJ_ROOT, "output/figures")
OUT_REP    <- file.path(PROJ_ROOT, "output/reports")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_REP, showWarnings = FALSE, recursive = TRUE)

# ── 16S rRNA データの検索候補パス ────────────────────────────────────────────
MICROBIOME_CANDIDATES <- c(
  file.path(PROJ_ROOT, "data/ref/microbiome_16S.csv"),
  file.path(PROJ_ROOT, "data/production/microbiome_16S.csv"),
  file.path(PROJ_ROOT, "data/production/raw/microbiome_16S.csv"),
  file.path(PROJ_ROOT, "data/ref/16S_abundance.csv"),
  file.path(PROJ_ROOT, "data/ref/cscindens_abundance.csv")
)

# ── DCA の Alignment_ID（過去解析で確認済み）──────────────────────────────────
DCA_ID          <- 10813L
DCA_ID_FALLBACK <- 3590L   # H1A 解析でも使われる別エントリ

# ── 16S rRNA データ存在確認 ────────────────────────────────────────────────────
micro_path <- NA_character_
for (p in MICROBIOME_CANDIDATES) {
  if (file.exists(p)) { micro_path <- p; break }
}

if (is.na(micro_path)) {
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # マイクロバイオームデータが存在しない場合のレポート
  # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  cat("
=======================================================================
  C. scindens × DCA 相関解析 — マイクロバイオームデータが見つかりません
=======================================================================

検索したパス:
", paste(" -", MICROBIOME_CANDIDATES, collapse = "\n"), "

【期待するファイル形式】
  ファイル名 : data/ref/microbiome_16S.csv（推奨）
  必須列     : Subject, cscindens_rel_abund
  任意列     : Timepoint（縦断データの場合），Sample_ID，total_reads

  例:
    Subject, Timepoint, cscindens_rel_abund, total_reads
    1,       1,         0.82,               85000
    1,       4,         1.24,               90000
    2,       1,         0.10,               78000
    ...

  * cscindens_rel_abund は C. scindens の相対存在量（%）
  * Timepoint が無い場合はベースライン値として Subject 列のみ使用

【代替戦略】
  16S rRNA データなしで解析を代替する方法:
  (1) DCA log2FC × TMAO log2FC の散布図（既存 H1A スクリプトで実施済み）
  (2) DCA の T4→T5 FC が高い被験者のグループが Both_up と一致するか確認
  (3) 公共データ（Human Microbiome Project 等）との比較
  (4) プロキシとして glycochenodeoxycholate (ID 3788) 等の
      C. scindens 関連代謝物でクラスタリング

【メモ】
  マイクロバイオームデータを入手した際は，上記パスに CSV を配置して
  このスクリプトを再実行してください。
=======================================================================
")
  message("マイクロバイオームデータが存在しないため，代替解析（DCA 時系列）のみ実行します。")

  # ── 代替解析: DCA × Both_up グループの時系列確認 ─────────────────────────
  if (!file.exists(LFC_CSV)) {
    message(sprintf("LFC CSV も見つかりません: %s\nスキプします。", LFC_CSV))
    quit(save = "no")
  }
  lfc    <- read_csv(LFC_CSV, show_col_types = FALSE)
  groups <- read_csv(GROUPS_CSV, show_col_types = FALSE)

  WEAK_GROUPS <- c("Sparse", "Both_down", "X_down_only", "Y_down_only", "Weak")
  GROUP_COLORS <- c(
    Both_up      = "#C00000", Y_only       = "#E06060",
    X_down_Y_up  = "#70AD47", X_up_Y_down  = "#ED7D31",
    X_only       = "#FFC000", Non_producer = "#4472C4",
    Weak         = "#A5A5A5"
  )

  grp_col <- if ("tmao_group" %in% colnames(groups)) "tmao_group" else colnames(groups)[2]
  groups <- groups %>%
    rename(tmao_group = all_of(grp_col)) %>%
    mutate(group = if_else(tmao_group %in% WEAK_GROUPS, "Weak", tmao_group),
           group = factor(group, levels = names(GROUP_COLORS)))

  # DCA ID を確認（優先 → フォールバック）
  avail_ids <- unique(lfc$Alignment_ID)
  use_dca   <- if (DCA_ID %in% avail_ids) DCA_ID else DCA_ID_FALLBACK
  if (!use_dca %in% avail_ids) {
    message(sprintf("DCA ID %d も %d も LFC CSV に存在しません。処理を中断します。",
                    DCA_ID, DCA_ID_FALLBACK))
    quit(save = "no")
  }
  cat(sprintf("使用する DCA Alignment_ID: %d\n", use_dca))

  dca_df <- lfc %>%
    filter(Alignment_ID == use_dca) %>%
    left_join(groups %>% select(Subject, group), by = "Subject") %>%
    filter(!is.na(group))

  tp_labels <- c("T1\n(X base)", "T2", "T3", "T4\n(Y base)", "T5", "T6")
  grp_sum <- dca_df %>%
    group_by(group, Timepoint) %>%
    summarise(mean_lfc = mean(intensity, na.rm = TRUE),
              se_lfc   = sd(intensity, na.rm = TRUE) / sqrt(n()),
              .groups  = "drop")

  p_alt <- ggplot() +
    annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
             fill = "#EEEEEE", alpha = 1) +
    geom_line(data = dca_df,
              aes(x = Timepoint, y = intensity, group = Subject),
              color = "#CCCCCC", linewidth = 0.4, alpha = 0.5) +
    geom_ribbon(data = grp_sum,
                aes(x = Timepoint, ymin = mean_lfc - se_lfc,
                    ymax = mean_lfc + se_lfc, fill = group),
                alpha = 0.25, inherit.aes = FALSE) +
    geom_line(data = grp_sum,
              aes(x = Timepoint, y = mean_lfc, color = group),
              linewidth = 1.5) +
    geom_point(data = grp_sum,
               aes(x = Timepoint, y = mean_lfc, color = group), size = 3) +
    geom_hline(yintercept = 0, color = "#AAAAAA", linetype = "dashed") +
    geom_vline(xintercept = 3.5, color = "#888888", linetype = "dashed") +
    scale_color_manual(values = GROUP_COLORS, name = "TMAO Group") +
    scale_fill_manual( values = GROUP_COLORS, name = "TMAO Group") +
    scale_x_continuous(breaks = 1:6, labels = tp_labels) +
    facet_wrap(~group, ncol = 4) +
    labs(
      title    = sprintf("DCA (ID %d) log₂FC — Individual Trajectories by TMAO Group", use_dca),
      subtitle = "16S rRNA データなし。DCA 時系列のみ表示。\nC. scindens データ入手後に再解析予定。",
      x = "Timepoint", y = "log₂FC (within-period)"
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = "#F0F0F0"),
          plot.title  = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "#666666"),
          panel.grid.minor = element_blank())

  out_png <- file.path(OUT_FIG, "cscindens_dca_nodata_timeseries.png")
  ggsave(out_png, p_alt, width = 14, height = 7, dpi = 150)
  message(sprintf("代替図（DCA 時系列）保存: %s", out_png))

  quit(save = "no")
}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# マイクロバイオームデータが存在する場合の本解析
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
cat(sprintf("16S rRNA データ読み込み: %s\n", micro_path))
micro <- read_csv(micro_path, show_col_types = FALSE)
cat(sprintf("  %d 行 × %d 列: %s\n", nrow(micro), ncol(micro),
            paste(colnames(micro), collapse = ", ")))

# ── C. scindens 列を自動検出 ─────────────────────────────────────────────────
cs_col <- grep("scindens|cscindens|c\\.scindens", colnames(micro),
               ignore.case = TRUE, value = TRUE)[1]
if (is.na(cs_col)) {
  stop("C. scindens 存在量列が見つかりません。列名を確認してください。\n",
       "  検索対象列名: 'scindens', 'cscindens', 'c.scindens'\n",
       "  実際の列名  : ", paste(colnames(micro), collapse = ", "))
}
cat(sprintf("C. scindens 列: '%s'\n", cs_col))
micro <- micro %>% rename(cscindens = all_of(cs_col))

# ── DCA の log2FC を取得 ───────────────────────────────────────────────────────
if (!file.exists(LFC_CSV)) stop(sprintf("LFC CSV が見つかりません: %s", LFC_CSV))
lfc <- read_csv(LFC_CSV, show_col_types = FALSE)

avail_ids <- unique(lfc$Alignment_ID)
use_dca   <- if (DCA_ID %in% avail_ids) DCA_ID else DCA_ID_FALLBACK
cat(sprintf("使用する DCA Alignment_ID: %d\n", use_dca))

dca_lfc <- lfc %>%
  filter(Alignment_ID == use_dca) %>%
  select(Subject, Timepoint, dca_lfc = intensity)

# ── TMAO グループ ────────────────────────────────────────────────────────────
groups <- read_csv(GROUPS_CSV, show_col_types = FALSE)
WEAK_GROUPS <- c("Sparse", "Both_down", "X_down_only", "Y_down_only", "Weak")
GROUP_COLORS <- c(
  Both_up = "#C00000", Y_only = "#E06060",
  X_down_Y_up = "#70AD47", X_up_Y_down = "#ED7D31",
  X_only = "#FFC000", Non_producer = "#4472C4", Weak = "#A5A5A5"
)
grp_col <- if ("tmao_group" %in% colnames(groups)) "tmao_group" else colnames(groups)[2]
groups <- groups %>%
  rename(tmao_group = all_of(grp_col)) %>%
  mutate(group = if_else(tmao_group %in% WEAK_GROUPS, "Weak", tmao_group),
         group = factor(group, levels = names(GROUP_COLORS)))

# ── マイクロバイオームデータのマージ ─────────────────────────────────────────
# Timepoint 列がある場合: 縦断的に DCA log2FC とマージ
# ない場合: ベースライン値として全期間の DCA 平均とマージ
has_tp <- "Timepoint" %in% colnames(micro)

if (has_tp) {
  cat("縦断的マイクロバイオームデータとして処理\n")
  merged <- dca_lfc %>%
    inner_join(micro %>% select(Subject, Timepoint, cscindens), by = c("Subject", "Timepoint")) %>%
    left_join(groups %>% select(Subject, group), by = "Subject")
} else {
  cat("ベースラインマイクロバイオームデータとして処理（全時点の DCA 平均と相関）\n")
  dca_mean_subj <- dca_lfc %>%
    group_by(Subject) %>%
    summarise(dca_mean = mean(dca_lfc, na.rm = TRUE), .groups = "drop")
  merged <- micro %>%
    select(Subject, cscindens) %>%
    inner_join(dca_mean_subj, by = "Subject") %>%
    left_join(groups %>% select(Subject, group), by = "Subject") %>%
    rename(dca_lfc = dca_mean)
}

cat(sprintf("マージ後: %d 行\n", nrow(merged)))

# ── Pearson 相関 ──────────────────────────────────────────────────────────────
cor_all <- tryCatch({
  ct <- cor.test(merged$cscindens, merged$dca_lfc, use = "complete.obs")
  list(r = round(ct$estimate, 3), p = round(ct$p.value, 4),
       label = sprintf("r = %.2f\np = %.3f%s",
                       ct$estimate, ct$p.value,
                       if (ct$p.value < 0.05) " *" else ""))
}, error = function(e) list(r = NA, p = NA, label = "n.a."))

cat(sprintf("\n全被験者 Pearson 相関: r = %s, p = %s\n",
            cor_all$r, cor_all$p))

# ── ggrepel があれば使用 ──────────────────────────────────────────────────────
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel) else
  message("ggrepel が未インストールのため geom_text を使用します。")

# ── 散布図 ───────────────────────────────────────────────────────────────────
y_lab <- if (has_tp) "DCA log₂FC (within-period)" else "DCA log₂FC (subject mean)"

p_scatter <- ggplot(merged %>% filter(!is.na(group)),
                    aes(x = cscindens, y = dca_lfc, color = group)) +
  geom_hline(yintercept = 0, color = "#CCCCCC") +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "#444444", fill = "#DDDDDD", linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.85) +
  {
    if (has_ggrepel)
      ggrepel::geom_text_repel(aes(label = Subject),
                               size = 2.8, max.overlaps = 15,
                               color = "#444444", segment.color = "#AAAAAA")
    else
      geom_text(aes(label = Subject), size = 2.8, vjust = -0.7, color = "#444444")
  } +
  annotate("label", x = Inf, y = Inf,
           label = cor_all$label,
           hjust = 1.1, vjust = 1.2, size = 3.8,
           color = if (!is.na(cor_all$p) && cor_all$p < 0.05) "#C00000" else "#555555",
           fontface = "bold", fill = "white", label.size = 0.4) +
  scale_color_manual(values = GROUP_COLORS, name = "TMAO Group") +
  labs(
    title    = sprintf("C. scindens Abundance × DCA (ID %d) log₂FC", use_dca),
    subtitle = sprintf("n = %d observations  |  Pearson correlation  |  Source: %s",
                       sum(!is.na(merged$cscindens) & !is.na(merged$dca_lfc)),
                       basename(micro_path)),
    x = "C. scindens relative abundance (%)",
    y = y_lab
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank())

# ── グループ別相関係数テーブル ────────────────────────────────────────────────
cor_by_group <- merged %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    n   = sum(!is.na(cscindens) & !is.na(dca_lfc)),
    r   = tryCatch(round(cor.test(cscindens, dca_lfc)$estimate, 3),
                   error = function(e) NA_real_),
    p   = tryCatch(round(cor.test(cscindens, dca_lfc)$p.value, 4),
                   error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                         p < 0.05  ~ "*",   TRUE ~ ""))

cat("\n=== グループ別 Pearson 相関 ===\n")
print(cor_by_group)

# ── 保存 ──────────────────────────────────────────────────────────────────────
out_pdf <- file.path(OUT_FIG, "cscindens_dca_scatter.pdf")
ggsave(out_pdf, p_scatter, width = 8, height = 6)
message(sprintf("PDF 保存: %s", out_pdf))

out_png <- file.path(OUT_FIG, "cscindens_dca_scatter.png")
ggsave(out_png, p_scatter, width = 8, height = 6, dpi = 150)
message(sprintf("PNG 保存: %s", out_png))

out_csv <- file.path(OUT_REP, "cscindens_dca_correlation.csv")
write_csv(bind_rows(
  merged %>%
    filter(!is.na(group)) %>%
    select(Subject, Timepoint = any_of("Timepoint"),
           cscindens, dca_lfc, group),
  .id = NULL
), out_csv)
message(sprintf("CSV 保存: %s", out_csv))

message("\n06_cscindens_dca_correlation.R 完了。")
