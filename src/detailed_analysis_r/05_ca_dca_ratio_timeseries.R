# ==============================================================================
# 05_ca_dca_ratio_timeseries.R
# Purpose : 個別被験者レベルで CA/DCA 比率を時系列（T1-T6）で可視化
#
#   - 既存解析はグループ中央値でまとめているため，被験者内変動を捉えにくい
#   - 本スクリプトは絶対強度（log2変換済み）を用いて CA/DCA 比率を算出し，
#     スパゲッティプロット（個別線 + グループ平均）として出力する
#   - DCA の Alignment_ID は feature_metadata から動的に検索する
#     （ID 10813 を優先するが，InChIKey でもフォールバック検索する）
#
# 入力:
#   data/production/processed/samplesheet.csv
#   data/production/processed/feature_metadata.csv
#   data/production/processed/feature_matrix.csv
#   data/ref/tmao_subject_groups.csv
#
# 出力:
#   output/figures/ca_dca_ratio_individual.pdf   -- スパゲッティ + グループ平均
#   output/figures/ca_dca_ratio_individual.png
#   output/reports/ca_dca_ratio_table.csv        -- 被験者×時点×比率テーブル
#
# 前提:
#   プロジェクトルートを作業ディレクトリに設定してから Source すること。
#   依存: install.packages(c("tidyverse", "patchwork", "scales"))
# ==============================================================================

library(tidyverse)
library(patchwork)
library(scales)

# ── パス設定 ──────────────────────────────────────────────────────────────────
PROJ_ROOT   <- "."
SS_CSV      <- file.path(PROJ_ROOT, "data/production/processed/samplesheet.csv")
FMETA_CSV   <- file.path(PROJ_ROOT, "data/production/processed/feature_metadata.csv")
FMAT_CSV    <- file.path(PROJ_ROOT, "data/production/processed/feature_matrix.csv")
GROUPS_CSV  <- file.path(PROJ_ROOT, "data/ref/tmao_subject_groups.csv")
OUT_FIG     <- file.path(PROJ_ROOT, "output/figures")
OUT_REP     <- file.path(PROJ_ROOT, "output/reports")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_REP, showWarnings = FALSE, recursive = TRUE)

# ── InChIKey プレフィックス（14文字） ─────────────────────────────────────────
# Cholic acid (CA)              : BHQCQFFYRZLCQQ
# Deoxycholic acid (DCA)        : KXGVEGMKQFWNSR
IK_CA  <- "BHQCQFFYRZLCQQ"
IK_DCA <- "KXGVEGMKQFWNSR"

# HMDB ID
HMDB_CA  <- "HMDB0000619"
HMDB_DCA <- "HMDB0000626"

# DCA は過去解析で ID 10813 が同定済み。
# 以下で feature_metadata から確認 → 見つからなければ名前/InChIKey で再検索。
PREFERRED_DCA_ID <- 10813L

# ── データ読み込み ────────────────────────────────────────────────────────────
for (f in c(SS_CSV, FMETA_CSV, FMAT_CSV)) {
  if (!file.exists(f)) stop(sprintf("ファイルが見つかりません: %s", f))
}

samplesheet <- read_csv(SS_CSV, show_col_types = FALSE)
feat_meta   <- read_csv(FMETA_CSV, show_col_types = FALSE)
feat_mat_raw <- read_csv(FMAT_CSV, show_col_types = FALSE)

# feat_mat: 1列目が Alignment_ID, 残りがサンプルラベル
feat_mat <- feat_mat_raw %>%
  rename(Alignment_ID = 1) %>%
  mutate(Alignment_ID = as.integer(Alignment_ID))

cat(sprintf("samplesheet : %d rows\n", nrow(samplesheet)))
cat(sprintf("feat_meta   : %d features\n", nrow(feat_meta)))
cat(sprintf("feat_mat    : %d features x %d samples\n",
            nrow(feat_mat), ncol(feat_mat) - 1))

# ── CA / DCA の Alignment_ID を検索 ──────────────────────────────────────────

find_ids <- function(meta, ik_prefix, hmdb_id, name_pattern,
                     name_exclude = NULL, preferred_id = NULL) {
  # NA を除去してから返す共通処理
  clean <- function(x) sort(unique(x[!is.na(x)]))

  # 1. Preferred ID が feat_meta に存在するか確認
  if (!is.null(preferred_id) && !is.na(preferred_id) &&
      preferred_id %in% meta$Alignment_ID) {
    cat(sprintf("  preferred ID %d が feat_meta に存在します\n", preferred_id))
    return(preferred_id)
  }
  # 2. InChIKey プレフィックスで検索
  ik_col <- grep("inchikey|INCHIKEY", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(ik_col)) {
    ik_vals <- meta[[ik_col]]
    hits <- clean(meta$Alignment_ID[
      !is.na(ik_vals) & startsWith(toupper(ik_vals), toupper(ik_prefix))
    ])
    if (length(hits) > 0) {
      cat(sprintf("  InChIKey (%s...) で %d 件ヒット: %s\n",
                  ik_prefix, length(hits), paste(hits, collapse = ", ")))
      return(hits)
    }
  }
  # 3. HMDB ID で検索
  hmdb_col <- grep("hmdb", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(hmdb_col)) {
    hm_vals <- meta[[hmdb_col]]
    hits <- clean(meta$Alignment_ID[
      !is.na(hm_vals) & grepl(hmdb_id, hm_vals, ignore.case = TRUE)
    ])
    if (length(hits) > 0) {
      cat(sprintf("  HMDB ID (%s) で %d 件ヒット: %s\n",
                  hmdb_id, length(hits), paste(hits, collapse = ", ")))
      return(hits)
    }
  }
  # 4. 名前パターンで検索（除外パターン付き）
  name_col <- grep("^(Metabolite.name|Metabolite_name|Name)", colnames(meta),
                   ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(name_col)) {
    nm_vals <- meta[[name_col]]
    keep <- !is.na(nm_vals) & grepl(name_pattern, nm_vals, ignore.case = TRUE)
    if (!is.null(name_exclude))
      keep <- keep & !grepl(name_exclude, nm_vals, ignore.case = TRUE)
    hits <- clean(meta$Alignment_ID[keep])
    if (length(hits) > 0) {
      cat(sprintf("  名前パターン ('%s', 除外: '%s') で %d 件ヒット: %s\n",
                  name_pattern,
                  if (is.null(name_exclude)) "" else name_exclude,
                  length(hits), paste(hits, collapse = ", ")))
      return(hits)
    }
  }
  integer(0)
}

cat("CA を検索中...\n")
# name_exclude で "deoxy", "tauro", "glyco" 等の修飾 CA を除く
ca_ids  <- find_ids(feat_meta, IK_CA, HMDB_CA,
                    name_pattern = "cholic acid|cholate",
                    name_exclude = "deoxy|tauro|glyco|sulfo|taurocheno|glycheno")

cat("DCA を検索中...\n")
dca_ids <- find_ids(feat_meta, IK_DCA, HMDB_DCA,
                    name_pattern = "deoxycholic acid|\\bDCA\\b",
                    name_exclude = "ursodeoxy|tauro|glyco|chenodeoxy|hyodeoxy",
                    preferred_id = PREFERRED_DCA_ID)

# ── CA が見つからない場合の診断メッセージ ─────────────────────────────────────
if (length(ca_ids) == 0) {
  name_col_diag <- grep("^(Metabolite.name|Metabolite_name|Name)", colnames(feat_meta),
                        ignore.case = TRUE, value = TRUE)[1]
  cat("\n[診断] feature_metadata 中の 'cholic' を含む代謝物名:\n")
  if (!is.na(name_col_diag)) {
    chol_rows <- feat_meta[grepl("cholic", feat_meta[[name_col_diag]], ignore.case = TRUE) &
                             !is.na(feat_meta[[name_col_diag]]), ]
    if (nrow(chol_rows) > 0) {
      print(chol_rows %>% select(Alignment_ID, all_of(name_col_diag)) %>% head(20))
    } else {
      cat("  'cholic' を含む名前が見つかりません。CA は未検出または Unknown の可能性があります。\n")
    }
  }
  stop("CA (Cholic acid) が feature_metadata に見つかりませんでした。\n",
       "  上記の診断結果を確認し、スクリプト先頭の IK_CA / HMDB_CA を修正するか、\n",
       "  CA_ID を手動で設定してください（例: CA_ID <- 1234L）。")
}
if (length(dca_ids) == 0) stop("DCA (Deoxycholic acid) が feature_metadata に見つかりませんでした。")

# 複数ヒットした場合は最も S/N が高いものを使う（または最初の1つ）
pick_best <- function(ids, meta) {
  ids <- ids[!is.na(ids)]          # 念のため NA を除去
  if (length(ids) == 0) return(NA_integer_)
  if (length(ids) == 1) return(ids)
  sn_col <- grep("S.N|signal.noise|SN_average|sn_average", colnames(meta),
                 ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sn_col)) {
    sub <- meta[meta$Alignment_ID %in% ids, c("Alignment_ID", sn_col)]
    sub[[sn_col]] <- as.numeric(sub[[sn_col]])
    if (nrow(sub) > 0 && any(!is.na(sub[[sn_col]]))) {
      best <- sub$Alignment_ID[which.max(sub[[sn_col]])]
      cat(sprintf("  複数ヒット → S/N 最大の ID %d を採用\n", best))
      return(best)
    }
  }
  cat(sprintf("  複数ヒット → 最初の ID %d を採用\n", ids[1]))
  ids[1]
}

CA_ID  <- pick_best(ca_ids,  feat_meta)
DCA_ID <- pick_best(dca_ids, feat_meta)

cat(sprintf("\n使用する Alignment_ID: CA = %d, DCA = %d\n", CA_ID, DCA_ID))

# ── 生物学的サンプルのラベル一覧 ─────────────────────────────────────────────
bio <- samplesheet %>% filter(type == "biological")

# ── 絶対強度を取得（log2 変換）────────────────────────────────────────────────
get_intensity_long <- function(aid) {
  row <- feat_mat %>% filter(Alignment_ID == aid)
  if (nrow(row) == 0) {
    warning(sprintf("Alignment_ID %d が feature_matrix に存在しません", aid))
    return(NULL)
  }
  row %>%
    select(-Alignment_ID) %>%
    pivot_longer(everything(), names_to = "label", values_to = "raw_intensity") %>%
    mutate(raw_intensity = as.numeric(raw_intensity))
}

ca_long  <- get_intensity_long(CA_ID)  %>% rename(ca_int  = raw_intensity)
dca_long <- get_intensity_long(DCA_ID) %>% rename(dca_int = raw_intensity)

if (is.null(ca_long) || is.null(dca_long)) stop("CA または DCA の強度データが取得できませんでした。")

# ── サンプルシートとマージ → Subject × Timepoint に平均化 ─────────────────────
intensity_wide <- bio %>%
  select(label, Subject, Timepoint, period) %>%
  left_join(ca_long,  by = "label") %>%
  left_join(dca_long, by = "label") %>%
  filter(!is.na(ca_int), !is.na(dca_int), ca_int > 0, dca_int > 0) %>%
  group_by(Subject, Timepoint, period) %>%
  summarise(
    ca_mean  = mean(ca_int,  na.rm = TRUE),
    dca_mean = mean(dca_int, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    log2_ca    = log2(ca_mean),
    log2_dca   = log2(dca_mean),
    ca_dca_ratio = ca_mean / dca_mean,
    log2_ratio   = log2_ca - log2_dca
  )

cat(sprintf("\n強度データ: %d 被験者 × %d 時点\n",
            n_distinct(intensity_wide$Subject),
            n_distinct(intensity_wide$Timepoint)))

# ── TMAO グループを付加 ───────────────────────────────────────────────────────
groups_raw <- read_csv(GROUPS_CSV, show_col_types = FALSE)

# グループ列の柔軟な取得
grp_col <- if ("tmao_group" %in% colnames(groups_raw)) "tmao_group" else colnames(groups_raw)[2]

WEAK_GROUPS <- c("Sparse", "Both_down", "X_down_only", "Y_down_only", "Weak")
GROUP_COLORS <- c(
  Both_up      = "#C00000", Y_only       = "#E06060",
  X_down_Y_up  = "#70AD47", X_up_Y_down  = "#ED7D31",
  X_only       = "#FFC000", Non_producer = "#4472C4",
  Weak         = "#A5A5A5"
)
GROUP_ORDER <- names(GROUP_COLORS)

groups <- groups_raw %>%
  rename(tmao_group = all_of(grp_col)) %>%
  mutate(group = if_else(tmao_group %in% WEAK_GROUPS, "Weak", tmao_group),
         group = factor(group, levels = GROUP_ORDER))

df_plot <- intensity_wide %>%
  left_join(groups %>% select(Subject, group), by = "Subject")

# ── グループ × 時点の平均・SD ─────────────────────────────────────────────────
group_summary <- df_plot %>%
  filter(!is.na(group)) %>%
  group_by(group, Timepoint) %>%
  summarise(
    mean_ratio = mean(log2_ratio, na.rm = TRUE),
    sd_ratio   = sd(log2_ratio,   na.rm = TRUE),
    n          = n(),
    se_ratio   = sd_ratio / sqrt(n),
    .groups    = "drop"
  )

# ── TP ラベル ─────────────────────────────────────────────────────────────────
tp_labels <- c("T1\n(X base)", "T2", "T3", "T4\n(Y base)", "T5", "T6")

# ══════════════════════════════════════════════════════════════════════════════
# 図 1: 全被験者スパゲッティ + グループ平均（4×2 ファセット）
# ══════════════════════════════════════════════════════════════════════════════
p_all <- ggplot() +
  # 背景：X 期
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = "#EEEEEE", alpha = 1) +
  # 個別線（薄いグレー）
  geom_line(data = df_plot %>% filter(!is.na(group)),
            aes(x = Timepoint, y = log2_ratio, group = Subject),
            color = "#BBBBBB", linewidth = 0.4, alpha = 0.6) +
  geom_point(data = df_plot %>% filter(!is.na(group)),
             aes(x = Timepoint, y = log2_ratio, group = Subject),
             color = "#BBBBBB", size = 0.8, alpha = 0.6) +
  # ± SD リボン
  geom_ribbon(data = group_summary,
              aes(x = Timepoint,
                  ymin = mean_ratio - se_ratio,
                  ymax = mean_ratio + se_ratio,
                  fill = group),
              alpha = 0.25, inherit.aes = FALSE) +
  # グループ平均線
  geom_line(data = group_summary,
            aes(x = Timepoint, y = mean_ratio, color = group),
            linewidth = 1.5) +
  geom_point(data = group_summary,
             aes(x = Timepoint, y = mean_ratio, color = group),
             size = 3) +
  geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 3.5, color = "#888888", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = GROUP_COLORS, breaks = GROUP_ORDER, name = "TMAO Group") +
  scale_fill_manual( values = GROUP_COLORS, breaks = GROUP_ORDER, name = "TMAO Group") +
  scale_x_continuous(breaks = 1:6, labels = tp_labels) +
  facet_wrap(~group, ncol = 4) +
  labs(
    title    = sprintf("CA/DCA Ratio (log₂) — Individual Trajectories  |  CA ID %d, DCA ID %d",
                       CA_ID, DCA_ID),
    subtitle = "Gray lines = individual subjects  |  Colored line = group mean ± SE  |  Gray shading = X period",
    x        = "Timepoint",
    y        = "log₂(CA intensity / DCA intensity)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "#F0F0F0"),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, color = "#555555"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ── 図 2: TMAO グループ別カラー（全被験者重ね書き，ファセットなし） ────────────
p_overlay <- ggplot() +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = "#F5F5F5", alpha = 1) +
  geom_line(data = df_plot %>% filter(!is.na(group)),
            aes(x = Timepoint, y = log2_ratio, group = Subject, color = group),
            linewidth = 0.6, alpha = 0.45) +
  geom_line(data = group_summary,
            aes(x = Timepoint, y = mean_ratio, color = group),
            linewidth = 2.0) +
  geom_point(data = group_summary,
             aes(x = Timepoint, y = mean_ratio, color = group),
             size = 3.5) +
  geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = 3.5, color = "#888888", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = GROUP_COLORS, breaks = GROUP_ORDER, name = "TMAO Group") +
  scale_x_continuous(breaks = 1:6, labels = tp_labels) +
  labs(
    title    = "CA/DCA Ratio — All Subjects Overlaid by TMAO Group",
    subtitle = "Thin lines = individual subjects  |  Bold line = group mean",
    x = "Timepoint", y = "log₂(CA / DCA)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank())

# ── PDF 保存 ──────────────────────────────────────────────────────────────────
out_pdf <- file.path(OUT_FIG, "ca_dca_ratio_individual.pdf")
pdf(out_pdf, width = 14, height = 8)
print(p_all)
print(p_overlay)
dev.off()
message(sprintf("PDF 保存: %s", out_pdf))

out_png <- file.path(OUT_FIG, "ca_dca_ratio_individual.png")
ggsave(out_png, p_all, width = 14, height = 8, dpi = 150)
message(sprintf("PNG 保存: %s", out_png))

# ── CSV 保存 ──────────────────────────────────────────────────────────────────
out_csv <- file.path(OUT_REP, "ca_dca_ratio_table.csv")
write_csv(df_plot %>%
            select(Subject, Timepoint, period,
                   ca_int = ca_mean, dca_int = dca_mean,
                   log2_ca, log2_dca, ca_dca_ratio, log2_ratio, group),
          out_csv)
message(sprintf("CSV 保存: %s", out_csv))

# ── サマリー表示 ──────────────────────────────────────────────────────────────
cat("\n=== CA/DCA log2比率 グループ×時点サマリー ===\n")
group_summary %>%
  mutate(across(c(mean_ratio, sd_ratio, se_ratio), ~round(.x, 3))) %>%
  print(n = Inf)

message("\n05_ca_dca_ratio_timeseries.R 完了。")
