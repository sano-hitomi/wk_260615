# ==============================================================================
# 08_dca_adduct_scan.R
# Purpose : feature_metadata 中の DCA 全エントリを確認
#           （異なるアダクト・MS モードによる重複エントリの有無を検証）
#
#   過去解析では DCA Alignment_ID = 10813 を使用してきた。
#   MS-DIAL では同一化合物が複数のアダクト（[M-H]⁻, [M+HCOO]⁻, [M+Cl]⁻ 等）
#   または正負イオンモードで別々に検出されることがある。
#   本スクリプトは全候補を洗い出し，最適なエントリを選定するための
#   定量的根拠（S/N, マッチスコア，強度 CV）を提供する。
#
# 出力:
#   output/reports/dca_adduct_scan.csv      -- 全 DCA エントリの詳細テーブル
#   output/figures/dca_adduct_intensities.pdf / .png
#   output/figures/dca_adduct_correlation.pdf / .png
#
# 依存: install.packages(c("tidyverse", "patchwork"))
# ==============================================================================

library(tidyverse)
library(patchwork)

# ── パス設定 ──────────────────────────────────────────────────────────────────
PROJ_ROOT  <- "."
FMETA_CSV  <- file.path(PROJ_ROOT, "data/production/processed/feature_metadata.csv")
FMAT_CSV   <- file.path(PROJ_ROOT, "data/production/processed/feature_matrix.csv")
DFCONF_CSV <- file.path(PROJ_ROOT, "data/production/raw/df_conf.csv")
GROUPS_CSV <- file.path(PROJ_ROOT, "data/ref/tmao_subject_groups.csv")
OUT_FIG    <- file.path(PROJ_ROOT, "output/figures")
OUT_REP    <- file.path(PROJ_ROOT, "output/reports")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_REP, showWarnings = FALSE, recursive = TRUE)

# ── DCA 検索パラメータ ─────────────────────────────────────────────────────────
PREFERRED_DCA_ID <- 10813L
IK_DCA_PREFIX    <- "KXGVEGMKQFWNSR"   # DCA InChIKey 前半 14 文字
HMDB_DCA         <- "HMDB0000626"
# DCA の正確な分子量: 392.2927 Da（C24H40O4）
DCA_EXACT_MZ     <- 392.2927
MZ_TOL_PPM       <- 10.0              # m/z 許容誤差 [ppm]

# ── feature_metadata 読み込み ─────────────────────────────────────────────────
if (!file.exists(FMETA_CSV)) stop(sprintf("feature_metadata が見つかりません: %s", FMETA_CSV))
feat_meta <- read_csv(FMETA_CSV, show_col_types = FALSE)
cat(sprintf("feature_metadata: %d features × %d 列\n", nrow(feat_meta), ncol(feat_meta)))
cat("列名:", paste(colnames(feat_meta), collapse = ", "), "\n\n")

# ── df_conf.csv（MS-DIAL 生出力）が利用可能か確認 ─────────────────────────────
has_dfconf <- file.exists(DFCONF_CSV)
if (has_dfconf) {
  cat("df_conf.csv を読み込み（アダクト情報取得）\n")
  # MS-DIAL 出力は Row1=ダミー, Row2=実ヘッダー
  raw <- read.csv(DFCONF_CSV, header = FALSE, stringsAsFactors = FALSE,
                  skip = 1, check.names = FALSE)
  colnames(raw) <- raw[1, ]
  raw <- raw[-1, ]

  # 列名の正規化
  colnames(raw)[colnames(raw) == "Alignment ID"]  <- "Alignment_ID"
  colnames(raw)[colnames(raw) == "Metabolite name"] <- "Metabolite_name"
  colnames(raw)[colnames(raw) == "Adduct type"]   <- "Adduct"
  colnames(raw)[colnames(raw) == "Average Mz"]    <- "Average_Mz"
  raw$Alignment_ID <- as.integer(raw$Alignment_ID)

  df_adduct <- raw %>%
    dplyr::select(Alignment_ID,
           Metabolite_name = any_of(c("Metabolite_name", "Metabolite name")),
           Adduct          = any_of(c("Adduct", "Adduct type")),
           Average_Mz      = any_of(c("Average_Mz", "Average Mz", "Precursor m/z")),
           INCHIKEY         = any_of(c("INCHIKEY", "InChIKey")),
           everything()) %>%
    mutate(Average_Mz = as.numeric(Average_Mz))

  cat(sprintf("df_conf.csv: %d features\n", nrow(df_adduct)))
} else {
  cat("df_conf.csv が見つかりません。feature_metadata のみで解析します。\n")
  df_adduct <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. DCA エントリの検索
# ══════════════════════════════════════════════════════════════════════════════

dca_ids_all <- integer(0)

# (a) known_id
if (PREFERRED_DCA_ID %in% feat_meta$Alignment_ID) {
  dca_ids_all <- c(dca_ids_all, PREFERRED_DCA_ID)
  cat(sprintf("(a) known ID %d: 存在\n", PREFERRED_DCA_ID))
}

# (b) InChIKey
ik_col <- grep("inchikey|INCHIKEY", colnames(feat_meta), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(ik_col)) {
  hits <- feat_meta$Alignment_ID[
    !is.na(feat_meta[[ik_col]]) &
      startsWith(toupper(feat_meta[[ik_col]]), IK_DCA_PREFIX)
  ]
  cat(sprintf("(b) InChIKey (%s...) ヒット: %d 件 → %s\n",
              IK_DCA_PREFIX, length(hits),
              if (length(hits) > 0) paste(hits, collapse = ", ") else "なし"))
  dca_ids_all <- union(dca_ids_all, hits)
}

# (c) HMDB
hmdb_col <- grep("hmdb", colnames(feat_meta), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(hmdb_col)) {
  hits <- feat_meta$Alignment_ID[
    !is.na(feat_meta[[hmdb_col]]) &
      trimws(feat_meta[[hmdb_col]]) %in% c(HMDB_DCA,
                                             sub("HMDB0+", "HMDB", HMDB_DCA))
  ]
  cat(sprintf("(c) HMDB (%s) ヒット: %d 件 → %s\n",
              HMDB_DCA, length(hits),
              if (length(hits) > 0) paste(hits, collapse = ", ") else "なし"))
  dca_ids_all <- union(dca_ids_all, hits)
}

# (d) 名前パターン
name_col <- grep("^(Metabolite.name|Metabolite_name|Name)", colnames(feat_meta),
                 ignore.case = TRUE, value = TRUE)[1]
if (!is.na(name_col)) {
  hits <- feat_meta$Alignment_ID[
    !is.na(feat_meta[[name_col]]) &
      grepl("deoxychol", feat_meta[[name_col]], ignore.case = TRUE) &
      !grepl("chenodeoxy|glycodeoxy|taurodeoxy|hyodeoxy", feat_meta[[name_col]], ignore.case = TRUE)
  ]
  cat(sprintf("(d) 名前パターン 'deoxychol' ヒット: %d 件 → %s\n",
              length(hits),
              if (length(hits) > 0) paste(hits, collapse = ", ") else "なし"))
  dca_ids_all <- union(dca_ids_all, hits)
}

# (e) m/z 検索（df_conf が利用可能な場合）
if (!is.null(df_adduct) && "Average_Mz" %in% colnames(df_adduct)) {
  # DCA の主要アダクト m/z:
  #   [M-H]⁻   : 391.2849  (392.2927 - 1.0073)
  #   [M+HCOO]⁻: 437.2955  (392.2927 + 44.9977 + 1.0073 - 1.0073 ... use formate)
  #   [M+Cl]⁻  : 427.2588
  #   [M-H2O-H]⁻: 373.2743
  dca_mz_candidates <- c(391.2849, 437.2955, 427.2588, 373.2743)
  for (mz in dca_mz_candidates) {
    ppm_err <- abs(df_adduct$Average_Mz - mz) / mz * 1e6
    hits <- df_adduct$Alignment_ID[!is.na(ppm_err) & ppm_err <= MZ_TOL_PPM]
    if (length(hits) > 0) {
      cat(sprintf("(e) m/z %.4f ± %g ppm ヒット: %d 件 → %s\n",
                  mz, MZ_TOL_PPM, length(hits), paste(hits, collapse = ", ")))
      dca_ids_all <- union(dca_ids_all, hits)
    }
  }
}

dca_ids_all <- sort(unique(dca_ids_all))
cat(sprintf("\n合計 DCA 候補エントリ数: %d\n", length(dca_ids_all)))
cat("IDs:", paste(dca_ids_all, collapse = ", "), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. 候補エントリの詳細テーブル作成
# ══════════════════════════════════════════════════════════════════════════════

# feature_metadata から基本情報
meta_sub <- feat_meta %>%
  filter(Alignment_ID %in% dca_ids_all)

# df_conf があればアダクト / m/z を取得してマージ
if (!is.null(df_adduct)) {
  adduct_sub <- df_adduct %>%
    filter(Alignment_ID %in% dca_ids_all) %>%
    dplyr::select(Alignment_ID,
           Adduct = any_of(c("Adduct", "Adduct type")),
           Average_Mz,
           `S/N average`   = any_of(c("S/N average", "SN_average")),
           `Total score`   = any_of(c("Total score", "Total_score")),
           `Matched peaks percentage` = any_of("Matched peaks percentage"))
  meta_sub <- meta_sub %>% left_join(adduct_sub, by = "Alignment_ID")
}

# feature_matrix から強度統計（生物学的サンプル）
if (file.exists(FMAT_CSV)) {
  fmat <- read_csv(FMAT_CSV, show_col_types = FALSE) %>%
    rename(Alignment_ID = 1) %>%
    mutate(Alignment_ID = as.integer(Alignment_ID)) %>%
    filter(Alignment_ID %in% dca_ids_all)

  # 強度の CV（変動係数）と平均を計算
  intensity_stats <- fmat %>%
    pivot_longer(-Alignment_ID, names_to = "label", values_to = "intensity") %>%
    mutate(intensity = as.numeric(intensity)) %>%
    filter(!is.na(intensity), intensity > 0) %>%
    group_by(Alignment_ID) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      median_intensity = median(intensity, na.rm = TRUE),
      cv_pct         = sd(intensity, na.rm = TRUE) / mean(intensity, na.rm = TRUE) * 100,
      n_detected     = sum(intensity > 0, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))

  meta_sub <- meta_sub %>% left_join(intensity_stats, by = "Alignment_ID")
}

# preferred_id フラグ
meta_sub <- meta_sub %>%
  mutate(is_current = Alignment_ID == PREFERRED_DCA_ID)

cat("=== DCA 全エントリ詳細 ===\n")
print(meta_sub %>% dplyr::select(Alignment_ID, is_current, any_of(
  c("Metabolite_name", "Adduct", "Average_Mz", "S/N average",
    "Total score", "Matched peaks percentage",
    "mean_intensity", "cv_pct", "n_detected")
)))

# CSV 保存
out_csv <- file.path(OUT_REP, "dca_adduct_scan.csv")
write_csv(meta_sub, out_csv)
message(sprintf("CSV 保存: %s", out_csv))

# ══════════════════════════════════════════════════════════════════════════════
# 3. 複数エントリがある場合: 強度の相関・比較図
# ══════════════════════════════════════════════════════════════════════════════
if (length(dca_ids_all) < 2 || !file.exists(FMAT_CSV)) {
  message(sprintf(
    "DCA エントリが %d 件（または feature_matrix なし）のため比較図をスキップします。",
    length(dca_ids_all)
  ))
} else {
  fmat_dca <- read_csv(FMAT_CSV, show_col_types = FALSE) %>%
    rename(Alignment_ID = 1) %>%
    mutate(Alignment_ID = as.integer(Alignment_ID)) %>%
    filter(Alignment_ID %in% dca_ids_all) %>%
    pivot_longer(-Alignment_ID, names_to = "label", values_to = "intensity") %>%
    mutate(intensity = log2(as.numeric(intensity) + 1),
           ID_label  = sprintf("ID %d", Alignment_ID))

  # ── 図 A: 各エントリの強度分布（バイオリン） ────────────────────────────────
  p_violin <- ggplot(fmat_dca %>% filter(is.finite(intensity) & intensity > 0),
                     aes(x = ID_label, y = intensity, fill = ID_label)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = "DCA Candidate Entries — Intensity Distribution",
      subtitle = sprintf("n_entries = %d  |  values = log₂(intensity + 1)  |  all samples",
                         length(dca_ids_all)),
      x = "Alignment ID", y = "log₂(intensity + 1)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          panel.grid.minor = element_blank())

  # ── 図 B: エントリ間の相関（全サンプルの散布図） ───────────────────────────
  dca_wide <- fmat_dca %>%
    pivot_wider(names_from = ID_label, values_from = intensity)

  id_cols <- colnames(dca_wide)[colnames(dca_wide) != "label"]
  n_ids   <- length(id_cols)

  scatter_plots <- list()
  if (n_ids >= 2) {
    for (i in 1:(n_ids - 1)) {
      for (j in (i + 1):n_ids) {
        ci <- id_cols[i]; cj <- id_cols[j]
        df_sc <- dca_wide %>%
          dplyr::select(label, xi = all_of(ci), xj = all_of(cj)) %>%
          filter(is.finite(xi) & is.finite(xj) & xi > 0 & xj > 0)

        ct <- tryCatch(cor.test(df_sc$xi, df_sc$xj), error = function(e) NULL)
        r_lab <- if (!is.null(ct))
          sprintf("r = %.3f\np = %.4f", ct$estimate, ct$p.value) else "n.a."

        p <- ggplot(df_sc, aes(x = xi, y = xj)) +
          geom_point(alpha = 0.4, size = 1.2, color = "#2E75B6") +
          geom_smooth(method = "lm", se = FALSE, color = "#C00000", linewidth = 0.8) +
          annotate("label", x = Inf, y = -Inf,
                   label = r_lab, hjust = 1.05, vjust = -0.1,
                   size = 3, fontface = "bold", fill = "white", label.size = 0.3) +
          labs(title = sprintf("%s vs %s", ci, cj),
               x = ci, y = cj) +
          theme_bw(base_size = 9) +
          theme(plot.title = element_text(size = 9, face = "bold"),
                panel.grid.minor = element_blank())
        scatter_plots[[sprintf("%s_vs_%s", ci, cj)]] <- p
      }
    }
  }

  # ── 図 C: 推奨エントリの選定根拠サマリー ─────────────────────────────────
  if (!is.null(meta_sub$`S/N average`) | !is.null(meta_sub$mean_intensity)) {
    ranking_cols <- intersect(
      c("Alignment_ID", "is_current", "Adduct", "Average_Mz",
        "S/N average", "Total score", "mean_intensity", "cv_pct", "n_detected"),
      colnames(meta_sub)
    )
    ranking_df <- meta_sub %>%
      dplyr::select(all_of(ranking_cols)) %>%
      arrange(desc(.data[[intersect(c("S/N average", "mean_intensity"),
                                     colnames(.))[1]]]))

    cat("\n=== DCA エントリ推奨ランキング ===\n")
    print(ranking_df)

    # テキストアノテーション用
    best_id <- ranking_df$Alignment_ID[1]
    msg <- sprintf(
      "推奨 ID: %d\n（現在使用 ID: %d  %s）",
      best_id, PREFERRED_DCA_ID,
      if (best_id == PREFERRED_DCA_ID) "✓ 一致" else "← 変更検討"
    )
    cat("\n", msg, "\n")
  }

  # ── 保存 ──────────────────────────────────────────────────────────────────
  out_pdf <- file.path(OUT_FIG, "dca_adduct_intensities.pdf")
  pdf(out_pdf, width = 10, height = 5)
  print(p_violin)
  if (length(scatter_plots) > 0) {
    layout_n <- min(length(scatter_plots), 4)
    print(wrap_plots(scatter_plots[1:layout_n], ncol = min(2, layout_n)))
  }
  dev.off()
  message(sprintf("PDF 保存: %s", out_pdf))

  out_png <- file.path(OUT_FIG, "dca_adduct_intensities.png")
  ggsave(out_png, p_violin, width = max(6, n_ids * 2), height = 5, dpi = 150)
  message(sprintf("PNG 保存: %s", out_png))
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. 最終サマリー出力
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("========================================================\n")
cat(sprintf("  DCA スキャン完了: %d エントリが見つかりました\n", length(dca_ids_all)))
cat("========================================================\n")
if (length(dca_ids_all) == 0) {
  cat("  → DCA エントリが見つかりません。\n")
  cat("    InChIKey / HMDB / 名前の列名を feat_meta で確認してください。\n")
} else if (length(dca_ids_all) == 1) {
  cat(sprintf("  → ID %d のみ。アダクト重複なし。\n", dca_ids_all))
  cat(sprintf("    現在使用の ID %d で問題ありません。\n", PREFERRED_DCA_ID))
} else {
  cat(sprintf("  → %d エントリ検出: %s\n",
              length(dca_ids_all), paste(dca_ids_all, collapse = ", ")))
  cat(sprintf("    詳細: %s\n", out_csv))
  cat("    S/N 最高・CV 最小のエントリを採用することを推奨します。\n")
}

message("\n08_dca_adduct_scan.R 完了。")
