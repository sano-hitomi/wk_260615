# =============================================================================
# classify_metabolites_metaboanalystr.R
#
# PLS-DA の VIP スコアとコンポーネント 1 ローディングに基づいて代謝物を
# 5 クラスに分類し、MetaboAnalystR へ流す CSV 群を出力する。
#
# 分類ロジック（Comp1 のみ使用）:
#   sig_early = VIP_early >= VIP_THRESH
#   sig_late  = VIP_late  >= VIP_THRESH
#   same_dir  = sign(loading_early) == sign(loading_late)
#
#   early_specific    : sig_early  AND !sig_late
#   late_specific     : !sig_early AND  sig_late
#   shared_concordant : sig_early  AND  sig_late  AND  same_dir
#                         → Early / Late 共通の機序（同方向）
#   shared_discordant : sig_early  AND  sig_late  AND !same_dir
#                         → 両期間で重要だが方向が逆（異なる機序）
#   ns                : !sig_early AND !sig_late
#
# アノテーション優先度:
#   1. MS-DIAL 同定済み (Metabolite_name != "Unknown")
#   2. MS-FINDER 構造提案 (MSFINDER_structure あり)
#   3. Unknown（除外 or 保留）
#
# Input:
#   plsda_vip_loadings.csv          (plsda_loadings_vip.R の出力)
#   feature_metadata_msfinder_2090.csv
#   pairwise_test_results_{TRANSFORM}.csv  (p 値付加 / 任意)
#
# Output (OUT_DIR に保存):
#   metabolite_classified.csv
#     → 全代謝物の VIP・ローディング・分類・アノテーション統合表
#   metaboanalystr_ranked_early.csv / _late.csv
#     → GSEA-style ranked list: Name, Score (signed VIP), InChIKey, Ontology
#   metaboanalystr_ora_early_specific.csv    など各クラスの ORA 用名称リスト
#   metaboanalystr_background.csv
#     → アノテーション付き全代謝物の background リスト
#
# Usage (RStudio):
#   # 1. Working directory をプロジェクトルートに設定
#   # 2. source("src/r/classify_metabolites_metaboanalystr.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Settings（必要に応じて変更）
# -----------------------------------------------------------------------------

# "production" or "dummy"
MODE <- "production"

# 本番データが揃っていない場合はダミーで検証
if (!dir.exists(file.path("data", MODE, "processed"))) {
  message(sprintf("'data/%s/processed' が見つかりません。dummy で代替します。", MODE))
  MODE <- "dummy"
}

DATA_DIR <- file.path("data", MODE, "processed")
OUT_DIR  <- DATA_DIR

# VIP 閾値（探索的検討のため甘めに設定）
VIP_THRESH_EARLY <- 0.8
VIP_THRESH_LATE  <- 0.8

# pairwise 結果を合わせて使う場合のトランスフォーム選択
PAIRWISE_TRANSFORM <- "log2FC"   # "none" | "log2" | "log2FC"

# アノテーションなし代謝物を出力に含めるか
#   TRUE  → ns も含め全代謝物を metabolite_classified.csv に出力
#   FALSE → 何らかのアノテーションがある代謝物のみ出力
INCLUDE_UNANNOTATED <- FALSE

# -----------------------------------------------------------------------------
# 0. パッケージ
# -----------------------------------------------------------------------------
pkgs <- c("dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# -----------------------------------------------------------------------------
# 1. VIP / ローディング CSV の読み込み
# -----------------------------------------------------------------------------
vip_csv <- file.path(DATA_DIR, "plsda_vip_loadings.csv")

if (!file.exists(vip_csv)) {
  stop(
    sprintf("plsda_vip_loadings.csv が見つかりません: %s\n", vip_csv),
    "plsda_loadings_vip.R を先に実行してください。"
  )
}

vip_df <- read.csv(vip_csv, stringsAsFactors = FALSE)
cat(sprintf("[1] VIP/loading テーブル読み込み: %d features\n", nrow(vip_df)))

# -----------------------------------------------------------------------------
# 2. アノテーション情報の読み込みとマージ
#    優先列: feat_meta の Metabolite_name / Ontology / INCHIKEY
#            + MS-FINDER rank-1 の MSFINDER_structure / MSFINDER_ontology /
#              MSFINDER_inchikey
# -----------------------------------------------------------------------------
msfinder_csv <- file.path(DATA_DIR, "feature_metadata_msfinder_2090.csv")

# --- 2a. featmeta（MS-DIAL アノテーション）
# plsda_loadings_vip.R はすでに feat_meta から Metabolite_name を持ってくる。
# ここでは Ontology, INCHIKEY, Formula を補完するためにあらためて読む。
# load_data.R / load_data_production.R が source 済みなら feat_meta を使う。
if (exists("feat_meta")) {
  meta_base <- feat_meta[, intersect(
    c("Alignment_ID", "Metabolite_name", "Ontology", "INCHIKEY",
      "Formula", "SMILES", "Adduct_type",
      "MSFINDER_structure", "MSFINDER_total_score",
      "MSFINDER_ontology",  "MSFINDER_inchikey",
      "MSFINDER_smiles",    "MSFINDER_annotated"),
    colnames(feat_meta)
  )]
  cat("[2] feat_meta をメモリから使用\n")
} else if (file.exists(msfinder_csv)) {
  msf <- read.csv(msfinder_csv, stringsAsFactors = FALSE)
  # rank-1 のみ
  msf_r1 <- msf[!is.na(msf$MSFINDER_rank) & msf$MSFINDER_rank == 1L, ]
  meta_base <- msf_r1[, intersect(
    c("Alignment_ID", "Metabolite_name", "Ontology", "INCHIKEY",
      "Formula", "SMILES", "Adduct_type",
      "MSFINDER_structure", "MSFINDER_total_score",
      "MSFINDER_ontology",  "MSFINDER_inchikey",
      "MSFINDER_smiles"),
    colnames(msf_r1)
  )]
  meta_base$MSFINDER_annotated <- !is.na(meta_base$MSFINDER_structure) &
                                    meta_base$MSFINDER_structure != ""
  cat(sprintf("[2] feature_metadata_msfinder_2090.csv から %d features 読み込み\n",
              nrow(meta_base)))
} else {
  stop("アノテーション情報が見つかりません。\n",
       "  load_data.R / load_data_production.R を先に source するか、\n",
       sprintf("  %s を確認してください。", msfinder_csv))
}

# Alignment_ID を整数に統一
meta_base$Alignment_ID <- as.integer(meta_base$Alignment_ID)
vip_df$Alignment_ID    <- as.integer(vip_df$Alignment_ID)

df <- left_join(vip_df, meta_base, by = "Alignment_ID")
cat(sprintf("[2] マージ後: %d features\n", nrow(df)))

# -----------------------------------------------------------------------------
# 2b. pairwise p 値の付加（任意）
# -----------------------------------------------------------------------------
# "dummy" モードは "_dummy" サフィックス、"production" モードは "_production"
suffix_map  <- c(dummy = "_dummy", production = "_production")
file_suffix <- suffix_map[[MODE]]

pairwise_csv <- file.path(
  DATA_DIR,
  sprintf("pairwise_test_results_%s%s.csv", PAIRWISE_TRANSFORM, file_suffix)
)

if (file.exists(pairwise_csv)) {
  pw <- read.csv(pairwise_csv, stringsAsFactors = FALSE,
                 check.names = FALSE)
  pw$Alignment_ID <- as.integer(pw$Alignment_ID)
  # 取得列: p 値・padj のみ（重複を避けるため Metabolite_name などは除外）
  pw_cols <- c("Alignment_ID",
               grep("^(T[0-9]v|padj_)", colnames(pw), value = TRUE))
  df <- left_join(df, pw[, pw_cols], by = "Alignment_ID")
  cat(sprintf("[2b] pairwise 結果 (%s) を付加\n", PAIRWISE_TRANSFORM))
} else {
  cat(sprintf("[2b] pairwise CSV なし（スキップ）: %s\n", pairwise_csv))
}

# -----------------------------------------------------------------------------
# 3. "best_name" と "best_inchikey" / "best_ontology" の決定
#    優先度: MS-DIAL 同定 > MS-FINDER 提案 > "Unknown"
# -----------------------------------------------------------------------------
is_msdial_id <- !is.na(df$Metabolite_name) &
                  df$Metabolite_name != "Unknown" &
                  !grepl("^low score:", df$Metabolite_name)

is_msf_id    <- !is_msdial_id &
                  !is.na(df$MSFINDER_structure) & df$MSFINDER_structure != ""

df$best_name <- ifelse(
  is_msdial_id, df$Metabolite_name,
  ifelse(is_msf_id, df$MSFINDER_structure, "Unknown")
)

df$best_inchikey <- ifelse(
  is_msdial_id & !is.na(df$INCHIKEY) & df$INCHIKEY != "",
  df$INCHIKEY,
  ifelse(is_msf_id, df$MSFINDER_inchikey, NA_character_)
)

df$best_ontology <- ifelse(
  is_msdial_id & !is.na(df$Ontology) & df$Ontology != "",
  df$Ontology,
  ifelse(is_msf_id, df$MSFINDER_ontology, NA_character_)
)

df$annotation_source <- ifelse(
  is_msdial_id, "MSDIAL",
  ifelse(is_msf_id, "MSFINDER", "Unknown")
)

cat(sprintf("[3] アノテーション内訳:\n"))
cat(sprintf("    MS-DIAL 同定 : %d\n",   sum(is_msdial_id)))
cat(sprintf("    MS-FINDER   : %d\n",   sum(is_msf_id)))
cat(sprintf("    Unknown     : %d\n",   sum(!is_msdial_id & !is_msf_id)))

# -----------------------------------------------------------------------------
# 4. VIP / ローディングによる分類
# -----------------------------------------------------------------------------
sig_early <- df$VIP_early >= VIP_THRESH_EARLY
sig_late  <- df$VIP_late  >= VIP_THRESH_LATE
same_dir  <- sign(df$loading_early) == sign(df$loading_late)

df$vip_class <- ifelse(
  sig_early & !sig_late,  "early_specific",
  ifelse(!sig_early & sig_late, "late_specific",
  ifelse(sig_early & sig_late & same_dir,  "shared_concordant",
  ifelse(sig_early & sig_late & !same_dir, "shared_discordant",
         "ns"))))

cat(sprintf("\n[4] 分類結果（VIP 閾値 early=%.2f, late=%.2f）:\n",
            VIP_THRESH_EARLY, VIP_THRESH_LATE))
print(table(df$vip_class))

# アノテーション × 分類のクロス集計
cat("\n    アノテーション × 分類:\n")
print(table(df$annotation_source, df$vip_class))

# -----------------------------------------------------------------------------
# 5. Signed VIP スコアの計算（GSEA スタイルのランク付け用）
#    loading の符号 × VIP を signed score とする。
#    positive = group-A（X_only / both）と正の相関
#    negative = group-B（Y_only / other）と正の相関
#    ※符号の解釈は plsda_loadings_vip.R の sign alignment に依存
# -----------------------------------------------------------------------------
df$signed_vip_early <- sign(df$loading_early) * df$VIP_early
df$signed_vip_late  <- sign(df$loading_late)  * df$VIP_late

# -----------------------------------------------------------------------------
# 6. 出力用データフレームの構築
# -----------------------------------------------------------------------------
# 全分類（アノテーションあり / なし両方）
all_cols <- c(
  "Alignment_ID", "best_name", "best_inchikey", "best_ontology",
  "annotation_source",
  "VIP_early", "VIP_late", "loading_early", "loading_late",
  "signed_vip_early", "signed_vip_late",
  "vip_class",
  # MS-DIAL 由来
  "Metabolite_name", "INCHIKEY", "Ontology", "Formula", "Adduct_type",
  # MS-FINDER 由来
  "MSFINDER_structure", "MSFINDER_total_score",
  "MSFINDER_ontology",  "MSFINDER_inchikey",
  # pairwise p 値（存在する列のみ）
  grep("^(T[0-9]v|padj_)", colnames(df), value = TRUE)
)
all_cols <- intersect(all_cols, colnames(df))
out_all  <- df[, all_cols]
out_all  <- out_all[order(-pmax(out_all$VIP_early, out_all$VIP_late)), ]

# アノテーションありのみ
out_annotated <- out_all[out_all$annotation_source != "Unknown", ]

# -----------------------------------------------------------------------------
# 7. CSV 出力
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 7-1. 全代謝物統合表
if (INCLUDE_UNANNOTATED) {
  save_df <- out_all
} else {
  save_df <- out_annotated
}
out_path_all <- file.path(OUT_DIR, "metabolite_classified.csv")
write.csv(save_df, out_path_all, row.names = FALSE)
cat(sprintf("\n[7] metabolite_classified.csv 保存: %d rows → %s\n",
            nrow(save_df), out_path_all))

# 7-2. GSEA スタイル ranked list（アノテーションあり代謝物のみ）
make_ranked <- function(df_in, score_col, period_label) {
  d <- df_in[df_in$annotation_source != "Unknown",
             c("best_name", score_col, "best_inchikey", "best_ontology",
               "Alignment_ID", "Adduct_type")]
  colnames(d)[colnames(d) == score_col] <- "Score"
  d <- d[order(-d$Score), ]
  fname <- file.path(OUT_DIR,
                     sprintf("metaboanalystr_ranked_%s.csv", period_label))
  write.csv(d, fname, row.names = FALSE)
  cat(sprintf("[7] metaboanalystr_ranked_%s.csv 保存: %d rows\n",
              period_label, nrow(d)))
}
make_ranked(out_all, "signed_vip_early", "early")
make_ranked(out_all, "signed_vip_late",  "late")

# 7-3. ORA 用クラス別リスト（アノテーションあり代謝物のみ）
ora_classes <- c("early_specific", "late_specific",
                 "shared_concordant", "shared_discordant")

for (cls in ora_classes) {
  d_cls <- out_annotated[out_annotated$vip_class == cls, ]
  # MetaboAnalystR ORA 用: 1 列目が化合物名、2 列目以降は参照用
  ora_out <- data.frame(
    Name         = d_cls$best_name,
    InChIKey     = d_cls$best_inchikey,
    Ontology     = d_cls$best_ontology,
    VIP_early    = d_cls$VIP_early,
    VIP_late     = d_cls$VIP_late,
    loading_early = d_cls$loading_early,
    loading_late  = d_cls$loading_late,
    Alignment_ID = d_cls$Alignment_ID,
    stringsAsFactors = FALSE
  )
  fname <- file.path(OUT_DIR,
                     sprintf("metaboanalystr_ora_%s.csv", cls))
  write.csv(ora_out, fname, row.names = FALSE)
  cat(sprintf("[7] metaboanalystr_ora_%s.csv 保存: %d features\n",
              cls, nrow(ora_out)))
}

# 7-4. ORA background（アノテーションあり全代謝物）
bg_out <- data.frame(
  Name         = out_annotated$best_name,
  InChIKey     = out_annotated$best_inchikey,
  Ontology     = out_annotated$best_ontology,
  Alignment_ID = out_annotated$Alignment_ID,
  stringsAsFactors = FALSE
)
bg_path <- file.path(OUT_DIR, "metaboanalystr_background.csv")
write.csv(bg_out, bg_path, row.names = FALSE)
cat(sprintf("[7] metaboanalystr_background.csv 保存: %d features\n", nrow(bg_out)))

# -----------------------------------------------------------------------------
# 8. サマリー出力
# -----------------------------------------------------------------------------
cat(sprintf("\n%s\n", strrep("=", 60)))
cat("  分類サマリー（アノテーションあり代謝物のみ）\n")
cat(sprintf("%s\n", strrep("=", 60)))
cat(sprintf("  VIP 閾値: early >= %.2f, late >= %.2f\n",
            VIP_THRESH_EARLY, VIP_THRESH_LATE))
cat(sprintf("  アノテーション代謝物数: %d / %d\n",
            nrow(out_annotated), nrow(out_all)))
cat("\n")

tbl <- table(out_annotated$vip_class)
for (cls in c("early_specific", "late_specific",
              "shared_concordant", "shared_discordant", "ns")) {
  n <- if (cls %in% names(tbl)) tbl[[cls]] else 0
  cat(sprintf("  %-22s : %d features\n", cls, n))
}

cat(sprintf("\n  Output directory: %s/\n", normalizePath(OUT_DIR)))
cat(sprintf("  Files:\n"))
for (f in c("metabolite_classified.csv",
            "metaboanalystr_ranked_early.csv",
            "metaboanalystr_ranked_late.csv",
            "metaboanalystr_ora_early_specific.csv",
            "metaboanalystr_ora_late_specific.csv",
            "metaboanalystr_ora_shared_concordant.csv",
            "metaboanalystr_ora_shared_discordant.csv",
            "metaboanalystr_background.csv")) {
  p <- file.path(OUT_DIR, f)
  status <- if (file.exists(p)) "✓" else "–"
  cat(sprintf("    %s %s\n", status, f))
}
cat(sprintf("%s\n", strrep("=", 60)))
