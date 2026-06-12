# =============================================================================
# enrich_pairwise_results.R
#
# pairwise_test_results_{transform}.csv に MS-DIAL / MS-FINDER 由来の
# メタ情報を付加して新しい CSV を出力する。
#
# Input:
#   data/production/processed/pairwise_test_results_{TRANSFORM}.csv
#     └ boxplot_pairwise_production.R の出力
#   data/production/raw/df_conf.csv
#     └ MS-DIAL 生出力（Row 1 ダミーヘッダー / Row 2 実ヘッダー）
#   data/production/processed/feature_metadata_msfinder_2090.csv  ← 任意
#     └ annotate_unknown_msfinder.py の出力（存在しない場合はスキップ）
#
# Output:
#   data/production/processed/pairwise_test_results_{TRANSFORM}_enriched.csv
#
# 付加される列:
#   INCHIKEY_final  : MS-DIAL の INCHIKEY（同定済み）、
#                     Unknown かつ MS-FINDER ヒットありの場合は MSFINDER_inchikey、
#                     それ以外は NA
#   INCHIKEY_source : "MSDIAL" | "MSFINDER" | NA
#   + EXTRA_MSDIAL_COLS に列挙した MS-DIAL 追加列（任意）
#
# 拡張方法:
#   EXTRA_MSDIAL_COLS に df_conf.csv のカラム名を追加するだけで
#   任意のメタ情報を出力 CSV に加えられる。
#
# Usage (RStudio):
#   1. TRANSFORM / パスを必要に応じて変更
#   2. source("src/r/enrich_pairwise_results.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
TRANSFORM    <- "none"   # "none" | "log2" | "log2FC"

IN_PAIRWISE  <- file.path("data/production/processed",
                           sprintf("pairwise_test_results_%s_production.csv", TRANSFORM))
IN_DFCONF    <- "data/production/raw/df_conf.csv"
IN_MSFINDER  <- "data/production/processed/feature_metadata_msfinder_2090.csv"
OUT_CSV      <- file.path("data/production/processed",
                           sprintf("pairwise_test_results_%s_enriched.csv", TRANSFORM))

# df_conf.csv から追加で取得したい MS-DIAL 列名（実ヘッダー行の列名そのまま）
# → 不要なら空ベクトル c() にする
EXTRA_MSDIAL_COLS <- c(
  "Total score",
  "S/N average",
  "RT similarity",
  "m/z similarity",
  "Matched peaks count",
  "Matched peaks percentage"
)

# -----------------------------------------------------------------------------
# 0. 入力ファイルの存在チェック
# -----------------------------------------------------------------------------
for (f in c(IN_PAIRWISE, IN_DFCONF)) {
  if (!file.exists(f)) stop(sprintf("File not found: %s", f))
}

# -----------------------------------------------------------------------------
# 1. pairwise 結果を読み込む
# -----------------------------------------------------------------------------
cat("Loading pairwise results ...\n")
pairwise <- read.csv(IN_PAIRWISE, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  %d features x %d columns\n", nrow(pairwise), ncol(pairwise)))

# -----------------------------------------------------------------------------
# 2. df_conf.csv から MS-DIAL メタ情報を取得
#    Row 1: ダミーヘッダー（スキップ）
#    Row 2: 実ヘッダー
# -----------------------------------------------------------------------------
cat("Loading df_conf.csv ...\n")
df_conf <- read.csv(IN_DFCONF, header = FALSE, stringsAsFactors = FALSE,
                    skip = 1)            # 1行スキップ → Row 2 が実ヘッダー
colnames(df_conf) <- df_conf[1, ]        # Row 2 を列名に
df_conf <- df_conf[-1, ]                 # ヘッダー行本体を除去
rownames(df_conf) <- NULL

cat(sprintf("  %d features x %d columns\n", nrow(df_conf), ncol(df_conf)))

# Alignment ID を整数に変換
df_conf[["Alignment ID"]] <- as.integer(df_conf[["Alignment ID"]])

# INCHIKEY と追加列を選択
msdial_cols_needed <- c("Alignment ID", "INCHIKEY", EXTRA_MSDIAL_COLS)
missing_cols <- setdiff(msdial_cols_needed, colnames(df_conf))
if (length(missing_cols) > 0) {
  warning(sprintf("df_conf.csv に存在しない列をスキップします: %s",
                  paste(missing_cols, collapse = ", ")))
  msdial_cols_needed <- intersect(msdial_cols_needed, colnames(df_conf))
}

df_msdial <- df_conf[, msdial_cols_needed, drop = FALSE]
colnames(df_msdial)[colnames(df_msdial) == "Alignment ID"] <- "Alignment_ID"

# 数値列を変換
for (col in EXTRA_MSDIAL_COLS) {
  if (col %in% colnames(df_msdial)) {
    df_msdial[[col]] <- suppressWarnings(as.numeric(df_msdial[[col]]))
  }
}

# -----------------------------------------------------------------------------
# 3. MS-FINDER メタ情報を取得（rank-1 候補のみ）
# -----------------------------------------------------------------------------
if (file.exists(IN_MSFINDER)) {
  cat("Loading MS-FINDER annotations ...\n")
  msfinder <- read.csv(IN_MSFINDER, stringsAsFactors = FALSE)
  msfinder_rank1 <- msfinder[!is.na(msfinder$MSFINDER_rank) &
                                msfinder$MSFINDER_rank == 1L,
                              c("Alignment_ID", "MSFINDER_inchikey")]
  cat(sprintf("  %d features with MS-FINDER rank-1 candidate\n",
              nrow(msfinder_rank1)))
} else {
  cat("MS-FINDER file not found — MSFINDER_inchikey will be NA.\n")
  msfinder_rank1 <- data.frame(Alignment_ID   = integer(0),
                                MSFINDER_inchikey = character(0),
                                stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 4. pairwise 結果に MS-DIAL / MS-FINDER メタ情報をマージ
# -----------------------------------------------------------------------------
result <- merge(pairwise,  df_msdial,      by = "Alignment_ID", all.x = TRUE)
result <- merge(result,    msfinder_rank1, by = "Alignment_ID", all.x = TRUE)

# -----------------------------------------------------------------------------
# 5. INCHIKEY_final / INCHIKEY_source を決定
#
#   優先順位:
#     1. MS-DIAL INCHIKEY が空でない → MSDIAL
#     2. Unknown かつ MSFINDER_inchikey が空でない → MSFINDER
#     3. それ以外 → NA
# -----------------------------------------------------------------------------
msdial_ik   <- result[["INCHIKEY"]]
msfinder_ik <- result[["MSFINDER_inchikey"]]

has_msdial   <- !is.na(msdial_ik)   & nchar(trimws(msdial_ik))   > 0
has_msfinder <- !is.na(msfinder_ik) & nchar(trimws(msfinder_ik)) > 0
is_unknown   <- result[["Metabolite_name"]] == "Unknown"

result[["INCHIKEY_final"]] <- ifelse(
  has_msdial,
    msdial_ik,
  ifelse(is_unknown & has_msfinder,
    msfinder_ik,
    NA_character_
  )
)

result[["INCHIKEY_source"]] <- ifelse(
  has_msdial,
    "MSDIAL",
  ifelse(is_unknown & has_msfinder,
    "MSFINDER",
    NA_character_
  )
)

# -----------------------------------------------------------------------------
# 6. 列順を整理（pairwise 元列 → INCHIKEY 列 → MS-DIAL 追加列 → MS-FINDER 列）
# -----------------------------------------------------------------------------
pairwise_orig_cols   <- colnames(pairwise)
inchikey_cols        <- c("INCHIKEY", "INCHIKEY_final", "INCHIKEY_source")
extra_msdial_present <- intersect(EXTRA_MSDIAL_COLS, colnames(result))
msfinder_extra_cols  <- c("MSFINDER_inchikey")

col_order <- unique(c(
  pairwise_orig_cols,
  inchikey_cols,
  extra_msdial_present,
  msfinder_extra_cols
))
col_order <- intersect(col_order, colnames(result))

result <- result[, col_order, drop = FALSE]
result <- result[order(result[["Alignment_ID"]]), ]

# -----------------------------------------------------------------------------
# 7. サマリー表示
# -----------------------------------------------------------------------------
cat(sprintf("\n--- INCHIKEY 付与サマリー ---\n"))
tbl <- table(result[["INCHIKEY_source"]], useNA = "ifany")
for (src in names(tbl)) {
  label <- if (is.na(src)) "NA (なし)" else src
  cat(sprintf("  %-12s : %d features\n", label, tbl[[src]]))
}
cat(sprintf("  合計        : %d features\n", nrow(result)))

# -----------------------------------------------------------------------------
# 8. 保存
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(result, OUT_CSV, row.names = FALSE)
cat(sprintf("\nSaved → %s\n", OUT_CSV))
cat(sprintf("Output: %d rows x %d columns\n", nrow(result), ncol(result)))
