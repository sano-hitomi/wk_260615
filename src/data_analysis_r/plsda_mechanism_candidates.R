# =============================================================================
# plsda_mechanism_candidates.R
#
# Purpose:
#   PLS-DA 結果（VIP スコア + Comp1 ローディング）を用いて代謝物を
#   機構クラスに分類し，MetaboAnalystR 向け CSV を出力する。
#
# Classification logic:
#   3 値（VIP_early, VIP_late, loading quadrant）の組み合わせで以下を判定:
#
#   loading_quad  |  Q1: Early(+) Late(+)   Q2: Early(-) Late(+)
#                 |  Q3: Early(-) Late(-)   Q4: Early(+) Late(-)
#
#   mechanism_class:
#     Reversed_Q2  : Q2, VIP pass いずれかの期間 → 機構差の最有力候補
#     Reversed_Q4  : Q4, 同上
#     Shared_Q1    : Q1, VIP pass 両期間 → 両摂動で正方向に共通
#     Shared_Q3    : Q3, VIP pass 両期間 → 両摂動で負方向に共通
#     Early_specific: VIP pass Early のみ
#     Late_specific : VIP pass Late のみ
#     Low_VIP       : いずれの期間も VIP 閾値未満（デフォルトでは除外）
#
# Annotation filter:
#   - annotated（MS-DIAL 同定 or MS-FINDER 暫定）: VIP_THRESH_ANNOTATED（甘め）
#   - unknown（識別子なし）              : VIP_THRESH_UNKNOWN（厳しめ / 除外可）
#
# Output (MetaboAnalystR 対応):
#   ① plsda_mechanism_candidates.csv   : 全候補（詳細）
#   ② plsda_metaboanalyst_ora.csv      : ORA 用（ID + mechanism_class のみ）
#   ③ plsda_metaboanalyst_msea.csv     : MSEA 用（ID + ranking score）
#
# Prerequisites:
#   1. load_data_production.R  (samplesheet, feat_meta, feat_mat)
#   2. plsda_loadings_vip.R    (generates plsda_vip_loadings.csv)
# =============================================================================

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
VIP_CSV               <- "data/production/processed/plsda_vip_loadings.csv"
OUT_DIR               <- "data/production/processed"

# VIP thresholds
# annotated 代謝物は甘め → 候補を広く残す
VIP_THRESH_ANNOTATED  <- 0.8
# unknown（識別子なし）は standard threshold
# Inf に設定すると unknown を完全除外
VIP_THRESH_UNKNOWN    <- Inf   # Inf = exclude unknowns entirely

# ローディングが「有意な方向性を持つ」と見なす絶対値の最低ライン
# 0 にすると全特徴量を象限分類する（デフォルト推奨）
LOADING_ABS_MIN       <- 0.0

# -----------------------------------------------------------------------------
# 0. Prerequisites
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing_obj <- required[!sapply(required, exists)]
if (length(missing_obj) > 0)
  stop("Missing: ", paste(missing_obj, collapse = ", "),
       "\nRun load_data_production.R first.")

if (!file.exists(VIP_CSV))
  stop("VIP/loading CSV not found: ", VIP_CSV,
       "\nRun plsda_loadings_vip.R first.")

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
for (p in c("dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Load VIP / loading table
# -----------------------------------------------------------------------------
vip_df <- read.csv(VIP_CSV, stringsAsFactors = FALSE)
cat(sprintf("VIP table loaded: %d features\n", nrow(vip_df)))

# -----------------------------------------------------------------------------
# 3. Merge feat_meta annotation columns
# -----------------------------------------------------------------------------
meta_cols <- c("Alignment_ID", "Rt_min", "Mz", "Adduct_type",
               "Formula", "Ontology", "INCHIKEY", "SMILES")
# Add MSFINDER columns if present
msfinder_cols <- c("MSFINDER_annotated", "MSFINDER_structure",
                   "MSFINDER_formula",   "MSFINDER_ontology",
                   "MSFINDER_inchikey",  "MSFINDER_smiles",
                   "MSFINDER_total_score", "MSFINDER_databases")
meta_cols <- c(meta_cols, intersect(msfinder_cols, colnames(feat_meta)))

meta_sub <- feat_meta[, intersect(meta_cols, colnames(feat_meta)), drop = FALSE]
df <- merge(vip_df, meta_sub, by = "Alignment_ID", all.x = TRUE)

# -----------------------------------------------------------------------------
# 4. Annotation status
# -----------------------------------------------------------------------------
is_msdial_id <- !is.na(df$INCHIKEY) & nchar(trimws(df$INCHIKEY)) > 0
is_msfinder  <- if ("MSFINDER_annotated" %in% colnames(df)) {
                  !is.na(df$MSFINDER_annotated) & df$MSFINDER_annotated
                } else {
                  rep(FALSE, nrow(df))
                }
df$is_annotated <- is_msdial_id | is_msfinder

# Best available identifiers for MetaboAnalystR
df$best_name    <- ifelse(is_msfinder,
                          df$MSFINDER_structure,
                          df$Metabolite_name)
df$best_name    <- ifelse(is.na(df$best_name) | df$best_name == "",
                          paste0("Unknown_", df$Alignment_ID),
                          df$best_name)

df$best_inchikey <- ifelse(is_msdial_id,
                            trimws(df$INCHIKEY),
                            ifelse(is_msfinder & !is.na(df$MSFINDER_inchikey),
                                   df$MSFINDER_inchikey, NA_character_))

df$best_formula  <- ifelse(!is.na(df$Formula) & nchar(trimws(df$Formula)) > 0,
                            trimws(df$Formula),
                            ifelse("MSFINDER_formula" %in% colnames(df),
                                   df$MSFINDER_formula, NA_character_))

df$inchikey_source <- case_when(
  is_msdial_id ~ "MSDIAL",
  is_msfinder  ~ "MSFINDER",
  TRUE         ~ NA_character_
)

# -----------------------------------------------------------------------------
# 5. Per-feature VIP threshold (depends on annotation status)
# -----------------------------------------------------------------------------
df$vip_thresh_used <- ifelse(df$is_annotated,
                              VIP_THRESH_ANNOTATED,
                              VIP_THRESH_UNKNOWN)

df$vip_early_pass  <- df$VIP_early >= df$vip_thresh_used
df$vip_late_pass   <- df$VIP_late  >= df$vip_thresh_used
df$vip_any_pass    <- df$vip_early_pass | df$vip_late_pass

cat(sprintf("\nVIP filter summary (thresh annotated=%.2f, unknown=%.2f):\n",
            VIP_THRESH_ANNOTATED, VIP_THRESH_UNKNOWN))
cat(sprintf("  VIP pass (annotated) : %d / %d\n",
            sum(df$vip_any_pass &  df$is_annotated),
            sum(df$is_annotated)))
cat(sprintf("  VIP pass (unknown)   : %d / %d\n",
            sum(df$vip_any_pass & !df$is_annotated),
            sum(!df$is_annotated)))

# -----------------------------------------------------------------------------
# 6. Loading quadrant
# -----------------------------------------------------------------------------
df$loading_quad <- case_when(
  abs(df$loading_early) < LOADING_ABS_MIN |
  abs(df$loading_late)  < LOADING_ABS_MIN  ~ "Near_zero",
  df$loading_early >= 0 & df$loading_late >= 0 ~ "Q1",
  df$loading_early <  0 & df$loading_late >= 0 ~ "Q2",
  df$loading_early <  0 & df$loading_late <  0 ~ "Q3",
  df$loading_early >= 0 & df$loading_late <  0 ~ "Q4"
)

# -----------------------------------------------------------------------------
# 7. Mechanism class
# -----------------------------------------------------------------------------
df$mechanism_class <- case_when(
  # 方向逆転 (Q2/Q4) → 機構差候補（VIP条件: いずれかの期間でpass）
  df$loading_quad == "Q2" & df$vip_any_pass ~ "Reversed_Q2",
  df$loading_quad == "Q4" & df$vip_any_pass ~ "Reversed_Q4",
  # 両期間 VIP pass → 共通機構
  df$loading_quad == "Q1" & df$vip_early_pass & df$vip_late_pass ~ "Shared_Q1",
  df$loading_quad == "Q3" & df$vip_early_pass & df$vip_late_pass ~ "Shared_Q3",
  # 期間特異的
  df$vip_early_pass & !df$vip_late_pass ~ "Early_specific",
  !df$vip_early_pass & df$vip_late_pass ~ "Late_specific",
  # 残り
  TRUE ~ "Low_VIP"
)

# -----------------------------------------------------------------------------
# 8. MetaboAnalystR 用スコア列
# -----------------------------------------------------------------------------
# mechanism_score = loading_early × loading_late
#   正値 → Q1/Q3（一致方向）
#   負値 → Q2/Q4（逆転 = 機構差候補）
df$mechanism_score   <- round(df$loading_early * df$loading_late, 8)

# vip_max: いずれかの期間の最大 VIP（MSEA ランキング用）
df$vip_max           <- pmax(df$VIP_early, df$VIP_late)

# signed_vip: 機構差候補を負, 一致を正に符号付け（MSEA 方向性スコア用）
df$signed_vip        <- df$vip_max * sign(df$mechanism_score)

# -----------------------------------------------------------------------------
# 9. Filter to candidates
# -----------------------------------------------------------------------------
candidates <- df[df$vip_any_pass, ]
candidates <- candidates[order(candidates$mechanism_class,
                               -candidates$vip_max), ]

cat(sprintf("\nCandidates after filter: %d features\n", nrow(candidates)))
print(table(candidates$mechanism_class))
cat(sprintf("\n  Annotated: %d\n", sum(candidates$is_annotated)))
cat(sprintf("  Unknown  : %d\n", sum(!candidates$is_annotated)))

# -----------------------------------------------------------------------------
# 10. Output ① — 全候補詳細 CSV
# -----------------------------------------------------------------------------
out_cols <- c(
  # 識別子
  "Alignment_ID", "best_name", "best_inchikey", "best_formula",
  "inchikey_source",
  # MS-DIAL metadata
  "Rt_min", "Mz", "Adduct_type", "Ontology",
  # MS-FINDER
  if ("MSFINDER_total_score" %in% colnames(candidates)) "MSFINDER_total_score",
  if ("MSFINDER_ontology"    %in% colnames(candidates)) "MSFINDER_ontology",
  # PLS-DA 結果
  "VIP_early", "VIP_late", "vip_max",
  "loading_early", "loading_late",
  "loading_quad", "mechanism_class",
  "mechanism_score", "signed_vip",
  # フィルタ補足
  "is_annotated", "vip_thresh_used"
)
out_cols <- out_cols[out_cols %in% colnames(candidates)]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_main <- file.path(OUT_DIR, "plsda_mechanism_candidates.csv")
write.csv(candidates[, out_cols], out_main, row.names = FALSE)
cat(sprintf("\n[①] Saved: %s\n", out_main))

# -----------------------------------------------------------------------------
# 11. Output ② — ORA 用 (MetaboAnalystR::PerformPathwayAnalysis)
#     1 列目 = compound identifier (InChIKey, available features only)
#     2 列目 = mechanism_class
# -----------------------------------------------------------------------------
ora_df <- candidates[!is.na(candidates$best_inchikey),
                     c("best_inchikey", "best_name", "mechanism_class", "vip_max")]
colnames(ora_df)[1] <- "InChIKey"
colnames(ora_df)[2] <- "Name"

out_ora <- file.path(OUT_DIR, "plsda_metaboanalyst_ora.csv")
write.csv(ora_df, out_ora, row.names = FALSE)
cat(sprintf("[②] ORA CSV saved: %s  (%d features with InChIKey)\n",
            out_ora, nrow(ora_df)))

# Class breakdown for ORA
cat("    Breakdown by mechanism_class:\n")
print(table(ora_df$mechanism_class))

# -----------------------------------------------------------------------------
# 12. Output ③ — MSEA 用 (MetaboAnalystR::PerformMSEA)
#     1 列目 = Name または InChIKey
#     2 列目 = ranking score (signed_vip)
#     ※ MetaboAnalystR の mSet$dataSet$cmpd と一致させる
# -----------------------------------------------------------------------------
msea_df <- candidates[!is.na(candidates$best_inchikey),
                      c("best_inchikey", "signed_vip")]
colnames(msea_df) <- c("InChIKey", "score")
msea_df <- msea_df[order(-abs(msea_df$score)), ]

out_msea <- file.path(OUT_DIR, "plsda_metaboanalyst_msea.csv")
write.csv(msea_df, out_msea, row.names = FALSE)
cat(sprintf("[③] MSEA CSV saved: %s  (%d features)\n",
            out_msea, nrow(msea_df)))

# -----------------------------------------------------------------------------
# 13. MetaboAnalystR 使用メモ（コンソール出力）
# -----------------------------------------------------------------------------
cat("
================================================================================
MetaboAnalystR 使用メモ
================================================================================

[ORA — pathway over-representation analysis]
  library(MetaboAnalystR)
  mSet <- InitDataObjects('conc', 'pathora', FALSE)
  mSet <- Read.TextData(mSet, 'plsda_metaboanalyst_ora.csv', 'rowu', 'disc')
  mSet <- CrossReferencing(mSet, 'inchikey')
  mSet <- CreateMappingResultTable(mSet)
  mSet <- SetKEGG.PathLib(mSet, 'hsa')
  mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg')
  # -> 各 mechanism_class でサブセットして別々に実行することを推奨

[MSEA — ranked list enrichment]
  library(MetaboAnalystR)
  mSet <- InitDataObjects('conc', 'msetqea', FALSE)
  mSet <- Read.TextData(mSet, 'plsda_metaboanalyst_msea.csv', 'rowu', 'cont')
  mSet <- CrossReferencing(mSet, 'inchikey')
  mSet <- PerformQEA(mSet, 'fisher', 'msea')
  # score: signed_vip (正 = 一致方向, 負 = 逆転 = 機構差候補)

[mechanism_class の優先順位]
  Reversed_Q2 / Reversed_Q4  → 機構差の最有力候補（最優先）
  Shared_Q1 / Shared_Q3      → 両摂動に共通する経路
  Early_specific              → 摂動X 特異的
  Late_specific               → 摂動Y 特異的
================================================================================
")
