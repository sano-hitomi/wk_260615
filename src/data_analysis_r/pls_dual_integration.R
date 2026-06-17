# =============================================================================
# pls_dual_integration.R
# Integrate PLS-DA results from two target metabolites (ID 414 and ID 4263)
# to identify shared pathway metabolites and directional relationships.
#
# Strategy:
#   Each feature has, from the two independent PLS-DAs:
#     loading_414_early,  loading_414_late   — Comp1 loadings w.r.t. ID 414 groups
#     loading_4263_early, loading_4263_late  — Comp1 loadings w.r.t. ID 4263 groups
#     VIP_414_early/late, VIP_4263_early/late
#
#   A 2×2 loading scatter (414 loading vs 4263 loading, Early / Late separately)
#   reveals the network topology:
#     Q1 (both+) / Q3 (both-) : feature co-varies with BOTH targets
#                               → likely same branch of the pathway
#     Q2 / Q4                  : feature moves in opposite directions relative
#                               to the two targets
#                               → metabolite may sit BETWEEN 414 and 4263,
#                                  or in a competing branch
#
# Additional integration metric:
#   dual_score = loading_414 × loading_4263  (positive = co-directional,
#                                              negative = opposing)
#
# Output:
#   data/production/processed/pls_dual_integration.csv   (feature table)
#   data/production/processed/pls_dual_msea.csv          (MSEA input: dual_score)
#   output/figures/pls_dual_loading_early.pdf / .png
#   output/figures/pls_dual_loading_late.pdf  / .png
#   output/figures/pls_dual_vip_overlap.pdf   / .png
#
# Prerequisites:
#   1. load_data_production.R
#   2. plsda_loadings_vip.R         → data/production/processed/plsda_vip_loadings.csv
#   3. plsda_loadings_vip_4263.R    → data/production/processed/plsda_vip_loadings_4263.csv
# =============================================================================

ID_A   <- 414
ID_B   <- 4263
CSV_A  <- "data/production/processed/plsda_vip_loadings.csv"
CSV_B  <- sprintf("data/production/processed/plsda_vip_loadings_%d.csv", ID_B)
OUT_DIR <- "output/figures"
OUT_CSV <- "data/production/processed/pls_dual_integration.csv"
OUT_MSEA <- "data/production/processed/pls_dual_msea.csv"

VIP_THRESH <- 1.0   # minimum VIP to highlight in plots

# -----------------------------------------------------------------------------
# 0. Prerequisites
# -----------------------------------------------------------------------------
required <- c("samplesheet", "feat_meta", "feat_mat")
missing  <- required[!sapply(required, exists)]
if (length(missing) > 0)
  stop("Missing: ", paste(missing, collapse = ", "),
       "\nRun load_data_production.R first.")
for (f in c(CSV_A, CSV_B)) {
  if (!file.exists(f))
    stop("CSV not found: ", f, "\nRun the corresponding plsda_loadings_vip script first.")
}

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
pkgs <- c("ggplot2", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Load VIP/loading tables
# -----------------------------------------------------------------------------
df_a <- read.csv(CSV_A, stringsAsFactors = FALSE)
df_b <- read.csv(CSV_B, stringsAsFactors = FALSE)

# Rename columns to distinguish the two analyses
rename_cols <- function(df, suffix) {
  colnames(df)[colnames(df) == "VIP_early"]     <- paste0("VIP_early_",     suffix)
  colnames(df)[colnames(df) == "VIP_late"]      <- paste0("VIP_late_",      suffix)
  colnames(df)[colnames(df) == "loading_early"] <- paste0("loading_early_", suffix)
  colnames(df)[colnames(df) == "loading_late"]  <- paste0("loading_late_",  suffix)
  df
}

df_a <- rename_cols(df_a, as.character(ID_A))
df_b <- rename_cols(df_b, as.character(ID_B))

# Drop duplicate Metabolite_name from df_b (keep from df_a)
df_b <- df_b[, !colnames(df_b) %in% "Metabolite_name", drop = FALSE]

# Merge on Alignment_ID (inner join → features present in both analyses)
merged <- merge(df_a, df_b, by = "Alignment_ID")
cat(sprintf("Features in both analyses: %d\n", nrow(merged)))

# -----------------------------------------------------------------------------
# 3. Merge annotation from feat_meta
# -----------------------------------------------------------------------------
meta_cols <- c("Alignment_ID", "Rt_min", "Mz", "Adduct_type",
               "Formula", "Ontology", "INCHIKEY", "SMILES")
msfinder_cols <- c("MSFINDER_annotated", "MSFINDER_structure",
                   "MSFINDER_formula", "MSFINDER_ontology",
                   "MSFINDER_inchikey", "MSFINDER_total_score")
keep_meta <- intersect(c(meta_cols, msfinder_cols), colnames(feat_meta))
meta_sub  <- feat_meta[, keep_meta, drop = FALSE]
merged    <- merge(merged, meta_sub, by = "Alignment_ID", all.x = TRUE)

# -----------------------------------------------------------------------------
# 4. Integration metrics
# -----------------------------------------------------------------------------
# co-directionality per period
merged$dual_score_early <- round(
  merged[[paste0("loading_early_", ID_A)]] *
  merged[[paste0("loading_early_", ID_B)]], 8)
merged$dual_score_late  <- round(
  merged[[paste0("loading_late_",  ID_A)]] *
  merged[[paste0("loading_late_",  ID_B)]], 8)
# Combined dual score (average of early and late)
merged$dual_score_avg   <- round((merged$dual_score_early + merged$dual_score_late) / 2, 8)

# VIP flags
merged$vip_A_any <- merged[[paste0("VIP_early_", ID_A)]] > VIP_THRESH |
                    merged[[paste0("VIP_late_",  ID_A)]] > VIP_THRESH
merged$vip_B_any <- merged[[paste0("VIP_early_", ID_B)]] > VIP_THRESH |
                    merged[[paste0("VIP_late_",  ID_B)]] > VIP_THRESH
merged$vip_both  <- merged$vip_A_any & merged$vip_B_any

# Quadrant classification (Early; both periods must agree for strong signal)
quad_class <- function(la, lb) {
  case_when(
    la >= 0 & lb >= 0 ~ "Q1_both_pos",
    la <  0 & lb <  0 ~ "Q3_both_neg",
    la >= 0 & lb <  0 ~ "Q4_A_pos_B_neg",
    la <  0 & lb >= 0 ~ "Q2_A_neg_B_pos",
    TRUE               ~ "NA"
  )
}
merged$quad_early <- quad_class(
  merged[[paste0("loading_early_", ID_A)]],
  merged[[paste0("loading_early_", ID_B)]])
merged$quad_late  <- quad_class(
  merged[[paste0("loading_late_",  ID_A)]],
  merged[[paste0("loading_late_",  ID_B)]])

# Pathway candidate class
merged$pathway_class <- case_when(
  merged$vip_both & merged$dual_score_avg > 0 ~ "Shared_pathway",    # co-directional, high VIP in both
  merged$vip_both & merged$dual_score_avg < 0 ~ "Opposing",          # opposing, high VIP in both → between 414 & 4263?
  merged$vip_A_any & !merged$vip_B_any        ~ sprintf("Specific_%d", ID_A),
  !merged$vip_A_any & merged$vip_B_any        ~ sprintf("Specific_%d", ID_B),
  TRUE                                         ~ "Low_VIP"
)

# max VIP across both analyses and both periods
merged$vip_max_overall <- pmax(
  merged[[paste0("VIP_early_", ID_A)]],
  merged[[paste0("VIP_late_",  ID_A)]],
  merged[[paste0("VIP_early_", ID_B)]],
  merged[[paste0("VIP_late_",  ID_B)]]
)

merged <- merged[order(merged$pathway_class, -merged$vip_max_overall), ]

cat(sprintf("\nPathway class summary:\n"))
print(table(merged$pathway_class))
cat(sprintf("\n  Shared_pathway (both high VIP, co-directional): %d features\n",
            sum(merged$pathway_class == "Shared_pathway")))
cat(sprintf("  Opposing       (both high VIP, opposing)       : %d features\n",
            sum(merged$pathway_class == "Opposing")))

# -----------------------------------------------------------------------------
# 5. Save main CSV
# -----------------------------------------------------------------------------
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(merged, OUT_CSV, row.names = FALSE)
cat(sprintf("\n[①] Integration table saved: %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 6. MSEA input (dual_score_avg as ranking score)
#    Use InChIKey if available
# -----------------------------------------------------------------------------
is_msdial  <- !is.na(merged$INCHIKEY) & nchar(trimws(merged$INCHIKEY)) > 0
is_msfinder <- if ("MSFINDER_annotated" %in% colnames(merged)) {
  !is.na(merged$MSFINDER_annotated) & merged$MSFINDER_annotated
} else rep(FALSE, nrow(merged))

merged$best_inchikey <- ifelse(is_msdial, trimws(merged$INCHIKEY),
                        ifelse(is_msfinder & "MSFINDER_inchikey" %in% colnames(merged),
                               merged$MSFINDER_inchikey, NA_character_))

msea_df <- merged[!is.na(merged$best_inchikey),
                  c("best_inchikey", "dual_score_avg")]
colnames(msea_df) <- c("InChIKey", "score")
msea_df <- msea_df[order(-abs(msea_df$score)), ]
write.csv(msea_df, OUT_MSEA, row.names = FALSE)
cat(sprintf("[②] MSEA input saved: %s  (%d features with InChIKey)\n",
            OUT_MSEA, nrow(msea_df)))

# -----------------------------------------------------------------------------
# 7. Plot helpers
# -----------------------------------------------------------------------------
trunc_name <- function(x, n = 30) {
  ifelse(!is.na(x) & nchar(x) > n, paste0(substr(x, 1, n), "…"), x)
}

make_dual_loading_plot <- function(df, period, id_a, id_b, vip_thresh) {
  xa <- paste0("loading_", period, "_", id_a)
  xb <- paste0("loading_", period, "_", id_b)
  va <- paste0("VIP_", period, "_", id_a)
  vb <- paste0("VIP_", period, "_", id_b)

  df$vip_flag <- case_when(
    df[[va]] > vip_thresh & df[[vb]] > vip_thresh ~ "Both",
    df[[va]] > vip_thresh                          ~ sprintf("ID %d only", id_a),
    df[[vb]] > vip_thresh                          ~ sprintf("ID %d only", id_b),
    TRUE                                           ~ "Low VIP"
  )

  # Label top features (VIP > threshold in either analysis)
  df_hi <- df[df[[va]] > vip_thresh | df[[vb]] > vip_thresh, ]
  label_col <- ifelse("MSFINDER_structure" %in% colnames(df) & !is.na(df_hi$MSFINDER_structure),
                      "MSFINDER_structure", "Metabolite_name")
  df_hi$label <- trunc_name(
    ifelse(!is.na(df_hi[[label_col]]) & df_hi[[label_col]] != "",
           paste0(df_hi$Alignment_ID, " ", df_hi[[label_col]]),
           as.character(df_hi$Alignment_ID)), n = 28)

  period_label <- if (period == "early") "Early (T1–T3)" else "Late (T4–T6)"

  ggplot(df, aes(x = .data[[xa]], y = .data[[xb]], color = vip_flag)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_point(size = 1.4, alpha = 0.5) +
    # Label high-VIP features
    ggplot2::geom_text(
      data = df_hi,
      aes(x = .data[[xa]], y = .data[[xb]], label = label),
      color = "black", size = 2.2, vjust = -0.7, hjust = 0.5,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = setNames(
        c("#E41A1C", "#377EB8", "#4DAF4A", "grey70"),
        c("Both", sprintf("ID %d only", id_a), sprintf("ID %d only", id_b), "Low VIP")
      ),
      name = sprintf("VIP > %g", vip_thresh)
    ) +
    labs(
      title    = sprintf("Dual loading scatter — %s", period_label),
      subtitle = sprintf(
        "X axis: Comp1 loading from PLS-DA (Y = ID %d groups)\nY axis: Comp1 loading from PLS-DA (Y = ID %d groups)\nQ1/Q3 = co-directional (same pathway branch) | Q2/Q4 = opposing (between-target metabolite?)",
        id_a, id_b),
      x = sprintf("Loading — Y = ID %d", id_a),
      y = sprintf("Loading — Y = ID %d", id_b)
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 7.5, color = "grey40"),
      legend.position  = "right"
    )
}

# -----------------------------------------------------------------------------
# 8. VIP overlap dot plot
#    X = max VIP from analysis A, Y = max VIP from analysis B
# -----------------------------------------------------------------------------
make_vip_overlap_plot <- function(df, id_a, id_b, vip_thresh) {
  va_max <- pmax(df[[paste0("VIP_early_", id_a)]], df[[paste0("VIP_late_", id_a)]])
  vb_max <- pmax(df[[paste0("VIP_early_", id_b)]], df[[paste0("VIP_late_", id_b)]])
  df2 <- data.frame(
    Alignment_ID    = df$Alignment_ID,
    Metabolite_name = df$Metabolite_name,
    vip_a           = va_max,
    vip_b           = vb_max,
    stringsAsFactors = FALSE
  )
  df2$flag <- case_when(
    df2$vip_a > vip_thresh & df2$vip_b > vip_thresh ~ "Both",
    df2$vip_a > vip_thresh                           ~ sprintf("ID %d only", id_a),
    df2$vip_b > vip_thresh                           ~ sprintf("ID %d only", id_b),
    TRUE                                             ~ "Low VIP"
  )
  df_hi <- df2[df2$vip_a > vip_thresh | df2$vip_b > vip_thresh, ]
  df_hi$label <- trunc_name(
    ifelse(!is.na(df_hi$Metabolite_name) & df_hi$Metabolite_name != "",
           paste0(df_hi$Alignment_ID, " ", df_hi$Metabolite_name),
           as.character(df_hi$Alignment_ID)), n = 28)

  ggplot(df2, aes(x = vip_a, y = vip_b, color = flag)) +
    geom_hline(yintercept = vip_thresh, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = vip_thresh, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_point(size = 1.4, alpha = 0.5) +
    ggplot2::geom_text(
      data = df_hi,
      aes(x = vip_a, y = vip_b, label = label),
      color = "black", size = 2.2, vjust = -0.7, hjust = 0.5,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = setNames(
        c("#E41A1C", "#377EB8", "#4DAF4A", "grey70"),
        c("Both", sprintf("ID %d only", id_a), sprintf("ID %d only", id_b), "Low VIP")
      ),
      name = sprintf("VIP > %g", vip_thresh)
    ) +
    labs(
      title    = "VIP overlap: ID 414 vs ID 4263 analyses",
      subtitle = sprintf(
        "X = max VIP (early/late) from PLS-DA Y=ID %d | Y = max VIP from PLS-DA Y=ID %d\nFeatures in upper-right quadrant are important for both targets → pathway candidates",
        id_a, id_b),
      x = sprintf("Max VIP — Y = ID %d", id_a),
      y = sprintf("Max VIP — Y = ID %d", id_b)
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 7.5, color = "grey40")
    )
}

# -----------------------------------------------------------------------------
# 9. Generate plots
# -----------------------------------------------------------------------------
p_dual_early  <- make_dual_loading_plot(merged, "early", ID_A, ID_B, VIP_THRESH)
p_dual_late   <- make_dual_loading_plot(merged, "late",  ID_A, ID_B, VIP_THRESH)
p_vip_overlap <- make_vip_overlap_plot(merged, ID_A, ID_B, VIP_THRESH)

# -----------------------------------------------------------------------------
# 10. Save plots
# -----------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, base_name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(base_name, ".pdf")),
         plot = p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(base_name, ".png")),
         plot = p, width = w, height = h, dpi = 150)
  cat(sprintf("  Saved: %s\n", base_name))
}

cat("\nSaving plots:\n")
save_plot(p_dual_early,  "pls_dual_loading_early",  7.5, 6)
save_plot(p_dual_late,   "pls_dual_loading_late",   7.5, 6)
save_plot(p_vip_overlap, "pls_dual_vip_overlap",    7,   5.5)

# -----------------------------------------------------------------------------
# 11. Interpretation memo
# -----------------------------------------------------------------------------
cat("
================================================================================
解釈ガイド (pls_dual_integration.R)
================================================================================

【dual loading scatter (Early / Late)】
  軸: X = PLS-DA(Y=ID 414)のComp1ローディング
      Y = PLS-DA(Y=ID 4263)のComp1ローディング

  Q1 (X+, Y+) / Q3 (X-, Y-):
    → 両ターゲットと同方向に動く代謝物
    → 414 と 4263 が同一パスウェイ上にある場合、
       その上流（または下流）の共通前駆体・産物候補

  Q2 (X-, Y+) / Q4 (X+, Y-):
    → 両ターゲットに対して逆方向
    → 414 と 4263 の間に位置する代謝物、または競合経路
    → 最も『パスウェイの分岐点』として興味深い

【VIP overlap plot】
  右上象限に入る特徴量 = 両方の PLS-DA で重要
  → パスウェイ解析の最優先入力候補

【pathway_class の優先順位】
  Shared_pathway : 両ターゲットと同方向・高VIP → 共通経路の主構成代謝物
  Opposing       : 両ターゲットと逆方向・高VIP → 414⇔4263 間の中間代謝物候補
  Specific_414   : 414にのみ関連
  Specific_4263  : 4263にのみ関連

【pls_dual_msea.csv の使い方 (MetaboAnalystR)】
  score列 = dual_score_avg (正=共方向, 負=逆方向)
  library(MetaboAnalystR)
  mSet <- InitDataObjects('conc', 'msetqea', FALSE)
  mSet <- Read.TextData(mSet, 'pls_dual_msea.csv', 'rowu', 'cont')
  mSet <- CrossReferencing(mSet, 'inchikey')
  mSet <- PerformQEA(mSet, 'fisher', 'msea')
================================================================================
")
