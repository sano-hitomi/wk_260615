# ==============================================================================
# 07_udca_detection.R
# Purpose : UDCA（ursodeoxycholic acid）の検出確認と DCA との経路比較
#
#   胆汁酸二次代謝経路の分岐を確認:
#     CDCA → (7α-脱水酸素化) → DCA       ← 腸内細菌（C. scindens 等）
#     CDCA → (7β-異性化) → UDCA           ← 腸内細菌（Ruminococcus torques 等）
#
#   本スクリプトの処理:
#   1. feature_metadata から UDCA エントリを InChIKey / HMDB / 名前で検索
#   2. 検出された場合: DCA・UDCA・CDCA の時系列を並べて比較
#   3. 検出されない場合: 類縁代謝物（タウロ/グリコ UDCA 等）で代替確認
#
# 出力:
#   output/figures/udca_detection_timeseries.pdf / .png
#   output/reports/udca_scan_result.csv
#
# 依存: install.packages(c("tidyverse", "patchwork"))
# ==============================================================================

library(tidyverse)
library(patchwork)

# ── パス設定 ──────────────────────────────────────────────────────────────────
PROJ_ROOT  <- "."
LFC_CSV    <- file.path(PROJ_ROOT, "data/production/processed/log2FC_values_production.csv")
FMETA_CSV  <- file.path(PROJ_ROOT, "data/production/processed/feature_metadata.csv")
GROUPS_CSV <- file.path(PROJ_ROOT, "data/ref/tmao_subject_groups.csv")
OUT_FIG    <- file.path(PROJ_ROOT, "output/figures")
OUT_REP    <- file.path(PROJ_ROOT, "output/reports")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_REP, showWarnings = FALSE, recursive = TRUE)

# ── InChIKey プレフィックス（14文字）──────────────────────────────────────────
#   DCA  (deoxycholic acid)       : KXGVEGMKQFWNSR  HMDB0000626
#   UDCA (ursodeoxycholic acid)   : RUDATHOIWHWSTL  HMDB0000491
#   CDCA (chenodeoxycholic acid)  : RUDATHOIWHWSTL は誤り; CDCA: CCBHSGSBSCKDSC HMDB0000518
#   ※ UDCA と CDCA は立体異性体なので InChIKey が近い点に注意
#     UDCA InChIKey : RUDATHOIWHWSTL-GNKTZXOSSA-N
#     CDCA InChIKey : RUDATHOIWHWSTL-HKUYNNGSSA-N  ← 後半のみ異なる

TARGETS <- tribble(
  ~short, ~full_name,                  ~ik_prefix,       ~hmdb_id,      ~known_id,
  "DCA",  "Deoxycholic acid",          "KXGVEGMKQFWNSR", "HMDB0000626", 10813L,
  "UDCA", "Ursodeoxycholic acid",      "RUDATHOIWHWSTL", "HMDB0000491", NA_integer_,
  "CDCA", "Chenodeoxycholic acid",     "CCBHSGSBSCKDSC", "HMDB0000518", NA_integer_,
  # 抱合体（存在した場合のサブ解析用）
  "TauroUDCA", "Tauro-UDCA",           "RMMJESDSFUWWGT", "HMDB0000850", 4394L,  # ID 4394 = tauroursoDCA（H1A より）
  "GlycoUDCA", "Glyco-UDCA",           NA_character_,    "HMDB0000946", NA_integer_
)

# ── feature_metadata 読み込み ─────────────────────────────────────────────────
if (!file.exists(FMETA_CSV)) stop(sprintf("feature_metadata が見つかりません: %s", FMETA_CSV))
feat_meta <- read_csv(FMETA_CSV, show_col_types = FALSE)
cat(sprintf("feature_metadata: %d features × %d columns\n",
            nrow(feat_meta), ncol(feat_meta)))
cat("列名:", paste(colnames(feat_meta), collapse = ", "), "\n\n")

# ── 胆汁酸を InChIKey / HMDB / 名前で検索するヘルパー ────────────────────────
search_metabolite <- function(meta, ik_prefix, hmdb_id, name_pattern, known_id = NULL) {
  result_ids <- integer(0)

  # 1. known_id が存在するか確認
  if (!is.null(known_id) && !is.na(known_id) && known_id %in% meta$Alignment_ID) {
    result_ids <- c(result_ids, known_id)
  }

  # 2. InChIKey 部分一致（full match: 14文字プレフィックス）
  if (!is.na(ik_prefix)) {
    ik_col <- grep("inchikey|INCHIKEY", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(ik_col)) {
      hits <- meta$Alignment_ID[
        !is.na(meta[[ik_col]]) &
          startsWith(toupper(meta[[ik_col]]), toupper(ik_prefix))
      ]
      result_ids <- union(result_ids, hits)
    }
  }

  # 3. HMDB ID
  hmdb_col <- grep("hmdb", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(hmdb_col)) {
    hits <- meta$Alignment_ID[
      !is.na(meta[[hmdb_col]]) &
        trimws(meta[[hmdb_col]]) == hmdb_id
    ]
    result_ids <- union(result_ids, hits)
  }

  # 4. 名前パターン
  if (!is.null(name_pattern)) {
    name_col <- grep("^(Metabolite.name|Name|metabolite_name)", colnames(meta),
                     ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(name_col)) {
      hits <- meta$Alignment_ID[
        !is.na(meta[[name_col]]) &
          grepl(name_pattern, meta[[name_col]], ignore.case = TRUE)
      ]
      result_ids <- union(result_ids, hits)
    }
  }

  sort(unique(result_ids))
}

# ── 各代謝物を検索 ────────────────────────────────────────────────────────────
name_patterns <- list(
  DCA      = "deoxychol",
  UDCA     = "ursodeoxy",
  CDCA     = "chenodeoxy",
  TauroUDCA = "tauroursodeo|taurourso",
  GlycoUDCA = "glycoursodeo|glycourso"
)

scan_result <- pmap_dfr(TARGETS, function(short, full_name, ik_prefix, hmdb_id, known_id) {
  ids <- search_metabolite(feat_meta,
                           ik_prefix    = ik_prefix,
                           hmdb_id      = hmdb_id,
                           name_pattern = name_patterns[[short]],
                           known_id     = known_id)
  tibble(
    metabolite     = short,
    full_name      = full_name,
    detected       = length(ids) > 0,
    n_entries      = length(ids),
    alignment_ids  = if (length(ids) > 0) paste(ids, collapse = ";") else NA_character_
  )
})

cat("=== UDCA スキャン結果 ===\n")
print(scan_result)

# CSV 保存
out_csv <- file.path(OUT_REP, "udca_scan_result.csv")
write_csv(scan_result, out_csv)
message(sprintf("スキャン結果 CSV 保存: %s", out_csv))

# ── LFC CSV 読み込み ──────────────────────────────────────────────────────────
if (!file.exists(LFC_CSV)) stop(sprintf("LFC CSV が見つかりません: %s", LFC_CSV))
lfc <- read_csv(LFC_CSV, show_col_types = FALSE)

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

tp_labels <- c("T1\n(X base)", "T2", "T3", "T4\n(Y base)", "T5", "T6")

# ── 時系列プロット関数 ────────────────────────────────────────────────────────
make_ts_panel <- function(aid, label, color) {
  df <- lfc %>%
    filter(Alignment_ID == aid) %>%
    left_join(groups %>% select(Subject, group), by = "Subject") %>%
    filter(!is.na(group))

  if (nrow(df) == 0) {
    return(ggplot() +
             labs(title = sprintf("%s (ID %d)", label, aid),
                  subtitle = "データなし") +
             theme_void(base_size = 9) +
             theme(plot.title = element_text(face = "bold")))
  }

  grp_sum <- df %>%
    group_by(Timepoint) %>%
    summarise(mean_lfc = mean(intensity, na.rm = TRUE),
              se_lfc   = sd(intensity, na.rm = TRUE) / sqrt(n()),
              .groups  = "drop")

  ggplot() +
    annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
             fill = "#EEEEEE", alpha = 1) +
    geom_line(data = df,
              aes(x = Timepoint, y = intensity, group = Subject),
              color = "#BBBBBB", linewidth = 0.4, alpha = 0.5) +
    geom_ribbon(data = grp_sum,
                aes(x = Timepoint,
                    ymin = mean_lfc - se_lfc,
                    ymax = mean_lfc + se_lfc),
                fill = color, alpha = 0.25, inherit.aes = FALSE) +
    geom_line(data = grp_sum,
              aes(x = Timepoint, y = mean_lfc),
              color = color, linewidth = 1.8) +
    geom_point(data = grp_sum,
               aes(x = Timepoint, y = mean_lfc),
               color = color, size = 3) +
    geom_hline(yintercept = 0, color = "#AAAAAA", linetype = "dashed") +
    geom_vline(xintercept = 3.5, color = "#888888", linetype = "dashed") +
    scale_x_continuous(breaks = 1:6, labels = tp_labels) +
    labs(title = sprintf("%s (ID %d)", label, aid),
         x = NULL, y = "log₂FC") +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          panel.grid.minor = element_blank())
}

# ── 主要代謝物のプロット ──────────────────────────────────────────────────────
# 描画対象: DCA, UDCA, CDCA（それぞれ検出された場合のみ）
PLOT_TARGETS <- list(
  list(short = "DCA",  color = "#C00000"),
  list(short = "UDCA", color = "#2E75B6"),
  list(short = "CDCA", color = "#70AD47"),
  list(short = "TauroUDCA", color = "#7030A0")
)

plots <- list()
for (pt in PLOT_TARGETS) {
  row <- scan_result %>% filter(metabolite == pt$short)
  if (!row$detected) {
    # 検出なし → 空パネルを追加
    p <- ggplot() +
      labs(title = sprintf("%s — Not detected", pt$short)) +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold", color = "#999999"),
            plot.background = element_rect(fill = "#F9F9F9", color = "#CCCCCC"))
    plots[[pt$short]] <- p
    cat(sprintf("  %s: 未検出 → 空パネル\n", pt$short))
    next
  }

  # 複数エントリがある場合は S/N 最大を使用
  ids <- as.integer(strsplit(row$alignment_ids, ";")[[1]])
  if (length(ids) > 1) {
    sn_col <- grep("S.N|sn_average|SN", colnames(feat_meta), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(sn_col)) {
      sub <- feat_meta %>%
        filter(Alignment_ID %in% ids) %>%
        mutate(sn = as.numeric(.data[[sn_col]])) %>%
        arrange(desc(sn))
      use_id <- sub$Alignment_ID[1]
    } else {
      use_id <- ids[1]
    }
  } else {
    use_id <- ids[1]
  }

  cat(sprintf("  %s (ID %d): 検出済み → プロット作成\n", pt$short, use_id))
  plots[[pt$short]] <- make_ts_panel(use_id, pt$short, pt$color)
}

# ── UDCA / DCA スキャッタープロット（両方検出された場合） ──────────────────────
p_scatter <- NULL
dca_row  <- scan_result %>% filter(metabolite == "DCA",  detected)
udca_row <- scan_result %>% filter(metabolite == "UDCA", detected)

if (nrow(dca_row) > 0 && nrow(udca_row) > 0) {
  dca_id  <- as.integer(strsplit(dca_row$alignment_ids,  ";")[[1]])[1]
  udca_id <- as.integer(strsplit(udca_row$alignment_ids, ";")[[1]])[1]

  # Y 期平均 FC（T4→T5, T4→T6 の平均）
  get_y_fc <- function(aid) {
    lfc %>%
      filter(Alignment_ID == aid) %>%
      pivot_wider(names_from = Timepoint, values_from = intensity,
                  names_prefix = "tp", values_fill = 0) %>%
      { d <- .; for (tp in 1:6) { c <- paste0("tp", tp); if (!c %in% names(d)) d[[c]] <- 0 }; d }() %>%
      mutate(fc = ((tp5 - tp4) + (tp6 - tp4)) / 2) %>%
      select(Subject, fc)
  }

  dca_fc  <- get_y_fc(dca_id)  %>% rename(dca_fc  = fc)
  udca_fc <- get_y_fc(udca_id) %>% rename(udca_fc = fc)

  sc_df <- inner_join(dca_fc, udca_fc, by = "Subject") %>%
    left_join(groups %>% select(Subject, group), by = "Subject") %>%
    filter(!is.na(group))

  ct <- tryCatch(
    cor.test(sc_df$dca_fc, sc_df$udca_fc),
    error = function(e) NULL
  )
  r_lab <- if (!is.null(ct))
    sprintf("r = %.2f\np = %.3f%s", ct$estimate, ct$p.value,
            if (ct$p.value < 0.05) " *" else "")
  else "n.a."

  p_scatter <- ggplot(sc_df, aes(x = dca_fc, y = udca_fc, color = group)) +
    geom_hline(yintercept = 0, color = "#CCCCCC") +
    geom_vline(xintercept = 0, color = "#CCCCCC") +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "#444444", linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.85) +
    annotate("label", x = Inf, y = Inf, label = r_lab,
             hjust = 1.1, vjust = 1.2, size = 3.5,
             color = "#333333", fontface = "bold",
             fill = "white", label.size = 0.3) +
    scale_color_manual(values = GROUP_COLORS, name = "TMAO Group") +
    labs(
      title    = sprintf("DCA (ID %d) vs UDCA (ID %d) — Y period mean FC",
                         dca_id, udca_id),
      subtitle = "両経路が競合するなら DCA↑ かつ UDCA↓ の反相関が期待される",
      x = "DCA log₂FC (Y mean)", y = "UDCA log₂FC (Y mean)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          panel.grid.minor = element_blank())
}

# ── レイアウト組立 ────────────────────────────────────────────────────────────
main_panels <- (plots$DCA | plots$UDCA | plots$CDCA | plots$TauroUDCA) +
  plot_annotation(
    title    = "Bile Acid Pathway Branches: DCA / UDCA / CDCA Trajectories",
    subtitle = "DCA path: CDCA → (7α-dehydroxylation by C. scindens) → DCA\nUDCA path: CDCA → (7β-epimerization) → UDCA\nGray = individual subjects  |  Colored = all-subjects mean ± SE",
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "#555555")
    )
  )

# ── 保存 ──────────────────────────────────────────────────────────────────────
out_pdf <- file.path(OUT_FIG, "udca_detection_timeseries.pdf")
pdf(out_pdf, width = 14, height = if (is.null(p_scatter)) 5 else 10)
print(main_panels)
if (!is.null(p_scatter)) print(p_scatter)
dev.off()
message(sprintf("PDF 保存: %s", out_pdf))

out_png <- file.path(OUT_FIG, "udca_detection_timeseries.png")
ggsave(out_png, main_panels, width = 14, height = 5, dpi = 150)
message(sprintf("PNG 保存: %s", out_png))

if (!is.null(p_scatter)) {
  out_sc_png <- file.path(OUT_FIG, "udca_vs_dca_scatter.png")
  ggsave(out_sc_png, p_scatter, width = 7, height = 6, dpi = 150)
  message(sprintf("散布図 PNG 保存: %s", out_sc_png))
}

message("\n07_udca_detection.R 完了。")
message(sprintf("UDCA 検出: %s",
                if (scan_result$detected[scan_result$metabolite == "UDCA"])
                  sprintf("YES (IDs: %s)", scan_result$alignment_ids[scan_result$metabolite == "UDCA"])
                else "NOT DETECTED"))
