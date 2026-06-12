# ==============================================================================
# 06_group_composite_scatter.R
# Purpose : グループ別コンポジット胆汁酸スコア vs TMAO 散布図
#   コンポジットスコア定義:
#     各期間（X or Y）の胆汁酸 16 種 × 全被験者 × 2 TP ペアの log2FC 値を
#     一括で z 標準化し，被験者×TP ペアごとに 16 種の z スコアを平均したもの。
#   TMAO: Alignment ID 414 + 415 の平均，Subject 18 除外。
#
#   Outputs:
#     output/scatter_xy/{group}_scatter.png  (7 ファイル)
#     output/scatter_xy/r_summary.csv
#
# 実行方法: プロジェクトルート（TMAO_pathway_analysis/）を作業ディレクトリに設定して Source
# 依存パッケージ: install.packages(c("tidyverse", "patchwork"))
# ==============================================================================

library(tidyverse)
library(patchwork)

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJ_ROOT <- "."

LFC_CSV    <- file.path(PROJ_ROOT, "data/log2FC_values_production.csv")
GROUPS_CSV <- file.path(PROJ_ROOT, "output/tmao_subject_groups.csv")
OUT_DIR    <- file.path(PROJ_ROOT, "output/scatter_xy")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 定数 ──────────────────────────────────────────────────────────────────────
BILE_IDS <- c(16774, 3590, 5279, 4394, 5534, 4863, 4657,
              17108, 17189, 4492, 3788, 4370, 3593, 14128, 4236, 5209)

WEAK_GROUPS <- c("Sparse", "Both_down", "X_down_only", "Y_down_only", "Weak")

GROUP_COLORS <- c(
  Both_up      = "#C00000", Y_only       = "#E06060",
  X_down_Y_up  = "#70AD47", X_up_Y_down  = "#ED7D31",
  X_only       = "#FFC000", Non_producer = "#4472C4",
  Weak         = "#A5A5A5"
)

GROUP_LABELS <- c(
  Both_up      = "Both↑  —  High TMAO producer",
  Y_only       = "Y↑ only  —  Carnitine-responsive",
  X_down_Y_up  = "X↓ Y↑  —  Primed suppression",
  X_up_Y_down  = "X↑ Y↓  —  X-period responder",
  X_only       = "X↑ only  —  X-period only",
  Non_producer = "Non-producer  —  No TMAO response",
  Weak         = "Weak / Other"
)

GROUP_ORDER <- names(GROUP_COLORS)

# ── データ読み込み ────────────────────────────────────────────────────────────
lfc <- read_csv(LFC_CSV, show_col_types = FALSE)

groups <- read_csv(GROUPS_CSV, show_col_types = FALSE) |>
  mutate(group = if_else(tmao_group %in% WEAK_GROUPS, "Weak", tmao_group))

# ── TMAO wide テーブル（ID414+415 平均，Subject18 除外）─────────────────────
tmao_wide <- lfc |>
  filter(Alignment_ID %in% c(414, 415), Subject != 18) |>
  group_by(Subject, Timepoint) |>
  summarise(val = mean(intensity, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = Timepoint, values_from = val,
              names_prefix = "tp", values_fill = 0)
for (tp in 1:6) {
  col <- paste0("tp", tp)
  if (!col %in% names(tmao_wide)) tmao_wide[[col]] <- 0
}

# ── 胆汁酸 log2FC 長形式テーブル（全 aid × 全 period pairs）──────────────────
get_bile_long <- function(period = c("Y", "X")) {
  period <- match.arg(period)
  if (period == "Y") {
    tp_pairs <- list(c(4,5), c(4,6))
    labels   <- c("T4vsT5", "T4vsT6")
  } else {
    tp_pairs <- list(c(1,2), c(1,3))
    labels   <- c("T1vsT2", "T1vsT3")
  }

  lfc |>
    filter(Alignment_ID %in% BILE_IDS, Subject != 18) |>
    select(Subject, Alignment_ID, Timepoint, intensity) |>
    pivot_wider(names_from = Timepoint, values_from = intensity,
                names_prefix = "tp", values_fill = 0) |>
    { d <- .; for (tp in 1:6) { c <- paste0("tp",tp); if(!c %in% names(d)) d[[c]] <- 0 }; d }() |>
    mutate(
      !!labels[1] := .data[[paste0("tp", tp_pairs[[1]][2])]] - .data[[paste0("tp", tp_pairs[[1]][1])]],
      !!labels[2] := .data[[paste0("tp", tp_pairs[[2]][2])]] - .data[[paste0("tp", tp_pairs[[2]][1])]]
    ) |>
    select(Subject, Alignment_ID, all_of(labels)) |>
    pivot_longer(all_of(labels), names_to = "tp_pair", values_to = "bile_fc")
}

# ── コンポジットスコア計算 ────────────────────────────────────────────────────
# 全 (Subject × tp_pair × Alignment_ID) の bile_fc を一括 z 標準化し，
# per (Subject, tp_pair) で 16 種平均 → composite
build_composite <- function(period = c("Y", "X")) {
  period <- match.arg(period)
  df <- get_bile_long(period)

  # 全値で z 標準化
  mu <- mean(df$bile_fc, na.rm = TRUE)
  sd_val <- sd(df$bile_fc, na.rm = TRUE)
  df <- mutate(df, z = (bile_fc - mu) / sd_val)

  # composite = 16 bile acid の z 平均 per (Subject, tp_pair)
  df |>
    group_by(Subject, tp_pair) |>
    summarise(composite = mean(z, na.rm = TRUE), .groups = "drop") |>
    mutate(period = period)
}

composite_Y <- build_composite("Y")
composite_X <- build_composite("X")

# ── TMAO log2FC per (Subject, tp_pair) ───────────────────────────────────────
tmao_Y <- bind_rows(
  transmute(tmao_wide, Subject, tp_pair = "T4vsT5", tmao_fc = tp5 - tp4),
  transmute(tmao_wide, Subject, tp_pair = "T4vsT6", tmao_fc = tp6 - tp4)
)
tmao_X <- bind_rows(
  transmute(tmao_wide, Subject, tp_pair = "T1vsT2", tmao_fc = tp2 - tp1),
  transmute(tmao_wide, Subject, tp_pair = "T1vsT3", tmao_fc = tp3 - tp1)
)

# composite + TMAO を結合
scatter_Y <- inner_join(composite_Y, tmao_Y, by = c("Subject", "tp_pair")) |>
  left_join(groups |> select(Subject, group), by = "Subject")
scatter_X <- inner_join(composite_X, tmao_X, by = c("Subject", "tp_pair")) |>
  left_join(groups |> select(Subject, group), by = "Subject")

# 全データの共通軸範囲（全グループで統一）
all_comp <- c(scatter_Y$composite, scatter_X$composite)
all_tmao <- c(scatter_Y$tmao_fc,   scatter_X$tmao_fc)
xlim_rng <- c(quantile(all_comp, 0.01, na.rm=TRUE) - 0.2,
               quantile(all_comp, 0.99, na.rm=TRUE) + 0.2)
ylim_rng <- c(quantile(all_tmao, 0.01, na.rm=TRUE) - 0.3,
               quantile(all_tmao, 0.99, na.rm=TRUE) + 0.3)

# ── ヘルパー ─────────────────────────────────────────────────────────────────
r_stat <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(label = "n.a.", r = NA_real_, p = NA_real_))
  ct  <- cor.test(x[ok], y[ok])
  r   <- unname(ct$estimate)
  p   <- ct$p.value
  sig <- if (p < 0.05) " *" else ""
  ps  <- if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
  list(label = sprintf("r = %.2f  %s%s", r, ps, sig), r = r, p = p)
}

# ── 1 グループ分の散布図パネル（Y 期または X 期）────────────────────────────
make_panel <- function(scatter_all, grp_subjs, period_lbl, tp_lbl, pt_col) {
  all_ok   <- filter(scatter_all, is.finite(composite), is.finite(tmao_fc))
  grp_data <- filter(all_ok, Subject %in% grp_subjs)
  n_subj   <- length(grp_subjs)
  n_pairs  <- nrow(grp_data)

  # 全体の回帰線（灰色破線）
  lm_all  <- lm(tmao_fc ~ composite, data = all_ok)
  line_all <- tibble(
    x = xlim_rng,
    y = predict(lm_all, newdata = data.frame(composite = xlim_rng))
  )

  # フォーカルグループの回帰線
  rs <- r_stat(grp_data$composite, grp_data$tmao_fc)
  line_grp <- NULL
  if (n_pairs >= 3) {
    lm_grp <- lm(tmao_fc ~ composite, data = grp_data)
    line_grp <- tibble(
      x = xlim_rng,
      y = predict(lm_grp, newdata = data.frame(composite = xlim_rng))
    )
  }

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
    geom_point(data = all_ok, aes(composite, tmao_fc),
               color = "#DDDDDD", size = 2.2) +
    geom_line(data = line_all, aes(x, y),
              color = "#AAAAAA", linewidth = 0.9, linetype = "dashed") +
    geom_point(data = grp_data, aes(composite, tmao_fc),
               color = pt_col, size = 3.2) +
    coord_cartesian(xlim = xlim_rng, ylim = ylim_rng) +
    annotate("text",
             x = xlim_rng[1] + diff(xlim_rng) * 0.04,
             y = ylim_rng[2] - diff(ylim_rng) * 0.04,
             label = rs$label, hjust = 0, vjust = 1,
             size = 3, color = pt_col, fontface = "bold") +
    annotate("text",
             x = xlim_rng[1] + diff(xlim_rng) * 0.04,
             y = ylim_rng[2] - diff(ylim_rng) * 0.16,
             label = sprintf("%d subj × 2 TP = %d pairs", n_subj, n_pairs),
             hjust = 0, vjust = 1, size = 2.6, color = "#777777") +
    labs(title = sprintf("%s  (%s)", period_lbl, tp_lbl),
         x = sprintf("Bile acid composite  (%s)", period_lbl),
         y = sprintf("TMAO log₂FC  (%s)", period_lbl)) +
    theme_classic(base_size = 9) +
    theme(plot.title  = element_text(size = 9, face = "bold", color = "#444444"),
          axis.title  = element_text(size = 8))

  if (!is.null(line_grp)) {
    p <- p + geom_line(data = line_grp, aes(x, y),
                        color = pt_col, linewidth = 1.8)
  }
  p
}

# ── グループごとに PNG 保存 ───────────────────────────────────────────────────
for (grp in GROUP_ORDER) {
  col       <- GROUP_COLORS[grp]
  grp_subjs <- groups |> filter(group == grp) |> pull(Subject)
  n_subj    <- length(grp_subjs)

  pY <- make_panel(scatter_Y, grp_subjs, "Y period", "T4→T5/T6", col)
  pX <- make_panel(scatter_X, grp_subjs, "X period", "T1→T2/T3", col)

  combined <- (pY | pX) +
    plot_annotation(
      title = sprintf("%s   (n = %d)", GROUP_LABELS[grp], n_subj),
      theme = theme(plot.title = element_text(size = 11, face = "bold", color = col,
                                               margin = margin(b = 4)))
    ) &
    theme(plot.background = element_rect(fill = "white", color = NA))

  out_path <- file.path(OUT_DIR, sprintf("%s_scatter.png", grp))
  ggsave(out_path, combined, width = 7.0, height = 3.2, dpi = 160, bg = "white")
  message("Saved: ", basename(out_path))
}

# ── r サマリー CSV ────────────────────────────────────────────────────────────
r_summary <- map_dfr(GROUP_ORDER, function(grp) {
  grp_subjs <- groups |> filter(group == grp) |> pull(Subject)
  map_dfr(c("Y", "X"), function(per) {
    sc <- if (per == "Y") scatter_Y else scatter_X
    gd <- filter(sc, Subject %in% grp_subjs)
    rs <- r_stat(gd$composite, gd$tmao_fc)
    n_pairs <- sum(is.finite(gd$composite) & is.finite(gd$tmao_fc))
    tibble(group = grp, period = per, n = length(grp_subjs),
           n_pairs = n_pairs, r = round(rs$r, 3), p = round(rs$p, 3))
  })
})

write_csv(r_summary, file.path(OUT_DIR, "r_summary.csv"))
message("Saved: r_summary.csv")
print(r_summary)
message("06_group_composite_scatter.R 完了。")
