# ==============================================================================
# 07_bothup_report.R
# Purpose : Both_up グループ詳細 PDF レポート（17 ページ，A4 縦）
#   Page 1      : TMAO 時系列（ID 414 と 415 を別々に，被験者ごと）
#   Pages 2-17  : 胆汁酸 16 種 × TMAO（ID 414 / 415）散布図
#                 期間別（T1vsT2 / T1vsT3 / T4vsT5 / T4vsT6）× 4 種/ページ
#                 左列 = ID 414，右列 = ID 415
#
# 実行方法: プロジェクトルートを作業ディレクトリにして Source
# 依存パッケージ: install.packages(c("tidyverse", "patchwork"))
# ==============================================================================

library(tidyverse)
library(patchwork)

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJ_ROOT <- "."

LFC_CSV    <- file.path(PROJ_ROOT, "data/log2FC_values_production.csv")
GROUPS_CSV <- file.path(PROJ_ROOT, "output/tmao_subject_groups.csv")
OUT_PDF    <- file.path(PROJ_ROOT, "output/Both_up_detailed_report.pdf")

# ── 定数 ──────────────────────────────────────────────────────────────────────
BILE_IDS <- c(16774, 3590, 5279, 4394, 5534, 4863, 4657,
              17108, 17189, 4492, 3788, 4370, 3593, 14128, 4236, 5209)

BILE_NAMES <- c(
  "16774" = "Thr-CDCA (ID 16774)",      "3590"  = "Deoxycholic Acid (ID 3590)",
  "5279"  = "His-CDCA (ID 5279)",       "4394"  = "tauroursoDCA (ID 4394)",
  "5534"  = "Arg-CDCA (ID 5534)",       "4863"  = "Asn-CDCA (ID 4863)",
  "4657"  = "Val-CDCA (ID 4657)",       "17108" = "Glu-CDCA (ID 17108)",
  "17189" = "Gln-CDCA (ID 17189)",      "4492"  = "Asp-CDCA (ID 4492)",
  "3788"  = "GLYCOCHENODCA (ID 3788)",  "4370"  = "Cys-CDCA (ID 4370)",
  "3593"  = "HyoDCA (ID 3593)",         "14128" = "glycolithocholic (ID 14128)",
  "4236"  = "Ala-CDCA (ID 4236)",       "5209"  = "Met-CDCA (ID 5209)"
)

PERIODS <- tribble(
  ~key,      ~base_tp, ~end_tp, ~label,            ~color,    ~X_period,
  "T1vsT2",  1,        2,       "X  T1 → T2", "#005A8E", TRUE,
  "T1vsT3",  1,        3,       "X  T1 → T3", "#2E75B6", TRUE,
  "T4vsT5",  4,        5,       "Y  T4 → T5", "#7B1A1A", FALSE,
  "T4vsT6",  4,        6,       "Y  T4 → T6", "#C00000", FALSE
)

TMAO_IDS    <- c(414L, 415L)
SUBJ_COLORS <- c("#C00000", "#E06030", "#2E75B6", "#70AD47", "#7030A0")
A4_W <- 8.27; A4_H <- 11.69

# ── データ読み込み ────────────────────────────────────────────────────────────
lfc    <- read_csv(LFC_CSV,    show_col_types = FALSE)
groups <- read_csv(GROUPS_CSV, show_col_types = FALSE)

BOTH_SUBJS <- sort(groups |> filter(tmao_group == "Both_up") |> pull(Subject))
message("Both_up subjects: ", paste(BOTH_SUBJS, collapse = ", "))

SUBJ_COLOR_MAP <- setNames(SUBJ_COLORS[seq_along(BOTH_SUBJS)], as.character(BOTH_SUBJS))

# ── ヘルパー ─────────────────────────────────────────────────────────────────
# 指定 aid × 被験者の wide pivot（tp1〜tp6 列）
pivot_aid <- function(aid, subjs = BOTH_SUBJS) {
  d <- lfc |>
    filter(Alignment_ID == aid, Subject %in% subjs) |>
    select(Subject, Timepoint, intensity) |>
    pivot_wider(names_from = Timepoint, values_from = intensity,
                names_prefix = "tp", values_fill = 0)
  for (tp in 1:6) {
    col <- paste0("tp", tp)
    if (!col %in% names(d)) d[[col]] <- 0
  }
  arrange(d, Subject)
}

# log2FC per subject
get_fc_df <- function(aid, base_tp, end_tp, subjs = BOTH_SUBJS) {
  pivot_aid(aid, subjs) |>
    transmute(Subject, fc = .data[[paste0("tp", end_tp)]] - .data[[paste0("tp", base_tp)]])
}

# Pearson r ラベル
r_stat <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(label = "n.a.", p = 1))
  ct  <- cor.test(x[ok], y[ok])
  r   <- unname(ct$estimate)
  p   <- ct$p.value
  sig <- if (p < 0.05) "*" else ""
  ps  <- if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
  list(label = sprintf("r=%.2f %s%s", r, ps, sig), p = p)
}

# ═════════════════════════════════════════════════════════════════════════════
# ページ 1: TMAO 時系列
# ═════════════════════════════════════════════════════════════════════════════
make_ts_page <- function() {
  tps       <- 1:6
  tp_labels <- c("T1\n(X base)", "T2", "T3", "T4\n(Y base)", "T5", "T6")

  # 統一 y 軸の範囲を求める
  all_vals <- unlist(lapply(TMAO_IDS, function(tid) {
    d <- pivot_aid(tid, BOTH_SUBJS)
    unlist(d[, paste0("tp", 1:6)])
  }))
  ymin <- min(all_vals, na.rm = TRUE) - 0.5
  ymax <- max(all_vals, na.rm = TRUE) + 0.5

  # 長形式データ（全 TMAO_ID × 全被験者）
  ts_long <- bind_rows(lapply(TMAO_IDS, function(tid) {
    pivot_aid(tid, BOTH_SUBJS) |>
      pivot_longer(starts_with("tp"), names_to = "tp_col", values_to = "fc") |>
      mutate(Timepoint = as.integer(sub("tp", "", tp_col)),
             TMAO_ID   = factor(tid))
  }))

  # 被験者ごとのサブプロット
  plots_ts <- lapply(seq_along(BOTH_SUBJS), function(si) {
    subj <- BOTH_SUBJS[si]
    d    <- filter(ts_long, Subject == subj)

    p <- ggplot(d, aes(x = Timepoint, y = fc,
                        color = TMAO_ID, linetype = TMAO_ID, shape = TMAO_ID)) +
      annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
               fill = "#EEEEEE", alpha = 1) +
      geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
      geom_vline(xintercept = 3.5, color = "#AAAAAA", linewidth = 0.8, linetype = "dashed") +
      geom_line(linewidth = 1.4) +
      geom_point(size = 3) +
      scale_color_manual(values    = c("414" = "#C00000", "415" = "#2E75B6"),
                         labels    = c("414" = "ID 414",  "415" = "ID 415")) +
      scale_linetype_manual(values = c("414" = "solid",   "415" = "dashed"),
                            labels = c("414" = "ID 414",  "415" = "ID 415")) +
      scale_shape_manual(values    = c("414" = 16,        "415" = 17),
                         labels    = c("414" = "ID 414",  "415" = "ID 415")) +
      scale_x_continuous(breaks = 1:6, labels = tp_labels) +
      coord_cartesian(xlim = c(0.5, 6.5), ylim = c(ymin, ymax)) +
      labs(title = sprintf("Subject %d", subj),
           x = NULL, y = "log₂FC",
           color = NULL, linetype = NULL, shape = NULL) +
      theme_classic(base_size = 9) +
      theme(plot.title             = element_text(size = 10, face = "bold",
                                                   color = "#C00000"),
            legend.position        = if (si == 1) c(0.02, 0.98) else "none",
            legend.justification   = c(0, 1),
            legend.background      = element_rect(fill = alpha("white", 0.7),
                                                   color = NA),
            legend.key.size        = unit(0.5, "cm"),
            legend.text            = element_text(size = 7.5),
            axis.text.x            = element_text(size = 7),
            panel.grid.major.y     = element_line(color = "#EEEEEE", linewidth = 0.4))
    p
  })

  # 5 被験者 → 3×2 グリッド（空白パディング）
  if (length(plots_ts) == 5) plots_ts[[6]] <- ggplot() + theme_void()

  wrap_plots(plots_ts, ncol = 3, nrow = 2) +
    plot_annotation(
      title    = "Both↑ Group: TMAO Time Series  (Alignment ID 414 & 415, not averaged)",
      subtitle = "n = 5 subjects  |  Gray shading = X period (T1=baseline)  |  Y period: T4=baseline  |  Vertical dashed line separates X and Y periods",
      caption  = "X period: T1–T3 (baseline=T1) | Y period: T4–T6 (baseline=T4) | Values = log₂FC vs period baseline",
      theme    = theme(plot.title    = element_text(size = 12, face = "bold",
                                                     color = "#8B0000"),
                       plot.subtitle = element_text(size = 8.5, color = "#555555"),
                       plot.caption  = element_text(size = 7,   color = "#888888"))
    )
}

# ═════════════════════════════════════════════════════════════════════════════
# ページ 2-17: 散布図（期間 × 胆汁酸チャンク × TMAO ID）
# ═════════════════════════════════════════════════════════════════════════════
make_scatter_page <- function(period_row, bile_chunk, chunk_idx,
                               page_num, total_pages) {
  key     <- period_row$key
  base_tp <- period_row$base_tp
  end_tp  <- period_row$end_tp
  lbl     <- period_row$label
  col     <- period_row$color
  n_subj  <- length(BOTH_SUBJS)

  # TMAO FC（各 ID）
  tmao_fcs <- setNames(
    lapply(TMAO_IDS, function(tid) get_fc_df(tid, base_tp, end_tp)),
    as.character(TMAO_IDS)
  )

  # 4 bile acids × 2 TMAO IDs = 最大 8 プロット，左列 ID414・右列 ID415 の順
  plots <- vector("list", length(bile_chunk) * 2)
  idx   <- 0
  for (bi in seq_along(bile_chunk)) {
    bile_id   <- bile_chunk[bi]
    bile_name <- BILE_NAMES[as.character(bile_id)]
    bile_fc   <- get_fc_df(bile_id, base_tp, end_tp) |> rename(bile_fc = fc)

    for (ci in seq_along(TMAO_IDS)) {
      idx      <- idx + 1
      tmao_id  <- TMAO_IDS[ci]
      tmao_fc  <- tmao_fcs[[as.character(tmao_id)]] |> rename(tmao_fc = fc)

      plot_df <- inner_join(bile_fc, tmao_fc, by = "Subject") |>
        mutate(color = SUBJ_COLOR_MAP[as.character(Subject)],
               label = paste0("S", Subject))

      ok   <- is.finite(plot_df$bile_fc) & is.finite(plot_df$tmao_fc)
      rs   <- r_stat(plot_df$bile_fc, plot_df$tmao_fc)
      r_col <- if (rs$p < 0.05) col else "#888888"

      # 軸パディング
      if (sum(ok) >= 2) {
        xpad <- diff(range(plot_df$bile_fc[ok])) * 0.20 + 0.3
        ypad <- diff(range(plot_df$tmao_fc[ok])) * 0.20 + 0.3
        xlim <- range(plot_df$bile_fc[ok]) + c(-xpad, xpad)
        ylim <- range(plot_df$tmao_fc[ok]) + c(-ypad, ypad)
      } else {
        xlim <- ylim <- c(-4, 4)
      }

      p <- ggplot(filter(plot_df, ok),
                  aes(x = bile_fc, y = tmao_fc)) +
        geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
        geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.5) +
        geom_smooth(method = "lm", se = FALSE,
                    color = col, linewidth = 1.0, alpha = 0.8) +
        geom_point(aes(color = factor(Subject)), size = 3.5,
                   show.legend = FALSE) +
        geom_text(aes(label = label, color = factor(Subject)),
                  hjust = -0.25, vjust = 0.5, size = 2.5,
                  show.legend = FALSE) +
        scale_color_manual(values = SUBJ_COLOR_MAP) +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        annotate("label",
                 x = xlim[2], y = ylim[2],
                 label = rs$label, hjust = 1, vjust = 1,
                 size = 2.6, color = r_col, fontface = "bold",
                 fill = "white", label.size = 0.3, alpha = 0.9) +
        labs(
          title = if (ci == 1)
                    sprintf("%s", bile_name)
                  else
                    sprintf("ID %d", tmao_id),
          x = sprintf("%s  log₂FC", sub(" \\(ID.*", "", bile_name)),
          y = sprintf("TMAO ID %d  log₂FC", tmao_id)
        ) +
        theme_classic(base_size = 8) +
        theme(plot.title  = element_text(size = 8.5, face = "bold",
                                          color = "#333333"),
              axis.title  = element_text(size = 7.5),
              axis.text   = element_text(size = 7))

      plots[[idx]] <- p
    }
  }

  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title    = sprintf("Both↑ Group  |  %s  |  Bile Acid vs TMAO (ID 414 & 415)", lbl),
      subtitle = sprintf("Bile acids %d–%d of 16  |  Each point = 1 subject (n=%d)  |  Page %d / %d",
                          (chunk_idx - 1) * 4 + 1,
                          (chunk_idx - 1) * 4 + length(bile_chunk),
                          n_subj, page_num, total_pages),
      theme    = theme(
        plot.title    = element_text(size = 10, face = "bold", color = col),
        plot.subtitle = element_text(size = 8, color = "#555555")
      )
    )
}

# ═════════════════════════════════════════════════════════════════════════════
# PDF 出力
# ═════════════════════════════════════════════════════════════════════════════
bile_chunks  <- split(BILE_IDS, ceiling(seq_along(BILE_IDS) / 4))
TOTAL_PAGES  <- 1 + nrow(PERIODS) * length(bile_chunks)
message(sprintf("Generating %d pages → %s", TOTAL_PAGES, OUT_PDF))

# cairo_pdf は Unicode 文字を確実に扱える
cairo_pdf(OUT_PDF, width = A4_W, height = A4_H, onefile = TRUE)

# Page 1: 時系列
print(make_ts_page())
message("Page 1 done: time series")

# Pages 2-17: 散布図
page_num <- 1L
for (pi in seq_len(nrow(PERIODS))) {
  period_row <- PERIODS[pi, ]
  for (ci in seq_along(bile_chunks)) {
    page_num <- page_num + 1L
    print(make_scatter_page(period_row, bile_chunks[[ci]],
                             chunk_idx  = ci,
                             page_num   = page_num,
                             total_pages = TOTAL_PAGES))
    message(sprintf("Page %d done: %s acids %d-%d",
                    page_num, period_row$key,
                    (ci-1)*4 + 1, (ci-1)*4 + length(bile_chunks[[ci]])))
  }
}

dev.off()
message(sprintf("PDF saved (%d pages): %s", TOTAL_PAGES, OUT_PDF))
