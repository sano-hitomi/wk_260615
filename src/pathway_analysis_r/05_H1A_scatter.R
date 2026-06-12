# ==============================================================================
# 05_H1A_scatter.R
# Purpose : H1-A 仮説：二次胆汁酸 16 種と TMAO の相関可視化
#   Outputs:
#     output/figures/H1A_boxplot.png          -- グループ別ボックスプロット
#     output/figures/H1A_heatmap.png          -- 被験者×胆汁酸ヒートマップ
#     output/figures/H1A_scatter_all.png      -- 16 種散布図（Y 期平均，グループ色分け）
#     output/figures/H1A_scatter_tp_split.png -- 16 種散布図（Y 期：T4vsT5 ○ / T4vsT6 △）
#     output/figures/H1A_scatter_tp_split_Xperiod.png -- 同 X 期
#
# 実行方法（RStudio）:
#   Session > Set Working Directory > To Source File Location の後、
#   PROJ_ROOT <- "../.." を確認してから Source する。
#   または RStudio プロジェクト（.Rproj）のルートで PROJ_ROOT <- "." で OK。
#
# 依存パッケージ:
#   install.packages(c("tidyverse", "patchwork", "pheatmap", "scales"))
# ==============================================================================

library(tidyverse)
library(patchwork)
library(pheatmap)
library(scales)

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJ_ROOT <- "."          # RStudio プロジェクトルートから実行する場合
# PROJ_ROOT <- "../.."   # src/r/ から直接 source する場合はこちら

LFC_CSV    <- file.path(PROJ_ROOT, "data/log2FC_values_production.csv")
GROUPS_CSV <- file.path(PROJ_ROOT, "output/tmao_subject_groups.csv")
OUT_FIG    <- file.path(PROJ_ROOT, "output/figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# ── 定数 ──────────────────────────────────────────────────────────────────────
BILE_IDS <- c(16774, 3590, 5279, 4394, 5534, 4863, 4657,
              17108, 17189, 4492, 3788, 4370, 3593, 14128, 4236, 5209)

BILE_NAMES <- c(
  "16774" = "Thr-CDCA (16774)",      "3590"  = "Deoxycholic Acid (3590)",
  "5279"  = "His-CDCA (5279)",       "4394"  = "tauroursoDCA (4394)",
  "5534"  = "Arg-CDCA (5534)",       "4863"  = "Asn-CDCA (4863)",
  "4657"  = "Val-CDCA (4657)",       "17108" = "Glu-CDCA (17108)",
  "17189" = "Gln-CDCA (17189)",      "4492"  = "Asp-CDCA (4492)",
  "3788"  = "GLYCOCHENODCA (3788)",  "4370"  = "Cys-CDCA (4370)",
  "3593"  = "HyoDCA (3593)",         "14128" = "glycolithocholic (14128)",
  "4236"  = "Ala-CDCA (4236)",       "5209"  = "Met-CDCA (5209)"
)

WEAK_GROUPS <- c("Sparse", "Both_down", "X_down_only", "Y_down_only", "Weak")

GROUP_COLORS <- c(
  Both_up      = "#C00000", Y_only       = "#E06060",
  X_down_Y_up  = "#70AD47", X_up_Y_down  = "#ED7D31",
  X_only       = "#FFC000", Non_producer = "#4472C4",
  Weak         = "#A5A5A5"
)
GROUP_ORDER <- names(GROUP_COLORS)

# ── データ読み込み ────────────────────────────────────────────────────────────
lfc <- read_csv(LFC_CSV, show_col_types = FALSE)

groups <- read_csv(GROUPS_CSV, show_col_types = FALSE) |>
  mutate(group = if_else(tmao_group %in% WEAK_GROUPS, "Weak", tmao_group),
         group = factor(group, levels = GROUP_ORDER))

# ── TMAO シグナル：ID 414+415 平均，Subject 18 除外 ───────────────────────────
tmao_wide <- lfc |>
  filter(Alignment_ID %in% c(414, 415), Subject != 18) |>
  group_by(Subject, Timepoint) |>
  summarise(val = mean(intensity, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = Timepoint, values_from = val,
              names_prefix = "tp", values_fill = 0)

# tp列を確実に確保
for (tp in 1:6) {
  col <- paste0("tp", tp)
  if (!col %in% names(tmao_wide)) tmao_wide[[col]] <- 0
}

# ── 汎用 log2FC 計算 ───────────────────────────────────────────────────────────
# aide: Alignment_ID / base_tp, end_tp: タイムポイント番号
get_fc <- function(aid, base_tp, end_tp, exclude_subj = 18) {
  lfc |>
    filter(Alignment_ID == aid, !Subject %in% exclude_subj) |>
    select(Subject, Timepoint, intensity) |>
    pivot_wider(names_from = Timepoint, values_from = intensity,
                names_prefix = "tp", values_fill = 0) |>
    { d <- .; for (tp in 1:6) { c <- paste0("tp",tp); if(!c %in% names(d)) d[[c]] <- 0 }; d }() |>
    mutate(fc = .data[[paste0("tp", end_tp)]] - .data[[paste0("tp", base_tp)]]) |>
    select(Subject, fc)
}

# ── 全胆汁酸 × 全期間の長形式データ作成 ──────────────────────────────────────
# period: "Y" or "X" / 返り値: Subject, Alignment_ID, bile_name, tp_pair, bile_fc, tmao_fc, group
build_long <- function(period = c("Y", "X")) {
  period <- match.arg(period)
  if (period == "Y") {
    pairs  <- list(c(4,5), c(4,6))
    labels <- c("T4vsT5", "T4vsT6")
  } else {
    pairs  <- list(c(1,2), c(1,3))
    labels <- c("T1vsT2", "T1vsT3")
  }

  map2_dfr(pairs, labels, function(pr, lbl) {
    base_tp <- pr[1]; end_tp <- pr[2]

    tmao_fc <- tmao_wide |>
      mutate(tmao_fc = .data[[paste0("tp", end_tp)]] - .data[[paste0("tp", base_tp)]]) |>
      select(Subject, tmao_fc)

    map_dfr(BILE_IDS, function(aid) {
      get_fc(aid, base_tp, end_tp) |>
        rename(bile_fc = fc) |>
        inner_join(tmao_fc, by = "Subject") |>
        mutate(Alignment_ID = aid,
               bile_name    = BILE_NAMES[as.character(aid)],
               tp_pair      = lbl)
    })
  }) |>
    left_join(groups |> select(Subject, group), by = "Subject") |>
    mutate(bile_name = factor(bile_name, levels = BILE_NAMES))
}

long_Y <- build_long("Y")
long_X <- build_long("X")

# ── ヘルパー：Pearson r ラベル文字列 ─────────────────────────────────────────
r_label <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return("n.a.")
  ct  <- cor.test(x[ok], y[ok])
  r   <- unname(ct$estimate)
  p   <- ct$p.value
  sig <- if (p < 0.05) "*" else ""
  ps  <- if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
  sprintf("r=%.2f\n%s%s", r, ps, sig)
}

# ──────────────────────────────────────────────────────────────────────────────
# 図 1: H1A_boxplot.png
# グループ別 Y 期 log2FC（胆汁酸別ファセット，各グループの中央値をボックス表示）
# ──────────────────────────────────────────────────────────────────────────────
# Y 期 mean FC per subject (T4vsT5, T4vsT6 平均)
box_df <- long_Y |>
  group_by(Subject, Alignment_ID, bile_name, group) |>
  summarise(fc = mean(bile_fc, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(group))

p_box <- ggplot(box_df, aes(x = group, y = fc, fill = group)) +
  geom_hline(yintercept = 0, color = "#CCCCCC", lwd = 0.5) +
  geom_boxplot(outlier.size = 0.8, lwd = 0.4, width = 0.65, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "#333333") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(limits = GROUP_ORDER) +
  facet_wrap(~bile_name, ncol = 4, scales = "free_y") +
  labs(
    title    = "H1-A: Secondary Bile Acids by TMAO Group (Y period log₂FC)",
    subtitle = "Y period: mean of T4vsT5 and T4vsT6 per subject",
    x = NULL, y = "log₂FC (Y period mean)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 6),
    strip.text       = element_text(size = 6.5, face = "bold"),
    strip.background = element_rect(fill = "#F0F0F0"),
    plot.title       = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_FIG, "H1A_boxplot.png"),
       p_box, width = 13, height = 11, dpi = 150)
message("Saved: H1A_boxplot.png")

# ──────────────────────────────────────────────────────────────────────────────
# 図 2: H1A_heatmap.png
# 被験者（行）× 胆汁酸（列）ヒートマップ，グループでソート
# ──────────────────────────────────────────────────────────────────────────────
heat_mat <- box_df |>
  pivot_wider(id_cols = c(Subject, group), names_from = bile_name,
              values_from = fc, values_fill = NA) |>
  arrange(match(as.character(group), GROUP_ORDER))

row_annot <- data.frame(Group = heat_mat$group, row.names = as.character(heat_mat$Subject))
annot_colors <- list(Group = setNames(GROUP_COLORS[levels(heat_mat$group)],
                                       levels(heat_mat$group)))

mat <- as.matrix(heat_mat |> select(-Subject, -group))
rownames(mat) <- as.character(heat_mat$Subject)

# クリップ（見やすさのため）
mat_clip <- pmax(pmin(mat, 5), -5)

png(file.path(OUT_FIG, "H1A_heatmap.png"), width = 1600, height = 1000, res = 120)
pheatmap(mat_clip,
         annotation_row    = row_annot,
         annotation_colors = annot_colors,
         cluster_rows      = FALSE,
         cluster_cols      = TRUE,
         color             = colorRampPalette(c("#4472C4","white","#C00000"))(100),
         breaks            = seq(-5, 5, length.out = 101),
         fontsize          = 8,
         fontsize_row      = 7,
         fontsize_col      = 8,
         main              = "H1-A: Bile Acid log₂FC (Y period) by Subject  [clipped at ±5]",
         border_color      = NA)
dev.off()
message("Saved: H1A_heatmap.png")

# ──────────────────────────────────────────────────────────────────────────────
# 図 3: H1A_scatter_all.png
# Y 期平均 log2FC，グループ色分け，16 種散布図（4×4）
# ──────────────────────────────────────────────────────────────────────────────
scatter_avg <- long_Y |>
  group_by(Subject, Alignment_ID, bile_name, group) |>
  summarise(bile_fc  = mean(bile_fc,  na.rm = TRUE),
            tmao_fc  = mean(tmao_fc,  na.rm = TRUE),
            .groups  = "drop")

r_df_all <- scatter_avg |>
  group_by(bile_name) |>
  summarise(rlabel = r_label(bile_fc, tmao_fc), .groups = "drop") |>
  left_join(
    scatter_avg |>
      group_by(bile_name) |>
      summarise(x_pos = quantile(bile_fc, 0.98, na.rm=TRUE),
                y_pos = quantile(tmao_fc, 0.98, na.rm=TRUE),
                .groups = "drop"),
    by = "bile_name"
  )

p_all <- ggplot(scatter_avg, aes(x = bile_fc, y = tmao_fc, color = group)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "#444444", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "#CCCCCC", lwd = 0.4) +
  geom_vline(xintercept = 0, color = "#CCCCCC", lwd = 0.4) +
  geom_text(data = r_df_all, aes(x = x_pos, y = y_pos, label = rlabel),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.5, color = "#333333") +
  scale_color_manual(values = GROUP_COLORS, breaks = GROUP_ORDER, name = "Group") +
  facet_wrap(~bile_name, ncol = 4, scales = "free") +
  labs(
    title    = "H1-A: Secondary Bile Acids vs TMAO (Y period average)",
    subtitle = "Mean of T4vsT5 and T4vsT6 per subject  |  Pooled Pearson r shown",
    x = "Bile acid log₂FC (Y mean)", y = "TMAO log₂FC (Y mean, ID 414+415)"
  ) +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 6.5, face = "bold"),
        strip.background = element_rect(fill = "#F0F0F0"),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 7),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_FIG, "H1A_scatter_all.png"),
       p_all, width = 12, height = 10, dpi = 150)
message("Saved: H1A_scatter_all.png")

# ──────────────────────────────────────────────────────────────────────────────
# 図 4/5: scatter_tp_split (Y 期 / X 期)
# T4vsT5 = 実線●, T4vsT6 = △ を重ねて 16 種散布図，全被験者 pooled で Pearson r
# ──────────────────────────────────────────────────────────────────────────────
make_tpsplit <- function(long_df, period_lbl, tp1_lbl, tp2_lbl, pt_color = "#C00000") {
  r_df <- long_df |>
    group_by(bile_name) |>
    summarise(rlabel = r_label(bile_fc, tmao_fc),
              x_pos  = quantile(bile_fc, 0.98, na.rm=TRUE),
              y_pos  = quantile(tmao_fc, 0.98, na.rm=TRUE),
              .groups = "drop")

  ggplot(long_df, aes(x = bile_fc, y = tmao_fc, shape = tp_pair)) +
    geom_point(color = pt_color, size = 1.8, alpha = 0.75) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "#444444", linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "#CCCCCC", lwd = 0.4) +
    geom_vline(xintercept = 0, color = "#CCCCCC", lwd = 0.4) +
    geom_text(data = r_df, aes(x = x_pos, y = y_pos, label = rlabel),
              inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.5, color = "#333333") +
    scale_shape_manual(values = c(16, 17),
                       labels = c(tp1_lbl, tp2_lbl),
                       name   = "TP pair") +
    facet_wrap(~bile_name, ncol = 4, scales = "free") +
    labs(
      title    = sprintf("H1-A: Secondary Bile Acids vs TMAO (%s)", period_lbl),
      subtitle = sprintf("%s = ● (filled circle)  |  %s = ▲ (triangle)  |  Pearson r pooled (all subjects × 2 TP)", tp1_lbl, tp2_lbl),
      x = "Bile acid log₂FC", y = "TMAO log₂FC (mean ID 414+415)"
    ) +
    theme_bw(base_size = 8) +
    theme(strip.text = element_text(size = 6.5, face = "bold"),
          strip.background = element_rect(fill = "#F0F0F0"),
          plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 7.5, color = "#555555"),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
}

p_tpY <- make_tpsplit(long_Y, "Y period", "T4vsT5", "T4vsT6", "#C00000")
ggsave(file.path(OUT_FIG, "H1A_scatter_tp_split.png"),
       p_tpY, width = 12, height = 10, dpi = 150)
message("Saved: H1A_scatter_tp_split.png")

p_tpX <- make_tpsplit(long_X, "X period", "T1vsT2", "T1vsT3", "#2E75B6")
ggsave(file.path(OUT_FIG, "H1A_scatter_tp_split_Xperiod.png"),
       p_tpX, width = 12, height = 10, dpi = 150)
message("Saved: H1A_scatter_tp_split_Xperiod.png")

# ── 相関係数テーブル保存 ───────────────────────────────────────────────────────
cor_table <- bind_rows(
  long_Y |> group_by(bile_name, Alignment_ID) |>
    group_modify(~{
      ok <- is.finite(.x$bile_fc) & is.finite(.x$tmao_fc)
      r45 <- cor.test(.x$bile_fc[.x$tp_pair=="T4vsT5" & ok], .x$tmao_fc[.x$tp_pair=="T4vsT5" & ok])
      r46 <- cor.test(.x$bile_fc[.x$tp_pair=="T4vsT6" & ok], .x$tmao_fc[.x$tp_pair=="T4vsT6" & ok])
      ry  <- cor.test(.x$bile_fc[ok], .x$tmao_fc[ok])
      data.frame(r_T4vsT5 = round(r45$estimate,3),
                 r_T4vsT6 = round(r46$estimate,3),
                 r_Y_pooled = round(ry$estimate,3),
                 p_Y_pooled = round(ry$p.value,4))
    })
) |>
  arrange(desc(r_Y_pooled))

write_csv(cor_table, file.path(PROJ_ROOT, "output/H1A_correlation_table.csv"))
message("Saved: H1A_correlation_table.csv")
message("05_H1A_scatter.R 完了。")
