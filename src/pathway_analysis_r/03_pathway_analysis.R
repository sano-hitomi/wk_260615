# ==============================================================================
# 03_pathway_analysis.R
# Purpose : KEGG pathway ORA using KEGGREST + base R (MetaboAnalystR-free)
#           Avoids MetaboAnalystR crashes on Apple Silicon macOS.
# Inputs  : filtered_features.csv  (columns: hmdb_id, padj_X, padj_Y)
# Outputs : pertX/pathway_results_pertX.csv
#           pertY/pathway_results_pertY.csv
#           pathway_comparison_XY.csv
#           bubble charts (PNG)
# ==============================================================================

suppressPackageStartupMessages({
  library(KEGGREST)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
PROCESSED_DIR     <- "."
OUT_BASE          <- "."
KEGG_ORG          <- "hsa"
PVAL_LABEL_CUTOFF <- 0.1
PADJ_THRESHOLD    <- 0.05
CACHE_DIR         <- file.path(OUT_BASE, ".kegg_cache")

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helper: KEGG API with simple caching
# ------------------------------------------------------------------------------
cached_kegg <- function(key, expr) {
  cache_file <- file.path(CACHE_DIR, paste0(key, ".rds"))
  if (file.exists(cache_file)) {
    cat("  [cache] Loading", key, "\n")
    return(readRDS(cache_file))
  }
  result <- expr
  saveRDS(result, cache_file)
  result
}

# ------------------------------------------------------------------------------
# Step 1: HMDB → KEGG compound ID mapping
# ------------------------------------------------------------------------------
get_hmdb_kegg_map <- function() {
  cat("Fetching HMDB→KEGG mapping...\n")
  cached_kegg("hmdb_kegg_map", {
    map <- tryCatch(
      keggConv("compound", "hmdb"),
      error = function(e) { cat("  Error:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(map) || length(map) == 0) return(data.frame(hmdb_id=character(), kegg_cid=character()))
    data.frame(
      hmdb_id  = sub("hmdb:", "", names(map)),
      kegg_cid = sub("cpd:",  "", unname(map)),
      stringsAsFactors = FALSE
    )
  })
}

# ------------------------------------------------------------------------------
# Step 2: KEGG pathway → compound mapping (human) using keggGet
# ------------------------------------------------------------------------------
get_pathway_compound_map <- function(org = "hsa") {
  cat("Fetching KEGG pathway list for", org, "...\n")

  pw_list <- cached_kegg(paste0("pathway_list_", org), {
    tryCatch(keggList("pathway", org), error = function(e) NULL)
  })
  if (is.null(pw_list)) stop("Could not fetch pathway list from KEGG")

  # pathway IDの形式: "path:hsa00010" → "hsa00010"
  pw_ids   <- sub("path:", "", names(pw_list))
  pw_names_vec <- setNames(unname(pw_list), pw_ids)
  cat("  ", length(pw_ids), "pathways found\n")

  cat("Fetching pathway details via keggGet (batched, ~", ceiling(length(pw_ids)/10), "calls)...\n")

  pw_compound_df <- cached_kegg(paste0("pathway_compound_keggget_", org), {
    batches <- split(pw_ids, ceiling(seq_along(pw_ids) / 10))
    rows <- list()

    for (i in seq_along(batches)) {
      if (i %% 10 == 1) cat("  batch", i, "/", length(batches), "\n")
      res <- tryCatch(keggGet(batches[[i]]), error = function(e) NULL)
      if (is.null(res)) { Sys.sleep(0.5); next }

      for (pw in res) {
        cpds <- pw$COMPOUND   # named char: name = KEGG cpd ID, value = compound name
        if (is.null(cpds) || length(cpds) == 0) next
        rows[[length(rows) + 1]] <- data.frame(
          pathway_id   = pw$ENTRY,
          pathway_name = if (!is.null(pw$NAME)) pw$NAME[1] else pw$ENTRY,
          kegg_cid     = names(cpds),
          stringsAsFactors = FALSE
        )
      }
      Sys.sleep(0.3)
    }

    if (length(rows) == 0) return(data.frame())
    do.call(rbind, rows)
  })

  if (is.null(pw_compound_df) || nrow(pw_compound_df) == 0) {
    stop("Could not fetch pathway-compound data from KEGG")
  }

  cat("  ", n_distinct(pw_compound_df$pathway_id), "pathways,",
      n_distinct(pw_compound_df$kegg_cid), "unique compounds\n")
  pw_compound_df
}

# ------------------------------------------------------------------------------
# Step 3: ORA (hypergeometric test)
# ------------------------------------------------------------------------------
run_ora <- function(hit_hmdb_ids, bg_hmdb_ids, hmdb_kegg_map, pathway_compound_map,
                    label = "pert") {
  cat("\n========== Pathway ORA:", label, "==========\n")

  # HMDB → KEGG
  hit_kegg <- hmdb_kegg_map %>% filter(hmdb_id %in% hit_hmdb_ids) %>% pull(kegg_cid) %>% unique()
  bg_kegg  <- hmdb_kegg_map %>% filter(hmdb_id %in% bg_hmdb_ids)  %>% pull(kegg_cid) %>% unique()

  cat("  Hit:", length(hit_hmdb_ids), "HMDB →", length(hit_kegg), "KEGG\n")
  cat("  BG :", length(bg_hmdb_ids),  "HMDB →", length(bg_kegg),  "KEGG\n")

  if (length(hit_kegg) == 0) {
    cat("  [skip] No KEGG IDs mapped for hit set\n")
    return(data.frame())
  }

  N <- length(bg_kegg)
  k <- length(hit_kegg)

  results <- pathway_compound_map %>%
    group_by(pathway_id, pathway_name) %>%
    summarise(
      Total    = n_distinct(kegg_cid),
      Expected = round(k * Total / N, 3),
      Hits     = sum(kegg_cid %in% hit_kegg),
      .groups  = "drop"
    ) %>%
    filter(Hits > 0) %>%
    mutate(
      Raw.p  = phyper(Hits - 1, Total, pmax(N - Total, 0), k, lower.tail = FALSE),
      Holm.p = p.adjust(Raw.p, method = "holm"),
      FDR    = p.adjust(Raw.p, method = "fdr"),
      Impact = Hits / Total   # hit ratio (proxy for pathway impact)
    ) %>%
    rename(Pathway = pathway_name) %>%
    select(Pathway, pathway_id, Total, Expected, Hits, Raw.p, Holm.p, FDR, Impact) %>%
    arrange(Raw.p)

  cat("  Pathways tested      :", nrow(results), "\n")
  cat("  Nominally sig (p<.05):", sum(results$Raw.p < 0.05, na.rm = TRUE), "\n")
  results
}

# ------------------------------------------------------------------------------
# Step 4: Bubble chart
# ------------------------------------------------------------------------------
plot_bubble <- function(results, label, out_dir, pval_cutoff = PVAL_LABEL_CUTOFF) {
  if (is.null(results) || nrow(results) == 0) {
    cat("  [skip] No results to plot for", label, "\n")
    return(invisible(NULL))
  }

  plot_df <- results %>%
    filter(!is.na(Raw.p), !is.na(Impact)) %>%
    mutate(
      neg_log10p = -log10(Raw.p + 1e-10),
      lab        = if_else(Raw.p < pval_cutoff, Pathway, NA_character_)
    )

  p <- ggplot(plot_df, aes(x = Impact, y = neg_log10p,
                            size = Hits, color = neg_log10p)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_text_repel(aes(label = lab), size = 2.8, max.overlaps = 20,
                    segment.color = "grey60", na.rm = TRUE) +
    scale_color_gradient(low = "#fee8c8", high = "#b30000", name = "-log10(p)") +
    scale_size_continuous(range = c(2, 10), name = "Hits") +
    labs(title    = paste("Pathway ORA —", label),
         subtitle = paste0("KEGG (", KEGG_ORG, "), hypergeometric | dashed: p = 0.05"),
         x = "Hit ratio (Hits / Total)", y = "-log10(p-value)") +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40", size = 9))

  out_path <- file.path(out_dir, paste0("pathway_bubble_", label, ".png"))
  ggsave(out_path, plot = p, width = 8, height = 6, dpi = 150)
  cat("  Saved:", out_path, "\n")
  invisible(p)
}

# ==============================================================================
# Main
# ==============================================================================

# --- Load inputs ---
cat("\nLoading input data...\n")
all_features <- read_csv(file.path(PROCESSED_DIR, "filtered_features.csv"),
                         show_col_types = FALSE)
bg_ids    <- all_features$hmdb_id
sig_X_ids <- all_features %>% filter(padj_X < PADJ_THRESHOLD) %>% pull(hmdb_id)
sig_Y_ids <- all_features %>% filter(padj_Y < PADJ_THRESHOLD) %>% pull(hmdb_id)

cat("Background:", length(bg_ids), "| Sig X:", length(sig_X_ids),
    "| Sig Y:", length(sig_Y_ids), "\n")

# --- Fetch KEGG data (cached after first run) ---
hmdb_kegg_map       <- get_hmdb_kegg_map()
pathway_compound_map <- get_pathway_compound_map(KEGG_ORG)

# --- Run ORA ---
out_X <- file.path(OUT_BASE, "pertX")
out_Y <- file.path(OUT_BASE, "pertY")
dir.create(out_X, showWarnings = FALSE, recursive = TRUE)
dir.create(out_Y, showWarnings = FALSE, recursive = TRUE)

res_X <- run_ora(sig_X_ids, bg_ids, hmdb_kegg_map, pathway_compound_map, "pertX")
res_Y <- run_ora(sig_Y_ids, bg_ids, hmdb_kegg_map, pathway_compound_map, "pertY")

write_csv(if (nrow(res_X) > 0) res_X else data.frame(),
          file.path(out_X, "pathway_results_pertX.csv"))
write_csv(if (nrow(res_Y) > 0) res_Y else data.frame(),
          file.path(out_Y, "pathway_results_pertY.csv"))

# --- Bubble charts ---
cat("\nGenerating bubble charts...\n")
plot_bubble(res_X, "pertX", out_X)
plot_bubble(res_Y, "pertY", out_Y)

# --- Comparison ---
cat("\n========== Pathway Comparison: X vs Y ==========\n")
if (nrow(res_X) == 0 && nrow(res_Y) == 0) {
  cat("Both pertX and pertY have no results. Skipping comparison.\n")
  quit(save = "no", status = 0)
}

comparison <- full_join(
  if (nrow(res_X) > 0)
    res_X %>% select(Pathway, Total_X=Total, Hits_X=Hits, pval_X=Raw.p, impact_X=Impact)
  else data.frame(),
  if (nrow(res_Y) > 0)
    res_Y %>% select(Pathway, Total_Y=Total, Hits_Y=Hits, pval_Y=Raw.p, impact_Y=Impact)
  else data.frame(),
  by = "Pathway"
) %>%
  mutate(
    sig_X    = replace_na(pval_X < 0.05, FALSE),
    sig_Y    = replace_na(pval_Y < 0.05, FALSE),
    category = case_when(
      sig_X &  sig_Y ~ "Both",
      sig_X & !sig_Y ~ "X only",
     !sig_X &  sig_Y ~ "Y only",
      TRUE            ~ "Neither"
    )
  ) %>%
  arrange(pmin(replace_na(pval_X, 1), replace_na(pval_Y, 1)))

write_csv(comparison, file.path(OUT_BASE, "pathway_comparison_XY.csv"))
cat("X only:", sum(comparison$category == "X only"),
    "| Y only:", sum(comparison$category == "Y only"),
    "| Both:", sum(comparison$category == "Both"), "\n")

# Scatter plot
scatter_df <- comparison %>%
  filter(!is.na(pval_X), !is.na(pval_Y)) %>%
  mutate(lab = if_else(category != "Neither" &
                         pmin(pval_X, pval_Y) < 0.1, Pathway, NA_character_))

if (nrow(scatter_df) > 0) {
  p_scatter <- ggplot(scatter_df,
                      aes(x = -log10(pval_X + 1e-10),
                          y = -log10(pval_Y + 1e-10),
                          color = category)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_vline(xintercept = -log10(0.05), linetype="dashed",
               color="grey50", linewidth=0.4) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed",
               color="grey50", linewidth=0.4) +
    geom_text_repel(aes(label = lab), size=2.5, max.overlaps=15,
                    segment.color="grey70", na.rm=TRUE) +
    scale_color_manual(
      values = c("Both"="red", "X only"="steelblue",
                 "Y only"="forestgreen", "Neither"="grey70"),
      name = "Significance") +
    labs(title = "Pathway enrichment: X vs Y",
         x = "-log10(p) Perturbation X", y = "-log10(p) Perturbation Y") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_BASE, "pathway_scatter_XvsY.png"),
         plot = p_scatter, width = 7, height = 6, dpi = 150)
  cat("Saved scatter →", file.path(OUT_BASE, "pathway_scatter_XvsY.png"), "\n")
}

cat("\nDone.\n")
