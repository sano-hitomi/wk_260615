# ==============================================================================
# 02_enrichment_ora.R
# Purpose : Metabolite Set Over-Representation Analysis (ORA) using
#           MetaboAnalystR, run separately for Perturbation X and Y.
# Inputs  : data/production/processed/input_ora_X.csv
#           data/production/processed/input_ora_Y.csv
#           data/production/processed/filtered_features.csv  (background)
# Outputs : output/enrichment/pertX/  — ORA results + plots
#           output/enrichment/pertY/  — ORA results + plots
#
# NOTE: Cairo must be loaded BEFORE MetaboAnalystR to avoid a recursive
#       default.dpi evaluation error on macOS (Apple Silicon).
# ==============================================================================

library(Cairo)          # must come first on macOS
library(MetaboAnalystR)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------------------------
# Patch MetaboAnalystR: fix recursive default argument reference for default.dpi
# Find any internal function whose formal argument "default.dpi" defaults to
# itself, and replace that default with 96L.
# ------------------------------------------------------------------------------
local({
  ns <- asNamespace("MetaboAnalystR")
  for (fn_name in ls(ns)) {
    fn <- tryCatch(get(fn_name, envir = ns), error = function(e) NULL)
    if (!is.function(fn)) next
    if (!"default.dpi" %in% names(formals(fn))) next
    tryCatch({
      formals(fn)$default.dpi <- 96L
      unlockBinding(fn_name, ns)
      assign(fn_name, fn, envir = ns)
      cat("[patch] Fixed default.dpi in:", fn_name, "\n")
    }, error = function(e) {
      cat("[patch] Could not patch", fn_name, ":", conditionMessage(e), "\n")
    })
  }
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
PROCESSED_DIR  <- "."
OUT_BASE       <- "."

# MetaboAnalyst metabolite set library options:
#   "smpdb_pathway"  — SMPDB pathway sets (recommended for metabolomics)
#   "kegg_pathway"   — KEGG pathway sets
#   "hmdb_disease"   — HMDB disease associations
MSET_LIBRARY   <- "smpdb_pathway"

# Minimum metabolite set size to include in testing
MIN_SET_SIZE   <- 2

# Use MetaboAnalystR's built-in human metabolome as background.
# Custom background (USE_CUSTOM_BG = TRUE) requires Setup.HMDBReferenceMetabolome
# which depends on file paths from MetaboAnalyst web server and can fail on macOS.
USE_CUSTOM_BG  <- FALSE

# ------------------------------------------------------------------------------
# Helper: run ORA for one perturbation
# ------------------------------------------------------------------------------
run_ora <- function(ora_ids,
                    bg_ids   = NULL,
                    label    = "pertX",
                    out_base = OUT_BASE,
                    lib      = MSET_LIBRARY) {

  out_dir <- file.path(out_base, label)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  orig_wd <- getwd()
  on.exit(setwd(orig_wd), add = TRUE)
  setwd(out_dir)

  cat("\n========== ORA:", label, "==========\n")
  cat("Input metabolites:", length(ora_ids), "\n")

  mSet <- InitDataObjects("conc", "msetora", FALSE, default.dpi = 96L)
  mSet <- Setup.MapData(mSet, ora_ids)
  mSet <- CrossReferencing(mSet, "hmdb")
  mSet <- CreateMappingResultTable(mSet)

  if (!is.null(bg_ids)) {
    mSet <- SetMetabolomeFilter(mSet, TRUE)
    mSet <- Setup.HMDBReferenceMetabolome(mSet, bg_ids)
  } else {
    mSet <- SetMetabolomeFilter(mSet, FALSE)
  }

  mSet <- SetCurrentMsetLib(mSet, lib, MIN_SET_SIZE)

  # Function name changed across MetaboAnalystR versions:
  #   ≤ 4.0.x : CalculateHypergeoScore
  #   ≥ 4.3.x : CalculateHyperScore
  if (exists("CalculateHyperScore", envir = asNamespace("MetaboAnalystR"), inherits = FALSE)) {
    mSet <- CalculateHyperScore(mSet)
  } else {
    mSet <- CalculateHypergeoScore(mSet)
  }

  # ---- MetaboAnalystR native plots (Cairo required) ----
  tryCatch({
    mSet <- PlotORA(mSet,
                    imgName  = paste0("ora_bar_", label),
                    plotType = "bar",
                    format   = "png",
                    dpi      = 150,
                    width    = NA)
    mSet <- PlotEnrichDotPlot(mSet,
                              imgName = paste0("ora_dot_", label),
                              format  = "png",
                              dpi     = 150,
                              width   = NA)
  }, error = function(e) {
    cat("  [note] MetaboAnalystR native plots skipped (Cairo unavailable):",
        conditionMessage(e), "\n")
  })

  # ---- Save results table ----
  results <- as.data.frame(mSet$analSet$ora.mat)
  results <- results %>%
    tibble::rownames_to_column("MetaboliteSet") %>%
    arrange(`Raw p`)

  write_csv(results, paste0("ora_results_", label, ".csv"))
  cat("Saved:", nrow(results), "tested sets →",
      sum(results$`Raw p` < 0.05, na.rm = TRUE), "nominally significant\n")

  mapping <- mSet$dataSet$map.table
  if (!is.null(mapping)) {
    write_csv(as.data.frame(mapping), paste0("id_mapping_", label, ".csv"))
  }

  return(invisible(mSet))
}

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------
ora_X <- read_csv(file.path(PROCESSED_DIR, "input_ora_X.csv"),
                  show_col_types = FALSE)
ora_Y <- read_csv(file.path(PROCESSED_DIR, "input_ora_Y.csv"),
                  show_col_types = FALSE)

if (USE_CUSTOM_BG) {
  all_features <- read_csv(file.path(PROCESSED_DIR, "filtered_features.csv"),
                            show_col_types = FALSE)
  bg_ids <- all_features$hmdb_id
} else {
  bg_ids <- NULL
}

cat("Background metabolome size:", length(bg_ids), "\n")

# ------------------------------------------------------------------------------
# Run ORA
# ------------------------------------------------------------------------------
mSet_ora_X <- run_ora(
  ora_ids  = ora_X$hmdb_id,
  bg_ids   = bg_ids,
  label    = "pertX",
  out_base = OUT_BASE,
  lib      = MSET_LIBRARY
)

mSet_ora_Y <- run_ora(
  ora_ids  = ora_Y$hmdb_id,
  bg_ids   = bg_ids,
  label    = "pertY",
  out_base = OUT_BASE,
  lib      = MSET_LIBRARY
)

# ------------------------------------------------------------------------------
# Comparative summary: X vs Y
# ------------------------------------------------------------------------------
cat("\n========== Comparative Summary ==========\n")

res_X <- read_csv(file.path(OUT_BASE, "pertX", "ora_results_pertX.csv"),
                  show_col_types = FALSE)
res_Y <- read_csv(file.path(OUT_BASE, "pertY", "ora_results_pertY.csv"),
                  show_col_types = FALSE)

comparison <- res_X %>%
  select(MetaboliteSet, pval_X = `Raw p`) %>%
  full_join(
    res_Y %>% select(MetaboliteSet, pval_Y = `Raw p`),
    by = "MetaboliteSet"
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
  arrange(pmin(pval_X, pval_Y, na.rm = TRUE))

write_csv(comparison, file.path(OUT_BASE, "ora_comparison_XY.csv"))
cat("Sets enriched in X only :", sum(comparison$category == "X only",  na.rm = TRUE), "\n")
cat("Sets enriched in Y only :", sum(comparison$category == "Y only",  na.rm = TRUE), "\n")
cat("Sets enriched in both   :", sum(comparison$category == "Both",    na.rm = TRUE), "\n")

# ------------------------------------------------------------------------------
# Comparison dot plot (ggplot2)
# ------------------------------------------------------------------------------
sig_both <- comparison %>%
  filter(category %in% c("Both", "X only", "Y only")) %>%
  slice_head(n = 30)

if (nrow(sig_both) > 0) {
  plot_data <- sig_both %>%
    tidyr::pivot_longer(cols = c(pval_X, pval_Y),
                        names_to  = "Perturbation",
                        values_to = "pval") %>%
    mutate(
      Perturbation  = recode(Perturbation,
                             pval_X = "Perturbation X",
                             pval_Y = "Perturbation Y"),
      neg_log10p    = -log10(pval + 1e-10),
      MetaboliteSet = factor(MetaboliteSet,
                             levels = rev(unique(sig_both$MetaboliteSet)))
    )

  p <- ggplot(plot_data,
              aes(x = Perturbation, y = MetaboliteSet,
                  size = neg_log10p, color = neg_log10p)) +
    geom_point() +
    scale_color_gradient(low = "#bdd7e7", high = "#08519c",
                         name = "-log10(p)") +
    scale_size_continuous(range = c(1, 8), name = "-log10(p)") +
    labs(title = "ORA: Enriched Metabolite Sets (X vs Y)",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(OUT_BASE, "ora_dotplot_comparison.png"),
         plot = p,
         width = 7,
         height = max(4, nrow(sig_both) * 0.35),
         dpi = 150)
  cat("Saved comparison dot plot →",
      file.path(OUT_BASE, "ora_dotplot_comparison.png"), "\n")
}
