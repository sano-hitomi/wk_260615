# ==============================================================================
# 04_cytoscape_export.R
# Purpose : Build Cytoscape-ready node and edge tables from pathway analysis
#           results and the filtered feature table.
#
# Output format
#   Node table (metabolites)
#     id, name, hmdb_id, log2FC_X, log2FC_Y, padj_X, padj_Y,
#     sig_X, sig_Y, classification, primary_pathway
#
#   Node table (pathways)
#     id, name, pval_X, pval_Y, impact_X, impact_Y, hits_X, hits_Y, category
#
#   Edge table (metabolite → pathway membership)
#     source (HMDB ID), target (Pathway name), interaction ("member of")
#
#   SIF file (minimal network for Cytoscape import)
#
# Cytoscape import workflow:
#   1. File → Import → Network from File → select .sif
#   2. File → Import → Table from File → select node/edge table CSVs
#   3. Map node attributes (log2FC_X / log2FC_Y) to node Fill Color via Style
# ==============================================================================

library(dplyr)
library(readr)
library(tidyr)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
PROCESSED_DIR  <- "data/production/processed"
PATHWAY_DIR    <- "output/pathway"
OUT_DIR        <- "output/cytoscape"

PADJ_THRESHOLD <- 0.05   # used to flag sig_X / sig_Y on metabolite nodes

# Which pathways to include in the network?
# Options: "sig_either" (p<0.05 in X or Y), "sig_both", "all_tested"
PATHWAY_FILTER <- "sig_either"

# Maximum number of pathways in the network (top by min p-value)
MAX_PATHWAYS   <- 50

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
features <- read_csv(file.path(PROCESSED_DIR, "filtered_features.csv"),
                     show_col_types = FALSE)

path_comp <- read_csv(file.path(PATHWAY_DIR, "pathway_comparison_XY.csv"),
                      show_col_types = FALSE)

# Also load per-perturbation hit lists from MetaboAnalystR if available
# MetaboAnalystR stores pathway-to-compound mapping in mSet$dataSet$path.hits
# Since we are working post-run, we reconstruct edges from KEGG data via
# the pathway results tables and the metabolite list.

# If MetaboAnalystR result .rda files were saved, load them here:
# load(file.path(PATHWAY_DIR, "pertX", "mSet_pertX.rda"))  # optional

# ------------------------------------------------------------------------------
# 1. Metabolite node table
# ------------------------------------------------------------------------------
cat("Building metabolite node table...\n")

met_nodes <- features %>%
  transmute(
    id              = hmdb_id,
    node_type       = "metabolite",
    name            = Metabolite_name,
    hmdb_id         = hmdb_id,
    classification  = classification,
    primary_pathway = primary_pathway,
    log2FC_X        = log2FC_X,
    log2FC_Y        = log2FC_Y,
    padj_X          = padj_X,
    padj_Y          = padj_Y,
    sig_X           = padj_X < PADJ_THRESHOLD,
    sig_Y           = padj_Y < PADJ_THRESHOLD,
    sig_category    = case_when(
      sig_X &  sig_Y ~ "Both",
      sig_X & !sig_Y ~ "X only",
     !sig_X &  sig_Y ~ "Y only",
      TRUE            ~ "Not significant"
    ),
    # Absolute mean fold change: useful for node size mapping
    mean_abs_log2FC = (abs(log2FC_X) + abs(log2FC_Y)) / 2
  )

# ------------------------------------------------------------------------------
# 2. Pathway node table
# ------------------------------------------------------------------------------
cat("Building pathway node table...\n")

# Filter pathways to include in network
path_nodes_all <- path_comp %>%
  filter(!is.na(pval_X) | !is.na(pval_Y)) %>%
  mutate(
    id        = paste0("pathway:", Pathway),
    node_type = "pathway",
    name      = Pathway,
    min_pval  = pmin(pval_X, pval_Y, na.rm = TRUE)
  )

path_nodes_filtered <- switch(PATHWAY_FILTER,
  "sig_either" = path_nodes_all %>% filter(sig_X | sig_Y),
  "sig_both"   = path_nodes_all %>% filter(sig_X & sig_Y),
  "all_tested" = path_nodes_all
)

path_nodes <- path_nodes_filtered %>%
  arrange(min_pval) %>%
  slice_head(n = MAX_PATHWAYS) %>%
  transmute(
    id           = id,
    node_type    = node_type,
    name         = name,
    pval_X       = pval_X,
    pval_Y       = pval_Y,
    impact_X     = impact_X,
    impact_Y     = impact_Y,
    hits_X       = Hits_X,
    hits_Y       = Hits_Y,
    sig_category = category,
    neg_log10p_X = -log10(pval_X + 1e-10),
    neg_log10p_Y = -log10(pval_Y + 1e-10)
  )

cat("Pathway nodes:", nrow(path_nodes), "\n")

# ------------------------------------------------------------------------------
# 3. Edge table: metabolite → pathway
#
# MetaboAnalystR stores per-pathway hit lists in mSet$analSet$hits (a named list
# mapping pathway name → vector of matched compound IDs). We attempt to load the
# saved MetaboAnalystR mSet objects; if not available, we use primary_pathway
# from the annotation as a fallback.
# ------------------------------------------------------------------------------
cat("Building edge table...\n")

build_edges_from_mset <- function(rda_path, label) {
  if (!file.exists(rda_path)) return(NULL)
  env <- new.env()
  load(rda_path, envir = env)
  mSet <- get(ls(env)[1], envir = env)
  if (is.null(mSet$analSet$hits)) return(NULL)

  lapply(names(mSet$analSet$hits), function(pw) {
    hits <- mSet$analSet$hits[[pw]]
    if (length(hits) == 0) return(NULL)
    data.frame(
      source      = hits,
      target      = paste0("pathway:", pw),
      interaction = "member of",
      perturbation = label,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()
}

edges_X <- build_edges_from_mset(
  file.path(PATHWAY_DIR, "pertX", "mSet_pertX.rda"), "X")
edges_Y <- build_edges_from_mset(
  file.path(PATHWAY_DIR, "pertY", "mSet_pertY.rda"), "Y")

if (!is.null(edges_X) || !is.null(edges_Y)) {
  edges <- bind_rows(edges_X, edges_Y) %>% distinct(source, target, .keep_all = TRUE)
  cat("Edges from MetaboAnalystR mSet:", nrow(edges), "\n")
} else {
  # Fallback: use primary_pathway annotation
  cat("mSet .rda not found — using primary_pathway annotation as fallback edges\n")
  edges <- met_nodes %>%
    filter(!is.na(primary_pathway), primary_pathway != "") %>%
    mutate(
      target      = paste0("pathway:", primary_pathway),
      interaction = "member of"
    ) %>%
    select(source = id, target, interaction) %>%
    # Keep only edges to pathways in our filtered set
    filter(target %in% path_nodes$id)
  cat("Fallback edges:", nrow(edges), "\n")
}

# Keep only metabolite nodes that appear in at least one edge
connected_met_ids <- unique(edges$source)
met_nodes_connected <- met_nodes %>%
  filter(id %in% connected_met_ids)

cat("Connected metabolite nodes:", nrow(met_nodes_connected), "\n")

# ------------------------------------------------------------------------------
# 4. Combined node table (metabolites + pathways, Cytoscape style)
# ------------------------------------------------------------------------------
# Cytoscape expects a single node table or separate tables.
# We write them separately for clarity, with a shared "id" column.

all_nodes <- bind_rows(
  met_nodes_connected %>%
    select(id, node_type, name, classification, primary_pathway,
           log2FC_X, log2FC_Y, padj_X, padj_Y,
           sig_X, sig_Y, sig_category, mean_abs_log2FC),
  path_nodes %>%
    select(id, node_type, name, pval_X, pval_Y,
           impact_X, impact_Y, sig_category,
           neg_log10p_X, neg_log10p_Y)
)

# ------------------------------------------------------------------------------
# 5. Write output files
# ------------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Node tables
write_csv(met_nodes_connected,
          file.path(OUT_DIR, "nodes_metabolites.csv"))
write_csv(path_nodes,
          file.path(OUT_DIR, "nodes_pathways.csv"))
write_csv(all_nodes,
          file.path(OUT_DIR, "nodes_combined.csv"))

# Edge table
write_csv(edges,
          file.path(OUT_DIR, "edges.csv"))

# SIF file (source <interaction> target, tab-separated)
# Cytoscape: File → Import → Network from File
sif_lines <- paste(edges$source, edges$interaction, edges$target, sep = "\t")
writeLines(sif_lines, file.path(OUT_DIR, "network.sif"))

cat("\n--- Cytoscape export complete ---\n")
cat("  Metabolite nodes :", nrow(met_nodes_connected), "\n")
cat("  Pathway nodes    :", nrow(path_nodes), "\n")
cat("  Edges            :", nrow(edges), "\n")
cat("  Output dir       :", OUT_DIR, "\n")
cat("\nImport into Cytoscape:\n")
cat("  1. File → Import → Network from File → network.sif\n")
cat("  2. File → Import → Table from File → nodes_combined.csv (Key: id)\n")
cat("  3. File → Import → Table from File → edges.csv\n")
cat("  4. Style panel: map log2FC_X / log2FC_Y to node Fill Color\n")
cat("  5. Style panel: map mean_abs_log2FC to node Size\n")
cat("  6. Style panel: map neg_log10p_X / neg_log10p_Y to pathway node Size\n")
