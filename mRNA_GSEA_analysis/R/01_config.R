# =============================================================================
# 01_config.R - Configuration and Parameters for mRNA GSEA Analysis
# =============================================================================
# This module defines all configuration parameters and global settings
# for the CSN subunit mRNA GSEA analysis pipeline.
# =============================================================================

# -----------------------------------------------------------------------------
# Random Seed
# -----------------------------------------------------------------------------
set.seed(1234)

# -----------------------------------------------------------------------------
# CSN Subunits Definition
# -----------------------------------------------------------------------------
csn_subunits <- c(
    "GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6",
    "COPS7A", "COPS7B", "COPS8", "COPS9"
)

# -----------------------------------------------------------------------------
# Analysis Parameters
# -----------------------------------------------------------------------------

# Data quality thresholds
min_frac_complete <- 0.75 # At least 75% non-NA per gene

# GSEA parameters
minSize <- 5 # Minimum gene set size
maxSize <- 1000 # Maximum gene set size
fgsea_eps <- 1e-10 # fgseaMultilevel() precision

# Limma parameters
min_per_group <- 8 # Minimum samples for each group (skip limma if insufficient)

# -----------------------------------------------------------------------------
# User Configuration - Analysis Passes
# -----------------------------------------------------------------------------
# Which passes to run:
#   - Only CovarAdj: RUN_PASSES <- c("CovarAdj")    (default)
#   - Only RAW:      RUN_PASSES <- "RAW"
#   - Run both:      RUN_PASSES <- c("RAW","CovarAdj")
RUN_PASSES <- c("CovarAdj")

# Write to options for getOption(...) to read elsewhere
options(csn.run_passes = RUN_PASSES)

# -----------------------------------------------------------------------------
# User Configuration - Heatmap Limits
# -----------------------------------------------------------------------------
# Heatmap filtering threshold for non-H collections
# - Single dataset: keep only top N and bottom M NES with padj<0.05 (based on CSN_SCORE)
# - PAN (meta-FDR): keep only top N and bottom M Z with padj_meta<0.05 (based on CSN_SCORE)
DATASET_HEATMAP_TOP_N <- 25 # Default top 25 (NES)
DATASET_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (NES)
PAN_HEATMAP_TOP_N <- 25 # Default top 25 (Z)
PAN_HEATMAP_BOTTOM_N <- 25 # Default bottom 25 (Z)

# -----------------------------------------------------------------------------
# User Configuration - Heatmap Collections
# -----------------------------------------------------------------------------
# Which collections to draw heatmaps for (string vector; e.g.: c("H", "C6", "C2:CP:BIOCARTA"))
# - Single dataset heatmap: leave empty or NULL means "draw all available collections"
# - PAN heatmap: leave empty or NULL means "draw all available collections"
PLOT_DATASET_COLLECTIONS <- NULL
PLOT_PAN_COLLECTIONS <- NULL

# -----------------------------------------------------------------------------
# Gene Set Groups to Run
# -----------------------------------------------------------------------------
GENESET_GROUPS_TO_RUN <- c("H")

# -----------------------------------------------------------------------------
# Output Prefix
# -----------------------------------------------------------------------------
OUTPUT_PREFIX <- "mRNA_GSEA"

# -----------------------------------------------------------------------------
# TP53 Strata Definition
# -----------------------------------------------------------------------------
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")

# -----------------------------------------------------------------------------
# Batch Cleanup Settings
# -----------------------------------------------------------------------------
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": values containing '|' set to NA or merged
BATCH_MIN_PER_LEVEL <- 2 # Levels below this threshold merged to b_small

# -----------------------------------------------------------------------------
# CSN Subunit Order for TPM Display
# -----------------------------------------------------------------------------
CSN_SUB_ORDER <- c(
    "GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
    "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9"
)
