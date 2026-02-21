# =============================================================================
# 01_config.R - Configuration and Global Parameters
# =============================================================================
# This file contains all configuration settings and global parameters for
# the phosphoproteomic gene-level GSEA analysis pipeline.
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

# For plotting order (if needed)
CSN_SUB_ORDER <- csn_subunits

# -----------------------------------------------------------------------------
# Analysis Parameters
# -----------------------------------------------------------------------------
min_frac_complete <- 0.75 # At least 75% non-NA per gene
minSize <- 5 # Minimum gene set size for GSEA
maxSize <- 1000 # Maximum gene set size for GSEA
fgsea_eps <- 1e-10 # fgseaMultilevel() precision

min_per_group <- 8 # Minimum samples per group (skip limma if insufficient)

# -----------------------------------------------------------------------------
# Pipeline Control Flags
# -----------------------------------------------------------------------------

# Which passes to run:
#   - BatchAdj only: RUN_PASSES <- c("BatchAdj")    (default)
#   - RAW only:      RUN_PASSES <- "RAW"
#   - Both:          RUN_PASSES <- c("RAW","BatchAdj")
RUN_PASSES <- c("BatchAdj")

# Write the above embedded settings to options for getOption(...) to read elsewhere
options(csn.run_passes = RUN_PASSES)

# Age missing indicator flag
USE_AGE_MISSING_INDICATOR <- FALSE

# Pipeline: limma_t only
.RUN_LIMMA <- TRUE

# -----------------------------------------------------------------------------
# Heatmap Configuration
# -----------------------------------------------------------------------------
# Heatmap filtering thresholds for non-H collections (single dataset uses NES, PAN uses Z; all filter by CSN_SCORE)
# - Single dataset: keep only top N and bottom M NES with padj<0.05 (by CSN_SCORE)
# - PAN (meta-FDR): keep only top N and bottom M Z with padj_meta<0.05 (by CSN_SCORE)
DATASET_HEATMAP_TOP_N <- 25 # default top 25 (NES)
DATASET_HEATMAP_BOTTOM_N <- 25 # default bottom 25 (NES)
PAN_HEATMAP_TOP_N <- 25 # default top 25 (Z)
PAN_HEATMAP_BOTTOM_N <- 25 # default bottom 25 (Z)

# Which collections to draw heatmaps for (string vector; e.g., c("H", "C6", "C2:CP:BIOCARTA"))
# - Single dataset heatmap: leave empty or NULL to draw all available collections
# - PAN heatmap: leave empty or NULL to draw all available collections
PLOT_DATASET_COLLECTIONS <- NULL
PLOT_PAN_COLLECTIONS <- NULL

# -----------------------------------------------------------------------------
# Gene Set Configuration
# -----------------------------------------------------------------------------
GENESET_GROUPS_TO_RUN <- c("H")

# -----------------------------------------------------------------------------
# Output Path Configuration
# -----------------------------------------------------------------------------
# Output prefix for phosphoproteomic analysis
OUTPUT_PREFIX <- "phosphoproteomic_gene_level_GSEA"

# TP53 stratification strata
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")
