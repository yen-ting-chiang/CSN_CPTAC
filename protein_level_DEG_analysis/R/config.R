## ============================================================================
## Configuration for Proteomic DEG Analysis
## CSN Subunits -> Differential Expression Gene (DEG) Analysis
## ============================================================================

# Package Dependencies ---------------------------------------------------
# Required packages (must be installed before running):
# tidyverse, data.table, janitor, glue, SummarizedExperiment,
# MultiAssayExperiment, S4Vectors, limma, imputeLCMD, msigdbr,
# fgsea, openxlsx, ComplexHeatmap, cowplot, matrixStats, yaml

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(janitor)
    library(glue)
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
    library(S4Vectors)
    library(limma)
    library(imputeLCMD)
    library(msigdbr)
    library(fgsea)
    library(openxlsx)
    library(ComplexHeatmap)
    library(cowplot)
    library(matrixStats)
    library(yaml)
})

# Global Settings --------------------------------------------------------

## Set random seed for reproducibility
set.seed(1234)

## Force single-threaded execution (improves reproducibility)
.force_serial_execution <- function() {
    options(mc.cores = 1L)
    if ("package:data.table" %in% search() &&
        exists("setDTthreads", asNamespace("data.table"))) {
        data.table::setDTthreads(1L)
    }
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::register(BiocParallel::SerialParam())
    }
    if (requireNamespace("foreach", quietly = TRUE)) {
        foreach::registerDoSEQ()
    }
    if (requireNamespace("future", quietly = TRUE)) {
        future::plan(future::sequential)
    }
    Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")
    options(.fgsea_nproc = 1L)
}
.force_serial_execution()

# Analysis Parameters ----------------------------------------------------

## CSN subunit genes
csn_subunits <- c(
    "GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
    "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9"
)

## Data quality thresholds
min_frac_complete <- 0.75 # Minimum fraction of non-NA values per gene
min_per_group <- 8 # Minimum samples per group for statistical tests

## CSN score calculation
csn_min_members <- 5L # Minimum CSN subunits required for PC1 calculation
csn_pca_min_samples <- 10L # Minimum samples required for PCA


# Dataset Configuration ---------------------------------------------------

## NOTE: Working directory should be set to CSN_CPTAC project root before
## sourcing main_analysis.R. See main_analysis.R for details.

## Dataset IDs for proteomic DEG analysis (7 CPTAC datasets)
dataset_ids <- c(
    "brca_cptac_2020", "luad_cptac_2020", "lusc_cptac_2021",
    "ucec_cptac_2020", "coad_cptac_2019", "gbm_cptac_2021",
    "paad_cptac_2021"
)

## TP53 stratification levels
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")

## Output directory prefix
OUTPUT_PREFIX <- "proteomic_DEG"

# TP53 Mutation Classification -------------------------------------------

## Variant classes considered as protein-altering mutations
## (only these will be counted as TP53-mutant)
TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION", "NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
    "IN_FRAME_DEL", "IN_FRAME_INS",
    "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

# Batch Effect Settings --------------------------------------------------

## Batch value cleaning policy
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": how to handle '|' in batch values
BATCH_MIN_PER_LEVEL <- 2 # Minimum samples per batch level

# Covariate Settings -----------------------------------------------------

## Age missing value indicator (experimental feature)
USE_AGE_MISSING_INDICATOR <- FALSE

## Covariate coverage thresholds (minimum fraction of non-NA values)
COVARIATE_COVERAGE_THRESHOLDS <- c(
    purity = 0.60,
    sex = 0.80,
    age = 0.80
)

# Gene Set Selection -----------------------------------------------------

## MSigDB gene set collections to use
GENESET_GROUPS_TO_RUN <- c("H")

# Heatmap Settings (for meta-analysis visualization) --------------------

## Number of top/bottom genes to display in heatmaps
DATASET_HEATMAP_TOP_N <- 25
DATASET_HEATMAP_BOTTOM_N <- 25
PAN_HEATMAP_TOP_N <- 10
PAN_HEATMAP_BOTTOM_N <- 5

# Helper Functions (must be global) -------------------------------------

## Safe option getter
opt <- function(nm, default) {
    if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
}

## Null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

## Safe filesystem name
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}

# Logging Function -------------------------------------------------------

log_msg <- function(text, ..., .envir = parent.frame()) {
    ts <- format(Sys.time(), "%H:%M:%S")
    msg <- tryCatch(
        {
            if (grepl("\\{[^}]+\\}", text)) {
                glue::glue(text, ..., .envir = .envir)
            } else if (grepl("%", text)) {
                do.call(sprintf, c(list(fmt = text), list(...)))
            } else {
                if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
            }
        },
        error = function(e) text
    )
    cat(sprintf("[%s] %s\n", ts, msg))
}

# Create Output Directories ----------------------------------------------

dir.create("run_info", recursive = TRUE, showWarnings = FALSE)

# Analysis Policy Documentation ------------------------------------------

cat(
    paste(
        "Multiple-testing policy:",
        " - Within each gene-set group and per statistic, we control FDR using fgsea padj (per-dataset/stratum).",
        " - Pan-cancer summaries aggregate directions and counts across datasets/strata (descriptive).",
        " - Additionally, we provide a simple meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
        "   followed by Benjamini-Hochberg to report meta-level FDR (see csn_gsea_pan_summary_TP53/meta_fdr).",
        sep = "\n"
    ),
    file = file.path("run_info", "analysis_notes.txt")
)

log_msg("Configuration loaded successfully")
