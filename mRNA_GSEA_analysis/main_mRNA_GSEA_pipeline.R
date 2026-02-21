# =============================================================================
# main_mRNA_GSEA_pipeline.R
# CSN Subunits -> mRNA-based GSEA Analysis Pipeline (CPTAC Datasets)
# =============================================================================
#
# Output is organized under: mRNA_GSEA/
#
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

# --- Working Directory Setup ---
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("mRNA_GSEA_analysis/main_mRNA_GSEA_pipeline.R")

message("====================================================================")
message("CSN CPTAC mRNA GSEA Analysis")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("mRNA_GSEA_analysis/main_mRNA_GSEA_pipeline.R")) {
    stop(
        "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
        "Current directory: ", getwd(), "\n",
        "Please setwd() to the CSN_CPTAC root directory."
    )
}

message("Working directory: ", getwd(), "\n")

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------
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
    library(readr)
})

# -----------------------------------------------------------------------------
# Source Utility Modules (helper functions)
# -----------------------------------------------------------------------------
script_dir <- "mRNA_GSEA_analysis/R"

# Core configuration and utility functions
source(file.path(script_dir, "01_config.R"))
source(file.path(script_dir, "02_utils_logging.R"))
source(file.path(script_dir, "03_utils_io.R"))
source(file.path(script_dir, "04_utils_geneset.R"))
source(file.path(script_dir, "05_utils_csn_score.R"))
source(file.path(script_dir, "06_utils_gsea.R"))
source(file.path(script_dir, "07_utils_limma.R"))
source(file.path(script_dir, "08_utils_tp53.R"))
source(file.path(script_dir, "09_utils_batch.R"))
source(file.path(script_dir, "10_utils_meta_analysis.R"))
source(file.path(script_dir, "11_utils_covariates.R"))
source(file.path(script_dir, "12_utils_correlation.R"))
source(file.path(script_dir, "13_utils_tpm.R"))

# Extended covariate utilities
source(file.path(script_dir, "14_utils_covariates_ext.R"))

# Core analysis functions (run_predictor_analyses, run_one_stratum)
source(file.path(script_dir, "15_core_analysis_functions.R"))

# Meta-FDR summary utilities
source(file.path(script_dir, "16_utils_meta_fdr_summary.R"))

log_msg("All utility modules loaded successfully")

# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

# Dataset root and IDs
datasets_root <- getwd()
dataset_ids <- c(
    "brca_cptac_2020", "coad_cptac_2019", "gbm_cptac_2021",
    "luad_cptac_2020", "lusc_cptac_2021", "paad_cptac_2021", "ucec_cptac_2020"
)

# Build dataset directories
dataset_dirs <- setNames(
    lapply(dataset_ids, function(x) file.path(datasets_root, x)),
    dataset_ids
)

# Strata for TP53 analysis
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")

# Which dataset to start from
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' is not in dataset_ids", start_from))
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

log_msg(
    "Datasets available this run (%d): %s",
    length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", ")
)

# =============================================================================
# GENE SET CONFIGURATION
# =============================================================================

# Prepare gene sets (from 04_utils_geneset.R)
genesets_by_group <- build_genesets_by_group(GENESET_GROUPS_TO_RUN)
if (is.null(genesets_by_group) || !length(genesets_by_group)) {
    stop("No gene sets prepared. Check GENESET_GROUPS_TO_RUN configuration.")
}
log_msg("Prepared %d gene set groups", length(genesets_by_group))

# Write gene set manifest
msigdbr_version <- write_geneset_manifest(genesets_by_group)

# Write run configuration
yaml::write_yaml(list(
    seed = 1234,
    min_frac_complete = min_frac_complete,
    minSize = minSize, maxSize = maxSize, fgsea_eps = fgsea_eps,
    min_per_group = min_per_group,
    datasets_root = datasets_root,
    dataset_ids = dataset_ids,
    strata = strata,
    geneset_groups_selected = GENESET_GROUPS_TO_RUN,
    msigdbr_version = msigdbr_version
), file.path("run_info", "run_manifest.yml"))

# =============================================================================
# MAIN ANALYSIS LOOP
# =============================================================================


options(csn.run_passes = c("CovarAdj"))

log_msg("Starting main analysis loop...")

# Helper function to check which strata to run
.should_run <- function(tag) {
    if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) {
        return(TRUE)
    }
    tag %in% only_strata
}

# Main loop: process each dataset
for (ds in names(dataset_dirs_run)) {
    ds_dir <- dataset_dirs_run[[ds]]
    log_msg("== Starting dataset: %s ==", ds)

    # All datasets use covariant-adjusted analysis
    options(csn.run_passes = c("CovarAdj"))

    # Load RNA matrix & TP53 status
    mat0_full <- load_matrix_from_dataset_dir(ds_dir)
    tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

    # Sample sets for three strata
    samples_ALL <- colnames(mat0_full)
    samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
    samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]

    # Output root directory for each stratum
    base_tp53_root <- file.path(OUTPUT_PREFIX, ds, "csn_gsea_results_TP53")

    # Run three strata in order
    if (.should_run("ALL")) {
        run_one_stratum(
            ds_id = ds, ds_dir = ds_dir,
            mat0_full = mat0_full,
            sample_keep = samples_ALL,
            out_root = file.path(base_tp53_root, "ALL"),
            genesets_by_group = genesets_by_group
        )
    }

    if (.should_run("TP53_mutant")) {
        run_one_stratum(
            ds_id = ds, ds_dir = ds_dir,
            mat0_full = mat0_full,
            sample_keep = samples_MUT,
            out_root = file.path(base_tp53_root, "TP53_mutant"),
            genesets_by_group = genesets_by_group
        )
    }

    if (.should_run("TP53_wild_type")) {
        run_one_stratum(
            ds_id = ds, ds_dir = ds_dir,
            mat0_full = mat0_full,
            sample_keep = samples_WT,
            out_root = file.path(base_tp53_root, "TP53_wild_type"),
            genesets_by_group = genesets_by_group
        )
    }

    log_msg("== Finished dataset: %s (TP53 stratified output -> %s) ==", ds, base_tp53_root)
}

# =============================================================================
# META-ANALYSIS
# =============================================================================

log_msg("Running cross-dataset meta-analysis (Stouffer's z)...")
meta_fdr_stouffer(
    dataset_dirs = dataset_dirs_run,
    strata = strata,
    stat_tags = c(
        "GSEA_limma_t_cont",
        "GSEA_limma_interaction"
    ),
    groups = names(genesets_by_group),
    out_root = file.path(OUTPUT_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")
)

# Summarize meta-FDR results across subunits
log_msg("Summarizing meta-FDR results across subunits...")
posthoc_summary_meta_fdr()
posthoc_summary_meta_fdr_interaction()

# =============================================================================
# CSN SUBUNITS PAIRWISE CORRELATION
# =============================================================================

log_msg("Running CSN subunits pairwise correlation analysis...")
for (ds in names(dataset_dirs_run)) {
    ds_dir <- dataset_dirs_run[[ds]]
    if (!dir.exists(ds_dir)) next
    tryCatch(
        csn_pairwise_correlation_one_ds(
            ds_id = ds,
            ds_dir = ds_dir,
            out_root = "RNA_CSN_subunits_correlation_coefficient",
            subunits = csn_subunits,
            min_pairs = 10L
        ),
        error = function(e) log_msg("[correlation] %s ERROR: %s", ds, conditionMessage(e))
    )
}

# =============================================================================
# TPM ANALYSIS FOR CSN SUBUNITS
# =============================================================================

log_msg("Running CSN subunits TPM analysis...")
tryCatch(
    run_csn_subunits_TPM(
        dataset_dirs_map = dataset_dirs_run,
        out_root = "TPM",
        subunits = CSN_SUB_ORDER,
        datasets_to_plot = dataset_ids,
        agg = "median",
        use_log2 = FALSE
    ),
    error = function(e) log_msg("[TPM] ERROR: %s", conditionMessage(e))
)


# =============================================================================
# COMPLETION
# =============================================================================

log_msg("=== Pipeline execution complete ===")
log_msg("Main GSEA output: %s/", OUTPUT_PREFIX)
