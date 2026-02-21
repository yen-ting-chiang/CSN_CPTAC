## =========================================================
## DPS_main.R
## Main Workflow for CSN Subunits Differential Phosphosite (DPS) Analysis
##
## Output strata: c("ALL", "TP53_mutant", "TP53_wild_type")
## =========================================================

# --- Working Directory Setup ---
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("phosphoprotein_DPS_analysis/DPS_main.R")

message("====================================================================")
message("CSN CPTAC Phosphoproteomic DPS Analysis")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("phosphoprotein_DPS_analysis/DPS_main.R")) {
    stop(
        "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
        "Current directory: ", getwd(), "\n",
        "Please setwd() to the CSN_CPTAC root directory."
    )
}

message("Working directory: ", getwd(), "\n")

## =========================================================
## SOURCE UTILITY SCRIPTS
## All utility scripts are sourced from the R/ folder
## =========================================================

## Get the directory where this main script is located
## This ensures sourcing from the correct R/ folder
SCRIPT_DIR <- file.path(getwd(), "phosphoprotein_DPS_analysis")
R_DIR <- file.path(SCRIPT_DIR, "R")

## Verify the R directory exists
if (!dir.exists(R_DIR)) {
    stop(
        "Cannot find R/ directory at: ", R_DIR,
        "\nPlease ensure phosphoprotein_DPS_analysis/R/ folder exists."
    )
}

cat("Sourcing utility scripts from:", R_DIR, "\n")

## Source all utility scripts in order
source(file.path(R_DIR, "00_packages.R"))
source(file.path(R_DIR, "01_utils_general.R"))
source(file.path(R_DIR, "02_data_loading.R"))
source(file.path(R_DIR, "03_batch_processing.R"))
source(file.path(R_DIR, "04_covariates.R"))
source(file.path(R_DIR, "05_tp53_status.R"))
source(file.path(R_DIR, "06_csn_score.R"))
source(file.path(R_DIR, "07_dps_analysis.R"))
source(file.path(R_DIR, "08_meta_analysis.R"))
source(file.path(R_DIR, "09_summary_output.R"))

log_msg("All utility scripts sourced successfully.")

## =========================================================
## CONFIGURATION
## =========================================================

## Force single-thread execution
.force_serial_execution()

## Create run_info directory
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)

## Write analysis description
cat(
    paste(
        "DPS Analysis (Differential Phosphosite):",
        " - Uses limma for differential phosphosite analysis.",
        " - Meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
        "   followed by Benjamini-Hochberg to report meta-level FDR.",
        sep = "\n"
    ),
    file = file.path("run_info", "analysis_notes.txt")
)

## ===== Analysis Parameters =====
set.seed(1234)

## CSN subunit genes
csn_subunits <- c(
    "GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
    "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9"
)

## Data quality thresholds
min_per_group <- 8 # Minimum samples for High/Low each

## Analysis version: BatchAdj (batch-adjusted)
RUN_PASSES <- c("BatchAdj")
options(csn.run_passes = RUN_PASSES)

## DPS-only mode
RUN_DPS_ONLY <- TRUE

## DPS output prefix (hardcoded for consistent output paths)
COMBO_PREFIX_DPS <- "phosphoproteomic_DPS"

## Create TP53 status directory
dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)

## =========================================================
## DATASET CONFIGURATION
## =========================================================

datasets_root <- getwd()

dataset_ids <- c(
    "brca_cptac_2020",
    "coad_cptac_2019",
    "gbm_cptac_2021",
    "luad_cptac_2020",
    "lusc_cptac_2021",
    "paad_cptac_2021",
    "ucec_cptac_2020"
)

dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")

log_msg("datasets_root = %s", datasets_root)

## Check for missing directories
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
    log_msg(
        "Detected %d missing folders, will skip: %s",
        length(missing_dirs), paste(missing_dirs, collapse = ", ")
    )
}

## Build list of available datasets
dataset_dirs_run <- dataset_dirs[
    dir.exists(dataset_dirs) &
        file.exists(file.path(dataset_dirs, "data_phosphoprotein_quantification.txt"))
]

if (!length(dataset_dirs_run)) {
    stop("dataset_dirs_run is empty, please verify folders and data_phosphoprotein_quantification.txt exist")
}

log_msg(
    "Available datasets this run (%d): %s",
    length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", ")
)

## ===== PTMsigDB Configuration =====
## PTMsigDB v2.0.0 (Krug et al., 2019) - required for DPS heatmap annotation
## Download from: https://proteogenomics.shinyapps.io/ptmsigdb/ptmsigdb.Rmd
PTMSIGDB_GMT_FP <- file.path(getwd(), "data_PTMsigDB_all_sites_v2.0.0.xlsx")

genesets_by_group_ptm <- list()
if (nzchar(PTMSIGDB_GMT_FP) && file.exists(PTMSIGDB_GMT_FP)) {
    genesets_by_group_ptm[["PTMsigDB"]] <- read_ptmsigdb_gmt(PTMSIGDB_GMT_FP)
    log_msg("Loaded PTMsigDB site collections: %d sets", length(genesets_by_group_ptm[["PTMsigDB"]]))
} else {
    log_msg("(Note) PTMSIGDB_GMT not set; site-level PTM-SEA will not be available.")
}

## Write run manifest
yaml::write_yaml(list(
    seed = 1234,
    min_per_group = min_per_group,
    datasets_root = datasets_root,
    dataset_ids = dataset_ids,
    strata = strata
), file.path("run_info", "run_manifest.yml"))

## ===== Start from specific dataset (optional) =====
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' not in dataset_ids", start_from))
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

## Filter to only available datasets
dataset_dirs_run <- dataset_dirs_run[
    dir.exists(dataset_dirs_run) &
        file.exists(file.path(dataset_dirs_run, "data_phosphoprotein_quantification.txt"))
]

## =========================================================
## MAIN DPS ANALYSIS
## =========================================================

if (isTRUE(RUN_DPS_ONLY)) {
    log_msg("[DPS] Will process datasets: %s", paste(names(dataset_dirs_run), collapse = ","))

    for (ds in names(dataset_dirs_run)) {
        ds_dir <- dataset_dirs_run[[ds]]
        ds_id <- ds

        ## Read site-level matrix
        M_site_full <- try(load_phosphosite_matrix_from_dataset_dir(ds_dir), silent = TRUE)
        if (inherits(M_site_full, "try-error") || is.null(dim(M_site_full))) {
            log_msg("[DPS] %s unable to read site matrix -> skip", ds)
            next
        }

        ## Read protein matrix (for building predictors)
        prot0_full <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
        if (inherits(prot0_full, "try-error") || is.null(dim(prot0_full))) {
            log_msg("[DPS] %s no protein matrix -> skip", ds)
            next
        }

        ## LUSC special case: skip site-level
        if (identical(ds, "lusc_cptac_2021")) {
            log_msg("[DPS] LUSC special case: skip site-level")
            next
        }

        ## Get TP53 stratification
        all_ids <- colnames(M_site_full)
        tp53_all <- get_tp53_status(ds_dir, all_ids)
        keep_all <- all_ids
        keep_mt <- names(tp53_all)[tp53_all == "TP53_mutant"]
        keep_wt <- names(tp53_all)[tp53_all == "TP53_wild_type"]

        ## Determine output root directory
        base_root <- if (!is.null(COMBO_PREFIX_DPS)) {
            file.path(COMBO_PREFIX_DPS, ds)
        } else {
            file.path(ds_dir, "phospho_DPS_results_TP53")
        }

        ## Execute each stratum
        run_dps_stratum(
            ds_id, ds_dir, M_site_full, keep_all,
            file.path(base_root, "ALL", "DPS"), prot0_full
        )

        if (length(keep_mt) >= (2L * min_per_group)) {
            run_dps_stratum(
                ds_id, ds_dir, M_site_full, keep_mt,
                file.path(base_root, "TP53_mutant", "DPS"), prot0_full
            )
        }

        if (length(keep_wt) >= (2L * min_per_group)) {
            run_dps_stratum(
                ds_id, ds_dir, M_site_full, keep_wt,
                file.path(base_root, "TP53_wild_type", "DPS"), prot0_full
            )
        }

        log_msg("[DPS] Completed dataset: %s -> %s", ds, base_root)
    }
}

## =========================================================
## META-ANALYSIS
## =========================================================

log_msg("Starting meta-analysis across datasets...")

## Build mapping table of DPS result root directories
dataset_dirs_run_dps_out <- setNames(
    if (!is.null(COMBO_PREFIX_DPS)) {
        file.path(COMBO_PREFIX_DPS, names(dataset_dirs_run))
    } else {
        file.path(dataset_dirs_run, "phospho_DPS_results_TP53")
    },
    names(dataset_dirs_run)
)

## Execute Stouffer meta-analysis
dps_meta_stouffer(
    dataset_dirs = dataset_dirs_run_dps_out,
    strata = strata,
    versions = c("BatchAdj"),
    out_root = if (is.null(COMBO_PREFIX_DPS)) {
        "phospho_DPS_meta"
    } else {
        file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
    }
)

## =========================================================
## SUMMARY GENERATION
## =========================================================

log_msg("Generating summary tables...")

## Run posthoc summary
posthoc_summary_meta_fdr(
    meta_root = if (is.null(COMBO_PREFIX_DPS)) {
        "phospho_DPS_meta"
    } else {
        file.path(COMBO_PREFIX_DPS, "phospho_DPS_meta")
    }
)

log_msg("DPS analysis pipeline completed successfully.")
log_msg(
    "Results saved to: %s",
    if (!is.null(COMBO_PREFIX_DPS)) COMBO_PREFIX_DPS else "phospho_DPS_meta"
)
