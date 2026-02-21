## ====================================================================
## main_analysis.R
##
## CSN CPTAC Proteomic GSEA Analysis
##
## ====================================================================

## ====================================================================
## PART 1: Environment Setup
## ====================================================================

# --- Working Directory Setup ---
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("protein_level_GSEA_analysis/main_analysis.R")

message("====================================================================")
message("CSN CPTAC Proteomic GSEA Analysis")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("protein_level_GSEA_analysis/main_analysis.R")) {
    stop(
        "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
        "Current directory: ", getwd(), "\n",
        "Please setwd() to the CSN_CPTAC root directory."
    )
}

message("Working directory: ", getwd(), "\n")

# --- Load Required Packages ---
message("[SETUP] Loading required packages...")
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(janitor)
    library(glue)
    library(limma)
    library(imputeLCMD)
    library(msigdbr)
    library(fgsea)
    library(openxlsx)
    library(cowplot)
    library(yaml)
})
message("[SETUP] Packages loaded successfully.\n")

# --- Load Configuration ---
message("[SETUP] Loading configuration...")
config_file <- file.path("protein_level_GSEA_analysis", "config", "analysis_parameters.yaml")
if (file.exists(config_file)) {
    config <- yaml::read_yaml(config_file)
    message("[SETUP] Configuration loaded from: ", config_file)
} else {
    stop("Configuration file not found: ", config_file)
}

# --- Source Modular Functions ---
message("[SETUP] Loading modular function libraries...")


script_dir <- getwd() # Should be CSN_CPTAC root
module_dir <- file.path(script_dir, "protein_level_GSEA_analysis", "R")

# Verify module directory exists
if (!dir.exists(module_dir)) {
    stop(
        "Module directory not found: ", module_dir, "\n",
        "Please ensure you are running from CSN_CPTAC root directory."
    )
}

modules <- c(
    "01_utility_functions.R",
    "02_data_loading.R",
    "03_csn_score.R",
    "04_batch_detection.R",
    "05_covariate_selection.R",
    "06_purity_extraction.R",
    "07_msigdb_utils.R",
    "08_meta_analysis.R",
    "09_meta_analysis_summary.R",
    "10_pairwise_correlation.R"
)

for (mod in modules) {
    mod_path <- file.path(module_dir, mod)

    # Verify module exists
    if (!file.exists(mod_path)) {
        stop("Module not found: ", mod_path)
    }

    source(mod_path, encoding = "UTF-8")
    message("  Loaded: ", mod)
}
message("[SETUP] All modules loaded successfully.\n")


## ====================================================================
## PART 2: Global Parameters from Configuration
## ====================================================================

message("[CONFIG] Setting global parameters...")

# Random seed for reproducibility
set.seed(config$random_seed)
message("  Random seed: ", config$random_seed)

# CSN subunits
csn_subunits <- config$csn_subunits
message("  CSN subunits: ", paste(csn_subunits, collapse = ", "))

# Analysis parameters
min_frac_complete <- config$min_frac_complete
minSize <- config$gsea$minSize
maxSize <- config$gsea$maxSize
fgsea_eps <- config$gsea$fgsea_eps
min_per_group <- config$min_per_group
MAKE_PLOTS <- config$make_plots
.RUN_LIMMA <- config$run_limma

# Batch parameters
BATCH_PIPE_POLICY <- config$batch$pipe_policy
BATCH_MIN_PER_LEVEL <- config$batch$min_per_level

# Heatmap parameters
DATASET_HEATMAP_TOP_N <- config$heatmap$dataset$top_n
DATASET_HEATMAP_BOTTOM_N <- config$heatmap$dataset$bottom_n
PAN_HEATMAP_TOP_N <- config$heatmap$pan$top_n
PAN_HEATMAP_BOTTOM_N <- config$heatmap$pan$bottom_n
PLOT_DATASET_COLLECTIONS <- config$heatmap$plot_dataset_collections
PLOT_PAN_COLLECTIONS <- config$heatmap$plot_pan_collections

# Dataset configuration
dataset_ids <- config$datasets$ids
GSEA_PREFIX <- config$output$gsea_prefix
strata <- config$strata

# Gene set groups
GENESET_GROUPS_TO_RUN <- config$geneset_groups_to_run

# TP53 variant classes
TP53_KEEP_CLASSES <- config$tp53_keep_classes

message("[CONFIG] Parameters set successfully.\n")

# Force serial execution
message("[SETUP] Enforcing serial execution for reproducibility...")
.force_serial_execution()
message("[SETUP] Serial execution enabled.\n")

# Create output directories
message("[SETUP] Creating output directories...")
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("run_info", "csn_score_audit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("run_info", "covars_audit"), recursive = TRUE, showWarnings = FALSE)
message("[SETUP] Output directories created.\n")

# Write analysis notes
cat(
    paste(
        "Multiple-testing policy:",
        " - Within each gene-set group and per statistic, we control FDR using fgsea padj (per-dataset/stratum).",
        " - Additionally, we provide a simple meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
        "   followed by Benjamini-Hochberg to report meta-level FDR.",
        sep = "\n"
    ),
    file = file.path("run_info", "analysis_notes.txt")
)

## ====================================================================
## PART 3: Load Gene Sets
## ====================================================================

message("====================================================================")
message("LOADING GENE SETS")
message("====================================================================\n")

genesets_by_group <- load_msigdb_genesets(GENESET_GROUPS_TO_RUN)

# Write geneset manifest
msigdb_version <- write_geneset_manifest(genesets_by_group)
if (!is.na(msigdb_version)) {
    log_msg("MSigDB version: %s", msigdb_version)
}

message(sprintf(
    "\n[GENESETS] Loaded %d gene set group(s): %s\n",
    length(genesets_by_group),
    paste(names(genesets_by_group), collapse = ", ")
))

## ====================================================================
## PART 4: Dataset Setup
## ====================================================================

message("====================================================================")
message("DATASET SETUP")
message("====================================================================\n")

datasets_root <- getwd()
dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)

message("Datasets root: ", datasets_root)
message("Expected datasets: ", paste(dataset_ids, collapse = ", "), "\n")

# Check for missing datasets
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
    log_msg(
        "Warning: %d dataset(s) not found, will skip: %s",
        length(missing_dirs), paste(missing_dirs, collapse = ", ")
    )
}

# Filter to available datasets with protein quantification data
dataset_dirs_run <- dataset_dirs[
    dir.exists(dataset_dirs) &
        file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
]

if (!length(dataset_dirs_run)) {
    stop("[ERROR] No valid datasets found. Check if data_protein_quantification.txt exists.")
}

log_msg(
    "Datasets available for analysis (%d): %s",
    length(dataset_dirs_run),
    paste(names(dataset_dirs_run), collapse = ", ")
)

message("")

## ====================================================================
## PART 5: Audit Clinical Data (Sex, Age, Batch, Purity)
## ====================================================================

message("====================================================================")
message("AUDITING CLINICAL DATA")
message("====================================================================\n")

audit_results <- audit_all_datasets_sa_batch(
    dataset_dirs_run,
    pipe_policy = BATCH_PIPE_POLICY,
    min_per_level = BATCH_MIN_PER_LEVEL
)

if (!is.null(audit_results)) {
    message("\n[AUDIT] Clinical data audit completed for %d datasets.\n", nrow(audit_results))
}

## ====================================================================
## Core GSEA Analysis Function
## ====================================================================
## This function contains the core Limma + fGSEA analysis logic.
## ====================================================================

message("====================================================================")
message("SETUP COMPLETED - Starting Analysis")
message("====================================================================\n")

## ====================================================================
## PART 6: Additional Helper Functions
## ====================================================================

# Get purity covariate (simplified wrapper for module function)
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
    pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)
    if (!is.null(pur$vec)) {
        return(pur$vec)
    }
    setNames(rep(NA_real_, length(sample_ids)), sample_ids)
}

# Age missing indicator flag
USE_AGE_MISSING_INDICATOR <- FALSE

## ====================================================================
## PART 7: Core GSEA Analysis Function
## ====================================================================


#' Core GSEA Ranking and Execution
#' @keywords internal
._run_fgsea_safe <- function(pw, stats, minSize, maxSize, eps = 1e-10) {
    if (is.null(pw) || !length(pw) || is.null(stats) || !length(stats)) {
        return(NULL)
    }

    pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
    if (!length(pw_use)) {
        return(NULL)
    }

    set.seed(1L)
    tryCatch(
        {
            suppressWarnings(fgsea::fgseaMultilevel(
                pathways = pw_use, stats = stats,
                minSize = minSize, maxSize = maxSize,
                eps = eps
            ))
        },
        error = function(e) {
            suppressWarnings(fgsea::fgseaSimple(
                pathways = pw_use, stats = stats,
                nperm = 10000, minSize = minSize, maxSize = maxSize
            ))
        }
    )
}

#' Run Predictor-Based GSEA Analysis
#'
#' @param predictor_name Name of the predictor
#' @param predictor_vec Named numeric vector of predictor values
#' @param exclude_genes Genes to exclude from ranking
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory
#' @param mat0 Original protein matrix
#' @param mat Imputed and filtered matrix
#' @param out_root Output directory
#' @param genesets_by_group List of gene sets by group
#' @param batch_all Batch factor or NULL
#' @param purity_all Purity vector or NULL
#' @param sa_all Sex/age data frame
#' @param tp53_num_all TP53 numeric status (0/1) or NULL
#' @param is_ALL Whether this is the ALL stratum
#' @param extra_covars_df Additional covariates
#' @return Invisible NULL
run_predictor_analyses <- function(
  predictor_name, predictor_vec, exclude_genes = NULL,
  ds_id, ds_dir, mat0, mat, out_root, genesets_by_group,
  batch_all = NULL, purity_all = NULL, sa_all = NULL,
  tp53_num_all = NULL, is_ALL = FALSE, extra_covars_df = NULL
) {
    if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
    keep <- intersect(colnames(mat), names(predictor_vec))
    pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
    names(pred_all) <- keep
    fin <- is.finite(pred_all)

    if (sum(fin) < (2L * min_per_group)) {
        log_msg("  [%s] Insufficient non-NA samples (%d < %d) → skip", predictor_name, sum(fin), 2L * min_per_group)
        return(invisible(NULL))
    }

    sample_order <- keep[fin]
    pred <- pred_all[fin]

    # 1) Build covariates table
    df_covars0 <- data.frame(row.names = sample_order, check.names = FALSE)
    if (!is.null(purity_all)) df_covars0$purity <- suppressWarnings(as.numeric(purity_all[sample_order]))
    if (!is.null(sa_all)) {
        if ("sex" %in% colnames(sa_all)) df_covars0$sex <- factor(sa_all[sample_order, "sex"])
        if ("age" %in% colnames(sa_all)) df_covars0$age <- suppressWarnings(as.numeric(sa_all[sample_order, "age"]))
    }
    if (!is.null(extra_covars_df) && !is.null(rownames(extra_covars_df))) {
        common_ex <- intersect(sample_order, rownames(extra_covars_df))
        if (length(common_ex) > 0) {
            df_covars0 <- cbind(df_covars0, extra_covars_df[sample_order, , drop = FALSE])
        }
    }

    # 2) Select covariates
    X_ba_cov <- tryCatch(
        select_covars_safely(df = df_covars0, sample_order = sample_order, label = "limma-cont:base", y = pred),
        error = function(e) NULL
    )
    if (!is.null(X_ba_cov)) X_ba_cov <- coerce_covariates_safely(X_ba_cov)

    # 3) Build design matrix
    DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
    if (!is.null(X_ba_cov)) for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
    if (!is.null(batch_all)) DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))

    all_na_col <- vapply(DF_ba, function(v) all(is.na(v)), logical(1))
    if (any(all_na_col)) DF_ba <- DF_ba[, !all_na_col, drop = FALSE]

    use_cols <- unique(c("predictor", intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba))))
    use_cols <- use_cols[use_cols %in% colnames(DF_ba)]
    ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])

    if (sum(ok_rows) < 16L) {
        log_msg("  [%s|BatchAdj] available samples not enough (%d < 16) → skip", predictor_name, sum(ok_rows))
        return(invisible(NULL))
    }

    DF_ba <- DF_ba[ok_rows, , drop = FALSE]
    DF_ba$predictor <- DF_ba$predictor - mean(DF_ba$predictor, na.rm = TRUE)
    M_ba <- as.matrix(mat[, rownames(DF_ba), drop = FALSE])

    des_ba <- stats::model.matrix(as.formula(paste("~ 1 +", paste(colnames(DF_ba), collapse = " + "))), data = DF_ba)

    # Audit alignment
    audit_covars_coverage(tag = sprintf("%s|BatchAdj", predictor_name), ds_id = ds_id, stratum = basename(out_root), su = predictor_name, sample_ids = rownames(DF_ba), batch = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_ba)])) else NULL, covars = DF_ba[, intersect(c("purity", "sex", "age", "batch"), colnames(DF_ba)), drop = FALSE])

    # 4) Regression + GSEA Execution
    if (.RUN_LIMMA) {
        fit2 <- limma::eBayes(limma::lmFit(M_ba, des_ba))
        pname <- if ("predictor" %in% colnames(des_ba)) "predictor" else colnames(des_ba)[ncol(des_ba)]
        t2 <- ._ensure_stats_names(fit2$t[, pname], rownames(M_ba))

        if (length(exclude_genes)) t2 <- t2[setdiff(names(t2), exclude_genes)]
        t2 <- ._finite_rank_stats(t2, label = paste0("A2-", predictor_name))

        if (!is.null(t2)) {
            for (grp_name in names(genesets_by_group)) {
                out_dir_A2 <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root), predictor_name)
                dir.create(out_dir_A2, recursive = TRUE, showWarnings = FALSE)

                res2 <- ._run_fgsea_safe(genesets_by_group[[grp_name]], t2, minSize = minSize, maxSize = maxSize, eps = fgsea_eps)
                if (!is.null(res2) && nrow(res2)) data.table::fwrite(res2, file.path(out_dir_A2, "GSEA_limma_t_cont.csv"))

                ## Interaction Analysis (ALL stratum only)
                if (isTRUE(is_ALL) && !is.null(tp53_num_all)) {
                    tp53_ba_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_ba)]))
                    tp53_ba_fac <- factor(ifelse(tp53_ba_num == 1, "MT", "WT"), levels = c("WT", "MT"))

                    if (nlevels(droplevels(tp53_ba_fac)) >= 2) {
                        df_int <- data.frame(pred = DF_ba$predictor, tp53 = tp53_ba_fac, DF_ba[, setdiff(colnames(DF_ba), c("predictor", "TP53_mutant")), drop = FALSE], check.names = FALSE)
                        des_int <- stats::model.matrix(~ pred * tp53 + ., data = df_int)

                        fit_int <- limma::eBayes(limma::lmFit(M_ba, des_int))
                        if ("pred:tp53MT" %in% colnames(fit_int)) {
                            t_int <- ._ensure_stats_names(fit_int$t[, "pred:tp53MT"], rownames(M_ba))
                            t_int <- ._finite_rank_stats(t_int, label = paste0("A2-", predictor_name, "-TP53_int"))
                            res_int <- ._run_fgsea_safe(genesets_by_group[[grp_name]], t_int, minSize = minSize, maxSize = maxSize, eps = fgsea_eps)
                            if (!is.null(res_int) && nrow(res_int)) data.table::fwrite(res_int, file.path(out_dir_A2, "GSEA_limma_interaction.csv"))
                        }
                    }
                }
            }
        }
    }
    invisible(NULL)
}

## ====================================================================
## PART 8: Per-Stratum Analysis Runner
## ====================================================================

#' Run complete analysis for one stratum
#'
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory
#' @param mat0_full Full protein matrix
#' @param sample_keep Samples to include in this stratum
#' @param out_root Output directory for this stratum
#' @param genesets_by_group Gene sets
#' @return Invisible NULL
run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
    log_msg("  -- stratum: %s | N(sample_keep)=%d", basename(out_root), length(sample_keep))

    ## Subset samples + size check + log scaling
    keep <- intersect(colnames(mat0_full), sample_keep)
    if (length(keep) < 4) {
        log_msg("  [skip] Too few samples: %d", length(keep))
        return(invisible(NULL))
    }
    mat0 <- mat0_full[, keep, drop = FALSE]
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)

    ## Impute & filter + CSN Score audit
    mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
    present_sub <- intersect(csn_subunits, rownames(mat0))
    if (!length(present_sub)) {
        log_msg("  [skip] no CSN subunit in this stratum")
        return(invisible(NULL))
    }
    audit_csn_score_feasibility_safe(
        ds_id = ds_id,
        stratum = basename(out_root),
        mat0 = mat0,
        present_sub = present_sub,
        min_members = 5L,
        pca_min_samples = 10L,
        min_per_group = min_per_group,
        out_dir = file.path("run_info", "csn_score_audit")
    )

    ## TP53 status (ALL stratum only)
    is_ALL <- identical(basename(out_root), "ALL")
    tp53_num_all <- NULL
    if (is_ALL) {
        tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
        tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
        names(tp53_num_all) <- colnames(mat0)
    }

    ## Covariates / Batch
    bi_all <- get_batch_factor(ds_dir, colnames(mat0))
    batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
    purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
    sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0))
    sa_all_limma <- sa_all
    if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
        keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
        sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
    }

    ## Analysis: per subunit, CSN_SCORE, and RESIDUAL_*
    {
        ## 4a) Each subunit (exclude itself)
        for (su in present_sub) {
            run_predictor_analyses(
                predictor_name = su,
                predictor_vec = mat0[su, ],
                exclude_genes = su,
                ds_id = ds_id, ds_dir = ds_dir,
                mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
                out_root = out_root,
                genesets_by_group = genesets_by_group,
                batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
                tp53_num_all = tp53_num_all, is_ALL = is_ALL
            )
        }

        ## 4b) CSN complex score
        csn_score <- build_csn_score_safe(
            mat0,
            subunits = present_sub, combine_7AB = TRUE,
            min_members = 5L, pca_min_samples = 10L
        )

        if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
            run_predictor_analyses(
                predictor_name = "CSN_SCORE",
                predictor_vec = csn_score,
                exclude_genes = present_sub,
                ds_id = ds_id, ds_dir = ds_dir,
                mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
                out_root = out_root,
                genesets_by_group = genesets_by_group,
                batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
                tp53_num_all = tp53_num_all, is_ALL = is_ALL
            )
        } else {
            log_msg("  [CSN_SCORE] non-NA samples not enough, skip")
        }

        ## 4c) RESIDUAL_<SU> (Gene ~ Subunit + CSN_SCORE + Covariates)
        min_n_resid <- min_per_group
        csn_nonNA <- sum(is.finite(csn_score))
        if (csn_nonNA < min_n_resid) {
            log_msg("  [RESIDUAL] CSN score non-NA samples not enough (%d < %d), skip residual_*", csn_nonNA, min_n_resid)
        } else {
            ex_df <- data.frame(CSN_SCORE = csn_score, row.names = names(csn_score))

            for (su in present_sub) {
                if (sum(is.finite(mat0[su, ])) < (2 * min_per_group)) {
                    log_msg("  [RESIDUAL_%s] non-NA samples not enough, skip", su)
                    next
                }


                run_predictor_analyses(
                    predictor_name = paste0("RESIDUAL_", su),
                    predictor_vec = mat0[su, ], # Use RAW subunit expression
                    exclude_genes = su,
                    ds_id = ds_id, ds_dir = ds_dir,
                    mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
                    out_root = out_root,
                    genesets_by_group = genesets_by_group,
                    batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
                    tp53_num_all = tp53_num_all, is_ALL = is_ALL,
                    extra_covars_df = ex_df # Pass CSN_SCORE as covariate
                )
            }
        }
    }

    ## 5) Summary outputs for each version
    present_sub <- intersect(csn_subunits, rownames(mat0))
    sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))

    for (grp_name in names(genesets_by_group)) {
        ver_root <- file.path(GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root))
        summarize_all_groups(
            out_root = ver_root,
            csn_subunits = sum_units,
            genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
            stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
        )
    }

    invisible(NULL)
}

## ====================================================================
## PART 9: Main Execution Loop
## ====================================================================

message("====================================================================")
message("MAIN ANALYSIS EXECUTION")
message("====================================================================\n")

## Run TP53 status auditing for all datasets
log_msg("Auditing TP53 mutation status for all datasets...")
summarize_tp53_counts_all_datasets <- function(dataset_dirs) {
    all_bin <- list()
    all_class <- list()
    k <- 1L
    j <- 1L
    for (ds in names(dataset_dirs)) {
        ds_dir <- dataset_dirs[[ds]]
        if (!dir.exists(ds_dir)) next
        log_msg("[TP53-audit] starting: %s", ds)

        M <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
        if (inherits(M, "try-error")) {
            log_msg("[TP53-audit] %s: cannot read matrix, skip", ds)
            next
        }
        sample_ids <- colnames(M)

        status <- get_tp53_status(ds_dir, sample_ids)
        tb_bin <- table(factor(status, levels = c("TP53_wild_type", "TP53_mutant")))
        bin_row <- data.frame(
            dataset = ds,
            in_matrix_n = length(sample_ids),
            WT_n = as.integer(tb_bin["TP53_wild_type"]),
            MUT_n = as.integer(tb_bin["TP53_mutant"]),
            stringsAsFactors = FALSE
        )
        all_bin[[k]] <- bin_row
        k <- k + 1L
    }

    if (length(all_bin)) {
        bin_df <- dplyr::bind_rows(all_bin) |>
            dplyr::mutate(Any_TP53_mutation_n = in_matrix_n - WT_n)
        data.table::fwrite(bin_df, file.path("run_info", "tp53_status", "tp53_binary_counts_by_dataset.csv"))
        log_msg("[TP53-audit] wrote: tp53_binary_counts_by_dataset.csv")
    }
}

summarize_tp53_counts_all_datasets(dataset_dirs_run)

## Main loop: iterate through datasets and strata
log_msg("\n====================================================================")
log_msg("Starting per-dataset, per-stratum GSEA analysis")
log_msg("====================================================================\n")

for (ds_id in names(dataset_dirs_run)) {
    ds_dir <- dataset_dirs_run[[ds_id]]
    log_msg("\n[DATASET] %s", ds_id)
    log_msg("================================================================\n")

    ## Load full matrix
    mat0_full <- load_matrix_from_dataset_dir(ds_dir)
    mx <- suppressWarnings(max(mat0_full, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) mat0_full <- log2(mat0_full + 1)

    ## Get TP53 status for all samples
    tp53_status_full <- get_tp53_status(ds_dir, colnames(mat0_full))

    ## Run analysis for each stratum
    for (stratum in strata) {
        log_msg(" [STRATUM] %s", stratum)

        if (stratum == "ALL") {
            sample_keep <- colnames(mat0_full)
        } else if (stratum == "TP53_mutant") {
            sample_keep <- names(tp53_status_full)[tp53_status_full == "TP53_mutant"]
        } else if (stratum == "TP53_wild_type") {
            sample_keep <- names(tp53_status_full)[tp53_status_full == "TP53_wild_type"]
        } else {
            log_msg("  Unknown stratum: %s, skip", stratum)
            next
        }

        if (length(sample_keep) < 4) {
            log_msg("  Too few samples (%d), skip", length(sample_keep))
            next
        }

        out_root <- file.path(ds_dir, stratum)

        ## Run analysis for this stratum
        run_one_stratum(
            ds_id = ds_id,
            ds_dir = ds_dir,
            mat0_full = mat0_full,
            sample_keep = sample_keep,
            out_root = out_root,
            genesets_by_group = genesets_by_group
        )
    }
}

## ====================================================================
## PART 10: Meta-Analysis
## ====================================================================

message("\n====================================================================")
message("META-ANALYSIS ACROSS DATASETS")
message("====================================================================\n")

log_msg("Running meta-analysis using Stouffer's Z-score method...")
meta_fdr_stouffer(
    dataset_dirs = dataset_dirs_run,
    strata = strata,
    stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction"),
    groups = names(genesets_by_group),
    out_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")
)


log_msg("Meta-analysis completed!")

## ====================================================================
## PART 11: Meta-Analysis Summary
## ====================================================================

message("\n====================================================================")
message("META-ANALYSIS SUMMARY ACROSS SUBUNITS")
message("====================================================================\n")

log_msg("Summarizing meta-FDR results across subunits...")

# Summary for limma continuous analysis
posthoc_summary_meta_fdr(
    meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata = strata,
    genesets_by_group = genesets_by_group
)

# Summary for limma interaction analysis
posthoc_summary_meta_fdr_interaction(
    meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata = strata,
    genesets_by_group = genesets_by_group
)

log_msg("Meta-analysis summary completed!")

## ====================================================================
## PART 12: Pairwise Correlation Analysis
## ====================================================================

message("\n====================================================================")
message("PAIRWISE CORRELATION ANALYSIS")
message("====================================================================\n")

# Control flag: Set to TRUE to run pairwise correlation analysis
RUN_PAIRWISE_CORRELATION <- get0("RUN_PAIRWISE_CORRELATION", ifnotfound = TRUE)

if (isTRUE(RUN_PAIRWISE_CORRELATION)) {
    log_msg("Running pairwise correlation analysis for CSN subunits...")

    run_csn_subunits_pairwise_correlations(
        dataset_dirs_map = dataset_dirs_run,
        out_root = "CSN_subunits_correlation_coefficient"
    )

    log_msg("Pairwise correlation analysis completed!")
    message("Results saved in: CSN_subunits_correlation_coefficient/")
} else {
    message("Pairwise correlation analysis skipped (RUN_PAIRWISE_CORRELATION = FALSE)")
    message("To enable, set RUN_PAIRWISE_CORRELATION <- TRUE before sourcing this script")
}

## ====================================================================
## ANALYSIS COMPLETED
## ====================================================================

message("\n====================================================================")
message("ANALYSIS COMPLETED SUCCESSFULLY")
message("====================================================================\n")
message("Results are saved in:\n")
message("  1. Per-dataset GSEA: ", GSEA_PREFIX %||% "csn_gsea_results_TP53_by_collection")
message("  2. Meta-analysis: ", if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"))
message("  3. Meta-analysis summary: ", if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr/summary" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr/summary"))
message("  4. Audit files: run_info/")
if (isTRUE(RUN_PAIRWISE_CORRELATION)) {
    message("  5. Pairwise correlations: CSN_subunits_correlation_coefficient/")
}
