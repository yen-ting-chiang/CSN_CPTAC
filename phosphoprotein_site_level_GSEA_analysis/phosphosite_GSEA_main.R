# =============================================================================
# Phosphosite-Level GSEA Analysis using PTMsigDB
# =============================================================================
#
# Analysis Pipeline:
#   1. Load phosphosite quantification matrices (protein-normalized)
#   2. Stratify samples by TP53 status (ALL, Mutant, Wild-type)
#   3. Perform fgsea enrichment using PTMsigDB
#   4. Aggregate results across datasets using Stouffer's method
#
# Output Structure:
#   phosphoproteomic_site_level_GSEA/
#     └── PTMsigDB/BatchAdj/{dataset}/{stratum}/{subunit}/
#           └── GSEA_limma_t_cont.csv
#     └── meta_fdr/
#           └── {stratum}/BatchAdj/{subunit}/PTMsigDB/
#                 └── GSEA_limma_t_cont_meta_fdr.csv
#
# =============================================================================

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

# --- Working Directory Setup ---
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("phosphoprotein_site_level_GSEA_analysis/phosphosite_GSEA_main.R")

message("====================================================================")
message("CSN CPTAC Phosphosite-Level GSEA Analysis (PTMsigDB)")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("phosphoprotein_site_level_GSEA_analysis/phosphosite_GSEA_main.R")) {
    stop(
        "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
        "Current directory: ", getwd(), "\n",
        "Please setwd() to the CSN_CPTAC root directory."
    )
}

message("Working directory: ", getwd(), "\n")

## Load required packages
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(glue)
    library(limma)
    library(imputeLCMD)
    library(fgsea)
    library(vroom)
})

## Source utility functions from this folder's R subdirectory
.script_dir <- "phosphoprotein_site_level_GSEA_analysis"
source(file.path(.script_dir, "R/utils_helpers.R"))
source(file.path(.script_dir, "R/utils_geneset_loading.R"))
source(file.path(.script_dir, "R/utils_data_loading.R"))
source(file.path(.script_dir, "R/utils_csn_score.R"))
source(file.path(.script_dir, "R/utils_gsea.R"))
source(file.path(.script_dir, "R/utils_limma_analysis.R"))
source(file.path(.script_dir, "R/utils_covariates.R"))
source(file.path(.script_dir, "R/utils_tp53.R"))
source(file.path(.script_dir, "R/utils_meta_analysis.R"))

## Force serial execution for reproducibility
.force_serial_execution()

## Create output directories
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)

## Set random seed
set.seed(1234)

# =============================================================================
# 2. ANALYSIS PARAMETERS
# =============================================================================

## CSN subunits to analyze
csn_subunits <- c(
    "GPS1", "COPS2", "COPS3", "COPS4", "COPS5",
    "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9"
)

## Data filtering parameters
min_frac_complete <- 0.30 # Minimum fraction of non-NA values per row (lower for sparse phosphosite data)
min_per_group <- 8 # Minimum samples per analysis group

## GSEA parameters
minSize <- 5 # Minimum gene set size
maxSize <- 1000 # Maximum gene set size
fgsea_eps <- 1e-10 # fgsea precision parameter

## Analysis configuration
RUN_PASSES <- c("BatchAdj")
PIPELINES_TO_RUN <- c("limma_t")
COMBO_PREFIX <- "phosphoproteomic_site_level_GSEA"

## TP53 stratification
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")

## Global flags
USE_AGE_MISSING_INDICATOR <- FALSE

# =============================================================================
# 3. DATASET CONFIGURATION
# =============================================================================

datasets_root <- getwd()

## CPTAC datasets to analyze
## Note: LUSC (lusc_cptac_2021) is excluded because its phosphoproteome data
##       format is incompatible with the other datasets.
dataset_ids <- c(
    "brca_cptac_2020",
    "luad_cptac_2020",
    "ucec_cptac_2020",
    "coad_cptac_2019",
    "gbm_cptac_2021",
    "paad_cptac_2021"
)

dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)

## Check for available datasets
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
    log_msg("Datasets not found: %s", paste(missing_dirs, collapse = ", "))
}

dataset_dirs_run <- dataset_dirs[
    dir.exists(dataset_dirs) &
        file.exists(file.path(dataset_dirs, "data_phosphoprotein_quantification.txt"))
]

if (!length(dataset_dirs_run)) {
    stop("No datasets available for analysis.")
}

log_msg(
    "Available datasets (%d): %s",
    length(dataset_dirs_run),
    paste(names(dataset_dirs_run), collapse = ", ")
)

# =============================================================================
# 4. LOAD GENE SETS (PTMsigDB)
# =============================================================================

## PTMsigDB gene set database file (relative to project root)
## Download from: https://proteogenomics.shinyapps.io/ptmsigdb/ptmsigdb.Rmd
PTMSIGDB_FILE <- "data_PTMsigDB_all_sites_v2.0.0.xlsx"

if (!file.exists(PTMSIGDB_FILE)) {
    stop(
        "PTMsigDB file not found: ", PTMSIGDB_FILE, "\n",
        "Please download from PTMsigDB repository and place in project root.\n",
        "See README.md for setup instructions."
    )
}

genesets_by_group_ptm <- list()
genesets_by_group_ptm[["PTMsigDB"]] <- read_ptmsigdb_gmt(PTMSIGDB_FILE)
log_msg("Loaded PTMsigDB: %d gene sets", length(genesets_by_group_ptm[["PTMsigDB"]]))

## Write run manifest
dir.create(file.path(.script_dir, "run_info"), recursive = TRUE, showWarnings = FALSE)
yaml::write_yaml(list(
    seed = 1234,
    min_frac_complete = min_frac_complete,
    minSize = minSize,
    maxSize = maxSize,
    fgsea_eps = fgsea_eps,
    min_per_group = min_per_group,
    datasets = names(dataset_dirs_run),
    strata = strata
), file.path(.script_dir, "run_info", "run_manifest.yml"))

# =============================================================================
# 5. COVARIATE EXTRACTION FUNCTIONS
# =============================================================================

#' Extract tumor purity covariate from clinical data
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param sample_ids Sample identifiers to extract
#' @return Named numeric vector of purity values
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    samp <- if (file.exists(samp_fp)) {
        suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#"))
    } else {
        NULL
    }

    if (!is.null(samp)) {
        samp <- as.data.frame(samp)
        names(samp) <- .norm_names(names(samp))
    }
    purity <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

    # Dataset-specific purity extraction
    if (ds_id == "brca_cptac_2020") {
        pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
        if (file.exists(pat_fp) && !is.null(samp)) {
            pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
            names(pat) <- .norm_names(names(pat))
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            pid_in_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
            pid_in_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
            if (!is.na(sid) && !is.na(pid_in_samp) && !is.na(pid_in_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat)) {
                map_pt <- setNames(as.character(samp[[pid_in_samp]]), samp[[sid]])
                pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_in_pat]]))
                purity[] <- unname(pt_pur[map_pt[sample_ids]])
            }
        }
    } else if (ds_id == "luad_cptac_2020") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp)) {
                purity[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "lusc_cptac_2021") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp)) {
                purity[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "paad_cptac_2021") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
                medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
                purity[] <- .to01(medv[match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "ucec_cptac_2020") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid)) {
                pc <- suppressWarnings(as.numeric(samp[["PURITY_CANCER"]]))
                pi <- suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]]))
                ps <- suppressWarnings(as.numeric(samp[["PURITY_STROMA"]]))
                idx <- match(sample_ids, samp[[sid]])
                purity_calc <- ifelse(is.finite(pc[idx]), pc[idx], 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0))
                purity[] <- pmin(pmax(purity_calc, 0), 1)
            }
        }
    }
    purity
}

#' Extract sex and age covariates from clinical data
#' @param ds_dir Dataset directory path
#' @param sample_ids Sample identifiers
#' @return Data frame with sex and age columns
get_sex_age_covariates <- function(ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")

    if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
        log_msg("[covars] Clinical files missing; returning NA for sex/age.")
        out <- cbind(
            sex = rep(NA_real_, length(sample_ids)),
            age = rep(NA_real_, length(sample_ids))
        )
        rownames(out) <- sample_ids
        return(out)
    }

    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .norm_names(names(samp))
    names(pat) <- .norm_names(names(pat))

    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]

    if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) {
        stop("[get_sex_age_covariates] Cannot establish sample-patient mapping.")
    }
    map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

    # SEX: Male=1, Female=0
    sex_col <- intersect(c("SEX", "GENDER"), names(pat))[1]
    sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    if (!is.na(sex_col)) {
        raw <- toupper(as.character(pat[[sex_col]]))
        val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
        names(val) <- as.character(pat[[pid_pat]])
        sex[] <- unname(val[map_pt[sample_ids]])
    }

    # AGE
    age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
    age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    if (!is.na(age_col)) {
        v <- suppressWarnings(as.numeric(pat[[age_col]]))
        names(v) <- as.character(pat[[pid_pat]])
        age_raw <- unname(v[map_pt[sample_ids]])
        age[] <- .z_no_impute(age_raw)
    }

    cov_sex <- mean(is.finite(sex)) * 100
    cov_age <- mean(is.finite(age)) * 100
    log_msg("  Covariate coverage: sex %.1f%%, age %.1f%%", cov_sex, cov_age)

    out <- data.frame(
        sex = as.numeric(sex),
        age = as.numeric(age),
        row.names = sample_ids,
        check.names = FALSE
    )
    out$age_z <- out$age
    out
}

# =============================================================================
# 6. COVARIATE SELECTION FUNCTION
# =============================================================================

#' Safely select and filter covariates for regression
#' @param df Data frame of covariates
#' @param sample_order Sample order to align to
#' @param label Label for logging
#' @param y Predictor vector for correlation filtering
#' @param min_cov_named Named vector of minimum coverage thresholds
#' @param max_abs_cor Maximum absolute correlation with predictor
#' @param min_pairs Minimum pairs for correlation
#' @return Filtered covariate data frame or NULL
select_covars_safely <- function(
  df,
  sample_order,
  label = "covars",
  y = NULL,
  min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
  max_abs_cor = 0.30,
  min_pairs = 20L
) {
    logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)

    if (is.null(df) || !nrow(df)) {
        return(NULL)
    }

    # Align to sample order
    if (is.null(rownames(df))) {
        logf("[covars-%s] df has no rownames, skip", label)
        return(NULL)
    }
    so <- as.character(sample_order)
    df <- df[so, , drop = FALSE]

    # Remove TP53 columns (policy)
    tp_cols_idx <- grep("^TP53($|_|)", colnames(df), ignore.case = FALSE)
    if (length(tp_cols_idx) > 0L) {
        tp_cols <- colnames(df)[tp_cols_idx]
        df <- df[, setdiff(colnames(df), tp_cols), drop = FALSE]
        logf("  [covars-%s] policy: drop TP53 columns: %s", label, paste(tp_cols, collapse = ","))
    }

    df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

    # Coverage filtering
    global_min <- if (length(min_cov_named)) min(min_cov_named, na.rm = TRUE) else 0
    keep_cov <- rep(TRUE, ncol(df))
    names(keep_cov) <- colnames(df)

    for (cn in colnames(df)) {
        v <- df[[cn]]
        cover <- if (is.numeric(v)) mean(is.finite(v)) else mean(!is.na(v))

        thr <- if (!is.null(names(min_cov_named)) && (cn %in% names(min_cov_named))) {
            min_cov_named[[cn]]
        } else {
            global_min
        }

        if (is.na(cover) || cover < thr) {
            logf("  [covars-%s] drop %s (coverage=%.0f%% < %.0f%%)", label, cn, 100 * ifelse(is.na(cover), 0, cover), 100 * thr)
            keep_cov[cn] <- FALSE
        }
    }
    df <- df[, keep_cov, drop = FALSE]
    if (!ncol(df)) {
        return(NULL)
    }

    # Correlation filtering with predictor
    keep_cov <- rep(TRUE, ncol(df))
    names(keep_cov) <- colnames(df)

    if (!is.null(y)) {
        y <- suppressWarnings(as.numeric(y))
        skip_rho_gate <- c("purity", "age", "sex")

        for (cn in colnames(df)) {
            v <- df[[cn]]
            if (is.numeric(v) && !(cn %in% skip_rho_gate)) {
                fin <- is.finite(v) & is.finite(y)
                if (sum(fin) >= min_pairs) {
                    r <- suppressWarnings(stats::cor(v[fin], y[fin], method = "spearman"))
                    if (is.finite(r) && abs(r) >= max_abs_cor) {
                        logf("  [covars-%s] drop %s (|rho|=%.3f >= %.2f vs predictor)", label, cn, abs(r), max_abs_cor)
                        keep_cov[cn] <- FALSE
                    }
                }
            }
        }
        df <- df[, keep_cov, drop = FALSE]
        if (!ncol(df)) {
            return(NULL)
        }
    }

    stopifnot(nrow(df) == length(so), identical(rownames(df), so))
    df
}

# =============================================================================
# 7. PREDICTOR ANALYSIS FUNCTION
# =============================================================================

#' Run limma-based GSEA analysis for a single predictor
#' @param predictor_name Name of the predictor (e.g., "COPS5", "CSN_SCORE")
#' @param predictor_vec Named numeric vector of predictor values
#' @param exclude_genes Genes to exclude from ranking (e.g., the predictor gene itself)
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param mat0 Raw expression matrix
#' @param mat Imputed/filtered expression matrix
#' @param out_root Output directory root
#' @param genesets_by_group List of gene set collections
#' @param batch_all Batch factor for all samples
#' @param purity_all Purity covariate for all samples
#' @param sa_all Sex/age covariates data frame
#' @param tp53_num_all TP53 mutation status (numeric)
#' @param is_ALL Whether this is the ALL stratum
run_predictor_analyses <- function(
  predictor_name,
  predictor_vec,
  exclude_genes = NULL,
  ds_id, ds_dir,
  mat0, mat,
  out_root,
  genesets_by_group,
  batch_all = NULL,
  purity_all = NULL,
  sa_all = NULL,
  tp53_num_all = NULL,
  is_ALL = FALSE,
  additional_covars = NULL
) {
    logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
    sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

    stopifnot(!is.null(genesets_by_group), length(genesets_by_group) > 0)

    # Align predictor to matrix samples
    if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
    keep <- intersect(colnames(mat), names(predictor_vec))
    pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
    names(pred_all) <- keep
    fin <- is.finite(pred_all)

    if (sum(fin) < (2L * min_per_group)) {
        logf("[%s] Insufficient non-NA samples (%d < %d), skip", predictor_name, sum(fin), 2L * min_per_group)
        return(invisible(NULL))
    }

    sample_order <- keep[fin]
    pred <- pred_all[fin]
    stopifnot(identical(names(pred), sample_order))

    # Build covariates data frame
    build_covars_df <- function(so) {
        df <- data.frame(row.names = so, check.names = FALSE)
        if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
        if (!is.null(sa_all)) {
            if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
            if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
            if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
                if ("age_missing" %in% colnames(sa_all)) df$age_missing <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
                if ("age_z_imputed" %in% colnames(sa_all)) df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
            }
        }
        # Handle additional_covars (e.g., CSN_SCORE for RESIDUAL_ analysis)
        if (!is.null(additional_covars)) {
            ac <- as.data.frame(additional_covars, check.names = FALSE)
            common_ac <- intersect(rownames(ac), so)
            if (length(common_ac) > 0) {
                ac <- ac[common_ac, , drop = FALSE]
                for (cn in colnames(ac)) {
                    df[common_ac, cn] <- ac[common_ac, cn]
                }
            }
        }
        df
    }
    df_covars0 <- build_covars_df(sample_order)

    # Pick and filter covariates
    pick_covars_df <- function(label) {
        if (is.null(df_covars0) || !nrow(df_covars0)) {
            return(NULL)
        }
        sel <- tryCatch(
            select_covars_safely(
                df = df_covars0,
                sample_order = sample_order,
                label = label,
                y = pred
            ),
            error = function(e) {
                logf("  [covars-%s] SELECT-FAIL: %s", label, conditionMessage(e))
                NULL
            }
        )
        if (is.null(sel)) {
            return(NULL)
        }
        X <- as.data.frame(sel, stringsAsFactors = FALSE)
        X <- X[sample_order, , drop = FALSE]
        if (!ncol(X)) {
            return(NULL)
        }
        X <- coerce_covariates_safely(X)
        good <- vapply(X, function(v) sum(!is.na(v)) >= 3, logical(1))
        if (any(!good)) X <- X[, good, drop = FALSE]
        if (!ncol(X)) {
            return(NULL)
        }

        # CRITICAL: Force-add additional_covars back even if filtered out by select_covars_safely
        # This ensures CSN_SCORE is included in RESIDUAL_ analysis
        if (!is.null(X) && !is.null(additional_covars)) {
            ac_names <- colnames(as.data.frame(additional_covars))
            for (acn in ac_names) {
                if (acn %in% colnames(df_covars0) && !(acn %in% colnames(X))) {
                    X[[acn]] <- df_covars0[rownames(X), acn]
                }
            }
        }

        logf("  [covars-picked:%s] kept = {%s}", label, paste(colnames(X), collapse = ","))
        X
    }

    X_ba_cov <- pick_covars_df("limma-cont:base")
    if (!is.null(X_ba_cov)) X_ba_cov <- coerce_covariates_safely(X_ba_cov)

    # Build design matrix for BatchAdj analysis
    DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
    if (!is.null(X_ba_cov)) {
        for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
    }
    if (!is.null(batch_all)) {
        DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))
    }

    # use_cols must include additional_covars columns (e.g., CSN_SCORE for RESIDUAL_ analysis)
    use_cols <- unique(c(
        "predictor",
        intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba)),
        grep("^SV\\d+$|^PC\\d+$", colnames(DF_ba), value = TRUE)
    ))
    # Add any additional_covars columns to use_cols
    if (!is.null(additional_covars)) {
        ac_names <- colnames(as.data.frame(additional_covars))
        use_cols <- unique(c(use_cols, intersect(ac_names, colnames(DF_ba))))
    }
    use_cols <- use_cols[use_cols %in% colnames(DF_ba)]

    ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])
    n_ok <- sum(ok_rows)

    if (n_ok < 16L) {
        logf("[%s|BatchAdj] Insufficient samples (%d < 16), skip", predictor_name, n_ok)
        return(invisible(NULL))
    }

    DF_ba <- DF_ba[ok_rows, , drop = FALSE]
    DF_ba$predictor <- DF_ba$predictor - mean(DF_ba$predictor, na.rm = TRUE)

    M_ba <- as.matrix(mat[, rownames(DF_ba), drop = FALSE])

    form_ba <- if (ncol(DF_ba) > 0) {
        as.formula(paste("~ 1 +", paste(colnames(DF_ba), collapse = " + ")))
    } else {
        ~1
    }
    des_ba <- stats::model.matrix(form_ba, data = DF_ba, na.action = stats::na.fail)
    stopifnot(nrow(des_ba) == ncol(M_ba))

    # Audit coverage
    audit_covars_coverage(
        tag = sprintf("%s|BatchAdj", predictor_name),
        ds_id = ds_id,
        stratum = basename(out_root),
        su = predictor_name,
        sample_ids = rownames(DF_ba),
        batch = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_ba)])) else NULL,
        covars = DF_ba[, intersect(c("purity", "sex", "age", "batch"), colnames(DF_ba)), drop = FALSE]
    )

    # Helper functions for GSEA
    coef_name <- function(des) if ("predictor" %in% colnames(des)) "predictor" else colnames(des)[ncol(des)]

    gsea_from <- function(pw, stats) {
        if (is.null(pw) || !length(pw) || is.null(stats) || !length(stats)) {
            return(NULL)
        }

        pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
        if (!length(pw_use)) {
            return(NULL)
        }

        set.seed(1L)
        res <- tryCatch(
            suppressWarnings(fgsea::fgseaMultilevel(
                pathways = pw_use, stats = stats,
                minSize = minSize, maxSize = maxSize,
                eps = fgsea_eps
            )),
            error = function(e) {
                suppressWarnings(fgsea::fgseaSimple(
                    pathways = pw_use, stats = stats,
                    nperm = 10000, minSize = minSize, maxSize = maxSize
                ))
            }
        )
        res
    }

    # ============================================================
    # TP53 Interaction Analysis (only for ALL stratum)
    # Generates GSEA_limma_interaction.csv
    # ============================================================
    if (isTRUE(is_ALL) && !is.null(tp53_num_all) &&
        length(genesets_by_group) &&
        any(vapply(genesets_by_group, function(x) length(x) > 0L, logical(1)))) {
        try(
            {
                tp53_ba_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_ba)]))
                tp53_ba_fac <- factor(ifelse(tp53_ba_num == 1, "MT", "WT"), levels = c("WT", "MT"))
                covars_only_ba <- DF_ba[, setdiff(
                    colnames(DF_ba),
                    c("predictor", "TP53_mutant", "TP53_status", "TP53")
                ), drop = FALSE]

                df_int_ba <- data.frame(
                    pred = DF_ba$predictor, tp53 = tp53_ba_fac,
                    covars_only_ba, row.names = rownames(DF_ba), check.names = FALSE
                )
                des_int_ba <- stats::model.matrix(~ pred * tp53 + ., data = df_int_ba)

                fit_int_ba <- limma::eBayes(limma::lmFit(M_ba, des_int_ba))
                coef_int <- "pred:tp53MT"

                if (coef_int %in% colnames(coef(fit_int_ba))) {
                    tt <- limma::topTable(fit_int_ba, coef = coef_int, number = nrow(M_ba), sort.by = "none")
                    tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
                    names(tvec) <- rownames(tt)

                    ranks <- ._ensure_stats_names(tvec, rownames(M_ba))
                    ranks <- ._finite_rank_stats(ranks, label = paste0("A2-", predictor_name, "-TP53_interaction"))

                    for (grp in names(genesets_by_group)) {
                        pw <- genesets_by_group[[grp]]
                        if (is.null(pw) || !length(pw)) next

                        out_dir_int_ba <- file.path(
                            COMBO_PREFIX %||% "phospho_site_csn_gsea_results_TP53_by_collection",
                            safe_fs_name(grp), "BatchAdj", ds_id, basename(out_root), predictor_name
                        )
                        dir.create(out_dir_int_ba, recursive = TRUE, showWarnings = FALSE)

                        # fgsea for interaction
                        res_fg <- .gsea_from_ranks(
                            pathways = pw,
                            stats = ranks,
                            minSize = minSize,
                            maxSize = maxSize,
                            gsea_eps = fgsea_eps,
                            label = paste0("A2-", predictor_name, "-TP53_interaction")
                        )
                        data.table::fwrite(as.data.frame(res_fg), file.path(out_dir_int_ba, "GSEA_limma_interaction.csv"))
                    }
                }
            },
            silent = TRUE
        )
    }

    # Run limma and GSEA for each gene set group (continuous analysis)
    for (grp_name in names(genesets_by_group)) {
        pw <- genesets_by_group[[grp_name]]

        out_dir_A2 <- file.path(
            COMBO_PREFIX %||% "phospho_site_csn_gsea_results_TP53_by_collection",
            safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root), predictor_name
        )
        dir.create(out_dir_A2, recursive = TRUE, showWarnings = FALSE)

        # Fit limma model
        fit2 <- limma::eBayes(limma::lmFit(M_ba, des_ba))
        t2 <- fit2$t[, coef_name(des_ba)]
        t2 <- ._ensure_stats_names(t2, rownames(M_ba))

        # Exclude predictor gene from ranking
        if (!is.null(exclude_genes) && length(exclude_genes)) {
            t2 <- t2[setdiff(names(t2), exclude_genes)]
        }
        t2 <- ._finite_rank_stats(t2, label = paste0("A2-", predictor_name))

        if (!is.null(t2)) {
            res2 <- gsea_from(pw, t2)
            if (!is.null(res2) && nrow(res2)) {
                data.table::fwrite(res2, file.path(out_dir_A2, "GSEA_limma_t_cont.csv"))
            }
        }
    }

    invisible(NULL)
}

# =============================================================================
# 8. STRATUM ANALYSIS FUNCTION
# =============================================================================

#' Run analysis for a single stratum (ALL, TP53_mutant, or TP53_wild_type)
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param mat0_full Full phosphosite expression matrix
#' @param sample_keep Samples to include in this stratum
#' @param out_root Output directory root for this stratum
#' @param genesets_by_group Gene set collections to analyze
run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
    log_msg("  -- stratum: %s | N(sample_keep)=%d", basename(out_root), length(sample_keep))

    # Subset to stratum samples
    keep <- intersect(colnames(mat0_full), sample_keep)
    if (length(keep) < 4) {
        log_msg("[Skip] Insufficient sample size: %d", length(keep))
        return(invisible(NULL))
    }

    mat0 <- mat0_full[, keep, drop = FALSE]
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)

    # Load protein matrix for CSN score calculation
    prot0_full <- load_matrix_from_dataset_dir(ds_dir)
    prot0 <- prot0_full[, intersect(colnames(prot0_full), keep), drop = FALSE]

    # Get CSN subunits present in protein matrix
    present_sub <- intersect(csn_subunits, rownames(prot0))
    if (!length(present_sub)) {
        log_msg("[Skip] No CSN subunits present in this stratum.")
        return(invisible(NULL))
    }

    # Audit CSN score feasibility
    audit_csn_score_feasibility_safe(
        ds_id = ds_id,
        stratum = basename(out_root),
        mat0 = mat0,
        prot0 = prot0,
        present_sub = present_sub,
        min_members = 5L,
        pca_min_samples = 10L,
        min_per_group = min_per_group,
        out_dir = file.path("run_info", "csn_score_audit")
    )

    # Determine if this is ALL stratum (for TP53 interaction analysis)
    is_ALL <- identical(basename(out_root), "ALL")
    tp53_num_all <- NULL
    if (is_ALL) {
        status_all <- get_tp53_status(ds_dir, colnames(mat0))
        tp53_num_all <- setNames(as.numeric(status_all == "TP53_mutant"), colnames(mat0))
    }

    # Get batch and covariates
    bi_all <- get_batch_factor_phospho(ds_dir, colnames(mat0))
    batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
    purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
    sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0))

    # Pre-compute imputed matrix (same for all predictors within this stratum)
    mat_imputed <- impute_and_filter(mat0, min_frac = min_frac_complete)

    # Run analysis for each CSN subunit
    for (su in present_sub) {
        run_predictor_analyses(
            predictor_name = su,
            predictor_vec = prot0[su, ],
            exclude_genes = su,
            ds_id = ds_id, ds_dir = ds_dir,
            mat0 = mat0, mat = mat_imputed,
            out_root = out_root,
            genesets_by_group = genesets_by_group,
            batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
            tp53_num_all = tp53_num_all, is_ALL = is_ALL
        )
    }

    # Run analysis for CSN composite score
    csn_score <- build_csn_score_safe(
        prot0,
        subunits = present_sub, combine_7AB = TRUE,
        min_members = 5L, pca_min_samples = 10L
    )

    if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
        run_predictor_analyses(
            predictor_name = "CSN_SCORE",
            predictor_vec = csn_score,
            exclude_genes = present_sub,
            ds_id = ds_id, ds_dir = ds_dir,
            mat0 = mat0, mat = mat_imputed,
            out_root = out_root,
            genesets_by_group = genesets_by_group,
            batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
            tp53_num_all = tp53_num_all, is_ALL = is_ALL
        )
    } else {
        log_msg("[CSN_SCORE] Insufficient non-NA samples, skipped.")
    }

    # ============================================================
    # RESIDUAL_ Analysis: Run subunit analysis with CSN_SCORE as covariate
    # This generates output in RESIDUAL_{subunit} directories
    # ============================================================
    min_n_resid <- min_per_group
    csn_nonNA <- sum(is.finite(csn_score))

    if (csn_nonNA < min_n_resid) {
        log_msg("[RESIDUAL] CSN score: Insufficient non-NA samples (%d < %d), skip residual_* in the entire batch.", csn_nonNA, min_n_resid)
    } else {
        # Prepare base covariates
        base_covars_all <- data.frame(
            purity = as.numeric(purity_all[colnames(mat0)]),
            sex = as.numeric(sa_all[colnames(mat0), "sex"]),
            age = as.numeric(sa_all[colnames(mat0), "age"]),
            row.names = colnames(mat0), check.names = FALSE
        )

        if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
            base_covars_all$age_missing <- as.numeric(sa_all[colnames(mat0), "age_missing"])
            base_covars_all$age_z_imputed <- as.numeric(sa_all[colnames(mat0), "age_z_imputed"])
        }
        if (is_ALL && !is.null(tp53_num_all)) {
            base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
        }

        for (su in present_sub) {
            # Check for sufficient samples based on CSN_SCORE availability
            if (sum(is.finite(csn_score)) < min_n_resid) {
                log_msg("[RESIDUAL_%s] Insufficient non-NA samples (CSN_SCORE), skipped.", su)
                next
            }

            # Prepare additional covariates (CSN_SCORE)
            add_cov <- data.frame(CSN_SCORE = csn_score, row.names = names(csn_score), check.names = FALSE)

            # Use the RAW subunit expression as the predictor vector
            # But name the predictor "RESIDUAL_su" to keep output filenames unchanged
            run_predictor_analyses(
                predictor_name = paste0("RESIDUAL_", su),
                predictor_vec = prot0[su, ],
                exclude_genes = su,
                ds_id = ds_id, ds_dir = ds_dir,
                mat0 = mat0, mat = mat_imputed,
                out_root = out_root,
                genesets_by_group = genesets_by_group,
                batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
                tp53_num_all = tp53_num_all, is_ALL = is_ALL,
                additional_covars = add_cov
            )
        }
    }

    # ============================================================
    # Summarize results for all predictors
    # ============================================================
    present_sub <- intersect(csn_subunits, rownames(prot0))
    sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))

    for (grp_name in names(genesets_by_group)) {
        ver_root <- file.path(
            COMBO_PREFIX %||% "phospho_site_csn_gsea_results_TP53_by_collection",
            safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root)
        )
        summarize_all_groups(
            out_root = ver_root,
            csn_subunits = sum_units,
            genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
            stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
        )
    }

    invisible(NULL)
}

# =============================================================================
# 9. MAIN ANALYSIS LOOP
# =============================================================================

log_msg("========== Starting Site-Level PTM-SEA Analysis ==========")

for (ds in names(dataset_dirs_run)) {
    ds_dir <- dataset_dirs_run[[ds]]
    log_msg("== Processing Dataset: %s ==", ds)

    # Load phosphosite matrix (protein-adjusted)
    mat0_full_site <- load_phosphosite_matrix_from_dataset_dir(ds_dir, protein_adjust = TRUE)

    # Get TP53 status for stratification
    tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full_site))
    samples_ALL <- colnames(mat0_full_site)
    samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
    samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]

    base_tp53_root_site <- file.path("phospho_site_combo_2", ds, "phospho_site_csn_gsea_results_TP53")

    # Check PTMsigDB availability
    if (length(genesets_by_group_ptm) == 0) {
        log_msg("[PTM-SEA] PTMsigDB not loaded, skipping dataset %s", ds)
        next
    }

    # Run analysis for each TP53 stratum
    run_one_stratum(
        ds_id = ds, ds_dir = ds_dir,
        mat0_full = mat0_full_site,
        sample_keep = samples_ALL,
        out_root = file.path(base_tp53_root_site, "ALL"),
        genesets_by_group = genesets_by_group_ptm
    )

    run_one_stratum(
        ds_id = ds, ds_dir = ds_dir,
        mat0_full = mat0_full_site,
        sample_keep = samples_MUT,
        out_root = file.path(base_tp53_root_site, "TP53_mutant"),
        genesets_by_group = genesets_by_group_ptm
    )

    run_one_stratum(
        ds_id = ds, ds_dir = ds_dir,
        mat0_full = mat0_full_site,
        sample_keep = samples_WT,
        out_root = file.path(base_tp53_root_site, "TP53_wild_type"),
        genesets_by_group = genesets_by_group_ptm
    )


    log_msg("== [PTM-SEA | site-level] Complete dataset: %s ==", ds)
}

# =============================================================================
# 10. META-ANALYSIS (Stouffer's Method)
# =============================================================================

log_msg("========== Running Meta-Analysis ==========")

meta_fdr_stouffer(
    dataset_dirs = dataset_dirs_run,
    strata = strata,
    stat_tags = c("GSEA_limma_t_cont"),
    groups = names(genesets_by_group_ptm),
    out_root = file.path(COMBO_PREFIX, "meta_fdr")
)

# Also run meta-analysis for interaction results
meta_fdr_stouffer(
    dataset_dirs = dataset_dirs_run,
    strata = strata,
    stat_tags = c("GSEA_limma_interaction"),
    groups = names(genesets_by_group_ptm),
    out_root = file.path(COMBO_PREFIX, "meta_fdr")
)

# =============================================================================
# 11. POST-HOC META-FDR SUMMARY
# =============================================================================

log_msg("========== Running Post-hoc Meta-FDR Summary ==========")

## Summarize meta-FDR for GSEA_limma_t_cont across subunits
posthoc_summary_meta_fdr(
    meta_root = file.path(COMBO_PREFIX, "meta_fdr"),
    genesets_by_group = genesets_by_group_ptm
)

## Summarize meta-FDR for GSEA_limma_interaction across subunits
posthoc_summary_meta_fdr_interaction(
    meta_root = file.path(COMBO_PREFIX, "meta_fdr"),
    genesets_by_group = genesets_by_group_ptm
)

log_msg("========== Analysis Complete ==========")
