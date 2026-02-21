# =============================================================================
# phosphoprotein_gene_level_main.R
# CSN Subunits -> Phosphoproteomic Gene-Level GSEA (CPTAC Datasets)
# With TP53 Stratification
# =============================================================================
# Output by strata = c("ALL","TP53_mutant","TP53_wild_type") to
#   <dataset>/csn_gsea_results_TP53/<STRATUM>/...
# Also generate pan-cancer summary and two summary tables to csn_gsea_pan_summary_TP53/
# =============================================================================

# --- Working Directory Setup ---
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("phosphoprotein_gene_level_GSEA_analysis/phosphoprotein_gene_level_main.R")

message("====================================================================")
message("CSN CPTAC Phosphoproteomic Gene-Level GSEA Analysis")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("phosphoprotein_gene_level_GSEA_analysis/phosphoprotein_gene_level_main.R")) {
  stop(
    "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
    "Current directory: ", getwd(), "\n",
    "Please setwd() to the CSN_CPTAC root directory."
  )
}

message("Working directory: ", getwd(), "\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================
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

# =============================================================================
# SOURCE UTILITY MODULES
# =============================================================================
# Source all utilities from phosphoprotein_gene_level_GSEA_analysis/R/
# These must be sourced in order due to dependencies

script_dir <- "phosphoprotein_gene_level_GSEA_analysis/R"

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
source(file.path(script_dir, "12_utils_audit.R"))

# =============================================================================
# INITIALIZATION
# =============================================================================

# Force serial execution for reproducibility
.force_serial_execution()

# Create run_info directory
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)

# Write analysis notes
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

# =============================================================================
# BUILD GENE SETS
# =============================================================================


set.seed(1234)

# Output available MSigDB collections for reference
output_msigdb_collections()

# Build gene sets based on configuration
genesets_by_group <- build_genesets_by_group(GENESET_GROUPS_TO_RUN)

# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

datasets_root <- getwd()
dataset_ids <- c(
  "brca_cptac_2020", "coad_cptac_2019", "gbm_cptac_2021",
  "luad_cptac_2020", "lusc_cptac_2021", "paad_cptac_2021", "ucec_cptac_2020"
)

dataset_dirs <- setNames(
  file.path(datasets_root, dataset_ids),
  dataset_ids
)

# Filter to only existing datasets with phosphoprotein data
dataset_dirs <- dataset_dirs[
  sapply(dataset_dirs, function(d) {
    dir.exists(d) && file.exists(file.path(d, "data_phosphoprotein_quantification.txt"))
  })
]

if (!length(dataset_dirs)) {
  stop("dataset_dirs_run is empty, please verify folder and data_phosphoprotein_quantification.txt exist")
}

log_msg("Available datasets: %s", paste(names(dataset_dirs), collapse = ", "))

# Write geneset manifest
msigdbr_version <- write_geneset_manifest(genesets_by_group)
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

# Decide which dataset to restart from (can be changed)
start_from <- "brca_cptac_2020"
ord <- names(dataset_dirs)
ix <- match(start_from, ord)
if (is.na(ix)) {
  # If start_from not found, start from first dataset
  ix <- 1
  log_msg("start_from '%s' not in dataset_dirs, starting from first dataset", start_from)
}
dataset_dirs_run <- dataset_dirs[ord[ix:length(ord)]]

log_msg("Datasets to run (%d): %s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))


# =============================================================================
# CORE ANALYSIS FUNCTIONS
# =============================================================================

## === robust run_predictor_analyses (center predictors for all limma models) ===
run_predictor_analyses <- function(
  predictor_name,
  predictor_vec, # named numeric by sample (may contain NA)
  exclude_genes = NULL, # genes to exclude before ranking (e.g. self/entire complex)
  ds_id, ds_dir,
  mat0, # RAW (this function handles limma continuous)
  mat, # impute+filter processed (for limma)
  out_root,
  genesets_by_group, # list(group -> pathways)
  batch_all = NULL, # factor by sample or NULL
  purity_all = NULL, # named numeric by sample or NULL
  sa_all = NULL, # data.frame(sex, age[, age_missing, age_z_imputed]); rownames=sample
  tp53_num_all = NULL, # named numeric (0/1) or NULL
  is_ALL = FALSE,
  extra_covars = NULL
) {
  ## -------- Safe parameters --------
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  min_per_group <- opt("min_per_group", 8L)
  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)
  fgsea_eps <- opt("fgsea_eps", 0)
  USE_AGE_MISSING_INDICATOR <- isTRUE(opt("USE_AGE_MISSING_INDICATOR", FALSE))

  logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

  stopifnot(!is.null(genesets_by_group), length(genesets_by_group) > 0)

  ## -------- Align samples (predictor-based) --------
  if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
  keep <- intersect(colnames(mat), names(predictor_vec))
  pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
  names(pred_all) <- keep
  fin <- is.finite(pred_all)
  if (sum(fin) < (2L * min_per_group)) {
    logf("  [%s] predictor non-NA samples insufficient (%d < %d) -> skip", predictor_name, sum(fin), 2L * min_per_group)
    return(invisible(NULL))
  }
  sample_order <- keep[fin]
  pred <- pred_all[fin] # <- original (uncentered) predictor
  stopifnot(identical(names(pred), sample_order))

  ## -------- Build covariate data.frame (no forced NA imputation; only align sample order) --------
  build_covars_df <- function(so) {
    df <- data.frame(row.names = so, check.names = FALSE)
    if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
    if (!is.null(sa_all)) {
      if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
      if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
      if (USE_AGE_MISSING_INDICATOR) {
        if ("age_missing" %in% colnames(sa_all)) df$age_missing <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
        if ("age_z_imputed" %in% colnames(sa_all)) df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
      }
    }
    # [POLICY] TP53 covariate disabled in ALL strata
    df
  }
  df_covars0 <- build_covars_df(sample_order)


  if (!is.null(extra_covars)) {
    # .align_to_samples expects vector or data.frame
    ext_aligned <- .align_to_samples(extra_covars, sample_order, what = "extra_covars")
    if (!is.null(ext_aligned)) {
      ext_df <- as.data.frame(ext_aligned, check.names = FALSE)
      # Merge columns into df_covars0
      for (cn in colnames(ext_df)) {
        # Check for name collision
        if (cn %in% colnames(df_covars0)) {
          logf("  [run_pred] Warning: extra_covars column '%s' overrides existing covariate", cn)
        }
        df_covars0[[cn]] <- ext_df[[cn]]
      }
    }
  }

  coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
      v <- df[[cn]]

      if (is.factor(v) || is.character(v) || is.logical(v)) {
        v <- factor(v) # Keep as factor
        # Only use non-NA samples to check actual level count
        lv <- levels(droplevels(v[!is.na(v)]))
        if (length(lv) <= 1) {
          # Single-level factor -> drop to avoid contrasts error
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        # Numeric stays numeric; don't forcibly convert to factor
        df[[cn]] <- suppressWarnings(as.numeric(v))
      }
    }

    df <- df[, keep, drop = FALSE]
    df
  }


  ## -------- Coverage threshold --------
  base_thr <- c(purity = 0.60, sex = 0.80, age = 0.80)
  present <- colnames(df_covars0)
  extra <- setdiff(present, names(base_thr))
  if (length(extra)) base_thr <- c(base_thr, stats::setNames(rep(min(base_thr), length(extra)), extra))

  ## -------- Select covariates (return numeric data.frame aligned with sample_order) --------
  pick_covars_df <- function(label) {
    cov_src <- df_covars0
    if (is.null(cov_src)) {
      logf("  [covars-%s] NO covariate table in scope ??skip", label)
      return(NULL)
    }
    sel <- tryCatch(
      select_covars_safely(
        df           = cov_src,
        sample_order = sample_order,
        label        = label,
        y            = pred
      ),
      error = function(e) {
        logf("  [covars-%s] SELECT-FAIL: %s", label, conditionMessage(e))
        ## Diagnostic info: also use cov_src to avoid error from non-existent object
        rns <- tryCatch(rownames(cov_src), error = function(e) NULL)
        nrowC <- tryCatch(NROW(cov_src), error = function(e) -1L)
        ncolC <- tryCatch(NCOL(cov_src), error = function(e) -1L)
        hasrn <- is.character(rns) || is.factor(rns)
        logf("  [covars-%s] diag: C_dim = %d x %d; has_rownames=%s", label, nrowC, ncolC, hasrn)
        miss <- tryCatch(setdiff(sample_order, rns), error = function(e) character(0))
        logf(
          "  [covars-%s] samples_not_in_covars = %d / %d (eg: %s)",
          label, length(miss), length(sample_order), paste(head(miss, 5), collapse = ",")
        )
        NULL
      }
    )
    if (is.null(sel)) {
      return(NULL)
    }
    # === Extract drop reasons from select_covars_safely() return (store in variable first) ===
    dl_sel <- tryCatch(attr(sel, "drop_log"), error = function(e) NULL)

    X <- NULL
    if (is.character(sel) && is.null(dim(sel))) {
      nm <- intersect(sel, colnames(cov_src))
      if (length(nm)) X <- cov_src[, nm, drop = FALSE]
    } else {
      X <- as.data.frame(sel, stringsAsFactors = FALSE)
      if (is.null(rownames(X))) {
        if (nrow(X) == length(sample_order)) {
          rownames(X) <- sample_order
        } else if (!is.null(colnames(X)) && ncol(X) == length(sample_order) && all(colnames(X) == sample_order)) {
          X <- as.data.frame(t(as.matrix(X)))
        } else {
          logf("  [covars-%s] Return shape/alignment unclear, discarding", label)
          return(NULL)
        }
      }
      X <- X[sample_order, , drop = FALSE]
    }
    if (is.null(X) || !ncol(X)) {
      return(NULL)
    }
    X <- coerce_covariates_safely(X) # Convert character/logical -> factor; numeric stays numeric
    # Factor-friendly minimum quality threshold: at least 3 non-NA
    good <- vapply(X, function(v) sum(!is.na(v)) >= 3, logical(1))
    if (any(!good)) {
      logf("  [covars-%s] drop low-coverage columns: %s", label, paste(colnames(X)[!good], collapse = ","))
      X <- X[, good, drop = FALSE]
    }
    if (!ncol(X)) {
      return(NULL)
    }
    logf("  [covars-picked:%s] kept = {%s}", label, if (!is.null(X)) paste(colnames(X), collapse = ",") else "NULL")

    if (!is.null(dl_sel) && nrow(dl_sel)) {
      dl_sel$dataset <- ds_id
      dl_sel$stratum <- basename(out_root)
      dl_sel$subunit <- predictor_name
      dl_sel$pass <- label # "limma-cont:RAW" "limma-cont:base"
      dl_sel$samples <- length(sample_order)

      out_dir <- file.path("run_info", "covars_audit")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fp <- file.path(out_dir, "covariate_drop_reasons.csv")

      data.table::fwrite(dl_sel, fp, append = file.exists(fp))
    }

    X
  }
  X_raw_cov <- pick_covars_df("limma-cont:RAW")
  if (!is.null(X_raw_cov)) X_raw_cov <- coerce_covariates_safely(X_raw_cov)
  X_ba_cov <- pick_covars_df("limma-cont:base")
  if (!is.null(X_ba_cov)) X_ba_cov <- coerce_covariates_safely(X_ba_cov)
  ## -------- RAW: build DF -> drop NA -> **center predictor** -> model.matrix --------
  DF_raw <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_raw_cov)) for (nm in colnames(X_raw_cov)) DF_raw[[nm]] <- X_raw_cov[[nm]]
  ok_raw <- stats::complete.cases(DF_raw)
  if (sum(ok_raw) < (2L * min_per_group)) {
    logf("  [%s|RAW] Insufficient usable samples (%d < %d) -> skip", predictor_name, sum(ok_raw), 2L * min_per_group)
    return(invisible(NULL))
  }
  DF_raw <- DF_raw[ok_raw, , drop = FALSE]
  ## -- Centering (only shift mean, no scaling): does not affect coefficients/tests, only intercept interpretation differs
  DF_raw$predictor <- DF_raw$predictor - mean(DF_raw$predictor, na.rm = TRUE)
  M_raw <- as.matrix(mat[, rownames(DF_raw), drop = FALSE])
  des_raw <- stats::model.matrix(~ 1 + ., data = DF_raw, na.action = stats::na.fail)
  stopifnot(nrow(des_raw) == ncol(M_raw))
  # Audit: coverage/covariate summary + alignment
  audit_covars_coverage(
    tag = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    ds_id = ds_id,
    stratum = basename(out_root),
    su = predictor_name,
    sample_ids = rownames(DF_raw),
    batch = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_raw)])) else NULL,
    covars = DF_raw[, setdiff(colnames(DF_raw), "predictor"), drop = FALSE]
  ) # -> will append to run_info/covars_audit/audit_rows.csv

  audit_design_alignment(
    tag = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    samples = colnames(M_raw),
    mod_interest = des_raw,
    mod_nuisance = NULL,
    out_dir = file.path("run_info", "covars_audit")
  ) # -> will write ALIGN_RAW.csv type filename

  ## -------- BatchAdj: same as above; batch enters design table -> **center predictor** --------
  DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_ba_cov)) for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
  if (!is.null(batch_all)) DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))


  all_na_col <- vapply(DF_ba, function(v) all(is.na(v)), logical(1))
  if (any(all_na_col)) {
    logf(
      "[BatchAdj] drop all-NA columns: %s",
      paste(names(all_na_col)[all_na_col], collapse = ",")
    )
    DF_ba <- DF_ba[, !all_na_col, drop = FALSE]
  }


  use_cols <- unique(c(
    "predictor",
    intersect(c("purity", "sex", "age", "batch", "TP53_mutant"), colnames(DF_ba))
  ))
  use_cols <- use_cols[use_cols %in% colnames(DF_ba)] # Guard

  ## Check complete rows using "columns to be used"
  ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])
  n_ok <- sum(ok_rows)

  if (n_ok < 16L) {
    logf("  [%s|BatchAdj] Insufficient usable samples (%d < 16) -> skip", predictor_name, n_ok)
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
  audit_covars_coverage(
    tag        = sprintf("%s|BatchAdj", predictor_name),
    ds_id      = ds_id,
    stratum    = basename(out_root),
    su         = predictor_name,
    sample_ids = rownames(DF_ba),
    batch      = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_ba)])) else NULL,
    covars     = DF_ba[, intersect(c("purity", "sex", "age", "batch"), colnames(DF_ba)), drop = FALSE]
  )


  tag_ba_align <- if (grepl("^RESIDUAL_", predictor_name)) {
    sprintf("%s_%s_%s|BatchAdj", ds_id, basename(out_root), predictor_name)
  } else {
    sprintf("%s|BatchAdj", predictor_name)
  }
  audit_design_alignment(
    tag = tag_ba_align,
    samples = colnames(M_ba),
    mod_interest = des_ba,
    mod_nuisance = NULL,
    out_dir = file.path("run_info", "covars_audit")
  )


  ## -------- (Optional) TP53: interaction and differential correlation (using centered DF_raw/DF_ba$predictor here) --------
  if (isTRUE(is_ALL) && !is.null(tp53_num_all) &&
    length(genesets_by_group) &&
    any(vapply(genesets_by_group, function(x) length(x) > 0L, logical(1)))) {
    ## ---------- RAW branch---------
    if (!("RAW" %in% getOption("csn.run_passes", c("BatchAdj")))) {
      log_msg("[RAW|interaction] skip: pass not selected")
    } else {
      try(
        {
          tp53_raw_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_raw)]))
          tp53_raw_fac <- factor(ifelse(tp53_raw_num == 1, "MT", "WT"), levels = c("WT", "MT"))
          covars_only_raw <- DF_raw[, setdiff(
            colnames(DF_raw),
            c("predictor", "TP53_mutant", "TP53_status", "TP53")
          ), drop = FALSE]
          df_int_raw <- data.frame(
            pred = DF_raw$predictor, tp53 = tp53_raw_fac,
            covars_only_raw, row.names = rownames(DF_raw), check.names = FALSE
          )

          if (nlevels(droplevels(tp53_raw_fac)) < 2) {
            log_msg("[RAW|interaction] skip: TP53 has a single level among usable samples")
          } else {
            des_int_raw <- stats::model.matrix(~ pred * tp53 + ., data = df_int_raw)
            ne <- tryCatch(limma::nonEstimable(des_int_raw), error = function(e) NULL)
            if (!is.null(ne) && any(grepl("^(tp53|pred:tp53)", ne))) {
              log_msg(
                "[RAW|interaction] skip: design has non-estimable TP53 terms (%s)",
                paste(ne[grepl("^(tp53|pred:tp53)", ne)], collapse = ",")
              )
            } else {
              if (.RUN_LIMMA) {
                fit_int_raw <- limma::eBayes(limma::lmFit(M_raw, des_int_raw))
                coef_int <- "pred:tp53MT"
                if (coef_int %in% colnames(coef(fit_int_raw))) {
                  tt <- limma::topTable(fit_int_raw, coef = coef_int, number = nrow(M_raw), sort.by = "none")
                  tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
                  names(tvec) <- rownames(tt)
                  # Consistent with main flow: add names, filter non-finite/duplicates and sort descending
                  ranks <- ._ensure_stats_names(tvec, rownames(M_raw))
                  ranks <- ._finite_rank_stats(ranks, label = paste0("A1-", predictor_name, "-TP53_interaction"))
                  for (grp in names(genesets_by_group)) {
                    pw <- genesets_by_group[[grp]]
                    if (is.null(pw) || !length(pw)) next
                    out_dir_int_raw <- file.path("phosphoproteomic_gene_level_GSEA", safe_fs_name(grp), "RAW", ds_id, basename(out_root), predictor_name)
                    dir.create(out_dir_int_raw, recursive = TRUE, showWarnings = FALSE)

                    # ---- Interaction: fgsea ----
                    res_fg <- .gsea_from_ranks(
                      pathways = pw,
                      stats    = ranks,
                      minSize  = opt("minSize", 15L),
                      maxSize  = opt("maxSize", 500L),
                      gsea_eps = 1e-10,
                      label    = paste0("A1-", predictor_name, "-TP53_interaction")
                    )
                    data.table::fwrite(as.data.frame(res_fg), file.path(out_dir_int_raw, "GSEA_limma_interaction.csv"))
                  }
                }
              }
            }
          }
        },
        silent = TRUE
      )
    }

    ## ---------- BatchAdj branch---------
    if (!("BatchAdj" %in% getOption("csn.run_passes", c("BatchAdj")))) {
      log_msg("[BatchAdj|interaction] skip: pass not selected")
    } else {
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
          if (.RUN_LIMMA) {
            fit_int_ba <- limma::eBayes(limma::lmFit(M_ba, des_int_ba))
            coef_int <- "pred:tp53MT"
            if (coef_int %in% colnames(coef(fit_int_ba))) {
              tt <- limma::topTable(fit_int_ba, coef = coef_int, number = nrow(M_ba), sort.by = "none")
              tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
              names(tvec) <- rownames(tt)
              # Consistent with main flow: add names, filter non-finite/duplicates and sort descending
              ranks <- ._ensure_stats_names(tvec, rownames(M_ba))
              ranks <- ._finite_rank_stats(ranks, label = paste0("A2-", predictor_name, "-TP53_interaction"))
              for (grp in names(genesets_by_group)) {
                pw <- genesets_by_group[[grp]]
                if (is.null(pw) || !length(pw)) next
                out_dir_int_ba <- file.path("phosphoproteomic_gene_level_GSEA", safe_fs_name(grp), "BatchAdj", ds_id, basename(out_root), predictor_name)
                dir.create(out_dir_int_ba, recursive = TRUE, showWarnings = FALSE)

                # ---- Interaction: fgsea ----
                res_fg <- .gsea_from_ranks(
                  pathways = pw,
                  stats    = ranks,
                  minSize  = opt("minSize", 15L),
                  maxSize  = opt("maxSize", 500L),
                  gsea_eps = 1e-10,
                  label    = paste0("A2-", predictor_name, "-TP53_interaction")
                )
                data.table::fwrite(as.data.frame(res_fg), file.path(out_dir_int_ba, "GSEA_limma_interaction.csv"))
              }
            }
          }
        },
        silent = TRUE
      )
    }
  }

  ## -------- Utility functions --------
  coef_name <- function(des) if ("predictor" %in% colnames(des)) "predictor" else colnames(des)[ncol(des)]
  rank_finite <- function(v, exclude = NULL, tag = NULL) {
    keep <- is.finite(v)
    if (any(!keep)) logf("[gsea-%s] drop non-finite: %d", ifelse(is.null(tag), "NA", tag), sum(!keep))
    v <- v[keep]
    if (!is.null(exclude) && length(exclude)) v <- v[setdiff(names(v), exclude)]
    if (!length(v)) {
      return(NULL)
    }
    v
  }
  gsea_from <- function(pw, stats) {
    if (is.null(pw) || !length(pw) || is.null(stats) || !length(stats)) {
      return(NULL)
    }
    # Read global thresholds
    minSize <- opt("minSize", 15L)
    maxSize <- opt("maxSize", 500L)
    # Explicit intersection + size filtering
    orig_n <- length(pw)
    pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf("[gsea] (gsea_from) |H| orig=%d -> use=%d", orig_n, length(pw_use)))
    if (!length(pw_use)) {
      return(NULL)
    }

    set.seed(1L)
    res <- tryCatch(
      {
        suppressWarnings(fgsea::fgseaMultilevel(
          pathways = pw_use, stats = stats,
          minSize = minSize, maxSize = maxSize,
          eps = 1e-10
        ))
      },
      error = function(e) {
        message(sprintf("[gsea] Multilevel failed: %s -> using fgseaSimple instead", conditionMessage(e)))
        suppressWarnings(fgsea::fgseaSimple(
          pathways = pw_use, stats = stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        ))
      }
    )
    res
  }

  ## -------- Per group: A1 (RAW) / A2 (BatchAdj) ----------
  if (.RUN_LIMMA) {
    for (grp_name in names(genesets_by_group)) {
      pw <- genesets_by_group[[grp_name]]

      if ("RAW" %in% getOption("csn.run_passes", c("BatchAdj"))) {
        out_dir_A1 <- file.path("phosphoproteomic_gene_level_GSEA", safe_fs_name(grp_name), "RAW", ds_id, basename(out_root), predictor_name)
        dir.create(out_dir_A1, recursive = TRUE, showWarnings = FALSE)
        fit1 <- limma::eBayes(limma::lmFit(M_raw, des_raw))
        t1 <- fit1$t[, coef_name(des_raw)]
        t1 <- ._ensure_stats_names(t1, rownames(M_raw))
        if (exists("exclude_genes") && length(exclude_genes)) {
          t1 <- t1[setdiff(names(t1), exclude_genes)]
        }
        t1 <- ._finite_rank_stats(t1, label = paste0("A1-", predictor_name))
        if (!is.null(t1)) {
          res1 <- gsea_from(pw, t1)
          if (!is.null(res1) && nrow(res1)) data.table::fwrite(res1, file.path(out_dir_A1, "GSEA_limma_t_cont.csv"))
        }
      }

      ## A2
      if ("BatchAdj" %in% getOption("csn.run_passes", c("BatchAdj"))) {
        out_dir_A2 <- file.path("phosphoproteomic_gene_level_GSEA", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root), predictor_name)
        dir.create(out_dir_A2, recursive = TRUE, showWarnings = FALSE)
        fit2 <- limma::eBayes(limma::lmFit(M_ba, des_ba))
        t2 <- fit2$t[, coef_name(des_ba)]
        t2 <- ._ensure_stats_names(t2, rownames(M_ba))
        if (exists("exclude_genes") && length(exclude_genes)) {
          t2 <- t2[setdiff(names(t2), exclude_genes)]
        }
        t2 <- ._finite_rank_stats(t2, label = paste0("A2-", predictor_name))
        if (!is.null(t2)) {
          res2 <- gsea_from(pw, t2)
          if (!is.null(res2) && nrow(res2)) data.table::fwrite(res2, file.path(out_dir_A2, "GSEA_limma_t_cont.csv"))
        }
      }
    }
  }
}


## ===== CSN subunits coverage & CSN_SCORE (PC1) feasibility audit (robust version) =====
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, prot0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L, # Consistent with threshold in build_csn_score
                                        min_per_group = 8L, # Your limma threshold
                                        out_dir = file.path("run_info", "csn_score_audit")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_samples <- ncol(prot0)
  if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
    log_msg("[CSN-audit] %s | %s: mat0 structure incomplete, skipping", ds_id, stratum)
    return(invisible(NULL))
  }
  present_sub <- intersect(present_sub, rownames(prot0))
  if (!length(present_sub) || n_samples == 0) {
    log_msg("[CSN-audit] %s | %s: No usable CSN subunits or samples = 0, skipping", ds_id, stratum)
    return(invisible(NULL))
  }

  ## 1) Coverage rate for each subunit
  sub_cov <- vapply(present_sub, function(g) mean(is.finite(prot0[g, ])) * 100, numeric(1))
  sub_tbl <- data.frame(
    dataset = ds_id,
    stratum = stratum,
    subunit = present_sub,
    nonNA_pct = round(sub_cov, 1),
    nonNA_n = vapply(present_sub, function(g) sum(is.finite(prot0[g, ])), integer(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  cov_min <- if (length(sub_cov)) round(min(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_med <- if (length(sub_cov)) round(stats::median(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_max <- if (length(sub_cov)) round(max(sub_cov, na.rm = TRUE), 1) else NA_real_

  ## 2) How many subunits are non-NA for each sample; number of samples satisfying >= min_members
  sample_counts <- colSums(is.finite(prot0[present_sub, , drop = FALSE]))
  enough <- sample_counts >= min_members
  n_enough <- sum(enough)

  ## 3) Try to compute CSN_SCORE (PC1) and record feasibility
  csn_score <- build_csn_score(prot0,
    subunits = present_sub,
    combine_7AB = TRUE, min_members = min_members
  )
  csn_nonNA <- sum(is.finite(csn_score))
  csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)

  ## Extra: PC1 explained variance (when feasible)
  pc1_var_pct <- NA_real_
  if (isTRUE(csn_can_pca)) {
    ok_sam <- names(sample_counts)[sample_counts >= min_members]
    get_z <- function(v) {
      v <- as.numeric(v)
      mu <- mean(v[is.finite(v)], na.rm = TRUE)
      sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[!is.finite(v)] <- mu
      (v - mu) / sdv
    }
    Z <- do.call(rbind, lapply(present_sub, function(g) get_z(prot0[g, ok_sam, drop = FALSE])))
    rownames(Z) <- present_sub
    if (all(c("COPS7A", "COPS7B") %in% rownames(Z))) {
      Z7 <- colMeans(Z[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
      Z <- rbind(Z[setdiff(rownames(Z), c("COPS7A", "COPS7B")), , drop = FALSE],
        "COPS7*" = Z7
      )
    }
    pc <- tryCatch(stats::prcomp(t(Z), center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (!is.null(pc)) {
      ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
      pc1_var_pct <- round(ve[1], 1)
    }
  }

  ## 4) Write output files
  # 4a) Summary table
  sum_row <- data.frame(
    dataset = ds_id, stratum = stratum,
    n_samples = n_samples,
    n_subunits_present = length(present_sub),
    subunit_cov_min = cov_min, subunit_cov_median = cov_med, subunit_cov_max = cov_max,
    min_members = min_members,
    samples_with_ge_min_members = n_enough,
    pca_min_samples = pca_min_samples,
    csn_score_nonNA = csn_nonNA,
    csn_pc1_feasible = csn_can_pca,
    csn_pc1_var_pct = pc1_var_pct,
    can_form_high_low_groups = (csn_nonNA >= 2 * min_per_group),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  fp_sum <- file.path(out_dir, "csn_score_feasibility_summary.csv")
  data.table::fwrite(sum_row, fp_sum, append = file.exists(fp_sum))

  # 4b) Subunit coverage (one file per stratum)
  tag <- paste(ds_id, stratum, sep = "_")
  fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
  data.table::fwrite(sub_tbl, fp_sub)

  # 4c) Non-NA subunit count per sample (one file per stratum)
  sample_tbl <- data.frame(
    dataset = ds_id, stratum = stratum,
    sample_id = colnames(prot0),
    nonNA_subunits_n = as.integer(sample_counts),
    ge_min_members = as.logical(enough),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  fp_sam <- file.path(out_dir, sprintf("%s_sample_subunit_counts.csv", tag))
  data.table::fwrite(sample_tbl, fp_sam)

  log_msg(
    "[CSN-audit] %s | %s: subunits=%d; min/median/max coverage=%.1f/%.1f/%.1f%%; eligible_samples=%d; CSN_SCORE nonNA=%d; PC1 feasible=%s (PC1%%=%.1f)",
    ds_id, stratum, length(present_sub), cov_min %||% NaN, cov_med %||% NaN, cov_max %||% NaN,
    n_enough, csn_nonNA, as.character(csn_can_pca), pc1_var_pct %||% NaN
  )

  invisible(list(summary = sum_row, per_subunit = sub_tbl, per_sample = sample_tbl))
}


## ===== Small function to run complete GSEA for "specified sample set" =====
## =========================================================
## Each stratum: produce both RAW and BatchAdj
## =========================================================


# =============================================================================
# run_one_stratum: Per-Stratum Analysis Orchestration
# =============================================================================

run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("  -- stratum: %s | N(sample_keep)=%d", basename(out_root), length(sample_keep))

  ## 1) Subset samples + size check + log scale
  keep <- intersect(colnames(mat0_full), sample_keep)
  if (length(keep) < 4) {
    log_msg("  [Skip] Too few samples: %d", length(keep))
    return(invisible(NULL))
  }
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)

  prot0_full <- load_matrix_from_dataset_dir(ds_dir)
  # Only take samples with intersection and same order as phospho
  sam_keep <- intersect(colnames(mat0), colnames(prot0_full))
  prot0 <- prot0_full[, sam_keep, drop = FALSE]
  mat0 <- mat0[, sam_keep, drop = FALSE]
  # If protein not yet log2, transform consistently with protein pipeline
  mxp <- suppressWarnings(max(prot0, na.rm = TRUE))
  if (is.finite(mxp) && mxp > 100) prot0 <- log2(prot0 + 1)

  ## 2) Imputed+filtered matrix for limma + CSN_SCORE feasibility audit
  mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
  present_sub <- intersect(csn_subunits, rownames(prot0))
  if (!length(present_sub)) {
    log_msg("  [Skip] This stratum has no CSN subunits")
    return(invisible(NULL))
  }
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

  ## === In ALL analysis: prepare TP53 status (0/1; mutant=1) ===
  is_ALL <- identical(basename(out_root), "ALL")
  tp53_num_all <- NULL
  if (is_ALL) {
    tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
    tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
    names(tp53_num_all) <- colnames(mat0)
  }

  ## 3) Covariates / batch (align samples) -- build limma covariates (optionally include missing-indicator)
  bi_all <- get_batch_factor_phospho(ds_dir, colnames(mat0))
  batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
  sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0)) # data.frame(sex, age [, age_missing, age_z_imputed])
  sa_all_limma <- sa_all
  if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
    keep_cols <- intersect(c("sex", "age", "age_missing", "age_z_imputed"), colnames(sa_all))
    sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
  }


  ## 4) Run analysis for each "predictor vector": subunit raw values, CSN_SCORE, RESIDUAL_<SU>
  {
    present_sub <- intersect(csn_subunits, rownames(prot0))
    if (!length(present_sub)) {
      log_msg("  [Skip] This stratum has no CSN subunits")
      return(invisible(NULL))
    }

    ## 4a) Per subunit
    for (su in present_sub) {
      run_predictor_analyses(
        predictor_name = su,
        predictor_vec = prot0[su, ],
        exclude_genes = su,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    }

    ## 4b) CSN complex score (PC1; direction corrected)
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
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    } else {
      log_msg("  [CSN_SCORE] Insufficient non-NA samples, skipping")
    }

    ## 4c) RESIDUAL_<SU> (with guard: if CSN score non-NA samples insufficient -> skip all)
    min_n_resid <- min_per_group
    csn_nonNA <- sum(is.finite(csn_score))
    if (csn_nonNA < min_n_resid) {
      log_msg("  [RESIDUAL] CSN score non-NA samples insufficient (%d < %d), skipping all residual_*", csn_nonNA, min_n_resid)
    } else {
      # Only build base covars when residuals needed (save computation)
      base_covars_all <- data.frame(
        purity = as.numeric(purity_all[colnames(mat0)]),
        sex = as.numeric(sa_all[colnames(mat0), "sex"]),
        age = as.numeric(sa_all[colnames(mat0), "age"]),
        row.names = colnames(mat0), check.names = FALSE
      )
      # (Optional) missing-indicator consistency: include in residualization
      if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing", "age_z_imputed") %in% colnames(sa_all))) {
        base_covars_all$age_missing <- as.numeric(sa_all[colnames(mat0), "age_missing"])
        base_covars_all$age_z_imputed <- as.numeric(sa_all[colnames(mat0), "age_z_imputed"])
      }
      if (is_ALL && !is.null(tp53_num_all)) {
        base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
      }

      for (su in present_sub) {
        # Use original protein vector
        orig_vec <- prot0[su, ]

        # Check sample sufficiency

        common_samples <- intersect(names(orig_vec), names(csn_score))
        valid_n <- sum(is.finite(orig_vec[common_samples]) & is.finite(csn_score[common_samples]))

        if (valid_n < (2 * min_per_group)) {
          log_msg("  [RESIDUAL_%s] Insufficient non-NA samples (One-Step check), skipping", su)
          next
        }

        # Run analysis with extra covariate

        run_predictor_analyses(
          predictor_name = paste0("RESIDUAL_", su),
          predictor_vec = orig_vec, # <--- Original vector
          exclude_genes = su, # Exclude self from GSEA ranking
          ds_id = ds_id, ds_dir = ds_dir,
          mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
          out_root = out_root,
          genesets_by_group = genesets_by_group,
          batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
          tp53_num_all = tp53_num_all, is_ALL = is_ALL,
          extra_covars = data.frame(CSN_SCORE = csn_score) # <--- Pass CSN_SCORE here
        )
      }
    }
  }

  ## Summarize each version separately
  present_sub <- intersect(csn_subunits, rownames(prot0))
  sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))

  # RAW: limma_t_cont + interaction (only summarize existing csv, no re-run)
  for (grp_name in names(genesets_by_group)) {
    ver_root <- file.path("phosphoproteomic_gene_level_GSEA", safe_fs_name(grp_name), "RAW", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
    )
  }


  # BatchAdj: limma_t_cont + interaction (only summarize existing csv, no re-run)
  for (grp_name in names(genesets_by_group)) {
    ver_root <- file.path("phosphoproteomic_gene_level_GSEA", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
    )
  }

  invisible(NULL)
}


## =========================================================
## Process datasets to run in order (complete main loop + define stratified sample sets)
## =========================================================
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  log_msg("== Starting dataset: %s ==", ds)

  ## Protein matrix & TP53 status
  mat0_full <- load_phospho_matrix_from_dataset_dir(ds_dir, protein_adjust = TRUE)
  if (exists("log_msg")) {
    log_msg(
      "[phospho] protein_adjusted = %s",
      if (isTRUE(attr(mat0_full, "protein_adjusted"))) "TRUE" else "FALSE"
    )
  }
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

  ## Sample sets for three strata
  samples_ALL <- colnames(mat0_full)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT <- names(tp53_status)[tp53_status == "TP53_wild_type"]

  ## Output root directory for each stratum
  base_tp53_root <- file.path("phosphoproteomic_gene_level_GSEA", ds, "phospho_csn_gsea_results_TP53")


  ## Run three strata in order
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

  log_msg("== Completed dataset: %s (TP53 stratified output -> %s) ==", ds, base_tp53_root)
}


## Generate cross-dataset meta-FDR (Stouffer -> BH)
## =========================================================
meta_fdr_stouffer(
  dataset_dirs = dataset_dirs_run,
  strata = strata,
  stat_tags = c(
    "GSEA_limma_t_cont",
    "GSEA_limma_interaction"
  ),
  groups = names(genesets_by_group),
  out_root = "phosphoproteomic_gene_level_GSEA/phospho_csn_gsea_pan_summary_TP53/meta_fdr"
)

## =========================================================
## Summarize meta-FDR across subunits (posthoc summary)
## =========================================================
posthoc_summary_meta_fdr()
posthoc_summary_meta_fdr_interaction()

