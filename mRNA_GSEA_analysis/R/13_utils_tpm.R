# =============================================================================
# 13_utils_tpm.R - TPM Calculation Utilities
# =============================================================================
# This module provides functions for converting FPKM/RPKM to TPM and
# generating TPM summary matrices for CSN subunits.
# =============================================================================

# -----------------------------------------------------------------------------
# Safe Filesystem Name Helper
# -----------------------------------------------------------------------------
.safe_fs_tpm <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# -----------------------------------------------------------------------------
# FPKM/RPKM to TPM Conversion
# -----------------------------------------------------------------------------
#' Convert FPKM/RPKM matrix to TPM
#' @param M Expression matrix (genes x samples)
#' @return TPM matrix where each column sums to 1e6
fpkm_rpkm_to_tpm <- function(M) {
    if (!is.matrix(M)) M <- as.matrix(M)
    M[!is.finite(M)] <- 0
    M[M < 0] <- 0
    denom <- colSums(M, na.rm = TRUE)
    denom[denom <= 0] <- NA_real_
    sweep(M, 2, denom / 1e6, "/")
}

# -----------------------------------------------------------------------------
# Read RNA Matrix (Helper)
# -----------------------------------------------------------------------------
#' Read RNA matrix from dataset directory
#' @param ds_dir Dataset directory path
#' @return Expression matrix
.read_rna_mat <- function(ds_dir) {
    stopifnot(dir.exists(ds_dir))
    if (!exists("load_matrix_from_dataset_dir", mode = "function")) {
        stop("Cannot find RNA file reader load_matrix_from_dataset_dir()")
    }
    load_matrix_from_dataset_dir(ds_dir)
}

# -----------------------------------------------------------------------------
# Build Covariates for TPM
# -----------------------------------------------------------------------------
#' Build covariates data frame for TPM analysis
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory
#' @param sample_ids Sample IDs
#' @return Data frame with purity, sex, age covariates
.build_covars_df_tpm <- function(ds_id, ds_dir, sample_ids) {
    pur <- if (exists("get_purity_covariate", mode = "function")) {
        get_purity_covariate(ds_id, ds_dir, sample_ids)
    } else {
        rep(NA_real_, length(sample_ids))
    }
    sa <- if (exists("get_sex_age_covariates", mode = "function")) {
        get_sex_age_covariates(ds_dir, sample_ids)
    } else {
        data.frame(sex = NA_real_, age = NA_real_)
    }
    df <- data.frame(
        purity = pur,
        sex = sa[, "sex"],
        age = sa[, "age"],
        row.names = sample_ids, check.names = FALSE
    )
    if (exists("coerce_covariates_safely", mode = "function")) {
        df <- coerce_covariates_safely(df)
    }
    df <- df[sample_ids, , drop = FALSE]
    df
}

# -----------------------------------------------------------------------------
# Residualize by Covariates
# -----------------------------------------------------------------------------
#' Residualize Y by covariates (sex/age/purity)
#' @param y Named numeric vector
#' @param covars Covariate data.frame
#' @return Adjusted values (residuals + original mean)
.residualize_y_by_covars <- function(y, covars) {
    if (exists("residualize_to_covars", mode = "function")) {
        r <- residualize_to_covars(y, batch = NULL, covars = covars)
    } else {
        DF <- as.data.frame(covars, check.names = FALSE)
        if (ncol(DF) > 0) {
            all_na <- vapply(DF, function(v) all(is.na(v)), logical(1))
            if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]
            is_const <- vapply(DF, function(v) {
                vv <- v[!is.na(v)]
                if (!length(vv)) TRUE else if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
            }, logical(1))
            if (any(is_const)) DF <- DF[, !is_const, drop = FALSE]
        }
        yv <- suppressWarnings(as.numeric(y))
        ok <- is.finite(yv) & (if (ncol(DF) == 0) TRUE else stats::complete.cases(DF))
        if (sum(ok) < 8L) {
            return(setNames(rep(NA_real_, length(y)), names(y)))
        }
        des <- if (ncol(DF) == 0) {
            matrix(1, nrow = sum(ok), ncol = 1)
        } else {
            stats::model.matrix(~ 1 + ., data = DF[ok, , drop = FALSE])
        }
        fit <- lm.fit(x = des, y = yv[ok])
        r <- setNames(rep(NA_real_, length(y)), names(y))
        r[ok] <- fit$residuals
    }
    mu <- suppressWarnings(mean(as.numeric(y[is.finite(y)]), na.rm = TRUE))
    r + ifelse(is.finite(mu), mu, 0)
}

# -----------------------------------------------------------------------------
# Build CSN Subunits TPM Matrix
# -----------------------------------------------------------------------------
#' Build TPM summary matrix for CSN subunits across datasets
#' @param dataset_dirs_map Named vector/list of dataset directories
#' @param subunits CSN subunit names
#' @param agg Aggregation method ("median" or "mean")
#' @param datasets_to_plot Optional subset of datasets
#' @param use_log2 Whether to log2-transform
#' @param adjust_covars Whether to adjust for covariates
#' @return List with mat, metric, order_ds
build_csn_subunits_tpm_matrix <- function(dataset_dirs_map,
                                          subunits = CSN_SUB_ORDER,
                                          agg = c("median", "mean"),
                                          datasets_to_plot = NULL,
                                          use_log2 = TRUE,
                                          adjust_covars = FALSE) {
    agg <- match.arg(agg)
    if (is.null(dataset_dirs_map)) {
        if (exists("dataset_dirs_run")) {
            dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
        } else if (exists("dataset_dirs")) {
            dataset_dirs_map <- get("dataset_dirs", inherits = TRUE)
        } else {
            stop("Cannot find dataset directory mapping")
        }
    }

    ds_all <- names(dataset_dirs_map)
    if (!is.null(datasets_to_plot)) {
        keep <- datasets_to_plot[datasets_to_plot %in% ds_all]
    } else {
        keep <- ds_all
    }
    if (!length(keep)) stop("No intersection between datasets_to_plot and available datasets")

    out_list <- list()
    for (ds in keep) {
        ds_dir <- dataset_dirs_map[[ds]]
        if (!dir.exists(ds_dir)) next
        mat0 <- .read_rna_mat(ds_dir)

        # Handle MYEOV2 -> COPS9 mapping for specific datasets
        if (ds %in% c("brca_cptac_2020", "gbm_cptac_2021")) {
            if ("MYEOV2" %in% rownames(mat0) && !("COPS9" %in% rownames(mat0))) {
                rownames(mat0)[rownames(mat0) == "MYEOV2"] <- "COPS9"
            }
        }

        present <- intersect(subunits, rownames(mat0))
        if (!length(present)) next

        tpm <- fpkm_rpkm_to_tpm(mat0)
        sub_mat <- tpm[present, , drop = FALSE]

        if (isTRUE(adjust_covars)) {
            sam_all <- colnames(sub_mat)
            cov_df <- .build_covars_df_tpm(ds, ds_dir, sam_all)
            adj_mat <- sub_mat
            for (g in rownames(sub_mat)) {
                y <- as.numeric(sub_mat[g, sam_all])
                names(y) <- sam_all
                adj_mat[g, sam_all] <- .residualize_y_by_covars(y, cov_df)
            }
            vec <- switch(agg,
                median = apply(adj_mat, 1, stats::median, na.rm = TRUE),
                mean   = apply(adj_mat, 1, mean, na.rm = TRUE)
            )
        } else {
            vec <- switch(agg,
                median = apply(sub_mat, 1, stats::median, na.rm = TRUE),
                mean   = apply(sub_mat, 1, mean, na.rm = TRUE)
            )
        }
        out_list[[ds]] <- vec
    }

    if (!length(out_list)) stop("No dataset produced TPM summary values")

    all_rows <- unique(unlist(lapply(out_list, names)))
    all_rows <- intersect(subunits, all_rows)
    MAT <- do.call(cbind, lapply(out_list, function(v) v[all_rows]))
    colnames(MAT) <- names(out_list)
    rownames(MAT) <- all_rows

    metric_tag <- if (isTRUE(use_log2)) "log2TPM" else "TPM"
    if (isTRUE(use_log2)) MAT <- log2(MAT + 1)

    list(mat = MAT, metric = metric_tag, order_ds = names(out_list))
}

# -----------------------------------------------------------------------------
# Run CSN Subunits TPM Analysis
# -----------------------------------------------------------------------------
#' Generate TPM matrix and summary CSV for CSN subunits
#' Outputs two versions: RAW (no covariate adjustment) and AdjCovars (sex/age/purity adjusted)
#' @param dataset_dirs_map Dataset directories mapping
#' @param out_root Output directory root
#' @param subunits CSN subunit names
#' @param datasets_to_plot Optional subset of datasets to process
#' @param agg Aggregation method ("median" or "mean")
#' @param use_log2 Whether to log2-transform (default FALSE to match original)
#' @return Invisible result list with both versions
run_csn_subunits_TPM <- function(dataset_dirs_map = NULL,
                                 out_root = "TPM",
                                 subunits = CSN_SUB_ORDER,
                                 datasets_to_plot = NULL,
                                 agg = "median",
                                 use_log2 = FALSE) {
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

    ## ---- RAW version (no covariate adjustment) ----
    res_raw <- build_csn_subunits_tpm_matrix(
        dataset_dirs_map = dataset_dirs_map,
        subunits = subunits,
        agg = agg,
        datasets_to_plot = datasets_to_plot,
        use_log2 = use_log2,
        adjust_covars = FALSE
    )
    MAT_raw <- res_raw$mat
    metric <- res_raw$metric

    out_csv_raw <- file.path(out_root, sprintf("CSN_subunits_TPM_%s_%s__RAW.csv", agg, metric))
    df_raw <- as.data.frame(MAT_raw)
    df_raw$subunit <- rownames(df_raw)
    df_raw <- df_raw[, c("subunit", setdiff(colnames(df_raw), "subunit"))]
    data.table::fwrite(df_raw, out_csv_raw)
    log_msg("[TPM] RAW matrix CSV: %s", out_csv_raw)

    ## ---- AdjCovars (sex/age/purity adjusted) version ----
    res_adj <- build_csn_subunits_tpm_matrix(
        dataset_dirs_map = dataset_dirs_map,
        subunits = subunits,
        agg = agg,
        datasets_to_plot = datasets_to_plot,
        use_log2 = use_log2,
        adjust_covars = TRUE
    )
    MAT_adj <- res_adj$mat

    out_csv_adj <- file.path(out_root, sprintf("CSN_subunits_TPM_%s_%s__AdjCovars.csv", agg, metric))
    df_adj <- as.data.frame(MAT_adj)
    df_adj$subunit <- rownames(df_adj)
    df_adj <- df_adj[, c("subunit", setdiff(colnames(df_adj), "subunit"))]
    data.table::fwrite(df_adj, out_csv_adj)
    log_msg("[TPM] AdjCovars matrix CSV: %s", out_csv_adj)

    invisible(list(
        csv_raw = out_csv_raw, matrix_raw = MAT_raw,
        csv_adj = out_csv_adj, matrix_adj = MAT_adj,
        metric = metric
    ))
}
