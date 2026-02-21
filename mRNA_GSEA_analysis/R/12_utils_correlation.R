# =============================================================================
# 12_utils_correlation.R - CSN Subunit Correlation Utilities
# =============================================================================
# This module provides functions for calculating pairwise correlations
# between CSN subunits, including covariate adjustment.
# =============================================================================

# -----------------------------------------------------------------------------
# Safe Filesystem Name
# -----------------------------------------------------------------------------
.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# -----------------------------------------------------------------------------
# Residualize to Covariates
# -----------------------------------------------------------------------------
#' Residualize expression to covariates
#' @param y Named numeric vector
#' @param batch Optional batch factor
#' @param covars Optional covariate data.frame
#' @param min_n Minimum samples required
#' @return Residual vector
residualize_to_covars <- function(y, batch = NULL, covars = NULL, min_n = 8L) {
    stopifnot(!is.null(names(y)))
    sam <- names(y)
    DF <- data.frame(row.names = sam, check.names = FALSE)
    if (!is.null(batch)) DF[["batch"]] <- droplevels(batch[sam])
    if (!is.null(covars)) {
        C <- as.data.frame(covars, check.names = FALSE)
        rn <- rownames(C)
        if (is.null(rn)) stop("[residualize_to_covars] covars must have rownames=samples")
        C <- C[sam, , drop = FALSE]
        if (exists("coerce_covariates_safely", mode = "function")) C <- coerce_covariates_safely(C)
        for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
    }
    yv <- suppressWarnings(as.numeric(y[sam]))

    # Remove all-NA/constant/single-level columns
    if (ncol(DF) > 0) {
        all_na <- vapply(DF, function(v) all(is.na(v)), logical(1))
        if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]
        is_const <- vapply(DF, function(v) {
            vv <- v[!is.na(v)]
            if (!length(vv)) TRUE else if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
        }, logical(1))
        if (any(is_const)) DF <- DF[, !is_const, drop = FALSE]
    }

    ok <- is.finite(yv) & (if (ncol(DF) == 0) TRUE else stats::complete.cases(DF))
    if (sum(ok) < min_n) {
        return(setNames(rep(NA_real_, length(y)), names(y)))
    }
    des <- if (ncol(DF) == 0) model.matrix(~1) else model.matrix(~ 1 + ., data = DF[ok, , drop = FALSE])
    fit <- lm.fit(x = des, y = yv[ok])
    out <- setNames(rep(NA_real_, length(y)), names(y))
    out[ok] <- fit$residuals
    out
}

# -----------------------------------------------------------------------------
# Pairwise Correlation for One Dataset
# -----------------------------------------------------------------------------
#' Calculate pairwise correlations between CSN subunits for one dataset
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory path
#' @param out_root Output root directory
#' @param subunits CSN subunit names
#' @param min_pairs Minimum sample pairs required
#' @param dpi Output DPI for plots
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param point_size Point size for scatter plots
#' @return Invisible results data.frame
csn_pairwise_correlation_one_ds <- function(
  ds_id, ds_dir,
  out_root = "RNA_CSN_subunits_correlation_coefficient",
  subunits = csn_subunits,
  min_pairs = 10L,
  dpi = 600,
  width = 4.5, height = 4.5,
  point_size = 1.6
) {
    # Prepare output directory
    dir.create(out_root, showWarnings = FALSE)
    ds_out <- file.path(out_root, ds_id)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

    # Read matrix
    mat0 <- load_matrix_from_dataset_dir(ds_dir)
    present <- intersect(subunits, rownames(mat0))
    if (length(present) < 2) {
        log_msg("[pairwise] %s: available CSN subunits < 2, skipping", ds_id)
        return(invisible(NULL))
    }
    sam_all <- colnames(mat0)

    # Build covariates
    build_covars_df_local <- function(ds_id, ds_dir, sample_ids) {
        pur <- if (exists("get_purity_covariate", mode = "function")) {
            get_purity_covariate(ds_id, ds_dir, sample_ids)
        } else {
            rep(NA_real_, length(sample_ids))
        }
        sa <- if (exists("get_sex_age_covariates", mode = "function")) {
            get_sex_age_covariates(ds_dir, sample_ids)
        } else {
            data.frame(sex = NA_real_, age = NA_real_, row.names = sample_ids)
        }
        df <- data.frame(
            purity = as.numeric(pur),
            sex = as.numeric(sa[, "sex"]),
            age = as.numeric(sa[, "age"]),
            row.names = sample_ids,
            check.names = FALSE
        )
        if (exists("coerce_covariates_safely", mode = "function")) {
            df <- coerce_covariates_safely(df)
        }
        df
    }
    cov0 <- build_covars_df_local(ds_id, ds_dir, sam_all)

    # Batch
    bi <- get_batch_factor(ds_dir, sam_all)
    batch <- NULL
    if (!is.null(bi) && !is.null(bi$fac)) {
        b <- bi$fac
        if (is.null(names(b))) {
            if (length(b) == length(sam_all)) {
                b <- stats::setNames(b, sam_all)
            } else {
                b <- NULL
            }
        }
        if (!is.null(b)) batch <- droplevels(b[sam_all])
    }

    # Pairwise combinations
    pairs <- utils::combn(present, 2, simplify = FALSE)
    all_rows <- list()


    # Main loop
    for (p in pairs) {
        gx <- p[1]
        gy <- p[2]
        x <- as.numeric(mat0[gx, sam_all])
        names(x) <- sam_all
        y <- as.numeric(mat0[gy, sam_all])
        names(y) <- sam_all
        ok <- stats::complete.cases(x, y)
        if (sum(ok) < min_pairs) next

        # Version 1: NoCovariate
        rr <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "pearson"))
        rs <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        fit <- stats::lm(y ~ x, data = data.frame(x = x[ok], y = y[ok]))
        all_rows[[length(all_rows) + 1L]] <- data.frame(
            dataset = ds_id, version = "NoCovariate",
            gene_x = gx, gene_y = gy, n = sum(ok),
            pearson_r = as.numeric(rr$estimate), pearson_p = rr$p.value,
            spearman_r = as.numeric(rs$estimate), spearman_p = rs$p.value,
            R2 = summary(fit)$r.squared,
            slope = unname(stats::coef(fit)[["x"]]),
            intercept = unname(stats::coef(fit)[["(Intercept)"]]),
            stringsAsFactors = FALSE, check.names = FALSE
        )


        # Version 2: RAW_covars
        cov_raw <- cov0
        xr <- residualize_to_covars(x, batch = NULL, covars = cov_raw)
        yr <- residualize_to_covars(y, batch = NULL, covars = cov_raw)
        ok2 <- stats::complete.cases(xr, yr)
        if (sum(ok2) >= min_pairs) {
            rr2 <- suppressWarnings(stats::cor.test(xr[ok2], yr[ok2], method = "pearson"))
            rs2 <- suppressWarnings(stats::cor.test(xr[ok2], yr[ok2], method = "spearman", exact = FALSE))
            fit2 <- stats::lm(yr ~ xr, data = data.frame(xr = xr[ok2], yr = yr[ok2]))
            all_rows[[length(all_rows) + 1L]] <- data.frame(
                dataset = ds_id, version = "RAW_covars",
                gene_x = gx, gene_y = gy, n = sum(ok2),
                pearson_r = as.numeric(rr2$estimate), pearson_p = rr2$p.value,
                spearman_r = as.numeric(rs2$estimate), spearman_p = rs2$p.value,
                R2 = summary(fit2)$r.squared,
                slope = unname(stats::coef(fit2)[["xr"]]),
                intercept = unname(stats::coef(fit2)[["(Intercept)"]]),
                stringsAsFactors = FALSE, check.names = FALSE
            )
        }
    }

    # Aggregate and output
    if (!length(all_rows)) {
        log_msg("[pairwise] %s: no available pairs", ds_id)
        return(invisible(NULL))
    }
    RES <- do.call(rbind, all_rows)

    # BH adjustment
    RES$pearson_padj <- ave(RES$pearson_p, interaction(RES$dataset, RES$version, drop = TRUE),
        FUN = function(p) stats::p.adjust(p, method = "BH")
    )
    RES$spearman_padj <- ave(RES$spearman_p,
        interaction(RES$dataset, RES$version, drop = TRUE),
        FUN = function(p) stats::p.adjust(p, method = "BH")
    )

    # Write outputs
    out_csv <- file.path(ds_out, "pairwise_correlations_all_versions.csv")
    data.table::fwrite(RES, out_csv)

    out_xlsx <- file.path(ds_out, "pairwise_correlations_all_versions.xlsx")
    wb <- openxlsx::createWorkbook()
    for (ver in unique(RES$version)) {
        openxlsx::addWorksheet(wb, ver)
        openxlsx::writeData(wb, ver, RES[RES$version == ver, ], withFilter = TRUE)
    }
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    log_msg("[pairwise] %s: completed. CSV: %s; XLSX: %s", ds_id, basename(out_csv), basename(out_xlsx))
    invisible(RES)
}

# -----------------------------------------------------------------------------
# Run Pairwise Correlations for All Datasets
# -----------------------------------------------------------------------------
#' Run CSN subunit pairwise correlations for all datasets
#' @param dataset_dirs_map Named list/vector of dataset directories
#' @param out_root Output root directory
#' @return Invisible TRUE
run_csn_subunits_pairwise_correlations <- function(
  dataset_dirs_map = NULL,
  out_root = "RNA_CSN_subunits_correlation_coefficient"
) {
    if (is.null(dataset_dirs_map)) {
        if (exists("dataset_dirs_run")) {
            dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
        } else if (exists("dataset_dirs")) {
            dd <- get("dataset_dirs", inherits = TRUE)
            dataset_dirs_map <- dd[dir.exists(dd) & file.exists(file.path(dd, "data_protein_quantification.txt"))]
        } else {
            stop("Cannot find dataset_dirs_run or dataset_dirs")
        }
    }
    dir.create(out_root, showWarnings = FALSE)
    for (ds in names(dataset_dirs_map)) {
        try(csn_pairwise_correlation_one_ds(ds_id = ds, ds_dir = dataset_dirs_map[[ds]], out_root = out_root), silent = FALSE)
    }
    invisible(TRUE)
}
