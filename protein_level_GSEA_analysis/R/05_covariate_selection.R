## ====================================================================
## 05_covariate_selection.R
##
## Purpose: Covariate selection and filtering functions
## Contains: Coverage-based selection, correlation filtering, TP53 policy
## ====================================================================

#' Select covariates safely with multiple filtering criteria
#'
#' @param df Data frame with candidate covariates (rownames = sample IDs)
#' @param sample_order Character vector of sample IDs (for alignment)
#' @param label Label for logging
#' @param y Optional numeric predictor for correlation filtering
#' @param min_cov_named Named vector of minimum coverage thresholds
#' @param max_abs_cor Maximum absolute correlation with predictor
#' @param min_pairs Minimum pairs for correlation calculation
#' @return Filtered data frame or NULL
#' @export
select_covars_safely <- function(
  df,
  sample_order,
  label = "covars",
  y = NULL,
  min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
  max_abs_cor = 0.30,
  min_pairs = 20L
) {
    ## Drop Cause Collector
    drop_log <- list()
    .append_drop <- function(step, col, reason, metric = NA_real_, threshold = NA_real_, extra = NA_character_) {
        drop_log[[length(drop_log) + 1L]] <<- data.frame(
            step = step,
            covariate = col,
            reason = reason,
            metric = metric,
            threshold = threshold,
            extra = extra,
            stringsAsFactors = FALSE
        )
    }

    if (is.null(df) || !nrow(df)) {
        return(NULL)
    }

    # 1) Align sample order
    if (is.null(rownames(df))) {
        log_msg("  [covars-%s] df no rownames → skip", label)
        return(NULL)
    }
    so <- as.character(sample_order)
    log_msg("  [covars-%s] before align: C_dim=%d x %d; has_rownames=%s", label, NROW(df), NCOL(df), !is.null(rownames(df)))
    df <- df[so, , drop = FALSE]

    ## [POLICY] Never include TP53 as a covariate (ALL / MT / WT strata)
    tp_cols_idx <- grep("^TP53($|_|)", colnames(df), ignore.case = FALSE)
    if (length(tp_cols_idx) > 0L) {
        tp_cols <- colnames(df)[tp_cols_idx]
        for (cc in tp_cols) {
            .append_drop(
                step = "policy", col = cc, reason = "exclude_TP53_as_covariate",
                metric = NA_real_, threshold = NA_real_, extra = "global policy"
            )
        }
        df <- df[, setdiff(colnames(df), tp_cols), drop = FALSE]
        log_msg(sprintf("  [covars-%s] policy: drop TP53 columns from covariates → %s", label, paste(tp_cols, collapse = ",")))
    }

    log_msg("  [covars-%s] after  align: C_dim=%d x %d", label, NROW(df), NCOL(df))

    # 2) Convert to data.frame, retain factors
    df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

    # 3) Coverage filtering
    `%||%` <- function(a, b) if (!is.null(a)) a else b
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
            log_msg(
                "  [covars-%s] drop %s (coverage=%s < %s)", label, cn,
                sprintf("%.0f%%", 100 * ifelse(is.na(cover), 0, cover)),
                sprintf("%.0f%%", 100 * thr)
            )
            keep_cov[cn] <- FALSE

            .append_drop("coverage", cn, "low_coverage",
                metric = ifelse(is.na(cover), 0, cover), threshold = thr
            )
        }
    }
    df <- df[, keep_cov, drop = FALSE]
    if (!ncol(df)) {
        return(NULL)
    }

    ## Reinitialize keep_cov
    keep_cov <- rep(TRUE, ncol(df))
    names(keep_cov) <- colnames(df)

    # 4) Correlation filtering (if predictor provided)
    if (!is.null(y)) {
        y <- suppressWarnings(as.numeric(y))
        ## Do not perform rho gate on these covariates
        skip_rho_gate <- c("purity", "age", "sex", "CSN_SCORE")

        for (cn in colnames(df)) {
            v <- df[[cn]]

            if (is.numeric(v) && !(cn %in% skip_rho_gate)) {
                fin <- is.finite(v) & is.finite(y)
                if (sum(fin) >= min_pairs) {
                    r <- suppressWarnings(stats::cor(v[fin], y[fin], method = "spearman"))
                    if (is.finite(r) && abs(r) >= max_abs_cor) {
                        log_msg("  [covars-%s] drop %s (|rho|=%.3f ≥ %.2f vs biology)", label, cn, abs(r), max_abs_cor)
                        keep_cov[cn] <- FALSE

                        .append_drop("correlation", cn, "high_correlation_with_predictor",
                            metric = abs(r), threshold = max_abs_cor
                        )
                    }
                }
            }
        }
        df <- df[, keep_cov, drop = FALSE]
        if (!ncol(df)) {
            return(NULL)
        }
    }

    # 5) Return: rownames = sample_order
    stopifnot(nrow(df) == length(so), identical(rownames(df), so))

    attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
    df
}


#' Audit covariate coverage for a dataset
#'
#' @param tag Audit tag
#' @param ds_id Dataset ID
#' @param stratum Stratum label
#' @param su Subunit name
#' @param sample_ids Sample IDs
#' @param batch Batch factor or NULL
#' @param covars Covariates data frame or NULL
#' @return NULL (side effect: writes audit CSV)
#' @export
audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL) {
    dir.create(file.path("run_info", "covars_audit"), recursive = TRUE, showWarnings = FALSE)

    cov_df <- if (is.null(covars)) NULL else as.data.frame(covars, check.names = FALSE)
    cov_nms <- if (!is.null(cov_df)) colnames(cov_df) else character(0)

    get_cov <- function(k) {
        if (is.null(cov_df) || !(k %in% names(cov_df))) {
            return(NA_real_)
        }
        v <- cov_df[[k]]
        if (is.numeric(v)) {
            mean(is.finite(v)) * 100
        } else {
            mean(!is.na(v)) * 100
        }
    }
    cov_purity <- get_cov("purity")
    cov_sex <- get_cov("sex")
    cov_age <- get_cov("age")

    line <- sprintf(
        "  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s",
        tag,
        if (length(cov_nms)) paste(cov_nms, collapse = ",") else "NULL",
        cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
        if (!is.null(batch)) nlevels(batch) else 0L
    )
    log_msg(line)

    df <- data.frame(
        dataset = ds_id, stratum = stratum, subunit = su, tag = tag,
        batch_levels = if (!is.null(batch)) nlevels(batch) else 0L,
        batch_sizes = if (!is.null(batch)) {
            paste(sprintf(
                "%s=%d", names(sort(table(batch), decreasing = TRUE)),
                as.integer(sort(table(batch), decreasing = TRUE))
            ), collapse = "; ")
        } else {
            NA_character_
        },
        covars_cols = if (length(cov_nms)) paste(cov_nms, collapse = ";") else "NULL",
        purity_cov = cov_purity, sex_cov = cov_sex, age_cov = cov_age,
        stringsAsFactors = FALSE
    )
    fp <- file.path("run_info", "covars_audit", "audit_rows.csv")
    data.table::fwrite(df, fp, append = file.exists(fp))

    invisible(df)
}


# ---- End of 05_covariate_selection.R ----
