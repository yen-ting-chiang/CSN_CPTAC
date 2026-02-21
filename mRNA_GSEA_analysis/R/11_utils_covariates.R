# =============================================================================
# 11_utils_covariates.R - Covariate Handling Utilities
# =============================================================================
# This module provides functions for covariate auditing and handling
# in the GSEA analysis pipeline.
# =============================================================================

# -----------------------------------------------------------------------------
# Audit Covariate Coverage
# -----------------------------------------------------------------------------
#' Audit and log covariate coverage for a specific analysis
#' @param tag Tag for logging
#' @param ds_id Dataset ID
#' @param stratum Stratum name
#' @param su Subunit name
#' @param sample_ids Sample IDs
#' @param batch Batch factor (optional)
#' @param covars Covariate data.frame (optional)
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

    # Append one row to summary table
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
}

