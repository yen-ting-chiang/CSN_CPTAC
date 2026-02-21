## ============================================================================
## Statistical Utilities and Covariate Processing
## Functions for covariate handling, model alignment, and statistical tools
## ============================================================================

#' Coerce Covariates Safely
#'
#' Converts covariates to appropriate types and removes problematic columns
#' that would cause issues in statistical models.
#'
#' @param df Data frame of covariates
#' @return Cleaned data frame with single-level factors removed
#'
#' @details
#' Processing steps:
#' 1. Converts character/logical to factor
#' 2. Keeps numeric as numeric
#' 3. Removes single-level factors (would cause contrasts error)
#' 4. Optionally removes constant numeric columns
#'
#' @examples
#' covars_clean <- coerce_covariates_safely(covars_df)
coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
        v <- df[[cn]]

        if (is.factor(v) || is.character(v) || is.logical(v)) {
            v <- factor(v)
            lv <- levels(droplevels(v[!is.na(v)]))

            if (length(lv) <= 1) {
                # Single-level factor -> drop, avoid contrasts error
                keep[cn] <- FALSE
                log_msg("  [covars] drop single-level factor: %s", cn)
            } else {
                df[[cn]] <- v
            }
        } else {
            # Keep numeric as numeric
            df[[cn]] <- suppressWarnings(as.numeric(v))
        }
    }

    df <- df[, keep, drop = FALSE]
    df
}

#' Audit Covariate Coverage
#'
#' Reports coverage statistics for covariates and batch variables, appending
#' summary to a CSV file for tracking.
#'
#' @param tag Label for this audit
#' @param ds_id Dataset ID
#' @param stratum Stratum name
#' @param su Subunit/predictor name
#' @param sample_ids Character vector of sample IDs
#' @param batch Batch factor (optional)
#' @param covars Data frame of covariates (optional)
#' @param sv Surrogate variables matrix (optional, legacy)
#' @param tech Technical PCs matrix (optional, legacy)
#'
#' @details
#' Reports:
#' - Covariate names and coverage percentages (purity, sex, age)
#' - Number of batch levels and sample counts per level
#' - SV and tech PC column names
#'
#' Appends one row to run_info/covars_audit/audit_rows.csv
#'
#' @examples
#' audit_covars_coverage(
#'     "brca_ALL_CSN", "brca_cptac_2020", "ALL", "Z_CSN_SCORE",
#'     sample_ids, batch, covars_df
#' )
audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL,
                                  sv = NULL,
                                  tech = NULL) {
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
            mean(!is.na(v)) * 100 # Factor/String -> Count by non-NA ratio
        }
    }

    cov_purity <- get_cov("purity")
    cov_sex <- get_cov("sex")
    cov_age <- get_cov("age")

    sv_nms <- if (!is.null(sv)) colnames(.ensure_mat_or_null(sv)) else NULL
    tech_nms <- if (!is.null(tech)) colnames(.ensure_mat_or_null(tech)) else NULL

    line <- sprintf(
        "  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s | SV={%s} | tech={%s}",
        tag,
        if (length(cov_nms)) paste(cov_nms, collapse = ",") else "NULL",
        cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
        if (!is.null(batch)) nlevels(batch) else 0L,
        if (length(sv_nms)) paste(sv_nms, collapse = ",") else "none",
        if (length(tech_nms)) paste(tech_nms, collapse = ",") else "none"
    )
    log_msg(line)

    # Append summary row
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
        sv_cols = if (length(sv_nms)) paste(sv_nms, collapse = ";") else "NULL",
        tech_cols = if (length(tech_nms)) paste(tech_nms, collapse = ";") else "NULL",
        purity_cov = cov_purity, sex_cov = cov_sex, age_cov = cov_age,
        stringsAsFactors = FALSE
    )

    fp <- file.path("run_info", "covars_audit", "audit_rows.csv")
    data.table::fwrite(df, fp, append = file.exists(fp))
}

#' Ensure Matrix or NULL
#'
#' Internal helper to safely convert object to matrix, returning NULL
#' if conversion fails or result is empty.
#'
#' @param x Object to convert
#' @return Matrix or NULL
#'
#' @keywords internal
.ensure_mat_or_null <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m) || length(m) == 0) {
        return(NULL)
    }
    if (is.null(dim(m))) m <- matrix(m, ncol = 1)
    if (nrow(m) == 0 || ncol(m) == 0) {
        return(NULL)
    }
    m[!is.finite(m)] <- NA_real_
    m
}
