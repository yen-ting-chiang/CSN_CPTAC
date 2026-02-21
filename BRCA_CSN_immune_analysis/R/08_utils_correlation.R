# =============================================================================
# 08_utils_correlation.R
# Correlation Analysis Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# Read Immune Scores from Patient Data
# -----------------------------------------------------------------------------

#' Read BRCA immune/stromal scores from patient clinical data
#'
#' Extracts ESTIMATE, xCell, and CIBERSORT scores from patient data,
#' aligned to provided sample IDs.
#'
#' @param ds_dir Dataset directory path
#' @param sample_ids Character vector of sample IDs
#' @return Data frame with immune score columns
.read_brca_scores_from_patient <- function(ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    if (!file.exists(pat_fp)) {
        message("[BRCA-immune] Cannot find data_clinical_patient.txt: ", ds_dir)
        out <- as.data.frame(matrix(NA_real_, nrow = length(sample_ids), ncol = 0),
            stringsAsFactors = FALSE
        )
        rownames(out) <- sample_ids
        return(out)
    }
    .NN <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))
    samp <- if (file.exists(samp_fp)) suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame() else NULL
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    if (!is.null(samp)) names(samp) <- .NN(names(samp))
    names(pat) <- .NN(names(pat))

    sid <- if (!is.null(samp)) intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1] else NA
    pid_samp <- if (!is.null(samp)) intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1] else NA
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
    map_pt <- NULL
    if (!is.na(sid) && !is.na(pid_samp)) map_pt <- setNames(as.character(samp[[pid_samp]]), as.character(samp[[sid]]))

    keep_cols <- c(
        "ESTIMATE_IMMUNE_SCORE", "ESTIMATE_STROMAL_SCORE",
        "XCELL_IMMUNE_SCORE", "XCELL_STROMAL_SCORE", "CIBERSORT_ABSOLUTE_SCORE"
    )
    out <- matrix(NA_real_,
        nrow = length(sample_ids), ncol = length(keep_cols),
        dimnames = list(sample_ids, keep_cols)
    )
    if (!is.null(pat) && !is.na(pid_pat)) {
        for (nm in intersect(keep_cols, names(pat))) {
            vec_pt <- suppressWarnings(as.numeric(pat[[nm]]))
            names(vec_pt) <- as.character(pat[[pid_pat]])
            if (!is.null(map_pt)) {
                out[, nm] <- vec_pt[map_pt[sample_ids]]
            } else {
                # Fallback: align using patient ID with suffix removed
                pids <- sub("([A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+).*", "\\1", toupper(sample_ids))
                out[, nm] <- vec_pt[pids]
            }
        }
    }
    as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# -----------------------------------------------------------------------------
# Residualize Pair
# -----------------------------------------------------------------------------

#' Residualize a predictor-outcome pair for partial correlation
#'
#' Removes effects of batch and covariates from both x and y.
#'
#' @param x Predictor numeric vector
#' @param y Outcome numeric vector
#' @param batch Optional batch factor
#' @param covars Optional covariate data frame
#' @param min_n Minimum samples required
#' @return List with residualized x and y
.residualize_pair <- function(x, y, batch = NULL, covars = NULL, min_n = 8L) {
    # If residualize_to_covars is defined globally, use it; otherwise use orthogonalize_to
    if (exists("residualize_to_covars", mode = "function")) {
        xr <- residualize_to_covars(x, batch = batch, covars = covars)
        yr <- residualize_to_covars(y, batch = batch, covars = covars)
    } else {
        X <- NULL
        if (!is.null(batch)) X <- cbind(X, stats::model.matrix(~ 0 + batch))
        if (!is.null(covars) && ncol(covars) > 0) X <- cbind(X, as.matrix(covars))
        xr <- orthogonalize_to(x, X)
        yr <- orthogonalize_to(y, X)
    }
    list(x = xr, y = yr)
}

# -----------------------------------------------------------------------------
# Partial Correlation Statistics
# -----------------------------------------------------------------------------

#' Extract partial correlation statistics from linear model fit
#'
#' @param fit Linear model fit object
#' @param term Term name to extract (default "x")
#' @return List with r (partial correlation) and p (p-value)
.extract_partial_stats <- function(fit, term = "x") {
    # If fit failed or singular
    if (inherits(fit, "try-error") || is.na(coef(fit)[term])) {
        return(list(r = NA_real_, p = NA_real_))
    }
    # Get t-statistic and df
    s <- summary(fit)
    coef_tbl <- s$coefficients
    if (!term %in% rownames(coef_tbl)) {
        return(list(r = NA_real_, p = NA_real_))
    }
    t_val <- coef_tbl[term, "t value"]
    p_val <- coef_tbl[term, "Pr(>|t|)"]
    df_resid <- s$df[2]

    # Convert t to partial r
    r <- t_val / sqrt(t_val^2 + df_resid)

    list(r = r, p = p_val)
}
