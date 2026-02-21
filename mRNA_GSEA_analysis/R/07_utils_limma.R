# =============================================================================
# 07_utils_limma.R - Limma Analysis Utilities
# =============================================================================
# This module provides utilities for limma differential expression analysis
# including design matrix alignment and covariate handling.
# =============================================================================

# -----------------------------------------------------------------------------
# Align Limma Inputs
# -----------------------------------------------------------------------------
#' Align expression matrix and design data frame for limma analysis
#' @param M Expression matrix (genes x samples)
#' @param design_df Design data frame with rownames as sample IDs
#' @param rhs_terms Right-hand side terms for model formula
#' @param label Label for error messages
#' @return List with aligned M, des (design matrix), and design_df
.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
    stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))

    design_df <- as.data.frame(design_df, check.names = FALSE)
    stopifnot(!is.null(rownames(design_df)))
    samp_M <- as.character(colnames(M))
    samp_des <- as.character(rownames(design_df))

    # Remove NA samples from design table
    na_mask <- rep(FALSE, nrow(design_df))
    rn <- rownames(design_df)

    # predictor must exist
    if (!("predictor" %in% colnames(design_df))) {
        stop(sprintf("[%s] design_df missing predictor column", label))
    }

    # predictor not finite
    v <- design_df$predictor
    na_mask <- na_mask | is.na(v) | !is.finite(v)

    # Other numeric columns
    num_cols <- setdiff(colnames(design_df), c("predictor"))
    for (cn in num_cols) {
        v <- design_df[[cn]]
        if (is.numeric(v)) {
            na_mask <- na_mask | is.na(v) | !is.finite(v)
        } else {
            na_mask <- na_mask | is.na(v)
        }
    }
    design_df <- design_df[!na_mask, , drop = FALSE]

    # Intersect and align order
    common <- intersect(as.character(colnames(M)), as.character(rownames(design_df)))
    if (!length(common)) {
        return(NULL)
    }
    common <- sort(common)

    M2 <- M[, common, drop = FALSE]
    design_df2 <- design_df[common, , drop = FALSE]

    # Build design matrix
    rhs <- unique(rhs_terms)
    des2 <- stats::model.matrix(stats::reformulate(rhs), data = design_df2)

    # Final check
    if (nrow(des2) != ncol(M2)) {
        msg <- sprintf(
            "[%s] nrow(des)=%d != ncol(M)=%d (still inconsistent after alignment)",
            label, nrow(des2), ncol(M2)
        )
        stop(msg)
    }
    list(M = M2, des = des2, design_df = design_df2)
}

# -----------------------------------------------------------------------------
# Audit Design Alignment
# -----------------------------------------------------------------------------
#' Audit and log design matrix alignment with expression matrix
#' @param tag Tag for logging
#' @param samples Sample names from matrix
#' @param mod_interest Interest model matrix
#' @param mod_nuisance Nuisance model matrix
#' @param out_dir Optional output directory for CSV
#' @return Invisible alignment data frame
audit_design_alignment <- function(tag, samples, mod_interest, mod_nuisance, out_dir = NULL) {
    safe_rn <- function(X) if (is.null(X)) rep(NA_character_, length(samples)) else rownames(as.matrix(X))
    df <- data.frame(
        pos = seq_along(samples),
        M_sample = samples,
        interest_rn = safe_rn(mod_interest),
        nuisance_rn = safe_rn(mod_nuisance),
        interest_ok = if (!is.null(mod_interest)) safe_rn(mod_interest) == samples else NA,
        nuisance_ok = if (!is.null(mod_nuisance)) safe_rn(mod_nuisance) == samples else NA,
        stringsAsFactors = FALSE
    )
    n_i <- sum(df$interest_ok, na.rm = TRUE)
    n_n <- sum(df$nuisance_ok, na.rm = TRUE)
    log_msg(
        "[align:%s] interest match = %d/%d, nuisance match = %d/%d, both = %d",
        tag, n_i, length(samples), n_n, length(samples), sum(df$interest_ok & df$nuisance_ok, na.rm = TRUE)
    )
    if (!is.null(out_dir)) {
        fn <- file.path(out_dir, sprintf("ALIGN_%s.csv", gsub("[^A-Za-z0-9]+", "_", tag)))
        utils::write.csv(df, fn, row.names = FALSE)
        log_msg("[align:%s] wrote %s", tag, fn)
    }
    head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
    log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
    invisible(df)
}

# -----------------------------------------------------------------------------
# Orthogonalize Matrix
# -----------------------------------------------------------------------------
#' Orthogonalize matrix Y to nuisance matrix X (remove effect of X from Y)
#' @param mat Matrix to orthogonalize
#' @param nuisance Nuisance matrix
#' @return Residual matrix
orthogonalize_to <- function(mat, nuisance) {
    if (is.null(mat) || is.null(nuisance)) {
        return(mat)
    }
    Y <- as.matrix(mat)
    X <- as.matrix(nuisance)

    # Add intercept
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)

    # Handle non-finite values + remove zero-variance columns
    for (j in seq_len(ncol(X))) {
        v <- X[, j]
        if (any(!is.finite(v))) {
            mu <- mean(v[is.finite(v)], na.rm = TRUE)
            if (!is.finite(mu)) mu <- 0
            v[!is.finite(v)] <- mu
            X[, j] <- v
        }
    }
    var_ok <- apply(X, 2, function(z) {
        z <- as.numeric(z)
        is.finite(var(z)) && var(z) > 0
    })
    if (!all(var_ok)) X <- X[, var_ok, drop = FALSE]

    beta <- tryCatch(qr.coef(qr(X), Y), error = function(e) NULL)
    if (is.null(beta)) {
        return(as.matrix(mat))
    }

    resid <- Y - X %*% beta
    rownames(resid) <- rownames(Y)
    colnames(resid) <- colnames(Y)
    return(resid)
}

# -----------------------------------------------------------------------------
# Coerce Covariates Safely
# -----------------------------------------------------------------------------
#' Safely coerce covariate data frame, dropping problematic columns
#' @param df Data frame of covariates
#' @return Cleaned data frame
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
                keep[cn] <- FALSE
                if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
            } else {
                df[[cn]] <- v
            }
        } else {
            df[[cn]] <- suppressWarnings(as.numeric(v))
        }
    }
    df <- df[, keep, drop = FALSE]
    df
}
