# =============================================================================
# 05_utils_csn_score.R - CSN Score Calculation Utilities
# =============================================================================
# This module provides functions for calculating CSN complex activity score
# based on PC1 of subunit expression, and related utilities.
# =============================================================================

# -----------------------------------------------------------------------------
# Clean Matrix for PCA
# -----------------------------------------------------------------------------
#' Prepare matrix for PCA by removing Inf and imputing NA
#' @param X Numeric matrix
#' @param min_samples Minimum non-NA values per row
#' @param min_genes Minimum genes after filtering
#' @return Cleaned matrix or NULL
.clean_for_pca <- function(X, min_samples = 10L, min_genes = 5L) {
    X <- as.matrix(X)
    X[!is.finite(X)] <- NA
    keep_rows <- rowSums(is.finite(X)) >= min_samples
    if (!any(keep_rows)) {
        return(NULL)
    }
    X <- X[keep_rows, , drop = FALSE]
    if (nrow(X) < min_genes) {
        return(NULL)
    }
    if (anyNA(X)) {
        med <- apply(X, 1, function(r) median(r[is.finite(r)], na.rm = TRUE))
        for (i in seq_len(nrow(X))) {
            xi <- X[i, ]
            xi[!is.finite(xi)] <- med[i]
            X[i, ] <- xi
        }
    }
    X
}

# -----------------------------------------------------------------------------
# Z-score Helper
# -----------------------------------------------------------------------------
#' Z-score a vector, preserving names
#' @param v Numeric vector with names
#' @return Z-scored vector with same names
.get_z <- function(v) {
    nm <- names(v)
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm
    out
}

# -----------------------------------------------------------------------------
# Build CSN Score (PC1-based)
# -----------------------------------------------------------------------------
#' Calculate CSN complex activity score using PC1 of subunit expression
#' @param mat0 Expression matrix (genes x samples)
#' @param subunits Character vector of CSN subunit gene names
#' @param combine_7AB Whether to average COPS7A and COPS7B
#' @param min_members Minimum subunits with data required
#' @return Named numeric vector of CSN scores
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
    present <- intersect(subunits, rownames(mat0))
    # Pre-create return skeleton (always has names)
    s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (!length(present)) {
        return(s)
    }

    # Z-score each subunit
    X <- do.call(rbind, lapply(present, function(g) .get_z(mat0[g, ])))
    rownames(X) <- present
    colnames(X) <- colnames(mat0)

    # Optional: merge COPS7A/7B
    if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
        Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
        X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
            "COPS7*" = Z7
        )
    }

    # Use "original mat0" non-NA count to determine coverage
    enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
    keep_sam <- names(s)[enough]

    if (length(keep_sam) >= 10) {
        pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
        sc <- pc$x[, 1]
        # Direction correction: same sign as subunit mean z
        mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
        if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
        s[keep_sam] <- sc
    }

    s
}

# -----------------------------------------------------------------------------
# Build CSN Score (Safe Version)
# -----------------------------------------------------------------------------
#' Safe wrapper for CSN score calculation with fallback
#' @param mat0 Expression matrix
#' @param subunits CSN subunit names
#' @param combine_7AB Whether to average COPS7A/7B
#' @param min_members Minimum subunits required
#' @param pca_min_samples Minimum samples for PCA
#' @return Named numeric vector of CSN scores
build_csn_score_safe <- function(mat0, subunits, combine_7AB = TRUE,
                                 min_members = 5L, pca_min_samples = 10L) {
    sub <- intersect(subunits, rownames(mat0))
    out_na <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (length(sub) < min_members) {
        return(out_na)
    }

    cs_try <- try(
        {
            build_csn_score(mat0, subunits = sub, combine_7AB = combine_7AB, min_members = min_members)
        },
        silent = TRUE
    )

    if (!inherits(cs_try, "try-error") && sum(is.finite(cs_try)) >= pca_min_samples) {
        return(cs_try)
    }

    X <- .clean_for_pca(mat0[sub, , drop = FALSE], min_samples = pca_min_samples, min_genes = min_members)
    if (is.null(X)) {
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] insufficient available subunits or samples -> all NA")
        return(out_na)
    }
    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
    if (inherits(pc, "try-error")) {
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed -> all NA")
        return(out_na)
    }
    sc <- pc$x[, 1]
    names(sc) <- rownames(pc$x)
    ref <- colMeans(X, na.rm = TRUE)
    rr <- suppressWarnings(stats::cor(sc, ref, use = "pairwise.complete.obs"))
    if (is.finite(rr) && rr < 0) sc <- -sc
    out <- out_na
    out[names(sc)] <- as.numeric(sc)
    if (exists("log_msg", mode = "function")) {
        varpc1 <- if (!is.null(pc$sdev)) (pc$sdev[1]^2) / sum(pc$sdev^2) else NA_real_
        log_msg(
            "[CSN_SCORE-safe] fallback: genes=%d; PC1%%=%.1f; nonNA=%d/%d",
            nrow(X), 100 * varpc1, sum(is.finite(out)), length(out)
        )
    }
    out
}

# -----------------------------------------------------------------------------
# Residualize Vector
# -----------------------------------------------------------------------------
#' Regress a subunit against CSN score and covariates, return residuals
#' @param y Named numeric vector (subunit expression)
#' @param csn_score Named numeric vector (CSN score)
#' @param batch Optional batch factor
#' @param covars Optional data.frame of covariates
#' @param min_n Minimum samples required
#' @return Named numeric vector of residuals (z-scored)
residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8L) {
    # 1) Align sample names
    if (is.null(names(y))) stop("[residualize_vector] y must be a numeric vector with sample names")
    sam <- names(y)
    common <- sam
    if (!is.null(csn_score)) common <- intersect(common, names(csn_score))
    if (!is.null(batch)) common <- intersect(common, names(batch))
    if (!is.null(covars)) {
        rn <- rownames(as.data.frame(covars, check.names = FALSE))
        if (is.null(rn)) stop("[residualize_vector] covars must have rownames=samples")
        common <- intersect(common, rn)
    }
    if (length(common) < min_n) {
        return(setNames(rep(NA_real_, length(y)), names(y)))
    }

    # 2) Build data frame
    DF <- data.frame(
        csn = as.numeric(csn_score[common]),
        row.names = common, check.names = FALSE
    )
    if (!is.null(batch)) {
        bf <- droplevels(batch[common])
        DF[["batch"]] <- bf
    }
    if (!is.null(covars)) {
        C <- as.data.frame(covars[common, , drop = FALSE], check.names = FALSE)
        C <- coerce_covariates_safely(C)
        for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
    }
    DF_y <- suppressWarnings(as.numeric(y[common]))

    # Drop all-NA columns
    all_na_col <- vapply(DF, function(v) all(is.na(v)), logical(1))
    if (any(all_na_col)) {
        if (exists("log_msg", mode = "function")) {
            try(log_msg(
                "  [resid] drop all-NA columns: %s",
                paste(names(all_na_col)[all_na_col], collapse = ",")
            ), silent = TRUE)
        }
        DF <- DF[, !all_na_col, drop = FALSE]
    }

    # Drop single-level factors or constant numeric columns
    is_const <- vapply(DF, function(v) {
        vv <- v[!is.na(v)]
        if (!length(vv)) {
            TRUE
        } else {
            if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
        }
    }, logical(1))
    if (any(is_const)) {
        if (exists("log_msg", mode = "function")) {
            try(log_msg(
                "  [resid] drop constant/single-level columns: %s",
                paste(names(is_const)[is_const], collapse = ",")
            ), silent = TRUE)
        }
        DF <- DF[, !is_const, drop = FALSE]
    }

    ok <- is.finite(DF_y) & stats::complete.cases(DF)
    if (sum(ok) < min_n) {
        out <- setNames(rep(NA_real_, length(y)), names(y))
        return(out)
    }
    DF <- DF[ok, , drop = FALSE]
    y_ok <- DF_y[ok]

    # 3) Expand design matrix
    des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
    stopifnot(nrow(des) == length(y_ok))

    # 4) Linear regression
    fit <- lm.fit(x = des, y = y_ok)
    res <- rep(NA_real_, length(DF_y))
    res[ok] <- fit$residuals

    # 5) Map back to full samples and z-score
    out <- setNames(rep(NA_real_, length(y)), names(y))
    out[common] <- res
    fin <- is.finite(out)
    if (sum(fin) >= 3L) {
        mu <- mean(out[fin])
        sdv <- stats::sd(out[fin])
        if (!is.finite(sdv) || sdv == 0) sdv <- 1
        out[fin] <- (out[fin] - mu) / sdv
    }
    out
}
