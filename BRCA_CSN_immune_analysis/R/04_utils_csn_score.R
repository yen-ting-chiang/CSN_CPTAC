# =============================================================================
# 04_utils_csn_score.R
# CSN Complex Score Building and Residualization Utilities
# =============================================================================

# -----------------------------------------------------------------------------
# Covariate Processing
# -----------------------------------------------------------------------------

#' Safely coerce covariates for linear modeling
#'
#' Converts character/logical to factor, drops single-level factors,
#' keeps numeric as numeric.
#'
#' @param df Data frame of covariates
#' @return Processed data frame safe for model.matrix
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
                # single-level factor -> drop to avoid contrasts error
                keep[cn] <- FALSE
                if (exists("log_msg", mode = "function")) {
                    try(log_msg("  [covars] drop single-level factor: %s", cn), silent = TRUE)
                }
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

# -----------------------------------------------------------------------------
# PCA Pre-cleaning
# -----------------------------------------------------------------------------

#' Pre-clean matrix for PCA
#'
#' Removes Inf values (converts to NA), filters rows with insufficient
#' finite values, and performs row-wise median imputation for NA.
#'
#' @param X Input matrix
#' @param min_samples Minimum finite values per row to keep
#' @param min_genes Minimum rows required
#' @return Cleaned matrix or NULL if insufficient data
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
# CSN Score Building
# -----------------------------------------------------------------------------

#' Build CSN complex score (PC1) from subunit expression
#'
#' Uses z-scores of subunits across samples for PCA, takes PC1 as CSN score,
#' and corrects direction to match mean z sign.
#'
#' @param mat0 Expression matrix (genes x samples)
#' @param subunits Character vector of subunit gene names
#' @param combine_7AB If TRUE, average COPS7A and COPS7B
#' @param min_members Minimum subunits with finite values per sample
#' @return Named numeric vector of CSN scores
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
    present <- intersect(subunits, rownames(mat0))
    # Pre-build return skeleton (must have names)
    s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (!length(present)) {
        return(s)
    }

    # z-score (preserve sample names)
    get_z <- function(v) {
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

    X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
    rownames(X) <- present
    colnames(X) <- colnames(mat0)

    # Optional: merge COPS7A/7B
    if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
        Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
        X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
            "COPS7*" = Z7
        )
    }

    # Use non-NA count from original mat0 to determine coverage
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

#' Safe CSN score building with fallback
#'
#' Tries build_csn_score first; if it fails, uses a simpler PCA approach
#' with scale. = TRUE and more robust handling.
#'
#' @param mat0 Expression matrix
#' @param subunits Character vector of subunit gene names
#' @param combine_7AB If TRUE, average COPS7A and COPS7B
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
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Insufficient subunits or samples -> all NA")
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
# Residualization
# -----------------------------------------------------------------------------

#' Regress vector on CSN score + batch + covariates, return residuals
#'
#' @param y Named numeric vector (e.g., a subunit expression)
#' @param csn_score Named numeric vector of CSN scores
#' @param batch Named factor of batch labels (or NULL)
#' @param covars Data frame of additional covariates (or NULL)
#' @param min_n Minimum samples required for regression
#' @return Named numeric vector of z-scored residuals
residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8L) {
    # 1) Align sample names
    if (is.null(names(y))) stop("[residualize_vector] y must be a named numeric vector")
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

    # 2) Build data frame -> single model.matrix expansion
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

    # 2a) Drop all-NA columns
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

    # 2b) Drop single-level factors or constant numeric columns
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

    # 4) Linear regression y ~ csn + batch + covars
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

#' Orthogonalize matrix columns to nuisance variables
#'
#' @param mat Data matrix or vector to orthogonalize
#' @param nuisance Matrix of nuisance variables
#' @return Residualized matrix
orthogonalize_to <- function(mat, nuisance) {
    if (is.null(mat) || is.null(nuisance)) {
        return(mat)
    }
    Y <- as.matrix(mat)
    X <- as.matrix(nuisance)

    # 1) Add intercept
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)

    # 2) Handle non-finite values + remove zero-variance columns
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
