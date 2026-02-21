## ====================================================================
## 03_csn_score.R
##
## Purpose: CSN complex score calculation functions
## Contains: PCA-based CSN score, safe fallback methods, feasibility auditing
## ====================================================================

#' Build CSN score via PCA on CSN subunit expression
#'
#' @param mat0 Protein expression matrix (genes × samples)
#' @param subunits Character vector of CSN subunit genes
#' @param combine_7AB Logical, whether to average COPS7A/COPS7B
#' @param min_members Minimum number of subunits required
#' @return Named numeric vector of CSN scores (PC1)
#' @export
build_csn_score <- function(
  mat0,
  subunits = c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9"),
  combine_7AB = TRUE,
  min_members = 5L
) {
    present <- intersect(subunits, rownames(mat0))
    s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (!length(present)) {
        return(s)
    }

    # Use global .zscore helper for consistency
    X <- do.call(rbind, lapply(present, function(g) .zscore(mat0[g, ])))
    rownames(X) <- present
    colnames(X) <- colnames(mat0)

    # Combined COPS7A/COPS7B if requested
    if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
        Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
        X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
            "COPS7*" = Z7
        )
    }

    # Check sample-level completeness
    enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
    keep_sam <- names(s)[enough]

    if (length(keep_sam) >= 10) {
        pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
        sc <- pc$x[, 1]
        # Direction correction: same sign as subunit average z
        mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
        if (suppressWarnings(stats::cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
        s[keep_sam] <- sc
    }

    s
}


#' Safe version of CSN score (cleanup + fallback)
#'
#' @param mat0 Protein expression matrix
#' @param subunits CSN subunit genes
#' @param combine_7AB Combine COPS7A/B
#' @param min_members Minimum subunits
#' @param pca_min_samples Minimum samples for PCA
#' @return CSN score vector
#' @export
build_csn_score_safe <- function(mat0, subunits, combine_7AB = TRUE,
                                 min_members = 5L, pca_min_samples = 10L) {
    sub <- intersect(subunits, rownames(mat0))
    out_na <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (length(sub) < min_members) {
        return(out_na)
    }

    # Try standard method first
    cs_try <- try(
        {
            build_csn_score(mat0, subunits = sub, combine_7AB = combine_7AB, min_members = min_members)
        },
        silent = TRUE
    )

    if (!inherits(cs_try, "try-error") && sum(is.finite(cs_try)) >= pca_min_samples) {
        return(cs_try)
    }

    # Fallback: clean data and retry
    X <- .clean_for_pca(mat0[sub, , drop = FALSE], min_samples = pca_min_samples, min_genes = min_members)
    if (is.null(X)) {
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Insufficient data → All NA")
        return(out_na)
    }
    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
    if (inherits(pc, "try-error")) {
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed → all NA")
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


#' Cleanup data for PCA (remove Inf→NA, require min finite values, median imputation)
#'
#' @param X Numeric matrix
#' @param min_samples Minimum finite samples per gene
#' @param min_genes Minimum genes retained
#' @return Cleaned matrix or NULL
#' @export
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


#' Audit CSN score feasibility (safe version)
#'
#' @param ds_id Dataset ID
#' @param stratum Stratum label
#' @param mat0 Protein matrix
#' @param present_sub Present CSN subunits
#' @param min_members Minimum subunits
#' @param pca_min_samples Minimum samples for PCA
#' @param min_per_group Minimum samples per group
#' @param out_dir Output directory for audit CSV
#' @return Logical (invisible)
#' @export
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, present_sub,
                                             min_members = 5L, pca_min_samples = 10L,
                                             min_per_group = get0("min_per_group", ifnotfound = 8L),
                                             out_dir = file.path("run_info", "csn_score_audit")) {
    X <- .clean_for_pca(mat0[intersect(present_sub, rownames(mat0)), , drop = FALSE],
        min_samples = pca_min_samples, min_genes = min_members
    )
    if (is.null(X)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[CSN-audit-safe] %s | %s: insufficient data, audit skipped", ds_id, stratum)
        }
        return(invisible(FALSE))
    }
    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
    if (inherits(pc, "try-error")) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[CSN-audit-safe] %s | %s: fallback prcomp failed, skip", ds_id, stratum)
        }
        return(invisible(FALSE))
    }
    varpc1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
    if (exists("log_msg", mode = "function")) {
        log_msg(
            "[CSN-audit-safe] %s | %s: fallback OK; genes=%d; PC1%%=%.1f",
            ds_id, stratum, nrow(X), 100 * varpc1
        )
    }
    invisible(TRUE)
}


# ---- End of 03_csn_score.R ----
