# =============================================================================
# CSN Score Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for building CSN scores.
# =============================================================================

# -----------------------------------------------------------------------------
# .clean_for_pca()
# -----------------------------------------------------------------------------
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
# build_csn_score_safe()
# -----------------------------------------------------------------------------
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
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] Insufficient available subcells or samples → All NA")
        return(out_na)
    }
    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
    if (inherits(pc, "try-error")) {
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] prcomp failed → All NA")
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
# audit_csn_score_feasibility_safe()
# -----------------------------------------------------------------------------
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, prot0, present_sub,
                                             min_members = 5L, pca_min_samples = 10L,
                                             min_per_group = get0("min_per_group", ifnotfound = 8L),
                                             out_dir = file.path("run_info", "csn_score_audit")) {
    ok <- try(
        {
            audit_csn_score_feasibility(
                ds_id = ds_id,
                stratum = stratum,
                mat0 = mat0,
                prot0 = prot0,
                present_sub = present_sub,
                min_members = min_members,
                pca_min_samples = pca_min_samples,
                min_per_group = min_per_group,
                out_dir = out_dir
            )
            TRUE
        },
        silent = TRUE
    )
    if (!inherits(ok, "try-error")) {
        return(invisible(TRUE))
    }

    X <- .clean_for_pca(prot0[intersect(present_sub, rownames(prot0)), , drop = FALSE],
        min_samples = pca_min_samples, min_genes = min_members
    )
    if (is.null(X)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[CSN-audit-safe] %s | %s: insufficient available subunits or samples, audit skipped.", ds_id, stratum)
        }
        return(invisible(FALSE))
    }
    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
    if (inherits(pc, "try-error")) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[CSN-audit-safe] %s | %s: fallback prcomp still failed, skipped (without aborting).", ds_id, stratum)
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
