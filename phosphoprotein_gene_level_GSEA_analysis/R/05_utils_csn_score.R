# =============================================================================
# 05_utils_csn_score.R - CSN Score Calculation Utilities
# =============================================================================
# This file contains functions for calculating CSN (COP9 Signalosome) scores.
# =============================================================================

# -----------------------------------------------------------------------------
# Build CSN Score from Subunit Expression Matrix
# -----------------------------------------------------------------------------
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
    present <- intersect(subunits, rownames(mat0))
    # Pre-build return skeleton (always has names)
    s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (!length(present)) {
        return(s)
    }

    # z-score (preserve sample names)
    get_z <- function(v) {
        nm <- names(v) # first save sample names
        v <- as.numeric(v)
        mu <- mean(v[is.finite(v)], na.rm = TRUE)
        sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
        if (!is.finite(sdv) || sdv == 0) sdv <- 1
        v[!is.finite(v)] <- mu
        out <- (v - mu) / sdv
        names(out) <- nm # put back sample names
        out
    }

    X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
    rownames(X) <- present
    colnames(X) <- colnames(mat0) # this line is important: ensure sample names exist

    # Optional: combine COPS7A/7B
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
        # Direction correction: match subunit average z sign
        mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
        if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
        s[keep_sam] <- sc # here can definitely align by sample names
    }

    s
}

# -----------------------------------------------------------------------------
# Safe Version CSN Score (Try Original Method First; Fallback on Failure)
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
        if (exists("log_msg", mode = "function")) log_msg("[CSN_SCORE-safe] insufficient subunits or samples -> all NA")
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
# Audit CSN Score Feasibility (Safe Version)
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
            log_msg("[CSN-audit-safe] %s | %s: insufficient subunits or samples, skipping audit (not stopping)", ds_id, stratum)
        }
        return(invisible(FALSE))
    }
    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
    if (inherits(pc, "try-error")) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[CSN-audit-safe] %s | %s: fallback prcomp still failed, skipping (not stopping)", ds_id, stratum)
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

# -----------------------------------------------------------------------------
# Audit CSN Score Feasibility
# -----------------------------------------------------------------------------
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, prot0, present_sub,
                                        min_members = 5L, pca_min_samples = 10L,
                                        min_per_group = 8L, out_dir = NULL) {
    # This function performs detailed CSN score feasibility audit
    
    X <- .clean_for_pca(prot0[intersect(present_sub, rownames(prot0)), , drop = FALSE],
        min_samples = pca_min_samples, min_genes = min_members
    )
    if (is.null(X)) {
        stop(sprintf("[audit] %s | %s: insufficient subunits or samples", ds_id, stratum))
    }
    pc <- stats::prcomp(t(X), center = TRUE, scale. = TRUE)
    varpc1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
    if (exists("log_msg", mode = "function")) {
        log_msg("[CSN-audit] %s | %s: genes=%d; PC1%%=%.1f", ds_id, stratum, nrow(X), 100 * varpc1)
    }
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        fn <- file.path(out_dir, sprintf("%s_%s.csv", ds_id, stratum))
        df <- data.frame(ds = ds_id, stratum = stratum, n_genes = nrow(X), pct_pc1 = 100 * varpc1)
        data.table::fwrite(df, fn, append = file.exists(fn))
    }
    invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Residualize Vector (Regress Subunit on CSN Score + Covariates)
# -----------------------------------------------------------------------------
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

    # Remove columns that would eliminate all complete.cases
    ok <- complete.cases(DF_y, DF)
    if (sum(ok) < min_n) {
        return(setNames(rep(NA_real_, length(y)), names(y)))
    }

    # Fit model and extract residuals
    formula_str <- paste("DF_y ~", paste(colnames(DF), collapse = " + "))
    fit <- try(stats::lm(as.formula(formula_str), data = DF, subset = ok), silent = TRUE)
    if (inherits(fit, "try-error")) {
        return(setNames(rep(NA_real_, length(y)), names(y)))
    }

    out <- setNames(rep(NA_real_, length(y)), names(y))
    out[common[ok]] <- stats::residuals(fit)
    out
}
