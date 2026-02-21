## =========================================================
## 06_csn_score.R
## CSN Complex Score Calculation for DPS Analysis
## =========================================================

# Safe version CSN SCORE (try original method first; fallback on failure)
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

    # z-score (keep sample names)
    get_z <- function(v) {
        nm <- names(v) # Store sample names first
        v <- as.numeric(v)
        mu <- mean(v[is.finite(v)], na.rm = TRUE)
        sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
        if (!is.finite(sdv) || sdv == 0) sdv <- 1
        v[!is.finite(v)] <- mu
        out <- (v - mu) / sdv
        names(out) <- nm # Put back sample names
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
        # Direction correction: same sign as average subunit z
        mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
        if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
        s[keep_sam] <- sc
    }

    s
}

## ===== Audit CSN score feasibility =====
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, prot0, present_sub,
                                        min_members = 5L,
                                        min_per_group = get0("min_per_group", ifnotfound = 8L),
                                        out_dir = file.path("run_info", "csn_score_audit")) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Check subunit coverage
    coverage <- rowSums(is.finite(mat0[present_sub, , drop = FALSE]))
    n_samples <- ncol(mat0)

    audit_df <- data.frame(
        ds_id = ds_id,
        stratum = stratum,
        n_samples = n_samples,
        n_subunits = length(present_sub),
        subunits = paste(present_sub, collapse = ";"),
        min_coverage = min(coverage),
        max_coverage = max(coverage),
        mean_coverage = mean(coverage),
        stringsAsFactors = FALSE
    )

    out_fp <- file.path(out_dir, sprintf("csn_audit_%s_%s.csv", ds_id, stratum))
    utils::write.csv(audit_df, out_fp, row.names = FALSE)

    invisible(audit_df)
}

# Safe version audit
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
                min_per_group = min_per_group,
                out_dir = out_dir
            )
        },
        silent = TRUE
    )

    if (inherits(ok, "try-error")) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[audit_csn_score_feasibility_safe] Error: %s", conditionMessage(attr(ok, "condition")))
        }
        return(invisible(NULL))
    }

    invisible(ok)
}
