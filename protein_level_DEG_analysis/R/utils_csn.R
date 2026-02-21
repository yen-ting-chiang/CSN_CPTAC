## ============================================================================
## CSN Complex Score Calculation Utilities
## Functions for computing CSN complex activity score using PCA
## ============================================================================

#' Build CSN Complex Score (PC1-based)
#'
#' Calculates a CSN complex activity score by performing PCA on z-scored
#' subunit expression levels and taking PC1 as the complex score.
#'
#' @param mat0 Numeric matrix (genes x samples) with row names as gene symbols
#' @param subunits Character vector of CSN subunit gene names
#' @param combine_7AB Logical; if TRUE, average COPS7A and COPS7B before PCA
#' @param min_members Minimum number of subunits required for PCA
#' @return Named numeric vector (samples) with CSN scores, or all NA if insufficient data
#'
#' @details
#' Algorithm:
#' 1. Z-score normalization of each subunit (cross-sample)
#' 2. Optional: Average COPS7A/COPS7B to single COPS7* feature
#' 3. PCA on z-scored subunit matrix (samples as observations)
#' 4. Take PC1 as CSN score
#' 5. Direction correction: Ensure positive correlation with mean subunit expression
#'
#' Requires:
#' - At least `min_members` subunits present in matrix
#' - At least 10 samples with sufficient coverage
#'
#' @examples
#' csn_score <- build_csn_score(mat, csn_subunits)
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
    present <- intersect(subunits, rownames(mat0))

    # Pre-build return skeleton (Must have names)
    s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))

    if (!length(present)) {
        return(s)
    }

    # Z-score (Keep sample names)
    get_z <- function(v) {
        nm <- names(v) # Save sample names first
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
    colnames(X) <- colnames(mat0) # Ensure sample names exist

    # Optional: Combine COPS7A/7B
    if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
        Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
        X <- rbind(
            X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
            "COPS7*" = Z7
        )
    }

    # Use non-NA count of "original mat0" to judge coverage
    enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
    keep_sam <- names(s)[enough]

    if (length(keep_sam) >= 10) {
        pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
        sc <- pc$x[, 1]

        # Direction correction: Ensure PC1 has same sign as subunit mean z-score
        mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)

        # Check if correlation is valid before using it
        rho <- suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs"))

        if (is.finite(rho) && rho < 0) {
            sc <- -sc
        } else if (!is.finite(rho)) {
            log_msg(
                "[CSN_SCORE] Warning: Cannot compute cor(PC1, mean_z) for direction correction (rho=%s). Using uncorrected PC1.",
                as.character(rho)
            )
        }

        s[keep_sam] <- sc # Alignment by sample name
    }

    s
}

#' Build CSN Score with Fallback (Safe Version)
#'
#' Attempts to build CSN score using standard method, with fallback to
#' cleaned matrix if initial attempt fails.
#'
#' @param mat0 Numeric matrix (genes x samples)
#' @param subunits Character vector of CSN subunit names
#' @param combine_7AB Logical; combine COPS7A/7B
#' @param min_members Minimum subunits required
#' @param pca_min_samples Minimum samples required for PCA
#' @return Named numeric vector with CSN scores
#'
#' @details
#' Fallback procedure:
#' 1. Try build_csn_score()
#' 2. If fails or too few valid samples, clean matrix (remove Inf, impute NA)
#' 3. Run PCA on cleaned matrix
#' 4. Return scores or all NA if still fails
#'
#' @examples
#' csn_score_safe <- build_csn_score_safe(mat, csn_subunits)
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
            build_csn_score(mat0,
                subunits = sub, combine_7AB = combine_7AB,
                min_members = min_members
            )
        },
        silent = TRUE
    )

    if (!inherits(cs_try, "try-error") && sum(is.finite(cs_try)) >= pca_min_samples) {
        return(cs_try)
    }

    # Fallback: Clean matrix for PCA
    X <- .clean_for_pca(mat0[sub, , drop = FALSE],
        min_samples = pca_min_samples,
        min_genes = min_members
    )

    if (is.null(X)) {
        log_msg("[CSN_SCORE-safe] Available subunits or samples insufficient -> All NA")
        return(out_na)
    }

    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)

    if (inherits(pc, "try-error")) {
        log_msg("[CSN_SCORE-safe] prcomp failed -> All NA")
        return(out_na)
    }

    sc <- pc$x[, 1]
    names(sc) <- rownames(pc$x)

    # Direction correction
    ref <- colMeans(X, na.rm = TRUE)
    rr <- suppressWarnings(stats::cor(sc, ref, use = "pairwise.complete.obs"))
    if (is.finite(rr) && rr < 0) sc <- -sc

    out <- out_na
    out[names(sc)] <- as.numeric(sc)

    varpc1 <- if (!is.null(pc$sdev)) (pc$sdev[1]^2) / sum(pc$sdev^2) else NA_real_
    log_msg(
        "[CSN_SCORE-safe] fallback: genes=%d; PC1%%=%.1f; nonNA=%d/%d",
        nrow(X), 100 * varpc1, sum(is.finite(out)), length(out)
    )

    out
}

#' Clean Matrix for PCA
#'
#' Internal helper to prepare matrix for PCA by removing Inf values,
#' filtering low-coverage rows, and imputing missing values.
#'
#' @param X Numeric matrix
#' @param min_samples Minimum finite values required per row
#' @param min_genes Minimum genes required after filtering
#' @return Cleaned matrix or NULL if insufficient data
#'
#' @keywords internal
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

    # Row-wise median imputation for NA
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

#' Audit CSN Score Feasibility
#'
#' Evaluates whether CSN complex score calculation is feasible for a dataset
#' and reports coverage statistics.
#'
#' @param ds_id Dataset identifier
#' @param stratum Stratum name (e.g., "ALL", "TP53_mutant")
#' @param mat0 Protein expression matrix
#' @param present_sub CSN subunits present in matrix
#' @param min_members Minimum subunits required
#' @param pca_min_samples Minimum samples required for PCA
#' @param min_per_group Minimum samples per group for statistical tests
#' @param out_dir Output directory for audit files
#' @return Invisible list with summary, per_subunit, and per_sample tables
#'
#' @details
#' Creates three output files:
#' 1. csn_score_feasibility_summary.csv - Overall feasibility metrics
#' 2. {dataset}_{stratum}_subunit_coverage.csv - Coverage per subunit
#' 3. {dataset}_{stratum}_sample_subunit_counts.csv - Subunits per sample
#'
#' @examples
#' audit_csn_score_feasibility("brca_cptac_2020", "ALL", mat, csn_subunits)
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L,
                                        min_per_group = 8L,
                                        out_dir = file.path("run_info", "csn_score_audit")) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    n_samples <- ncol(mat0)
    if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
        log_msg("[CSN-audit] %s | %s: mat0 structure incomplete, skip", ds_id, stratum)
        return(invisible(NULL))
    }

    present_sub <- intersect(present_sub, rownames(mat0))
    if (!length(present_sub) || n_samples == 0) {
        log_msg(
            "[CSN-audit] %s | %s: No available CSN subunit or samples is 0, skip",
            ds_id, stratum
        )
        return(invisible(NULL))
    }

    # 1) Coverage of each subunit
    sub_cov <- vapply(present_sub, function(g) {
        mean(is.finite(mat0[g, ])) * 100
    }, numeric(1))

    sub_tbl <- data.frame(
        dataset = ds_id,
        stratum = stratum,
        subunit = present_sub,
        nonNA_pct = round(sub_cov, 1),
        nonNA_n = vapply(present_sub, function(g) sum(is.finite(mat0[g, ])), integer(1)),
        stringsAsFactors = FALSE, check.names = FALSE
    )

    cov_min <- if (length(sub_cov)) round(min(sub_cov, na.rm = TRUE), 1) else NA_real_
    cov_med <- if (length(sub_cov)) round(stats::median(sub_cov, na.rm = TRUE), 1) else NA_real_
    cov_max <- if (length(sub_cov)) round(max(sub_cov, na.rm = TRUE), 1) else NA_real_

    # 2) Subunits per sample
    sample_counts <- colSums(is.finite(mat0[present_sub, , drop = FALSE]))
    enough <- sample_counts >= min_members
    n_enough <- sum(enough)

    # 3) Try calculating CSN_SCORE
    csn_score <- build_csn_score(mat0,
        subunits = present_sub,
        combine_7AB = TRUE,
        min_members = min_members
    )
    csn_nonNA <- sum(is.finite(csn_score))
    csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)

    # PC1 explained variance
    pc1_var_pct <- NA_real_
    if (isTRUE(csn_can_pca)) {
        ok_sam <- names(sample_counts)[sample_counts >= min_members]

        get_z <- function(v) {
            v <- as.numeric(v)
            mu <- mean(v[is.finite(v)], na.rm = TRUE)
            sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
            if (!is.finite(sdv) || sdv == 0) sdv <- 1
            v[!is.finite(v)] <- mu
            (v - mu) / sdv
        }

        Z <- do.call(rbind, lapply(present_sub, function(g) {
            get_z(mat0[g, ok_sam, drop = FALSE])
        }))
        rownames(Z) <- present_sub

        if (all(c("COPS7A", "COPS7B") %in% rownames(Z))) {
            Z7 <- colMeans(Z[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
            Z <- rbind(
                Z[setdiff(rownames(Z), c("COPS7A", "COPS7B")), , drop = FALSE],
                "COPS7*" = Z7
            )
        }

        pc <- tryCatch(stats::prcomp(t(Z), center = TRUE, scale. = FALSE),
            error = function(e) NULL
        )

        if (!is.null(pc)) {
            ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
            pc1_var_pct <- round(ve[1], 1)
        }
    }

    # 4) Write files
    sum_row <- data.frame(
        dataset = ds_id, stratum = stratum,
        n_samples = n_samples,
        n_subunits_present = length(present_sub),
        subunit_cov_min = cov_min,
        subunit_cov_median = cov_med,
        subunit_cov_max = cov_max,
        min_members = min_members,
        samples_with_ge_min_members = n_enough,
        pca_min_samples = pca_min_samples,
        csn_score_nonNA = csn_nonNA,
        csn_pc1_feasible = csn_can_pca,
        csn_pc1_var_pct = pc1_var_pct,
        can_form_high_low_groups = (csn_nonNA >= 2 * min_per_group),
        stringsAsFactors = FALSE, check.names = FALSE
    )

    fp_sum <- file.path(out_dir, "csn_score_feasibility_summary.csv")
    data.table::fwrite(sum_row, fp_sum, append = file.exists(fp_sum))

    # Per-subunit coverage
    tag <- paste(ds_id, stratum, sep = "_")
    fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
    data.table::fwrite(sub_tbl, fp_sub)

    # Per-sample subunit counts
    sample_tbl <- data.frame(
        dataset = ds_id, stratum = stratum,
        sample_id = colnames(mat0),
        nonNA_subunits_n = as.integer(sample_counts),
        ge_min_members = as.logical(enough),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    fp_sam <- file.path(out_dir, sprintf("%s_sample_subunit_counts.csv", tag))
    data.table::fwrite(sample_tbl, fp_sam)

    log_msg(
        "[CSN-audit] %s | %s: subunits=%d; min/median/max coverage=%.1f/%.1f/%.1f%%; eligible_samples=%d; CSN_SCORE nonNA=%d; PC1 feasible=%s (PC1%%=%.1f)",
        ds_id, stratum, length(present_sub), cov_min %||% NaN, cov_med %||% NaN, cov_max %||% NaN,
        n_enough, csn_nonNA, as.character(csn_can_pca), pc1_var_pct %||% NaN
    )

    invisible(list(summary = sum_row, per_subunit = sub_tbl, per_sample = sample_tbl))
}

#' Audit CSN Score Feasibility (Safe Version)
#'
#' Wrapper around audit_csn_score_feasibility that catches errors and
#' uses fallback methods if standard audit fails.
#'
#' @inheritParams audit_csn_score_feasibility
#' @return Invisible TRUE/FALSE indicating success
#'
#' @details
#' If standard audit fails, attempts fallback using cleaned matrix.
#' Does not stop execution on failure.
#'
#' @keywords internal
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, present_sub,
                                             min_members = 5L,
                                             pca_min_samples = 10L,
                                             min_per_group = NULL,
                                             out_dir = file.path("run_info", "csn_score_audit")) {
    # Resolve min_per_group safely
    if (is.null(min_per_group)) {
        if (exists("min_per_group", envir = .GlobalEnv, inherits = FALSE)) {
            min_per_group <- get("min_per_group", envir = .GlobalEnv, inherits = FALSE)
        } else {
            min_per_group <- 8L
        }
    }

    ok <- try(
        {
            audit_csn_score_feasibility(
                ds_id = ds_id,
                stratum = stratum,
                mat0 = mat0,
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

    # Fallback attempt
    X <- .clean_for_pca(mat0[intersect(present_sub, rownames(mat0)), , drop = FALSE],
        min_samples = pca_min_samples, min_genes = min_members
    )

    if (is.null(X)) {
        log_msg(
            "[CSN-audit-safe] %s | %s: Available subunits or samples insufficient, skip audit (do not stop)",
            ds_id, stratum
        )
        return(invisible(FALSE))
    }

    pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)

    if (inherits(pc, "try-error")) {
        log_msg(
            "[CSN-audit-safe] %s | %s: fallback prcomp still failed, skip (do not stop)",
            ds_id, stratum
        )
        return(invisible(FALSE))
    }

    varpc1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
    log_msg(
        "[CSN-audit-safe] %s | %s: fallback OK; genes=%d; PC1%%=%.1f",
        ds_id, stratum, nrow(X), 100 * varpc1
    )

    invisible(TRUE)
}
