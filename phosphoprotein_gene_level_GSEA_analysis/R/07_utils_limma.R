# =============================================================================
# 07_utils_limma.R - Limma Analysis Utilities
# =============================================================================
# This file contains functions for limma model fitting and alignment.
# =============================================================================

# -----------------------------------------------------------------------------
# Align Limma Inputs (Matrix, Design DataFrame, Create Design Matrix)
# -----------------------------------------------------------------------------
.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
    stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))
    # First convert design table to data.frame, ensure rownames are sample IDs
    design_df <- as.data.frame(design_df, check.names = FALSE)
    stopifnot(!is.null(rownames(design_df)))
    samp_M <- as.character(colnames(M))
    samp_des <- as.character(rownames(design_df))

    # First remove NA samples from design table (predictor / numeric covariate / factor NA are considered unusable)
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
            # Factor/string/logical: remove if NA
            na_mask <- na_mask | is.na(v)
        }
    }
    design_df <- design_df[!na_mask, , drop = FALSE]

    # Take intersection and sort
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

    # Final check again
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
# Audit Design Matrix Alignment
# -----------------------------------------------------------------------------
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
    # Extra: list first 6 entries for quick log review
    head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
    log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
    invisible(df)
}
