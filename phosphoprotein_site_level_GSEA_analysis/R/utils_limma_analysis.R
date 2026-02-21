# =============================================================================
# Limma Analysis Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for limma input alignment and analysis.
# =============================================================================


# -----------------------------------------------------------------------------
# .ensure_mat_or_null()
# -----------------------------------------------------------------------------
if (!exists(".ensure_mat_or_null", mode = "function")) {
    .ensure_mat_or_null <- function(x) {
        if (is.null(x)) {
            return(NULL)
        }
        m <- tryCatch(as.matrix(x), error = function(e) NULL)
        if (is.null(m)) {
            return(NULL)
        }
        if (is.null(dim(m))) { # vector -> single-column matrix
            m <- matrix(m, ncol = 1)
        }
        if (ncol(m) == 0) {
            return(NULL)
        }
        m
    }
}


# -----------------------------------------------------------------------------
# .take_first()
# -----------------------------------------------------------------------------
.take_first <- function(x) sub("\\|.*$", "", as.character(x))

# -----------------------------------------------------------------------------
# .sanitize_levels()
# -----------------------------------------------------------------------------
.sanitize_levels <- function(f) {
    if (is.null(f)) {
        return(NULL)
    }
    f <- droplevels(factor(f))
    lv <- levels(f)
    lv2 <- make.names(lv)
    lv2 <- paste0("b_", lv2)
    levels(f) <- lv2
    f
}

# -----------------------------------------------------------------------------
# .align_by_colnames()
# -----------------------------------------------------------------------------
.align_by_colnames <- function(vec_or_df, target_names) {
    if (is.null(vec_or_df)) {
        return(NULL)
    }
    if (is.vector(vec_or_df) || is.factor(vec_or_df)) {
        nm <- names(vec_or_df)
        if (!is.null(nm) && length(nm) == length(vec_or_df)) {
            return(vec_or_df[match(target_names, nm)])
        } else {
            return(vec_or_df)
        }
    } else {
        rn <- rownames(vec_or_df)
        if (!is.null(rn)) {
            return(vec_or_df[match(target_names, rn), , drop = FALSE])
        }
        return(vec_or_df)
    }
}
