# =============================================================================
# 09_utils_batch.R - Batch Effect Handling Utilities
# =============================================================================
# This module provides functions for detecting and handling batch effects
# in proteomics/transcriptomics data.
# =============================================================================

# -----------------------------------------------------------------------------
# Batch Cleanup Settings (defaults, can be overridden)
# -----------------------------------------------------------------------------
if (!exists("BATCH_PIPE_POLICY")) {
    BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": values containing '|' set to NA or merged
}
if (!exists("BATCH_MIN_PER_LEVEL")) {
    BATCH_MIN_PER_LEVEL <- 2 # Levels below this threshold merged to b_small
}

# -----------------------------------------------------------------------------
# Helper: Take First Value Before Pipe
# -----------------------------------------------------------------------------
#' Extract first value before pipe separator
#' @param x Character vector
#' @return Character vector with values before first pipe
.take_first <- function(x) sub("\\|.*$", "", as.character(x))

# -----------------------------------------------------------------------------
# Sanitize Factor Levels
# -----------------------------------------------------------------------------
#' Sanitize factor levels for model compatibility
#' @param f Factor to sanitize
#' @return Factor with cleaned levels
.sanitize_levels <- function(f) {
    if (is.null(f)) {
        return(NULL)
    }
    f <- droplevels(factor(f))
    lv <- levels(f)
    lv2 <- make.names(lv)
    lv2 <- paste0("b_", lv2) # Avoid starting with number
    levels(f) <- lv2
    f
}

# -----------------------------------------------------------------------------
# Align by Column Names
# -----------------------------------------------------------------------------
#' Align vector or data.frame to target sample names
#' @param vec_or_df Vector, factor, or data.frame
#' @param target_names Target sample names
#' @return Aligned vector/data.frame
.align_by_colnames <- function(vec_or_df, target_names) {
    if (is.null(vec_or_df)) {
        return(NULL)
    }
    if (is.vector(vec_or_df) || is.factor(vec_or_df)) {
        if (!is.null(names(vec_or_df))) {
            return(vec_or_df[target_names])
        } else {
            if (length(vec_or_df) == length(target_names)) {
                names(vec_or_df) <- target_names
                return(vec_or_df)
            }
        }
    }
    if (is.data.frame(vec_or_df)) {
        if (!is.null(rownames(vec_or_df))) {
            return(vec_or_df[target_names, , drop = FALSE])
        }
    }
    vec_or_df
}

# -----------------------------------------------------------------------------
# Sanitize Batch Levels
# -----------------------------------------------------------------------------
#' Clean batch factor: handle pipe values, make names valid, merge sparse levels
#' @param x Vector of batch values
#' @param pipe_policy How to handle pipe characters ("NA" or "b_small")
#' @param min_per_level Minimum samples per level before merging
#' @return Cleaned factor
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")
    if (any(has_pipe)) {
        log_msg("  [batch] Detected {sum(has_pipe)} values containing '|' -> processing according to policy {pipe_policy}")
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }
    # Make names valid
    fac <- factor(make.names(x0))
    fac <- droplevels(fac)

    # Merge sparse levels
    if (!is.null(min_per_level) && min_per_level > 1) {
        tab <- table(fac, useNA = "no")
        small <- names(tab)[tab < min_per_level]
        if (length(small)) {
            log_msg(
                "  [batch] Merging sparse levels to 'b_small': %s",
                paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
            )
            fac_chr <- as.character(fac)
            fac_chr[fac_chr %in% small] <- "b_small"
            fac <- droplevels(factor(fac_chr))
        }
    }
    fac
}

# -----------------------------------------------------------------------------
# Detect Batch Column
# -----------------------------------------------------------------------------
#' Auto-detect batch column from clinical metadata
#' @param meta Clinical metadata data.frame
#' @param pipe_policy How to handle pipe characters
#' @param min_per_level Minimum samples per level
#' @return List with name and fac, or NULL
detect_batch_column <- function(meta,
                                pipe_policy = BATCH_PIPE_POLICY,
                                min_per_level = BATCH_MIN_PER_LEVEL) {
    cand <- c("TMT_PLEX", "EXPERIMENT", "PROTEOMICS_TMT_BATCH")
    hit <- intersect(cand, colnames(meta))
    if (!length(hit)) {
        return(NULL)
    }

    for (cn in hit) {
        fac <- sanitize_batch_levels(meta[[cn]],
            pipe_policy   = pipe_policy,
            min_per_level = min_per_level
        )
        # Need at least 2 valid levels and valid sample count >= 3
        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            return(list(name = cn, fac = fac))
        }
    }
    NULL
}

# -----------------------------------------------------------------------------
# Get Batch Factor
# -----------------------------------------------------------------------------
#' Get batch factor for samples from clinical metadata
#' @param ds_dir Dataset directory path
#' @param sample_ids Sample IDs to align
#' @param pipe_policy How to handle pipe characters
#' @param min_per_level Minimum samples per level
#' @return List with name and fac, or NULL
get_batch_factor <- function(ds_dir, sample_ids,
                             pipe_policy = BATCH_PIPE_POLICY,
                             min_per_level = BATCH_MIN_PER_LEVEL) {
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    if (!file.exists(meta_fp)) {
        return(NULL)
    }

    meta <- suppressMessages(
        readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
    ) |> as.data.frame()

    id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
    if (!length(id_cols)) {
        return(NULL)
    }

    meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
    meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)

    # Align by sample_ids
    meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
    rownames(meta) <- sample_ids

    det <- detect_batch_column(meta,
        pipe_policy   = pipe_policy,
        min_per_level = min_per_level
    )
    if (is.null(det)) {
        return(NULL)
    }

    fac <- det$fac
    names(fac) <- sample_ids
    list(name = det$name, fac = fac)
}
