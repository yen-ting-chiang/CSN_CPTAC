# =============================================================================
# 05_utils_batch.R
# Batch Effect Handling Utilities
# =============================================================================

# -----------------------------------------------------------------------------
# Batch Level Sanitization
# -----------------------------------------------------------------------------

#' Sanitize batch level values
#'
#' Handles pipe-separated values, sanitizes names for model matrix,
#' and merges sparse levels to avoid non-estimable coefficients.
#'
#' @param x Character vector of batch labels
#' @param pipe_policy How to handle pipe '|' characters: "NA" or "b_small"
#' @param min_per_level Minimum samples per level (smaller merged to b_small)
#' @return Factor with cleaned batch levels
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")
    if (any(has_pipe)) {
        log_msg("  [batch] Detected {sum(has_pipe)} values containing '|' -> processing by policy {pipe_policy}")
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }
    # Sanitize names (to avoid eBayes/design matrix column name issues)
    fac <- factor(make.names(x0))
    fac <- droplevels(fac)

    # Merge sparse levels (e.g. those with only 1 sample)
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
# Batch Column Detection
# -----------------------------------------------------------------------------

#' Auto-detect batch column in clinical metadata
#'
#' Searches for common batch column names (TMT_PLEX, EXPERIMENT) and
#' validates that they have sufficient levels.
#'
#' @param meta Data frame of clinical metadata
#' @param pipe_policy How to handle pipe characters
#' @param min_per_level Minimum samples per level
#' @return List with name and fac, or NULL if not found
detect_batch_column <- function(meta,
                                pipe_policy = BATCH_PIPE_POLICY,
                                min_per_level = BATCH_MIN_PER_LEVEL) {
    cand <- c("TMT_PLEX", "EXPERIMENT")
    hit <- intersect(cand, colnames(meta))
    if (!length(hit)) {
        return(NULL)
    }

    for (cn in hit) {
        fac <- sanitize_batch_levels(meta[[cn]],
            pipe_policy   = pipe_policy,
            min_per_level = min_per_level
        )
        # Must have at least 2 valid levels and valid sample count >= 3
        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            return(list(name = cn, fac = fac))
        }
    }
    NULL
}

# -----------------------------------------------------------------------------
# Get Batch Factor
# -----------------------------------------------------------------------------

#' Get batch factor aligned to sample order
#'
#' Reads clinical metadata, detects batch column, and aligns to provided
#' sample IDs. Falls back to TMT_protein.csv if main detection fails.
#'
#' @param ds_dir Dataset directory path
#' @param sample_ids Character vector of sample IDs
#' @param pipe_policy How to handle pipe characters
#' @param min_per_level Minimum samples per level
#' @return List with name and fac, or NULL if no batch detected
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

    # First align by sample_ids, ensure detection and return vector length consistent
    meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
    rownames(meta) <- sample_ids

    det <- detect_batch_column(meta,
        pipe_policy   = pipe_policy,
        min_per_level = min_per_level
    )
    if (is.null(det)) {
        ## Fallback: derive TMT-plex as batch from TMT_protein.csv
        tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
        if (file.exists(tmt_fp)) {
            tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()

            ## Normalize column names (remove leading/trailing spaces)
            cn <- names(tmt)
            cn_trim <- trimws(cn)
            names(tmt) <- cn_trim
            run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
            tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

            if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
                run_col <- run_hits[1]

                ## Build sample_id -> plex mapping
                plex_by_sample <- list()
                nR <- nrow(tmt)
                for (i in seq_len(nR)) {
                    run_id <- as.character(tmt[[run_col]][i])
                    if (!nzchar(run_id) || is.na(run_id)) next
                    plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id) # take first two digits
                    if (!nzchar(plex2) || is.na(plex2)) next

                    for (tc in tmt_cols) {
                        cell <- tmt[[tc]][i]
                        if (is.na(cell)) next
                        cell <- as.character(cell)
                        if (!nzchar(cell)) next
                        sid <- sub("\\r?\\n.*$", "", cell) # take content before newline
                        sid <- trimws(sid)
                        if (!nzchar(sid)) next
                        if (is.null(plex_by_sample[[sid]])) plex_by_sample[[sid]] <- plex2
                    }
                }

                if (length(plex_by_sample)) {
                    v <- rep(NA_character_, length(sample_ids))
                    names(v) <- sample_ids
                    mm <- intersect(names(plex_by_sample), sample_ids)
                    if (length(mm)) v[mm] <- unlist(plex_by_sample[mm], use.names = FALSE)

                    fac2 <- sanitize_batch_levels(v,
                        pipe_policy   = pipe_policy,
                        min_per_level = min_per_level
                    )
                    ## Success condition: at least 2 valid levels and valid sample count >= 3
                    if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
                        names(fac2) <- sample_ids
                        return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
                    }
                }
            }
        }
        ## fallback also failed -> return NULL
        return(NULL)
    }

    fac <- det$fac # already aligned with sample_ids
    names(fac) <- sample_ids
    list(name = det$name, fac = fac)
}
