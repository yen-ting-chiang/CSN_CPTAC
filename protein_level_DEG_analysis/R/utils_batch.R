## ============================================================================
## Batch Effect Detection and Handling Utilities
## Functions for detecting and processing batch information (e.g., TMT plex)
## ============================================================================

#' Sanitize Batch Factor Levels
#'
#' Cleans batch variable values by handling special characters and merging
#' sparse levels to avoid estimation issues in statistical models.
#'
#' @param x Character or factor vector of batch values
#' @param pipe_policy How to handle values containing '|': "NA" or "b_small"
#' @param min_per_level Minimum samples per level; smaller levels merged to "b_small"
#' @return Factor with cleaned and merged levels
#'
#' @details
#' Cleaning steps:
#' 1. Handle '|' characters according to policy
#' 2. Legalize names (make.names) to avoid model.matrix issues
#' 3. Merge sparse levels (< min_per_level) to "b_small"
#'
#' @examples
#' sanitize_batch_levels(c("01", "02", "03|04", "05"), pipe_policy = "NA", min_per_level = 2)
sanitize_batch_levels <- function(x,
                                  pipe_policy = "NA",
                                  min_per_level = 2) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")

    if (any(has_pipe)) {
        log_msg("  [batch] Detected {sum(has_pipe)} values containing '|' -> Handle by policy {pipe_policy}")
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }

    # Legalize names (avoid eBayes/design matrix column name issues)
    fac <- factor(make.names(x0))
    fac <- droplevels(fac)

    # Merge sparse levels (e.g., only 1 sample)
    if (!is.null(min_per_level) && min_per_level > 1) {
        tab <- table(fac, useNA = "no")
        small <- names(tab)[tab < min_per_level]

        if (length(small)) {
            log_msg(
                "  [batch] Merge sparse level to 'b_small': %s",
                paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
            )
            fac_chr <- as.character(fac)
            fac_chr[fac_chr %in% small] <- "b_small"
            fac <- droplevels(factor(fac_chr))
        }
    }

    fac
}

#' Detect Batch Column from Metadata
#'
#' Searches for batch information in sample metadata, preferring TMT_PLEX
#' over EXPERIMENT column.
#'
#' @param meta Data frame of sample metadata
#' @param pipe_policy How to handle '|' in batch values
#' @param min_per_level Minimum samples per batch level
#' @return List with $name (column name) and $fac (factor), or NULL if not found
#'
#' @details
#' Looks for columns in order: TMT_PLEX, EXPERIMENT
#' Returns first valid batch column with:
#' - At least 2 levels
#' - At least 3 non-NA samples
#'
#' @examples
#' detect_batch_column(meta_df)
detect_batch_column <- function(meta,
                                pipe_policy = "NA",
                                min_per_level = 2) {
    cand <- c("TMT_PLEX", "EXPERIMENT")
    hit <- intersect(cand, colnames(meta))

    if (!length(hit)) {
        return(NULL)
    }

    for (cn in hit) {
        fac <- sanitize_batch_levels(meta[[cn]],
            pipe_policy = pipe_policy,
            min_per_level = min_per_level
        )

        # Must have at least 2 valid levels and valid sample count >= 3
        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            return(list(name = cn, fac = fac))
        }
    }

    NULL
}

#' Get Batch Factor by Sample Order
#'
#' Retrieves batch information from sample metadata file and aligns it
#' to the specified sample order. Includes fallback to TMT_protein.csv
#' if metadata batch column not found.
#'
#' @param ds_dir Path to dataset directory
#' @param sample_ids Character vector of sample IDs (in desired order)
#' @param pipe_policy How to handle '|' in batch values
#' @param min_per_level Minimum samples per batch level
#' @return List with $name (source) and $fac (factor aligned to sample_ids), or NULL
#'
#' @details
#' Search order:
#' 1. data_clinical_sample.txt (TMT_PLEX or EXPERIMENT column)
#' 2. TMT_protein.csv (Run Metadata ID column)
#'
#' The returned factor is aligned to sample_ids order with names set.
#'
#' @examples
#' batch_info <- get_batch_factor("brca_cptac_2020", colnames(mat))
get_batch_factor <- function(ds_dir, sample_ids,
                             pipe_policy = "NA",
                             min_per_level = 2) {
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")

    if (!file.exists(meta_fp)) {
        return(NULL)
    }

    meta <- suppressMessages(
        readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
    ) |> as.data.frame()

    # Identify sample ID column
    id_cols <- intersect(
        c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"),
        names(meta)
    )
    if (!length(id_cols)) {
        return(NULL)
    }

    meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
    meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)

    # Align by sample_ids first, ensure detected and returned vector lengths are consistent
    meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
    rownames(meta) <- sample_ids

    # Try to detect batch column from metadata
    det <- detect_batch_column(meta,
        pipe_policy = pipe_policy,
        min_per_level = min_per_level
    )

    if (is.null(det)) {
        ## Fallback: Infer TMT-plex from TMT_protein.csv as batch
        tmt_fp <- file.path(ds_dir, "TMT_protein.csv")

        if (file.exists(tmt_fp)) {
            tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |>
                as.data.frame()

            # Normalize column names (trim whitespace; handle " Run Metadata ID")
            cn <- names(tmt)
            cn_trim <- trimws(cn)
            names(tmt) <- cn_trim

            run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
            tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

            if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
                run_col <- run_hits[1]

                # Build sample_id -> plex mapping
                plex_by_sample <- list()
                nR <- nrow(tmt)

                for (i in seq_len(nR)) {
                    run_id <- as.character(tmt[[run_col]][i])
                    if (!nzchar(run_id) || is.na(run_id)) next

                    # Take first two digits as plex ID
                    plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id)
                    if (!nzchar(plex2) || is.na(plex2)) next

                    for (tc in tmt_cols) {
                        cell <- tmt[[tc]][i]
                        if (is.na(cell)) next

                        cell <- as.character(cell)
                        if (!nzchar(cell)) next

                        # Take before newline
                        sid <- sub("\\r?\\n.*$", "", cell)
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
                        pipe_policy = pipe_policy,
                        min_per_level = min_per_level
                    )

                    # Success condition: At least two valid levels and valid sample count >= 3
                    if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
                        names(fac2) <- sample_ids
                        return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
                    }
                }
            }
        }

        # Fallback also failed -> Return NULL
        return(NULL)
    }

    # Return detected batch from metadata
    fac <- det$fac
    names(fac) <- sample_ids
    list(name = det$name, fac = fac)
}

#' Screen Dataset for Batch Effect Necessity
#'
#' Evaluates whether batch adjustment is recommended for a dataset based on
#' PCA variance explained and gene-level F-tests.
#'
#' @param ds_dir Path to dataset directory
#' @param min_frac_complete Minimum fraction of non-NA values per gene
#' @return Invisible list with batch evaluation metrics
#'
#' @details
#' Recommendation criteria:
#' - PC1-5 R² >= 10% and p < 0.01, OR
#' - Gene-level F-test FDR < 0.05 proportion >= 5%
#'
#' Writes evaluation results to log and returns metrics invisibly.
#'