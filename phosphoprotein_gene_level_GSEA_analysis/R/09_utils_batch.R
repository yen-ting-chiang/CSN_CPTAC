# =============================================================================
# 09_utils_batch.R - Batch Correction Utilities
# =============================================================================
# This file contains functions for batch detection and handling.
# =============================================================================

# -----------------------------------------------------------------------------
# Batch Cleanup Settings
# -----------------------------------------------------------------------------
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": values containing '|' should be set to NA or merged
BATCH_MIN_PER_LEVEL <- 2 # levels below this threshold will be merged to b_small

# -----------------------------------------------------------------------------
# Force Serial Execution for Reproducibility
# -----------------------------------------------------------------------------
.force_serial_execution <- function() {
    # base R / data.table
    options(mc.cores = 1L)
    if ("package:data.table" %in% search() &&
        exists("setDTthreads", asNamespace("data.table"))) {
        data.table::setDTthreads(1L)
    }

    # BiocParallel
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::register(BiocParallel::SerialParam())
    }

    # foreach / doParallel
    if (requireNamespace("foreach", quietly = TRUE)) {
        foreach::registerDoSEQ()
    }

    # future
    if (requireNamespace("future", quietly = TRUE)) {
        future::plan(future::sequential)
    }
    Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")

    # fgsea: some versions support nproc, here we just set an option for wrapper to read
    options(.fgsea_nproc = 1L)
}

# -----------------------------------------------------------------------------
# Sanitize Batch Levels
# -----------------------------------------------------------------------------
sanitize_batch_levels <- function(x, pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")
    if (any(has_pipe)) {
        log_msg("  [batch] Detected {sum(has_pipe)} values containing '|' -> processing with policy {pipe_policy}")
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }
    # Legalize names (avoid eBayes/design matrix column name issues)
    fac <- factor(make.names(x0))
    fac <- droplevels(fac)

    # Merge sparse levels (e.g. only 1 sample)
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
# Detect Batch Column from Clinical Metadata
# -----------------------------------------------------------------------------
detect_batch_column <- function(meta, pipe_policy = BATCH_PIPE_POLICY,
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
# Get Batch Factor for Phosphoprotein Data
# -----------------------------------------------------------------------------
get_batch_factor_phospho <- function(ds_dir, sample_ids,
                                     pipe_policy = c("strict", "lenient", "NA"),
                                     min_per_level = 3L) {
    pipe_policy <- match.arg(pipe_policy)
    # First use clinical detection (consistent with protein version)
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    if (file.exists(meta_fp)) {
        meta <- suppressMessages(
            readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
        ) |> as.data.frame()
        id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
        if (length(id_cols)) {
            meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
            meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)
            # Consistent with protein version: use match to align and preserve sample_ids order and length
            meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
            rownames(meta) <- sample_ids
            det <- detect_batch_column(meta,
                pipe_policy   = pipe_policy,
                min_per_level = min_per_level
            )
            if (!is.null(det)) {
                fac <- det$fac
                names(fac) <- sample_ids
                return(list(name = det$name, fac = fac))
            }
        }
    }

    # FALLBACK: Read TMT_phos.csv
    tmt_fp <- file.path(ds_dir, "TMT_phos.csv")
    if (!file.exists(tmt_fp)) {
        return(NULL)
    }

    tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()
    cn <- names(tmt)
    cn_trim <- trimws(cn)
    names(tmt) <- cn_trim
    run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
    tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)
    if (!length(run_hits) || !length(tmt_cols)) {
        return(NULL)
    }

    run_col <- run_hits[1]
    plex_by_sample <- list()
    for (i in seq_len(nrow(tmt))) {
        run_id <- as.character(tmt[[run_col]][i])
        if (!nzchar(run_id) || is.na(run_id)) next
        plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id) # take first two digits
        if (!nzchar(plex2) || is.na(plex2)) next

        for (tc in tmt_cols) {
            cell <- tmt[[tc]][i]
            if (is.na(cell)) next
            cell <- as.character(cell)
            if (!nzchar(cell)) next
            sid <- sub("\\r?\\n.*$", "", cell) # take before newline
            sid <- trimws(sid)
            if (!nzchar(sid)) next
            if (is.null(plex_by_sample[[sid]])) plex_by_sample[[sid]] <- plex2
        }
    }

    if (!length(plex_by_sample)) {
        return(NULL)
    }
    v <- rep(NA_character_, length(sample_ids))
    names(v) <- sample_ids
    mm <- intersect(names(plex_by_sample), sample_ids)
    if (length(mm)) v[mm] <- unlist(plex_by_sample[mm], use.names = FALSE)

    fac2 <- sanitize_batch_levels(v,
        pipe_policy   = pipe_policy,
        min_per_level = min_per_level
    )
    if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
        names(fac2) <- sample_ids
        return(list(name = "TMT_phos.csv:RunMetadataID", fac = fac2))
    }
    return(NULL)
}

# -----------------------------------------------------------------------------
# Remove Confounding Effects from Matrix
# -----------------------------------------------------------------------------
.remove_confounding <- function(mat, nuisance) {
    if (is.null(nuisance)) {
        return(as.matrix(mat))
    }
    Y <- as.matrix(mat)
    X <- as.matrix(nuisance)

    # Add intercept
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)

    # Handle non-finite values + remove zero-variance columns
    for (j in seq_len(ncol(X))) {
        v <- X[, j]
        if (any(!is.finite(v))) {
            mu <- mean(v[is.finite(v)], na.rm = TRUE)
            if (!is.finite(mu)) mu <- 0
            v[!is.finite(v)] <- mu
            X[, j] <- v
        }
    }
    var_ok <- apply(X, 2, function(z) {
        z <- as.numeric(z)
        is.finite(var(z)) && var(z) > 0
    })
    if (!all(var_ok)) X <- X[, var_ok, drop = FALSE]

    beta <- tryCatch(qr.coef(qr(X), Y), error = function(e) NULL)
    if (is.null(beta)) {
        return(as.matrix(mat))
    }

    resid <- Y - X %*% beta
    rownames(resid) <- rownames(Y)
    colnames(resid) <- colnames(Y)
    return(resid)
}
