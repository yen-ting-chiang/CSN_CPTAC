## ====================================================================
## 04_batch_detection.R
##
## Purpose: Batch effect detection and correction functions
## Contains: TMT plex detection, batch sanitization, batch need assessment
## ====================================================================

#' Sanitize batch levels (handle pipe characters, merge sparse levels)
#'
#' @param x Character/factor vector
#' @param pipe_policy "NA" or "b_small" for values containing "|"
#' @param min_per_level Levels below this threshold merged to "b_small"
#' @return Factor with sanitized levels
#' @export
sanitize_batch_levels <- function(x,
                                  pipe_policy = get0("BATCH_PIPE_POLICY", ifnotfound = "NA"),
                                  min_per_level = get0("BATCH_MIN_PER_LEVEL", ifnotfound = 2)) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")
    if (any(has_pipe)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("  [batch] Detects {sum(has_pipe)} values containing '|' → policy {pipe_policy}")
        }
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }

    fac <- factor(make.names(x0))
    fac <- droplevels(fac)

    # Merge sparse levels to "b_small"
    if (!is.null(min_per_level) && min_per_level > 1) {
        tab <- table(fac, useNA = "no")
        small <- names(tab)[tab < min_per_level]
        if (length(small)) {
            if (exists("log_msg", mode = "function")) {
                log_msg(
                    "  [batch] Merge sparse levels to 'b_small': %s",
                    paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse = ", ")
                )
            }
            fac_chr <- as.character(fac)
            fac_chr[fac_chr %in% small] <- "b_small"
            fac <- droplevels(factor(fac_chr))
        }
    }
    fac
}


#' Detect batch column from clinical metadata
#'
#' @param meta Data frame with clinical metadata
#' @param pipe_policy Pipe handling policy
#' @param min_per_level Minimum samples per level
#' @return List with name and fac, or NULL
#' @export
detect_batch_column <- function(meta,
                                pipe_policy = get0("BATCH_PIPE_POLICY", ifnotfound = "NA"),
                                min_per_level = get0("BATCH_MIN_PER_LEVEL", ifnotfound = 2)) {
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
        # Must have >= 2 levels and >= 3 valid samples
        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            return(list(name = cn, fac = fac))
        }
    }
    NULL
}


#' Get batch factor for samples (from clinical data or TMT_protein.csv)
#'
#' @param ds_dir Dataset directory
#' @param sample_ids Character vector of sample IDs
#' @param pipe_policy Pipe handling policy
#' @param min_per_level Minimum samples per level
#' @return List with name and fac, or NULL
#' @export
get_batch_factor <- function(ds_dir, sample_ids,
                             pipe_policy = get0("BATCH_PIPE_POLICY", ifnotfound = "NA"),
                             min_per_level = get0("BATCH_MIN_PER_LEVEL", ifnotfound = 2)) {
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

    # Align to sample_ids
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

    ## === Fallback: Extract TMT plex from TMT_protein.csv ===
    tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
    if (file.exists(tmt_fp)) {
        tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()

        # Normalize column names
        cn <- names(tmt)
        cn_trim <- trimws(cn)
        names(tmt) <- cn_trim
        run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
        tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

        if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
            run_col <- run_hits[1]

            # Create sample_id → plex mapping
            plex_by_sample <- list()
            nR <- nrow(tmt)
            for (i in seq_len(nR)) {
                run_id <- as.character(tmt[[run_col]][i])
                if (!nzchar(run_id) || is.na(run_id)) next
                plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id) # First two digits
                if (!nzchar(plex2) || is.na(plex2)) next

                for (tc in tmt_cols) {
                    cell <- tmt[[tc]][i]
                    if (is.na(cell)) next
                    cell <- as.character(cell)
                    if (!nzchar(cell)) next
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
                    pipe_policy   = pipe_policy,
                    min_per_level = min_per_level
                )
                if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
                    names(fac2) <- sample_ids
                    return(list(name = "TMT_protein.csv:RunMetadataID", fac = fac2))
                }
            }
        }
    }

    # No batch found
    NULL
}


#' Screen whether batch correction is needed (PCA + limma F-test)
#'
#' @param ds_dir Dataset directory
#' @param min_frac_complete Minimum gene completeness
#' @return Invisible list with recommendation
#' @export
screen_batch_need <- function(ds_dir, min_frac_complete = 0.75) {
    if (exists("log_msg", mode = "function")) {
        log_msg("== Batch check: %s ==", basename(ds_dir))
    }

    mat0 <- load_matrix_from_dataset_dir(ds_dir)
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
    m <- impute_and_filter(mat0, min_frac = min_frac_complete)

    bi <- get_batch_factor(ds_dir, colnames(m))
    if (is.null(bi)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("  [batch] No batch column found → No correction")
        }
        return(invisible(list(
            dataset = basename(ds_dir),
            batch_col = NA, n_levels = 0, sizes = NA,
            pc_R2 = NA, pc_p = NA, frac_FDR05 = NA, frac_FDR25 = NA,
            recommend = FALSE
        )))
    }

    batch <- droplevels(bi$fac[colnames(m)])
    tab <- sort(table(batch), decreasing = TRUE)
    if (exists("log_msg", mode = "function")) {
        log_msg(
            "  [batch] column: %s; levels=%d; sizes: %s",
            bi$name, nlevels(batch),
            paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
        )
    }

    # PCA analysis
    X <- t(scale(t(m), center = TRUE, scale = TRUE))
    pc <- prcomp(t(X), scale. = FALSE)
    K <- min(5, ncol(pc$x))
    r2 <- vapply(seq_len(K), function(i) summary(lm(pc$x[, i] ~ batch))$r.squared, numeric(1))
    pv <- vapply(seq_len(K), function(i) {
        a <- anova(lm(pc$x[, i] ~ batch))
        as.numeric(a$`Pr(>F)`[1])
    }, numeric(1))
    if (exists("log_msg", mode = "function")) {
        log_msg(
            "  PCA (by batch) R²: %s ; p: %s",
            paste(round(r2, 3), collapse = ", "),
            paste(signif(pv, 3), collapse = ", ")
        )
    }

    # Limma F-test
    design <- model.matrix(~batch)
    fit <- limma::lmFit(m, design)
    fit <- limma::eBayes(fit)
    padj <- p.adjust(fit$F.p.value, "BH")
    prop05 <- mean(padj < 0.05, na.rm = TRUE)
    prop25 <- mean(padj < 0.25, na.rm = TRUE)
    if (exists("log_msg", mode = "function")) {
        log_msg(
            "  gene-level F test: FDR<0.05 = %.1f%%; FDR<0.25 = %.1f%%",
            100 * prop05, 100 * prop25
        )
    }

    recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
    if (recommend) {
        if (exists("log_msg", mode = "function")) {
            log_msg("  **Correction recommended**: R² or gene ratio exceeds threshold")
        }
    } else {
        if (exists("log_msg", mode = "function")) {
            log_msg("  **No correction needed**: No significant batch effect")
        }
    }

    invisible(list(
        dataset = basename(ds_dir),
        batch_col = bi$name, n_levels = nlevels(batch),
        sizes = tab, pc_R2 = r2, pc_p = pv,
        frac_FDR05 = prop05, frac_FDR25 = prop25,
        recommend = recommend
    ))
}


# ---- End of 04_batch_detection.R ----
