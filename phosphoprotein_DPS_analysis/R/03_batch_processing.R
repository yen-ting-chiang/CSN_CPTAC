## =========================================================
## 03_batch_processing.R
## Batch Detection and Processing for DPS Analysis
## =========================================================

## ===== Batch cleanup settings =====
BATCH_PIPE_POLICY <- "NA" # "NA" or "b_small": values with '|' should be set to NA or merged
BATCH_MIN_PER_LEVEL <- 2 # Levels below this threshold will be merged to b_small

## ===== Batch value cleanup: handle '|', legalize names, merge sparse levels =====
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")
    if (any(has_pipe)) {
        log_msg("  [batch] Detected {sum(has_pipe)} values with '|' -> processing per policy {pipe_policy}")
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }
    # Legalize names
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

## ===== Auto-detect batch column =====
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
        # Need at least 2 valid levels and valid sample count >= 3
        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            return(list(name = cn, fac = fac))
        }
    }
    NULL
}

## ===== Get batch factor for protein data =====
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
        # Fallback: infer TMT-plex from TMT_protein.csv
        tmt_fp <- file.path(ds_dir, "TMT_protein.csv")
        if (file.exists(tmt_fp)) {
            tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()

            cn <- names(tmt)
            cn_trim <- trimws(cn)
            names(tmt) <- cn_trim
            run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn_trim, ignore.case = TRUE, value = TRUE)
            tmt_cols <- grep("^tmt_", cn_trim, ignore.case = TRUE, value = TRUE)

            if (length(run_hits) >= 1 && length(tmt_cols) >= 1) {
                run_col <- run_hits[1]

                plex_by_sample <- list()
                nR <- nrow(tmt)
                for (i in seq_len(nR)) {
                    run_id <- as.character(tmt[[run_col]][i])
                    if (!nzchar(run_id) || is.na(run_id)) next
                    plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id)
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
        return(NULL)
    }

    fac <- det$fac
    names(fac) <- sample_ids
    list(name = det$name, fac = fac)
}

## ===== Get batch factor for phospho data =====
get_batch_factor_phospho <- function(ds_dir, sample_ids,
                                     pipe_policy = c("strict", "lenient", "NA"),
                                     min_per_level = 3L) {
    pipe_policy <- match.arg(pipe_policy)

    # First use clinical detection
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    if (file.exists(meta_fp)) {
        meta <- suppressMessages(
            readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
        ) |> as.data.frame()
        id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
        if (length(id_cols)) {
            meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
            meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)
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
        plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id)
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

## ===== Screen batch need =====
screen_batch_need <- function(ds_dir, min_frac_complete = 0.75) {
    log_msg("== Batch check: %s ==", basename(ds_dir))
    mat0 <- load_phospho_matrix_from_dataset_dir(ds_dir)
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
    m <- impute_and_filter(mat0, min_frac = min_frac_complete)

    bi <- get_batch_factor_phospho(ds_dir, colnames(m))
    if (is.null(bi)) {
        log_msg("  [batch] Cannot find clear batch column -> no correction for now")
        return(invisible(list(
            dataset = basename(ds_dir),
            batch_col = NA, n_levels = 0, sizes = NA,
            pc_R2 = NA, pc_p = NA, frac_FDR05 = NA, frac_FDR25 = NA,
            recommend = FALSE
        )))
    }

    batch <- droplevels(bi$fac[colnames(m)])
    tab <- sort(table(batch), decreasing = TRUE)
    log_msg(
        "  [batch] Column: %s; levels=%d; n per level: %s",
        bi$name, nlevels(batch),
        paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
    )

    # Safe PCA
    X <- .clean_for_pca(m, min_samples = 10L, min_genes = 5L)
    if (is.null(X)) {
        r2 <- pv <- numeric(0)
        log_msg("  PCA (by batch) skipped: insufficient genes/samples or still contains non-finite values")
    } else {
        Xs <- t(scale(t(X), center = TRUE, scale = TRUE))
        pc <- stats::prcomp(t(Xs), scale. = FALSE)
        K <- min(5L, ncol(pc$x))
        r2 <- vapply(seq_len(K), function(i) summary(lm(pc$x[, i] ~ batch))$r.squared, numeric(1))
        pv <- vapply(seq_len(K), function(i) {
            a <- anova(lm(pc$x[, i] ~ batch))
            as.numeric(a$`Pr(>F)`[1])
        }, numeric(1))
        log_msg(
            "  PCA (by batch) R2: %s ; p: %s", paste(round(r2, 3), collapse = ", "),
            paste(signif(pv, 3), collapse = ", ")
        )
    }

    design <- model.matrix(~batch)
    fit <- limma::lmFit(m, design)
    fit <- limma::eBayes(fit)
    padj <- p.adjust(fit$F.p.value, "BH")
    prop05 <- mean(padj < 0.05, na.rm = TRUE)
    prop25 <- mean(padj < 0.25, na.rm = TRUE)
    log_msg(
        "  Gene-level F test: FDR<0.05 proportion = %.1f%%; FDR<0.25 = %.1f%%",
        100 * prop05, 100 * prop25
    )

    recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
    if (recommend) {
        log_msg("  **Recommend correction**: R2 or gene proportion reached threshold (>=10%% R2 and p<0.01, or FDR<0.05 genes >=5%%)")
    } else {
        log_msg("  **No correction needed for now**: No obvious batch effect detected (record kept)")
    }

    invisible(list(
        dataset = basename(ds_dir),
        batch_col = bi$name, n_levels = nlevels(batch),
        sizes = tab, pc_R2 = r2, pc_p = pv,
        frac_FDR05 = prop05, frac_FDR25 = prop25,
        recommend = recommend
    ))
}
