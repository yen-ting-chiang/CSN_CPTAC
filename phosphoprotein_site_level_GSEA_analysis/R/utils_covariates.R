# =============================================================================
# Covariate Handling Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for batch factor and covariate processing.
# =============================================================================


build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
    present <- intersect(subunits, rownames(mat0))

    s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
    if (!length(present)) {
        return(s)
    }


    get_z <- function(v) {
        nm <- names(v)
        v <- as.numeric(v)
        mu <- mean(v[is.finite(v)], na.rm = TRUE)
        sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
        if (!is.finite(sdv) || sdv == 0) sdv <- 1
        v[!is.finite(v)] <- mu
        out <- (v - mu) / sdv
        names(out) <- nm
        out
    }

    X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
    rownames(X) <- present
    colnames(X) <- colnames(mat0)

    # combine COPS7A/7B
    if (combine_7AB && all(c("COPS7A", "COPS7B") %in% rownames(X))) {
        Z7 <- colMeans(X[c("COPS7A", "COPS7B"), , drop = FALSE], na.rm = TRUE)
        X <- rbind(X[setdiff(rownames(X), c("COPS7A", "COPS7B")), , drop = FALSE],
            "COPS7*" = Z7
        )
    }


    enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
    keep_sam <- names(s)[enough]

    if (length(keep_sam) >= 10) {
        pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
        sc <- pc$x[, 1]

        mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
        if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
        s[keep_sam] <- sc
    }

    s
}


# -----------------------------------------------------------------------------
# audit_covars_coverage()
# -----------------------------------------------------------------------------
audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL,
                                  sv = NULL,
                                  tech = NULL) {
    dir.create(file.path("run_info", "covars_audit"), recursive = TRUE, showWarnings = FALSE)

    cov_df <- if (is.null(covars)) NULL else as.data.frame(covars, check.names = FALSE)
    cov_nms <- if (!is.null(cov_df)) colnames(cov_df) else character(0)

    get_cov <- function(k) {
        if (is.null(cov_df) || !(k %in% names(cov_df))) {
            return(NA_real_)
        }
        v <- cov_df[[k]]
        if (is.numeric(v)) {
            mean(is.finite(v)) * 100
        } else {
            mean(!is.na(v)) * 100
        }
    }
    cov_purity <- get_cov("purity")
    cov_sex <- get_cov("sex")
    cov_age <- get_cov("age")

    sv_nms <- if (!is.null(sv)) colnames(.ensure_mat_or_null(sv)) else NULL
    tech_nms <- if (!is.null(tech)) colnames(.ensure_mat_or_null(tech)) else NULL

    line <- sprintf(
        "  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s | SV={%s} | tech={%s}",
        tag,
        if (length(cov_nms)) paste(cov_nms, collapse = ",") else "NULL",
        cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
        if (!is.null(batch)) nlevels(batch) else 0L,
        if (length(sv_nms)) paste(sv_nms, collapse = ",") else "none",
        if (length(tech_nms)) paste(tech_nms, collapse = ",") else "none"
    )
    log_msg(line)


    df <- data.frame(
        dataset = ds_id, stratum = stratum, subunit = su, tag = tag,
        batch_levels = if (!is.null(batch)) nlevels(batch) else 0L,
        batch_sizes = if (!is.null(batch)) {
            paste(sprintf(
                "%s=%d", names(sort(table(batch), decreasing = TRUE)),
                as.integer(sort(table(batch), decreasing = TRUE))
            ), collapse = "; ")
        } else {
            NA_character_
        },
        covars_cols = if (length(cov_nms)) paste(cov_nms, collapse = ";") else "NULL",
        sv_cols = if (length(sv_nms)) paste(sv_nms, collapse = ";") else "NULL",
        tech_cols = if (length(tech_nms)) paste(tech_nms, collapse = ";") else "NULL",
        purity_cov = cov_purity, sex_cov = cov_sex, age_cov = cov_age,
        stringsAsFactors = FALSE
    )
    fp <- file.path("run_info", "covars_audit", "audit_rows.csv")
    data.table::fwrite(df, fp, append = file.exists(fp))
}

# -----------------------------------------------------------------------------
# .mk_batch_factor()
# -----------------------------------------------------------------------------
.mk_batch_factor <- function(batch, sample_order, min_count_collapse = 1L) {
    if (is.null(batch)) {
        return(NULL)
    }
    b <- .align_by_colnames(batch, sample_order)
    b <- .take_first(b)
    b[!nzchar(b) | is.na(b)] <- "unknown"
    b <- .sanitize_levels(b)

    if (min_count_collapse > 0L) {
        tb <- table(b)
        small <- names(tb)[tb <= min_count_collapse]
        if (length(small)) {
            b[b %in% small] <- factor("b_small")
            b <- droplevels(factor(b))
        }
    }
    b
}

# -----------------------------------------------------------------------------
# BATCH_PIPE_POLICY and BATCH_MIN_PER_LEVEL
# -----------------------------------------------------------------------------
BATCH_PIPE_POLICY <- "NA"
BATCH_MIN_PER_LEVEL <- 2

# -----------------------------------------------------------------------------
# sanitize_batch_levels()
# -----------------------------------------------------------------------------
sanitize_batch_levels <- function(x,
                                  pipe_policy = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0 %||% "")
    if (any(has_pipe)) {
        log_msg("  [batch] Detected {sum(has_pipe)} values containing '|'; applying policy '{pipe_policy}'")
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }

    fac <- factor(make.names(x0))
    fac <- droplevels(fac)


    if (!is.null(min_per_level) && min_per_level > 1) {
        tab <- table(fac, useNA = "no")
        small <- names(tab)[tab < min_per_level]
        if (length(small)) {
            log_msg(
                "  [batch] Merge sparse levels to 'b_small': %s",
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
# detect_batch_column()
# -----------------------------------------------------------------------------
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

        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            return(list(name = cn, fac = fac))
        }
    }
    NULL
}

# -----------------------------------------------------------------------------
# get_batch_factor()
# -----------------------------------------------------------------------------
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


    meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
    rownames(meta) <- sample_ids

    det <- detect_batch_column(meta,
        pipe_policy   = pipe_policy,
        min_per_level = min_per_level
    )
    if (is.null(det)) {
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

# -----------------------------------------------------------------------------
# get_batch_factor_phospho()
# -----------------------------------------------------------------------------
get_batch_factor_phospho <- function(ds_dir, sample_ids,
                                     pipe_policy = c("strict", "lenient", "NA"),
                                     min_per_level = 3L) {
    pipe_policy <- match.arg(pipe_policy)

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
