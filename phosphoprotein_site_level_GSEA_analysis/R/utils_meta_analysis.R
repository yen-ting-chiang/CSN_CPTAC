# =============================================================================
# Meta-Analysis Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for meta-analysis via Stouffer's method.
# =============================================================================

# -----------------------------------------------------------------------------
# safe_read_gsea()
# -----------------------------------------------------------------------------
safe_read_gsea <- function(fp) {
    if (!file.exists(fp)) {
        return(NULL)
    }
    dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) {
        return(NULL)
    }
    need <- c("pathway", "NES", "pval")
    if (!all(need %in% names(dt))) {
        return(NULL)
    }
    as.data.frame(dt[, need])
}

# -----------------------------------------------------------------------------
# meta_fdr_stouffer()
# -----------------------------------------------------------------------------
meta_fdr_stouffer <- function(dataset_dirs,
                              strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                              stat_tags = c("GSEA_limma_t_cont"),
                              groups = names(genesets_by_group),
                              out_root = if (is.null(COMBO_PREFIX)) "phospho_csn_gsea_pan_summary_TP53/meta_fdr" else file.path(COMBO_PREFIX, "phospho_csn_gsea_pan_summary_TP53/meta_fdr")) {
    stopifnot(length(dataset_dirs) > 0)

    if (!is.null(names(dataset_dirs))) {
        keep_idx <- names(dataset_dirs) != "lusc_cptac_2021"
    } else {
        keep_idx <- basename(dataset_dirs) != "lusc_cptac_2021"
    }
    dataset_dirs <- dataset_dirs[keep_idx]
    if (!length(dataset_dirs)) {
        message("[meta] All datasets were excluded or no datasets were available (lusc_cptac_2021 was excluded); omitted.")
        return(invisible(TRUE))
    }
    if (!requireNamespace("data.table", quietly = TRUE)) stop("install data.table first.")
    if (!requireNamespace("stats", quietly = TRUE)) stop("Missing stats suite")
    data.table::setDTthreads(1L)

    # ---- helpers --------------------------------------------------------------
    normalize_cols <- function(DT) {
        nm <- names(DT)

        low <- tolower(nm)
        map <- c(
            "pathway" = "pathway",
            "term" = "pathway",
            "pathway_name" = "pathway",
            "nes" = "NES",
            "enrichment_score" = "NES",
            "pval" = "pval",
            "p.value" = "pval",
            "pvalue" = "pval",
            "p" = "pval",
            "padj" = "padj",
            "fdr" = "padj",
            "qval" = "padj",
            "size" = "size",
            "setsize" = "size",
            "n" = "size"
        )

        for (i in seq_along(nm)) {
            key <- low[i]
            if (key %in% names(map)) data.table::setnames(DT, nm[i], map[[key]])
        }
        DT
    }

    pick_cols <- function(DT, need) {
        need <- unique(need)
        if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)

        miss <- setdiff(need, names(DT))
        for (m in miss) DT[, (m) := NA]

        DT[, ..need]
    }

    stouffer_df <- function(D) {
        D <- D[is.finite(pval) & is.finite(NES)]
        if (nrow(D) == 0) {
            return(NULL)
        }

        D[, z := sign(NES) * stats::qnorm(p = pmax(pmin(pval / 2, 1 - 1e-100), 1e-100), lower.tail = FALSE)]

        Out <- D[, .(k = .N, Z = sum(z) / sqrt(.N)), by = .(pathway)]
        Out[, p_meta := 2 * stats::pnorm(abs(Z), lower.tail = FALSE)]
        Out[, padj_meta := p.adjust(p_meta, method = "BH")]
        data.table::setorder(Out, p_meta)
        Out[]
    }

    # ---- main loop ------------------------------------------------------------
    passes <- c("BatchAdj")
    base_out <- out_root
    dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

    for (stratum in strata) {
        for (pass_label in passes) {
            for (stat_tag in stat_tags) {
                for (grp in groups) {
                    subunits <- unique(unlist(lapply(dataset_dirs, function(dsdir) {
                        base <- file.path(
                            COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection,
                            safe_fs_name(grp), pass_label, basename(dsdir), stratum
                        )
                        if (!dir.exists(base)) {
                            return(character(0))
                        }

                        subs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
                        subs[subs != ""]
                    })))

                    if (length(subunits) == 0) {
                        message(sprintf(
                            "[meta] %s/%s/%s/%sNo subunit available, skip this step.",
                            stratum, pass_label, stat_tag, grp
                        ))
                        next
                    }

                    # 2) For each subunit, collect CSVs across datasets, run Stouffer meta-analysis, and write output
                    for (su in subunits) {
                        parts <- list()
                        for (dsdir in dataset_dirs) {
                            f <- file.path(
                                COMBO_PREFIX %||% phospho_site_csn_gsea_results_TP53_by_collection,
                                safe_fs_name(grp), pass_label, basename(dsdir), stratum, su,
                                sprintf("%s.csv", stat_tag)
                            )
                            if (!file.exists(f)) next
                            dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
                            if (is.null(dt) || nrow(dt) == 0) next
                            dt <- normalize_cols(dt)
                            dt <- pick_cols(dt, need = c("pathway", "NES", "pval", "padj", "size"))
                            dt[, dataset := basename(dsdir)]
                            parts[[length(parts) + 1L]] <- dt
                        }

                        if (length(parts) == 0) {
                            next
                        }

                        D <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

                        D <- D[!is.na(pathway)]
                        res <- stouffer_df(D)
                        if (is.null(res)) next


                        grp_safe <- safe_fs_name(grp)
                        out_dir <- file.path(base_out, stratum, pass_label, su, grp_safe)
                        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                        out_file <- file.path(out_dir, paste0(stat_tag, "_meta_fdr.csv"))
                        data.table::fwrite(res, out_file)
                        message(sprintf(
                            "[meta] %s/%s/%s/%s -> %s (n=%d pathways, k>=1 datasets)",
                            stratum, pass_label, su, stat_tag, out_file, nrow(res)
                        ))
                    } # su
                } # grp
            } # stat_tag
        } # pass
    } # stratum

    invisible(TRUE)
}

# =============================================================================
# Meta-FDR Summary Functions
# Summarize meta-FDR (Stouffer → BH) across subunits
# =============================================================================

# -----------------------------------------------------------------------------
# .read_meta_fdr_table() - Read meta-FDR table for a single subunit
# -----------------------------------------------------------------------------
.read_meta_fdr_table <- function(meta_root, stratum, version, subunit,
                                 group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first.")
    grp <- safe_fs_name(group_name)
    fp <- file.path(meta_root, stratum, version, subunit, grp, paste0(stat_tag, ".csv"))
    if (!file.exists(fp)) {
        return(NULL)
    }
    dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
    if (is.null(dt)) {
        return(NULL)
    }
    need <- c("pathway", "Z", "padj_meta")

    nms <- tolower(names(dt))
    if (!"pathway" %in% names(dt) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
    if (!"Z" %in% names(dt) && "z" %in% nms) names(dt)[match("z", nms)] <- "Z"
    if (!"padj_meta" %in% names(dt) && "padj_meta" %in% nms) names(dt)[match("padj_meta", nms)] <- "padj_meta"
    if (!all(need %in% names(dt))) {
        return(NULL)
    }
    out <- as.data.frame(dt[, ..need])
    names(out)[names(out) == "Z"] <- paste0("Z_", subunit)
    names(out)[names(out) == "padj_meta"] <- paste0("padj_meta_", subunit)
    out
}

# -----------------------------------------------------------------------------
# .merge_subunit_tables_meta() - Merge tables from multiple subunits
# -----------------------------------------------------------------------------
.merge_subunit_tables_meta <- function(tbl_list) {
    keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
    if (!length(keep)) {
        return(NULL)
    }
    out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
    num_cols <- setdiff(names(out), "pathway")
    out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
    out
}

# -----------------------------------------------------------------------------
# .add_sig_counts_meta() - Add significance count columns
# -----------------------------------------------------------------------------
.add_sig_counts_meta <- function(df, alphas = c(0.05, 0.25)) {
    if (is.null(df) || !nrow(df)) {
        return(df)
    }
    padj_cols <- grep("^padj_meta_", names(df), value = TRUE)
    subs <- sub("^padj_meta_", "", padj_cols)
    z_cols <- paste0("Z_", subs)
    keep_idx <- z_cols %in% names(df)
    padj_cols <- padj_cols[keep_idx]
    z_cols <- z_cols[keep_idx]
    if (!length(padj_cols)) {
        return(df)
    }

    for (a in alphas) {
        a_tag <- gsub("\\.", "_", sprintf("%.2f", a))
        sig_mat <- sapply(padj_cols, function(p) {
            pv <- df[[p]]
            as.integer(is.finite(pv) & pv < a)
        })
        if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
        df[[sprintf("sig_n_padj_meta_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)

        pos_mat <- mapply(
            function(p, zc) {
                pv <- df[[p]]
                zv <- df[[zc]]
                as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv > 0)
            },
            padj_cols, z_cols,
            SIMPLIFY = TRUE
        )
        if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
        df[[sprintf("pos_n_padj_meta_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)

        neg_mat <- mapply(
            function(p, zc) {
                pv <- df[[p]]
                zv <- df[[zc]]
                as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv < 0)
            },
            padj_cols, z_cols,
            SIMPLIFY = TRUE
        )
        if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
        df[[sprintf("neg_n_padj_meta_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
    }
    df
}

# -----------------------------------------------------------------------------
# .write_summary_outputs_meta_csv() - Write summary outputs (CSV + XLSX)
# -----------------------------------------------------------------------------
.write_summary_outputs_meta_csv <- function(df, out_dir, group_name,
                                            stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first.")
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx")
    if (is.null(df) || !nrow(df)) {
        if (exists("log_msg", mode = "function")) try(log_msg("[Output omitted] {group_name} | {stat_tag} No results available"), silent = TRUE)
        return(invisible(NULL))
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))


    ord_keys <- c("sig_n_padj_meta_0_05", "pos_n_padj_meta_0_05", "neg_n_padj_meta_0_05")
    ord_keys <- intersect(ord_keys, names(df))
    if (length(ord_keys)) {
        df <- df |>
            dplyr::arrange(
                dplyr::desc(.data[[ord_keys[1]]]),
                dplyr::desc(ifelse(length(ord_keys) > 1, .data[[ord_keys[2]]], 0)),
                dplyr::desc(ifelse(length(ord_keys) > 2, .data[[ord_keys[3]]], 0)),
                .data[["pathway"]]
            )
    } else {
        df <- df |> dplyr::arrange(.data[["pathway"]])
    }


    data.table::fwrite(df, paste0(base, "_ALL.csv"))

    df005 <- NULL
    df025 <- NULL
    if ("sig_n_padj_meta_0_05" %in% names(df)) {
        df005 <- df |> dplyr::filter(.data[["sig_n_padj_meta_0_05"]] > 0)
        dir.create(dirname(paste0(base, "_padjLT0.05.csv")), recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(df005, paste0(base, "_padjLT0.05.csv"))
    }
    if ("sig_n_padj_meta_0_25" %in% names(df)) {
        df025 <- df |> dplyr::filter(.data[["sig_n_padj_meta_0_25"]] > 0)
        dir.create(dirname(paste0(base, "_padjLT0.25.csv")), recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(df025, paste0(base, "_padjLT0.25.csv"))
    }

    ## ---- XLSX ----
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "ALL")
    openxlsx::writeData(wb, "ALL", df)
    openxlsx::freezePane(wb, "ALL", firstRow = TRUE)
    openxlsx::setColWidths(wb, "ALL", cols = 1:ncol(df), widths = "auto")

    openxlsx::addWorksheet(wb, "padjLT0.05")
    if (!is.null(df005)) {
        openxlsx::writeData(wb, "padjLT0.05", df005)
        openxlsx::freezePane(wb, "padjLT0.05", firstRow = TRUE)
        openxlsx::setColWidths(wb, "padjLT0.05", cols = 1:ncol(df), widths = "auto")
    }

    openxlsx::addWorksheet(wb, "padjLT0.25")
    if (!is.null(df025)) {
        openxlsx::writeData(wb, "padjLT0.25", df025)
        openxlsx::freezePane(wb, "padjLT0.25", firstRow = TRUE)
        openxlsx::setColWidths(wb, "padjLT0.25", cols = 1:ncol(df), widths = "auto")
    }

    openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)

    if (exists("log_msg", mode = "function")) {
        try(log_msg("[Completed Output] {group_name} | {stat_tag} -> {dirname(base)} (.csv + .xlsx)"), silent = TRUE)
    }
    invisible(NULL)
}

# -----------------------------------------------------------------------------
# summarize_meta_fdr_across_subunits() - Main summary function
# -----------------------------------------------------------------------------
summarize_meta_fdr_across_subunits <- function(
  meta_root,
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  genesets_by_group,
  stat_tag = "GSEA_limma_t_cont_meta_fdr"
) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first.")
    for (st in strata) {
        for (ver in versions) {
            base_dir <- file.path(meta_root, st, ver)
            if (!dir.exists(base_dir)) next

            subs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
            subs <- subs[nzchar(subs)]
            if (!length(subs)) next

            for (grp in names(genesets_by_group)) {
                if (exists("log_msg", mode = "function")) try(log_msg("== [meta-summary] stratum=%s | version=%s | group=%s ==", st, ver, grp), silent = TRUE)
                lst <- setNames(vector("list", length(subs)), subs)
                for (su in subs) lst[[su]] <- .read_meta_fdr_table(meta_root, st, ver, su, grp, stat_tag)
                wide <- .merge_subunit_tables_meta(lst)
                wide <- .add_sig_counts_meta(wide, alphas = c(0.05, 0.25))
                out_dir <- file.path(meta_root, "summary", st, ver, safe_fs_name(grp), stat_tag)
                .write_summary_outputs_meta_csv(wide, out_dir, grp, stat_tag)
            }
        }
    }
    invisible(TRUE)
}

# -----------------------------------------------------------------------------
# posthoc_summary_meta_fdr() - Post-hoc summary for GSEA_limma_t_cont
# -----------------------------------------------------------------------------
posthoc_summary_meta_fdr <- function(
  meta_root = NULL,
  genesets_by_group = NULL
) {
    # Get global variables if not provided
    if (is.null(meta_root)) {
        combo_prefix <- tryCatch(get("COMBO_PREFIX", envir = globalenv()), error = function(e) NULL)
        meta_root <- if (is.null(combo_prefix)) "PTMsigDB/meta_fdr" else file.path(combo_prefix, "meta_fdr")
    }
    if (is.null(genesets_by_group)) {
        genesets_by_group <- tryCatch(get("genesets_by_group_ptm", envir = globalenv()), error = function(e) list())
    }

    summarize_meta_fdr_across_subunits(
        meta_root = meta_root,
        strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
        versions = c("BatchAdj"),
        genesets_by_group = genesets_by_group,
        stat_tag = "GSEA_limma_t_cont_meta_fdr"
    )
    invisible(TRUE)
}

# -----------------------------------------------------------------------------
# posthoc_summary_meta_fdr_interaction() - Post-hoc summary for interaction
# -----------------------------------------------------------------------------
posthoc_summary_meta_fdr_interaction <- function(
  meta_root = NULL,
  genesets_by_group = NULL
) {
    # Get global variables if not provided
    if (is.null(meta_root)) {
        combo_prefix <- tryCatch(get("COMBO_PREFIX", envir = globalenv()), error = function(e) NULL)
        meta_root <- if (is.null(combo_prefix)) "PTMsigDB/meta_fdr" else file.path(combo_prefix, "meta_fdr")
    }
    if (is.null(genesets_by_group)) {
        genesets_by_group <- tryCatch(get("genesets_by_group_ptm", envir = globalenv()), error = function(e) list())
    }

    summarize_meta_fdr_across_subunits(
        meta_root = meta_root,
        strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
        versions = c("BatchAdj"),
        genesets_by_group = genesets_by_group,
        stat_tag = "GSEA_limma_interaction_meta_fdr"
    )
    invisible(TRUE)
}
