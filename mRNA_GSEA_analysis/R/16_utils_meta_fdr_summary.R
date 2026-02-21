# =============================================================================
# 16_utils_meta_fdr_summary.R - Meta-FDR Summary Utilities
# =============================================================================
# This module summarizes meta-FDR (Stouffer to BH) results across subunits.
#   - Input:  csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/<RAW|CovarAdj>/<SUBUNIT>/<GROUP>/GSEA_limma_t_cont_meta_fdr.csv
#   - Output: .../summary/<STRATUM>/<RAW|CovarAdj>/<GROUP>/GSEA_limma_t_cont_meta_fdr/
#              - Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_ALL.csv
#              - Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_padjLT0.05.csv
#              - Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_padjLT0.25.csv
#   - Columns merged per subunit: Z_<SUBUNIT>, padj_meta_<SUBUNIT>
#   - Counting uses padj_meta thresholds; direction by sign(Z).
# =============================================================================

.read_meta_fdr_table <- function(meta_root, stratum, version, subunit,
                                 group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
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
    # Lenient matching
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

.write_summary_outputs_meta_csv <- function(df, out_dir, group_name,
                                            stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx first")
    if (is.null(df) || !nrow(df)) {
        if (exists("log_msg", mode = "function")) try(log_msg("  [Skip output] {group_name} | {stat_tag} no available results"), silent = TRUE)
        return(invisible(NULL))
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))

    ## ---- Sorting: prioritize significant count (0.05), then positive/negative significant count ----
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

    ## ---- CSV output (ALL / padj<0.05 / padj<0.25) ----
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

    ## ---- XLSX output (ALL / padjLT0.05 / padjLT0.25 three worksheets) ----
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
        try(log_msg("  [Output completed] {group_name} | {stat_tag} -> {dirname(base)} (.csv + .xlsx)"), silent = TRUE)
    }
    invisible(NULL)
}


summarize_meta_fdr_across_subunits <- function(
  meta_root = file.path(OUTPUT_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("RAW", "CovarAdj"),
  genesets_by_group = genesets_by_group,
  stat_tag = "GSEA_limma_t_cont_meta_fdr"
) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
    for (st in strata) {
        for (ver in versions) {
            base_dir <- file.path(meta_root, st, ver)
            if (!dir.exists(base_dir)) next
            # Auto-detect subunit list: get the first-level subdirectory names at this level
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

## ---- One-click execute for limma_t_cont (do not re-run GSEA; only read existing *_meta_fdr.csv and summarize) ----
posthoc_summary_meta_fdr <- function() {
    summarize_meta_fdr_across_subunits(
        meta_root = file.path(OUTPUT_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
        strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
        versions = c("RAW", "CovarAdj"),
        genesets_by_group = genesets_by_group,
        stat_tag = "GSEA_limma_t_cont_meta_fdr"
    )
    invisible(TRUE)
}

## ---- One-click execute for limma interaction ----
posthoc_summary_meta_fdr_interaction <- function() {
    summarize_meta_fdr_across_subunits(
        meta_root = file.path(OUTPUT_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
        strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
        versions = c("RAW", "CovarAdj"),
        genesets_by_group = genesets_by_group,
        stat_tag = "GSEA_limma_interaction_meta_fdr"
    )
    invisible(TRUE)
}
