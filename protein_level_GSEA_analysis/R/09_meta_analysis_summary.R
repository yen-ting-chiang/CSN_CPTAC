## ============================================================================
## 09_meta_analysis_summary.R
##
## Summarize meta-FDR (Stouffer → BH) across subunits
## - Input: csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/BatchAdj/<SUBUNIT>/<GROUP>/GSEA_limma_t_cont_meta_fdr.csv
## ============================================================================

#' Read meta-FDR table for a specific subunit and group
#'
#' @param meta_root Root directory for meta-analysis results
#' @param stratum Stratum name (e.g., "ALL", "TP53_mutant", "TP53_wild_type")
#' @param version Version name (e.g., "BatchAdj")
#' @param subunit Subunit name
#' @param group_name Gene set group name
#' @param stat_tag Statistics tag (e.g., "GSEA_limma_t_cont_meta_fdr")
#' @return Data frame with pathway, Z, and padj_meta columns, or NULL if not found
.read_meta_fdr_table <- function(meta_root, stratum, version, subunit,
                                 group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("please install data.table")
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

#' Merge subunit tables into a wide format
#'
#' @param tbl_list List of data frames from .read_meta_fdr_table
#' @return Merged data frame with all subunits, or NULL if empty
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

#' Add significance counts across subunits
#'
#' @param df Data frame from .merge_subunit_tables_meta
#' @param alphas Significance thresholds
#' @return Data frame with added sig_n, pos_n, neg_n columns
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

#' Write summary outputs in CSV and XLSX formats
#'
#' @param df Summary data frame
#' @param out_dir Output directory
#' @param group_name Gene set group name
#' @param stat_tag Statistics tag
#' @return Invisible NULL
.write_summary_outputs_meta_csv <- function(df, out_dir, group_name,
                                            stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("please install data.table")
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("please install openxlsx")
    if (is.null(df) || !nrow(df)) {
        log_msg("  [skip output] %s | %s No results available", group_name, stat_tag)
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

    ## ---- XLSX output----
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

    log_msg("  [output completed] %s | %s -> %s (.csv + .xlsx)", group_name, stat_tag, dirname(base))
    invisible(NULL)
}


#' Summarize meta-FDR across all subunits
#'
#' Aggregates meta-analysis results from all subunits for each stratum/version/group
#' combination and generates summary CSV and XLSX files.
#'
#' @param meta_root Root directory for meta-analysis results
#' @param strata Character vector of strata to process
#' @param versions Character vector of versions to process
#' @param genesets_by_group Named list of gene sets by group (used for group names)
#' @param stat_tag Statistics tag (e.g., "GSEA_limma_t_cont_meta_fdr")
#' @return Invisible TRUE
summarize_meta_fdr_across_subunits <- function(
  meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  genesets_by_group = genesets_by_group,
  stat_tag = "GSEA_limma_t_cont_meta_fdr"
) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("please install data.table")
    for (st in strata) {
        for (ver in versions) {
            base_dir <- file.path(meta_root, st, ver)
            if (!dir.exists(base_dir)) next

            subs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
            subs <- subs[nzchar(subs)]
            if (!length(subs)) next

            for (grp in names(genesets_by_group)) {
                log_msg("== [meta-summary] stratum=%s | version=%s | group=%s ==", st, ver, grp)
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

#' Run meta-FDR summary for limma continuous analysis
#'
#' Wrapper function to run summary for GSEA_limma_t_cont_meta_fdr results
#'
#' @param meta_root Root directory for meta-analysis results
#' @param strata Character vector of strata to process
#' @param genesets_by_group Named list of gene sets by group
#' @return Invisible TRUE
posthoc_summary_meta_fdr <- function(
  meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  genesets_by_group = genesets_by_group
) {
    summarize_meta_fdr_across_subunits(
        meta_root = meta_root,
        strata = strata,
        versions = c("BatchAdj"),
        genesets_by_group = genesets_by_group,
        stat_tag = "GSEA_limma_t_cont_meta_fdr"
    )
    invisible(TRUE)
}

#' Run meta-FDR summary for limma interaction analysis
#'
#' Wrapper function to run summary for GSEA_limma_interaction_meta_fdr results
#'
#' @param meta_root Root directory for meta-analysis results
#' @param strata Character vector of strata to process
#' @param genesets_by_group Named list of gene sets by group
#' @return Invisible TRUE
posthoc_summary_meta_fdr_interaction <- function(
  meta_root = if (is.null(GSEA_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(GSEA_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  genesets_by_group = genesets_by_group
) {
    summarize_meta_fdr_across_subunits(
        meta_root = meta_root,
        strata = strata,
        versions = c("BatchAdj"),
        genesets_by_group = genesets_by_group,
        stat_tag = "GSEA_limma_interaction_meta_fdr"
    )
    invisible(TRUE)
}
