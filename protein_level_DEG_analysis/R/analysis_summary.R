## ============================================================================
## Summary Tools for DEG Meta-Analysis Results
## Aggregates meta-FDR results across subunits
## ============================================================================

#' Read Meta-FDR Table for DEG
#'
#' Reads a single meta-FDR CSV file for a given subunit
#'
#' @param meta_root Root directory for meta-FDR results
#' @param stratum Stratum name (ALL, TP53_mutant, TP53_wild_type)
#' @param version Version (BatchAdj, RAW)
#' @param subunit Subunit/predictor name
#' @param group_name Group name (always "GENE" for DEG)
#' @param stat_tag Stat tag (DEG_meta_fdr)
#'
#' @return Data frame with pathway, Z_{subunit}, padj_meta_{subunit}
.read_meta_fdr_table_DEG <- function(meta_root, stratum, version, subunit,
                                     group_name = "GENE", stat_tag = "DEG_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")

    grp <- safe_fs_name(group_name)
    fp <- file.path(meta_root, stratum, version, subunit, grp, paste0(stat_tag, ".csv"))

    # If no group subdirectory, also tolerate: .../<subunit>/DEG_meta_fdr.csv
    if (!file.exists(fp)) {
        fp <- file.path(meta_root, stratum, version, subunit, paste0(stat_tag, ".csv"))
    }

    if (!file.exists(fp)) {
        return(NULL)
    }

    dt <- tryCatch(
        data.table::fread(fp, na.strings = c("NA", "NaN", "")),
        error = function(e) NULL
    )

    if (is.null(dt)) {
        return(NULL)
    }

    nms <- tolower(names(dt))
    names(dt) <- nms

    # Check required columns
    if (!("pathway" %in% tolower(names(dt)))) {
        return(NULL)
    }
    if (!("z" %in% tolower(names(dt)))) {
        return(NULL)
    }
    if (!("padj_meta" %in% tolower(names(dt)))) {
        return(NULL)
    }

    # Unify column names
    if (!"pathway" %in% names(dt) && "pathway" %in% nms) {
        names(dt)[match("pathway", nms)] <- "pathway"
    }
    if (!"Z" %in% names(dt) && "z" %in% nms) {
        names(dt)[match("z", nms)] <- "Z"
    }
    if (!"padj_meta" %in% names(dt) && "padj_meta" %in% nms) {
        names(dt)[match("padj_meta", nms)] <- "padj_meta"
    }

    out <- as.data.frame(dt[, c("pathway", "Z", "padj_meta")])
    names(out)[names(out) == "Z"] <- paste0("Z_", subunit)
    names(out)[names(out) == "padj_meta"] <- paste0("padj_meta_", subunit)

    out
}

#' Merge Subunit Tables for Meta DEG
#'
#' Merges multiple subunit tables into wide format
#'
#' @param tbl_list List of data frames from .read_meta_fdr_table_DEG
#'
#' @return Wide-format data frame with all subunits
.merge_subunit_tables_meta_DEG <- function(tbl_list) {
    keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]

    if (!length(keep)) {
        return(NULL)
    }

    out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
    num_cols <- setdiff(names(out), "pathway")
    out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))

    out
}

#' Add Significance Counts for Meta DEG
#'
#' Adds columns counting significant genes at various thresholds
#'
#' @param df Wide-format data frame
#' @param alphas FDR thresholds
#'
#' @return Data frame with added count columns
.add_sig_counts_meta_DEG <- function(df, alphas = c(0.05, 0.25)) {
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

        # Total significant count
        sig_mat <- sapply(padj_cols, function(p) {
            pv <- df[[p]]
            as.integer(is.finite(pv) & pv < a)
        })
        if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
        df[[sprintf("sig_n_padj_meta_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)

        # Positive Z count (significant & Z > 0)
        pos_mat <- mapply(function(p, zc) {
            pv <- df[[p]]
            zv <- df[[zc]]
            as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv > 0)
        }, padj_cols, z_cols, SIMPLIFY = TRUE)
        if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
        df[[sprintf("pos_n_padj_meta_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)

        # Negative Z count (significant & Z < 0)
        neg_mat <- mapply(function(p, zc) {
            pv <- df[[p]]
            zv <- df[[zc]]
            as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv < 0)
        }, padj_cols, z_cols, SIMPLIFY = TRUE)
        if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
        df[[sprintf("neg_n_padj_meta_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
    }

    df
}

#' Write Summary Outputs for Meta DEG
#'
#' Writes CSV and Excel files with summary tables
#'
#' @param df Wide-format summary data frame
#' @param out_dir Output directory
#' @param tag File tag (DEG_meta_fdr)
#'
#' @return Invisible NULL
.write_summary_outputs_meta_csv_DEG <- function(df, out_dir, tag = "DEG_meta_fdr") {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
    if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx first")

    if (is.null(df) || !nrow(df)) {
        return(invisible(NULL))
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(out_dir, paste0("Summary_GENE_", tag))

    # Sort: Significant count (0.05) -> Positive/Negative count -> pathway
    ord_keys <- c("sig_n_padj_meta_0_05", "pos_n_padj_meta_0_05", "neg_n_padj_meta_0_05")
    ord_keys <- intersect(ord_keys, names(df))

    if (length(ord_keys)) {
        df <- df |>
            dplyr::arrange(
                dplyr::desc(.data[[ord_keys[1]]]),
                dplyr::desc(if (length(ord_keys) > 1) .data[[ord_keys[2]]] else 0),
                dplyr::desc(if (length(ord_keys) > 2) .data[[ord_keys[3]]] else 0),
                .data[["pathway"]]
            )
    } else {
        df <- df |> dplyr::arrange(.data[["pathway"]])
    }

    # Write CSV files
    data.table::fwrite(df, paste0(base, "_ALL.csv"))

    if ("sig_n_padj_meta_0_05" %in% names(df)) {
        data.table::fwrite(
            dplyr::filter(df, .data[["sig_n_padj_meta_0_05"]] > 0),
            paste0(base, "_padjLT0.05.csv")
        )
    }

    if ("sig_n_padj_meta_0_25" %in% names(df)) {
        data.table::fwrite(
            dplyr::filter(df, .data[["sig_n_padj_meta_0_25"]] > 0),
            paste0(base, "_padjLT0.25.csv")
        )
    }

    # Write Excel file
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "ALL")
    openxlsx::writeData(wb, "ALL", df)
    openxlsx::freezePane(wb, "ALL", firstRow = TRUE)
    openxlsx::setColWidths(wb, "ALL", cols = 1:ncol(df), widths = "auto")

    if ("sig_n_padj_meta_0_05" %in% names(df)) {
        openxlsx::addWorksheet(wb, "padjLT0.05")
        openxlsx::writeData(wb, "padjLT0.05", dplyr::filter(df, .data[["sig_n_padj_meta_0_05"]] > 0))
        openxlsx::freezePane(wb, "padjLT0.05", firstRow = TRUE)
        openxlsx::setColWidths(wb, "padjLT0.05", cols = 1:ncol(df), widths = "auto")
    }

    if ("sig_n_padj_meta_0_25" %in% names(df)) {
        openxlsx::addWorksheet(wb, "padjLT0.25")
        openxlsx::writeData(wb, "padjLT0.25", dplyr::filter(df, .data[["sig_n_padj_meta_0_25"]] > 0))
        openxlsx::freezePane(wb, "padjLT0.25", firstRow = TRUE)
        openxlsx::setColWidths(wb, "padjLT0.25", cols = 1:ncol(df), widths = "auto")
    }

    openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)

    invisible(NULL)
}

#' Summarize Meta-FDR Across Subunits for DEG
#'
#' Main function to aggregate meta-FDR results across all subunits
#'
#' @param meta_root Root directory for meta-FDR results
#' @param strata Strata to process
#' @param versions Versions to process
#'
#' @return Invisible TRUE
#'
#' @details
#' Creates summary tables in:
#' {meta_root}/summary/{STRATUM}/{VERSION}/GENE/DEG_meta_fdr/
#'
#' @examples
#' summarize_meta_fdr_across_subunits_DEG(
#'     meta_root = "proteomic_DEG/csn_deg_pan_summary_TP53/meta_fdr",
#'     strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
#'     versions = c("BatchAdj")
#' )
summarize_meta_fdr_across_subunits_DEG <- function(
  meta_root = file.path(OUTPUT_PREFIX, "csn_deg_pan_summary_TP53/meta_fdr"),
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj")
) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Please install data.table first")
    }

    for (st in strata) {
        for (ver in versions) {
            base_dir <- file.path(meta_root, st, ver)
            if (!dir.exists(base_dir)) next

            subs <- basename(list.dirs(base_dir, full.names = TRUE, recursive = FALSE))
            subs <- subs[nzchar(subs) & subs != "."]
            if (!length(subs)) next

            # DEG has no geneset grouping, group is fixed to "GENE"
            for (grp in "GENE") {
                log_msg("== [meta-summary-DEG] stratum=%s | version=%s ==", st, ver)

                lst <- setNames(vector("list", length(subs)), subs)
                for (su in subs) {
                    lst[[su]] <- .read_meta_fdr_table_DEG(meta_root, st, ver, su, grp, "DEG_meta_fdr")
                }

                wide <- .merge_subunit_tables_meta_DEG(lst)
                wide <- .add_sig_counts_meta_DEG(wide, alphas = c(0.05, 0.25))

                out_dir <- file.path(meta_root, "summary", st, ver, "GENE", "DEG_meta_fdr")
                .write_summary_outputs_meta_csv_DEG(wide, out_dir, "DEG_meta_fdr")

                log_msg("  [Saved] Summary tables to: %s", out_dir)
            }
        }
    }

    invisible(TRUE)
}
