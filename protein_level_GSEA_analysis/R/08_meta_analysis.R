## ====================================================================
## 08_meta_analysis.R
##
## Purpose: Meta-analysis across datasets using Stouffer's method
## Contains: Pan-cancer meta-FDR calculation, result aggregation
## ====================================================================

#' Perform meta-analysis via Stouffer's Z-score method
#'
#' @param dataset_dirs Named vector of dataset directories
#' @param strata Character vector of strata (e.g., "ALL", "TP53_mutant")
#' @param stat_tags Character vector of statistic types
#' @param groups Character vector of gene set groups
#' @param out_root Output root directory
#' @return Invisible TRUE
#' @export
meta_fdr_stouffer <- function(dataset_dirs,
                              strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                              stat_tags = c("GSEA_limma_t_cont"),
                              groups = get0("genesets_by_group", ifnotfound = list()) |> names(),
                              out_root = if (is.null(get0("GSEA_PREFIX"))) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(get0("GSEA_PREFIX"), "csn_gsea_pan_summary_TP53/meta_fdr")) {
    stopifnot(length(dataset_dirs) > 0)

    ## [EXCLUDE-FROM-META] exclude lusc_cptac_2021
    if (!is.null(names(dataset_dirs))) {
        keep_idx <- names(dataset_dirs) != "lusc_cptac_2021"
    } else {
        keep_idx <- basename(dataset_dirs) != "lusc_cptac_2021"
    }
    dataset_dirs <- dataset_dirs[keep_idx]

    if (!length(dataset_dirs)) {
        message("[meta] All datasets excluded (lusc_cptac_2021 removed); skipped.")
        return(invisible(TRUE))
    }

    if (!requireNamespace("data.table", quietly = TRUE)) stop("install data.table first.")
    data.table::setDTthreads(1L)

    # ---- helper functions ------------------------------------------------------

    # Normalize column names to standard names
    normalize_cols <- function(DT) {
        nm <- names(DT)
        low <- tolower(nm)
        map <- c(
            "pathway" = "pathway", "term" = "pathway", "pathway_name" = "pathway",
            "nes" = "NES", "enrichment_score" = "NES",
            "pval" = "pval", "p.value" = "pval", "pvalue" = "pval", "p" = "pval",
            "padj" = "padj", "fdr" = "padj", "qval" = "padj",
            "size" = "size", "setsize" = "size", "n" = "size"
        )

        for (i in seq_along(nm)) {
            key <- low[i]
            if (key %in% names(map)) data.table::setnames(DT, nm[i], map[[key]])
        }
        DT
    }

    # Select required columns (add NA if missing)
    pick_cols <- function(DT, need) {
        need <- unique(need)
        if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)

        miss <- setdiff(need, names(DT))
        for (m in miss) DT[, (m) := NA]

        DT[, ..need]
    }

    # Stouffer's method: combine p-values using Z-scores
    stouffer_df <- function(D) {
        D <- D[is.finite(pval) & is.finite(NES)]
        if (nrow(D) == 0) {
            return(NULL)
        }
        D[, pval := pmin(pmax(pval, 1e-100), 1 - 1e-100)]
        D[, z := sign(NES) * stats::qnorm(p = pval / 2, lower.tail = FALSE)]

        # Unweighted Stouffer
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
                    # 1) Find all subunits (by unioning folders from each dataset)
                    subunits <- unique(unlist(lapply(dataset_dirs, function(dsdir) {
                        base <- file.path(
                            get0("GSEA_PREFIX") %||% "csn_gsea_results_TP53_by_collection",
                            safe_fs_name(grp), pass_label, basename(dsdir), stratum
                        )
                        if (!dir.exists(base)) {
                            return(character(0))
                        }
                        # List top-level subdirectories (each subunit)
                        subs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
                        subs[subs != ""]
                    })))

                    if (length(subunits) == 0) {
                        log_msg("[meta] %s/%s/%s/%s: no available subunits, skip", stratum, pass_label, stat_tag, grp)
                        next
                    }

                    # 2) For each subunit, collect CSV files from datasets and perform Stouffer
                    for (su in subunits) {
                        parts <- list()
                        for (dsdir in dataset_dirs) {
                            f <- file.path(
                                get0("GSEA_PREFIX") %||% "csn_gsea_results_TP53_by_collection",
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

                        # 3) Write result file
                        grp_safe <- safe_fs_name(grp)
                        out_dir <- file.path(base_out, stratum, pass_label, su, grp_safe)
                        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                        out_file <- file.path(out_dir, paste0(stat_tag, "_meta_fdr.csv"))
                        data.table::fwrite(res, out_file)
                        log_msg("[meta] %s/%s/%s/%s -> %s (n=%d pathways, k=%d datasets)", stratum, pass_label, su, stat_tag, out_file, nrow(res), nrow(D[, .N, by = dataset]))
                    } # su
                } # grp
            } # stat_tag
        } # pass
    } # stratum

    invisible(TRUE)
}


#' Read GSEA results table from output structure
#'
#' @param out_root Output root directory
#' @param subunit CSN subunit name
#' @param group_name Gene set group name
#' @param stat_tag Statistic tag (e.g., "GSEA_limma_t")
#' @return Data frame with pathway, NES, padj columns
#' @export
read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
    grp <- safe_fs_name(group_name)
    su <- subunit
    f <- paste0(stat_tag, ".csv")

    candidates <- c(
        file.path(out_root, su, f), # collection-first: ./<subunit>/STAT.csv
        file.path(out_root, grp, su, "H", f), # ./<group>/<subunit>/H/STAT.csv
        file.path(out_root, grp, su, f), # no H folder
        file.path(out_root, su, grp, "H", f), # ./<subunit>/<group>/H/STAT.csv
        file.path(out_root, su, grp, f) # no H folder
    )

    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) {
        log_msg(" file not found: %s", paste(candidates, collapse = " | "))
        return(NULL)
    }
    fp <- hit[1L]

    dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
    if (is.null(dt)) {
        return(NULL)
    }

    need_cols <- c("pathway", "NES", "padj")
    if (!all(need_cols %in% names(dt))) {
        nms <- tolower(names(dt))
        if (!("pathway" %in% names(dt)) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
        if (!("NES" %in% names(dt)) && "nes" %in% nms) names(dt)[match("nes", nms)] <- "NES"
        if (!("padj" %in% names(dt))) dt$padj <- NA_real_
    }
    dt <- as.data.frame(dt[, c("pathway", "NES", "padj")])
    names(dt)[names(dt) == "NES"] <- paste0("NES_", subunit)
    names(dt)[names(dt) == "padj"] <- paste0("padj_", subunit)
    dt
}


#' Merge multiple subunit GSEA tables
#'
#' @param tbl_list List of data frames from read_gsea_table
#' @return Merged data frame with all subunits
#' @export
merge_subunit_tables <- function(tbl_list) {
    keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
    if (!length(keep)) {
        return(NULL)
    }
    out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
    num_cols <- setdiff(names(out), "pathway")
    out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
    out
}


#' Add significance count columns to summary table
#'
#' @param df Data frame with pathway, NES_*, padj_* columns
#' @param alphas FDR thresholds (default: c(0.05, 0.25))
#' @return Data frame with additional sig_n, pos_n, neg_n columns
#' @export
add_sig_counts <- function(df, alphas = c(0.05, 0.25)) {
    if (is.null(df) || !nrow(df)) {
        return(df)
    }
    padj_cols <- grep("^padj_", names(df), value = TRUE)
    subunits <- sub("^padj_", "", padj_cols)
    nes_cols <- paste0("NES_", subunits)
    keep_idx <- nes_cols %in% names(df)
    padj_cols <- padj_cols[keep_idx]
    nes_cols <- nes_cols[keep_idx]
    if (!length(padj_cols)) {
        return(df)
    }
    for (a in alphas) {
        a_tag <- gsub("\\.", "_", as.character(a))
        sig_mat <- mapply(function(p, n) {
            pv <- df[[p]]
            as.integer(!is.na(pv) & pv < a)
        }, padj_cols, nes_cols, SIMPLIFY = TRUE)
        if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
        df[[sprintf("sig_n_padj_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)

        pos_mat <- mapply(function(p, n) {
            pv <- df[[p]]
            nv <- df[[n]]
            as.integer(!is.na(pv) & pv < a & !is.na(nv) & nv > 0)
        }, padj_cols, nes_cols, SIMPLIFY = TRUE)
        if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
        df[[sprintf("pos_n_padj_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)

        neg_mat <- mapply(function(p, n) {
            pv <- df[[p]]
            nv <- df[[n]]
            as.integer(!is.na(pv) & pv < a & !is.na(nv) & nv < 0)
        }, padj_cols, nes_cols, SIMPLIFY = TRUE)
        if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
        df[[sprintf("neg_n_padj_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
    }
    df
}


#' Write summary output files (CSV + Excel)
#'
#' @param df Summary data frame
#' @param out_dir Output directory
#' @param group_name Gene set group name
#' @param stat_tag Statistic tag
#' @return Invisible NULL
#' @export
write_summary_outputs <- function(df, out_dir, group_name, stat_tag) {
    if (is.null(df) || !nrow(df)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("  [Skip output] {group_name} | {stat_tag}: No results available")
        }
        return(invisible(NULL))
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))

    ord_keys <- c("sig_n_padj_0_05", "pos_n_padj_0_05", "neg_n_padj_0_05")
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
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "ALL")
    openxlsx::writeData(wb, "ALL", df)

    if ("sig_n_padj_0_05" %in% names(df)) {
        df005 <- df |> dplyr::filter(.data[["sig_n_padj_0_05"]] > 0)
        openxlsx::addWorksheet(wb, "padj_lt_0.05")
        openxlsx::writeData(wb, "padj_lt_0.05", df005)
        data.table::fwrite(df005, paste0(base, "_padjLT0.05.csv"))
    }
    if ("sig_n_padj_0_25" %in% names(df)) {
        df025 <- df |> dplyr::filter(.data[["sig_n_padj_0_25"]] > 0)
        openxlsx::addWorksheet(wb, "padj_lt_0.25")
        openxlsx::writeData(wb, "padj_lt_0.25", df025)
        data.table::fwrite(df025, paste0(base, "_padjLT0.25.csv"))
    }
    openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)

    if (exists("log_msg", mode = "function")) {
        log_msg("  [Complete output] {group_name} | {stat_tag} -> {dirname(base)}")
    }
    invisible(NULL)
}


#' Summarize all gene set groups across subunits
#'
#' @param out_root Output root directory
#' @param csn_subunits CSN subunit names
#' @param genesets_by_group List of gene sets by group
#' @param stat_tags Statistic tags to summarize
#' @return Invisible NULL
#' @export
summarize_all_groups <- function(out_root, csn_subunits, genesets_by_group,
                                 stat_tags = c("GSEA_limma_t")) {
    sum_root <- file.path(out_root, "summary")
    dir.create(sum_root, recursive = TRUE, showWarnings = FALSE)
    for (grp_name in names(genesets_by_group)) {
        grp_safe <- safe_fs_name(grp_name)
        for (stat_tag in stat_tags) {
            if (exists("log_msg", mode = "function")) {
                log_msg("== Summary: group={grp_name} | stat={stat_tag} ==")
            }
            tbl_list <- setNames(vector("list", length(csn_subunits)), csn_subunits)
            for (su in csn_subunits) tbl_list[[su]] <- read_gsea_table(out_root, su, grp_name, stat_tag)
            wide <- merge_subunit_tables(tbl_list)
            wide <- add_sig_counts(wide, alphas = c(0.05, 0.25))
            out_dir <- file.path(sum_root, stat_tag)
            write_summary_outputs(wide, out_dir, grp_name, stat_tag)
        }
    }
    invisible(NULL)
}


# ---- End of 08_meta_analysis.R ----
