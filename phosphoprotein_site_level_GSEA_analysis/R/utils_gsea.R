# =============================================================================
# GSEA Analysis Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for running fgsea and processing GSEA results.
# =============================================================================

# -----------------------------------------------------------------------------
# ._ensure_stats_names()
# -----------------------------------------------------------------------------
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
    if (is.null(stats)) {
        return(NULL)
    }
    v <- suppressWarnings(as.numeric(stats))
    if (length(v) != length(gene_names)) {
        stop(sprintf(
            "The lengths of `[ensure-names%s] stats(%d)` and `gene_names(%d)` do not match.",
            if (!is.null(label)) paste0("-", label) else "",
            length(v), length(gene_names)
        ))
    }
    names(v) <- as.character(gene_names)

    v
}

# -----------------------------------------------------------------------------
# .gsea_from_ranks()
# -----------------------------------------------------------------------------
.gsea_from_ranks <- function(pathways, stats, minSize, maxSize, gsea_eps = 0, label = NULL) {
    if (is.null(pathways) || !length(pathways) || is.null(stats) || !length(stats)) {
        return(NULL)
    }
    orig_n <- length(pathways)
    pw_use <- ._intersect_and_filter_pathways(pathways, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf(
        "[gsea:%s] original=%d → use=%d",
        ifelse(is.null(label), "NA", label),
        orig_n, length(pw_use)
    ))
    if (!length(pw_use)) {
        return(NULL)
    }

    nproc <- getOption(".fgsea_nproc", 1L)
    set.seed(1L)
    res <- tryCatch(
        {
            suppressWarnings(fgsea::fgseaMultilevel(
                pathways = pw_use, stats = stats,
                minSize = minSize, maxSize = maxSize,
                eps = gsea_eps, nproc = nproc
            ))
        },
        error = function(e) {
            message(sprintf("[gsea_from_ranks] fgseaMultilevel failed: %s; fallback to fgseaSimple", conditionMessage(e)))
            suppressWarnings(fgsea::fgseaSimple(
                pathways = pw_use, stats = stats,
                nperm = 10000, minSize = minSize, maxSize = maxSize
            ))
        }
    )
    if (!is.null(label)) attr(res, "label") <- label
    res
}

# -----------------------------------------------------------------------------
# ._finite_rank_stats()
# -----------------------------------------------------------------------------
._finite_rank_stats <- function(stats, gene_names = NULL, label = "", min_n = 20L) {
    nm <- names(stats)
    stats <- suppressWarnings(as.numeric(stats))
    if (is.null(nm) && !is.null(gene_names) && length(stats) == length(gene_names)) {
        nm <- gene_names
    }
    names(stats) <- nm

    ok <- is.finite(stats) & !is.na(stats)
    if (sum(!ok) > 0L) {
        if (!is.null(label) && nzchar(label)) {
            log_msg("[gsea-%s] drop non-finite stats: %d", label, sum(!ok))
        }
    }
    stats <- stats[ok]
    nm <- nm[ok]
    names(stats) <- nm

    if (is.null(nm)) {
        log_msg("[gsea-%s] stats has no names after filtering → skip", label)
        return(NULL)
    }

    dup <- duplicated(nm)
    if (any(dup)) {
        if (!is.null(label) && nzchar(label)) {
            log_msg("[gsea-%s] drop duplicated gene names in stats: %d", label, sum(dup))
        }
        stats <- stats[!dup]
        nm <- nm[!dup]
    }
    names(stats) <- nm
    if (length(stats) < min_n) {
        log_msg("[gsea-%s] too few finite stats after filtering: %d < %d → skip", label, length(stats), min_n)
        return(NULL)
    }

    stats[order(stats, decreasing = TRUE)]
}

# -----------------------------------------------------------------------------
# ._intersect_and_filter_pathways()
# -----------------------------------------------------------------------------
._intersect_and_filter_pathways <- function(pathways, universe, minSize, maxSize, label = "") {
    pw <- lapply(pathways, function(gs) unique(intersect(gs, universe)))
    lens <- vapply(pw, length, integer(1))
    keep <- which(lens >= minSize & lens <= maxSize)
    if (length(keep) == 0L) {
        log_msg("[gsea-%s] no pathways within size bounds after intersect (min=%d, max=%d)", label, minSize, maxSize)
        return(NULL)
    }
    pw[keep]
}


# -----------------------------------------------------------------------------
# read_gsea_table()
# -----------------------------------------------------------------------------
read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
    grp <- safe_fs_name(group_name)
    su <- subunit
    f <- paste0(stat_tag, ".csv")


    candidates <- c(
        file.path(out_root, su, f),
        file.path(out_root, grp, su, "H", f),
        file.path(out_root, grp, su, f),
        file.path(out_root, su, grp, "H", f),
        file.path(out_root, su, grp, f)
    )

    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) {
        if (exists("log_msg", mode = "function")) try(log_msg("file not found: %s", paste(candidates, collapse = " | ")), silent = TRUE)
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
        if (!"pathway" %in% names(dt) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
        if (!"NES" %in% names(dt) && "nes" %in% nms) names(dt)[match("nes", nms)] <- "NES"
        if (!"padj" %in% names(dt)) dt$padj <- NA_real_
    }
    dt <- as.data.frame(dt[, c("pathway", "NES", "padj")])
    names(dt)[names(dt) == "NES"] <- paste0("NES_", subunit)
    names(dt)[names(dt) == "padj"] <- paste0("padj_", subunit)
    dt
}

# -----------------------------------------------------------------------------
# merge_subunit_tables()
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# add_sig_counts()
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# write_summary_outputs()
# -----------------------------------------------------------------------------
write_summary_outputs <- function(df, out_dir, group_name, stat_tag) {
    if (is.null(df) || !nrow(df)) {
        log_msg("[Output omitted] {group_name} | {stat_tag} No results available")
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
    log_msg("  [Output complete] {group_name} | {stat_tag} -> {dirname(base)}")
}

# -----------------------------------------------------------------------------
# summarize_all_groups()
# -----------------------------------------------------------------------------
summarize_all_groups <- function(out_root, csn_subunits, genesets_by_group,
                                 stat_tags = c("GSEA_limma_t")) {
    sum_root <- file.path(out_root, "summary")
    dir.create(sum_root, recursive = TRUE, showWarnings = FALSE)
    for (grp_name in names(genesets_by_group)) {
        grp_safe <- safe_fs_name(grp_name)
        for (stat_tag in stat_tags) {
            log_msg("==Summary:group={grp_name} | stat={stat_tag} ==")
            tbl_list <- setNames(vector("list", length(csn_subunits)), csn_subunits)
            for (su in csn_subunits) tbl_list[[su]] <- read_gsea_table(out_root, su, grp_name, stat_tag)
            wide <- merge_subunit_tables(tbl_list)
            wide <- add_sig_counts(wide, alphas = c(0.05, 0.25))
            out_dir <- file.path(sum_root, stat_tag)
            write_summary_outputs(wide, out_dir, grp_name, stat_tag)
        }
    }
}
