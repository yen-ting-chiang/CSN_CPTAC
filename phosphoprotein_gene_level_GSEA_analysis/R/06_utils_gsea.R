# =============================================================================
# 06_utils_gsea.R - GSEA Helper Functions
# =============================================================================
# This file contains functions for running and processing GSEA results.
# =============================================================================

# -----------------------------------------------------------------------------
# Ensure Stats Have Correct Gene Names
# -----------------------------------------------------------------------------
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
    # stats: numeric vector (limma t etc.)
    # gene_names: rownames of corresponding matrix (gene IDs)
    if (is.null(stats)) {
        return(NULL)
    }
    v <- suppressWarnings(as.numeric(stats))
    if (length(v) != length(gene_names)) {
        stop(sprintf(
            "[ensure-names%s] stats(%d) and gene_names(%d) length mismatch",
            if (!is.null(label)) paste0("-", label) else "",
            length(v), length(gene_names)
        ))
    }
    names(v) <- as.character(gene_names)
    # Do not filter yet, leave to next step for cleanup (will log)
    v
}

# -----------------------------------------------------------------------------
# Clean and Sort GSEA Stats (Keep Finite Values, Deduplicate, Sort)
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
    # Names must exist
    if (is.null(nm)) {
        log_msg("[gsea-%s] stats has no names after filtering -> skip", label)
        return(NULL)
    }
    # Remove duplicate names (keep first occurrence)
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
        log_msg("[gsea-%s] too few finite stats after filtering: %d < %d -> skip", label, length(stats), min_n)
        return(NULL)
    }
    # Sort by score from high to low (preranked)
    stats[order(stats, decreasing = TRUE)]
}

# -----------------------------------------------------------------------------
# Map Pathways to Stats Universe and Filter by Size
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
# Run GSEA from Sorted Statistics (fgseaMultilevel with fallback)
# -----------------------------------------------------------------------------
.gsea_from_ranks <- function(pathways, stats, minSize, maxSize, gsea_eps = 0, label = NULL) {
    if (is.null(pathways) || !length(pathways) || is.null(stats) || !length(stats)) {
        return(NULL)
    }
    orig_n <- length(pathways)
    pw_use <- ._intersect_and_filter_pathways(pathways, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf(
        "[gsea:%s] orig=%d -> used=%d",
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
# GSEA from Differential Correlation Z-scores
# -----------------------------------------------------------------------------
gsea_from_diffcorr <- function(pathways, z_stats, label = "diff-corr", out_prefix = NULL) {
    stopifnot(is.numeric(z_stats), !is.null(names(z_stats)))
    minSize <- opt("minSize", 15L)
    maxSize <- opt("maxSize", 500L)
    orig_n <- length(pathways)
    pw_use <- ._intersect_and_filter_pathways(pathways, names(z_stats), minSize = minSize, maxSize = maxSize)
    if (!is.null(out_prefix)) {
        message(sprintf("[gsea:%s] %s orig=%d -> used=%d", label, out_prefix, orig_n, length(pw_use)))
    } else {
        message(sprintf("[gsea:%s] orig=%d -> used=%d", label, orig_n, length(pw_use)))
    }
    if (!length(pw_use)) {
        return(invisible(NULL))
    }

    set.seed(1L)
    res <- tryCatch(
        {
            suppressWarnings(
                fgsea::fgseaMultilevel(
                    pathways = pw_use, stats = z_stats,
                    minSize = minSize, maxSize = maxSize,
                    eps = 1e-10, nproc = getOption(".fgsea_nproc", 1L)
                )
            )
        },
        error = function(e) {
            message(sprintf("[gsea_from_diffcorr] Multilevel failed: %s -> fallback to fgseaSimple", conditionMessage(e)))
            suppressWarnings(
                fgsea::fgseaSimple(
                    pathways = pw_use, stats = z_stats,
                    nperm = 10000, minSize = minSize, maxSize = maxSize
                )
            )
        }
    )
    invisible(as.data.frame(res))
}

# -----------------------------------------------------------------------------
# Run FGSEA and Save Results
# -----------------------------------------------------------------------------
run_fgsea_save <- function(stats, genesets, out_prefix, top_plot_n = 0, plot_title = "") {
    set.seed(1234)
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

    minSize <- opt("minSize", 15L)
    maxSize <- opt("maxSize", 500L)

    orig_n <- length(genesets)
    pw_use <- ._intersect_and_filter_pathways(genesets, names(stats), minSize = minSize, maxSize = maxSize)
    if (!length(pw_use)) {
        log_msg("[run_fgsea_save] no pathways to run")
        return(invisible(NULL))
    }

    bp <- BiocParallel::SerialParam()
    res <- tryCatch(
        {
            suppressWarnings(
                fgsea::fgseaMultilevel(
                    pathways = pw_use,
                    stats    = stats,
                    minSize  = minSize,
                    maxSize  = maxSize,
                    eps      = 1e-10,
                    nproc    = 1L,
                    BPPARAM  = bp
                )
            )
        },
        error = function(e) {
            message(sprintf("[run_fgsea_save] fgseaMultilevel failed: %s; fallback", conditionMessage(e)))
            suppressWarnings(
                fgsea::fgseaSimple(
                    pathways = pw_use,
                    stats    = stats,
                    minSize  = minSize,
                    maxSize  = maxSize,
                    nperm    = 1000,
                    BPPARAM  = bp
                )
            )
        }
    )
    res <- as.data.frame(res)
    data.table::fwrite(res, paste0(out_prefix, ".csv"))
    if (top_plot_n > 0) {
        
    }
    invisible(res)
}

# -----------------------------------------------------------------------------
# Read GSEA Result Table
# -----------------------------------------------------------------------------
read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
    grp <- safe_fs_name(group_name)
    su <- subunit
    f <- paste0(stat_tag, ".csv")

    # Based on actual file structure (RAW/BatchAdj) and old assumptions, list candidate paths
    candidates <- c(
        file.path(out_root, su, f), # collection-first: ./<subunit>/STAT.csv
        file.path(out_root, grp, su, "H", f), # old: ./<group>/<subunit>/H/STAT.csv
        file.path(out_root, grp, su, f), # old: no H subfolder
        file.path(out_root, su, grp, "H", f), # old: ./<subunit>/<group>/H/STAT.csv
        file.path(out_root, su, grp, f) # old: no H subfolder
    )

    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) {
        if (exists("log_msg", mode = "function")) try(log_msg("    File not found: %s", paste(candidates, collapse = " | ")), silent = TRUE)
        return(NULL)
    }
    fp <- hit[1L]

    dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
    if (is.null(dt)) {
        return(NULL)
    }

    need_cols <- c("pathway", "NES", "padj")
    if (!all(need_cols %in% names(dt))) {
        # Flexible column name mapping
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
# Merge Subunit Tables
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
# Add Significance Counts to Summary Table
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
# Write Summary Outputs (Excel and CSV)
# -----------------------------------------------------------------------------
write_summary_outputs <- function(df, out_dir, group_name, stat_tag) {
    if (is.null(df) || !nrow(df)) {
        log_msg("  [Skip output] {group_name} | {stat_tag} no available results")
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
