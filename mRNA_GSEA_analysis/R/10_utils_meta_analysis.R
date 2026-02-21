# =============================================================================
# 10_utils_meta_analysis.R - Meta-Analysis Utilities
# =============================================================================
# This module provides functions for cross-dataset meta-analysis using
# Stouffer's z-score method with FDR correction.
# =============================================================================

# -----------------------------------------------------------------------------
# Normalize Column Names
# -----------------------------------------------------------------------------
#' Normalize GSEA result column names to canonical names
#' @param DT data.table with GSEA results
#' @return data.table with standardized column names
.normalize_cols <- function(DT) {
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

# -----------------------------------------------------------------------------
# Pick Columns
# -----------------------------------------------------------------------------
#' Select specified columns, adding NA for missing
#' @param DT data.table
#' @param need Character vector of required column names
#' @return data.table with specified columns
.pick_cols <- function(DT, need) {
    need <- unique(need)
    if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
    miss <- setdiff(need, names(DT))
    for (m in miss) DT[, (m) := NA]
    DT[, ..need]
}

# -----------------------------------------------------------------------------
# Stouffer Meta-Analysis
# -----------------------------------------------------------------------------
#' Perform Stouffer's z-score meta-analysis on GSEA results
#' @param D data.table with pathway, NES, pval columns
#' @return data.table with meta-analysis results
.stouffer_df <- function(D) {
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

# -----------------------------------------------------------------------------
# Read Single Meta Table (for Summary)
# -----------------------------------------------------------------------------
#' Read meta-FDR table for a single subunit
#' @param meta_root Root directory for meta-analysis output
#' @param stratum Stratum name
#' @param version Analysis version (RAW/CovarAdj)
#' @param subunit Subunit name
#' @param group_name Gene set group name
#' @param stat_tag Statistics tag
#' @return Data frame with pathway, Z, padj_meta columns
.read_meta_table <- function(meta_root, stratum, version, subunit,
                             group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
    grp_safe <- safe_fs_name(group_name)
    f <- file.path(meta_root, stratum, version, subunit, grp_safe, paste0(stat_tag, ".csv"))
    if (!file.exists(f)) {
        return(NULL)
    }
    dt <- tryCatch(data.table::fread(f, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) {
        return(NULL)
    }
    need <- c("pathway", "Z", "padj_meta")
    if (!all(need %in% names(dt))) {
        return(NULL)
    }
    as.data.frame(dt[, need, with = FALSE])
}

# -----------------------------------------------------------------------------
# Meta-FDR Stouffer Analysis (Main Function)
# -----------------------------------------------------------------------------
#' Perform cross-dataset meta-analysis using Stouffer's z-method
#' @param dataset_dirs Named vector/list of dataset directories
#' @param strata Vector of strata to analyze
#' @param stat_tags Vector of statistics tags
#' @param groups Vector of gene set groups
#' @param out_root Output root directory
#' @return Invisible TRUE
meta_fdr_stouffer <- function(dataset_dirs,
                              strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                              stat_tags = c("GSEA_limma_t_cont"),
                              groups = names(genesets_by_group),
                              out_root = file.path(OUTPUT_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")) {
    stopifnot(length(dataset_dirs) > 0)

    # Exclude problematic datasets
    if (!is.null(names(dataset_dirs))) {
        keep_idx <- names(dataset_dirs) != "lusc_cptac_2021"
    } else {
        keep_idx <- basename(dataset_dirs) != "lusc_cptac_2021"
    }
    dataset_dirs <- dataset_dirs[keep_idx]
    if (!length(dataset_dirs)) {
        message("[meta] All datasets excluded or no available datasets; skipping.")
        return(invisible(TRUE))
    }

    if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table first")
    data.table::setDTthreads(1L)

    passes <- c("RAW", "CovarAdj")
    base_out <- out_root
    dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

    for (stratum in strata) {
        for (pass_label in passes) {
            for (stat_tag in stat_tags) {
                for (grp in groups) {
                    # Find all subunits
                    subunits <- unique(unlist(lapply(dataset_dirs, function(dsdir) {
                        base <- file.path(
                            OUTPUT_PREFIX,
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
                            "[meta] %s/%s/%s/%s no available subunit, skipping",
                            stratum, pass_label, stat_tag, grp
                        ))
                        next
                    }

                    # For each subunit, collect CSV from each dataset
                    for (su in subunits) {
                        parts <- list()
                        for (dsdir in dataset_dirs) {
                            f <- file.path(
                                OUTPUT_PREFIX,
                                safe_fs_name(grp), pass_label, basename(dsdir), stratum, su,
                                sprintf("%s.csv", stat_tag)
                            )
                            if (!file.exists(f)) next
                            dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
                            if (is.null(dt) || nrow(dt) == 0) next
                            dt <- .normalize_cols(dt)
                            dt <- .pick_cols(dt, need = c("pathway", "NES", "pval", "padj", "size"))
                            dt[, dataset := basename(dsdir)]
                            parts[[length(parts) + 1L]] <- dt
                        }

                        if (length(parts) == 0) {
                            next
                        }

                        D <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
                        D <- D[!is.na(pathway)]
                        res <- .stouffer_df(D)
                        if (is.null(res)) next

                        # Write file
                        grp_safe <- safe_fs_name(grp)
                        out_dir <- file.path(base_out, stratum, pass_label, su, grp_safe)
                        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                        out_file <- file.path(out_dir, paste0(stat_tag, "_meta_fdr.csv"))
                        data.table::fwrite(res, out_file)
                        message(sprintf(
                            "[meta] %s/%s/%s/%s -> %s (n=%d pathways)",
                            stratum, pass_label, su, stat_tag, out_file, nrow(res)
                        ))
                    }
                }
            }
        }
    }

    invisible(TRUE)
}
