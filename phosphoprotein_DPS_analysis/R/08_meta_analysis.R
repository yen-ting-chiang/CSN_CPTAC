## =========================================================
## 08_meta_analysis.R
## Stouffer Meta-Analysis Functions for DPS Analysis
## =========================================================

## ===== Safe DPS result file reading =====
.safe_read_dps <- function(fp) {
    if (!file.exists(fp)) {
        return(NULL)
    }
    dt <- tryCatch(
        data.table::fread(fp, na.strings = c("NA", "NaN", "")),
        error = function(e) NULL
    )
    if (is.null(dt) || !nrow(dt)) {
        return(NULL)
    }

    
    if (!"site_id" %in% names(dt) && "ID" %in% names(dt)) {
        dt$site_id <- dt$ID
    }
    dt
}

## ===== Combine results for one predictor =====
.combine_one_predictor <- function(tbl_list) {
    keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
    if (!length(keep)) {
        return(NULL)
    }

    long <- data.table::rbindlist(
        lapply(names(keep), function(ds_nm) {
            x <- keep[[ds_nm]]
            if (!all(c("site_id", "logFC", "P.Value", "adj.P.Val") %in% names(x))) {
                return(NULL)
            }
            data.table::data.table(
                dataset = ds_nm,
                site_id = x$site_id,
                logFC   = as.numeric(x$logFC),
                pval    = as.numeric(x$P.Value),
                padj    = as.numeric(x$adj.P.Val)
            )
        }),
        fill = TRUE
    )
    if (is.null(long) || !nrow(long)) {
        return(NULL)
    }

    ## Direction = sign of logFC, magnitude = two-tailed p-value converted to z
    long[, z_i := {
        safe_p <- pmax(pval, .Machine$double.xmin)
        dirn <- sign(logFC)
        zz <- stats::qnorm(1 - safe_p / 2)
        dirn * zz
    }]

    ## Aggregate by site_id to get meta
    summ <- long[
        , .(
            n_ds       = sum(is.finite(z_i)),
            Z          = if (sum(is.finite(z_i)) > 0) sum(z_i, na.rm = TRUE) / sqrt(sum(is.finite(z_i))) else NA_real_,
            logFC_mean = mean(logFC, na.rm = TRUE),
            pval_min   = suppressWarnings(min(pval, na.rm = TRUE)),
            padj_min   = suppressWarnings(min(padj, na.rm = TRUE))
        ),
        by = .(site_id)
    ]

    ## Z -> two-tailed p, then BH to get meta FDR
    summ[, p_meta := 2 * stats::pnorm(-abs(Z))]
    summ[, padj_meta := stats::p.adjust(p_meta, method = "BH")]

    ## Sorting
    summ[, absZ := abs(Z)]
    data.table::setorder(summ, padj_meta, p_meta, -absZ, site_id)
    summ[, absZ := NULL]

    summ[]
}

## ===== Stouffer's method meta-analysis across datasets =====
dps_meta_stouffer <- function(
  dataset_dirs,
  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
  versions = c("BatchAdj"),
  out_root = "phospho_DPS_meta"
) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Please install data.table first")
    }

    for (st in strata) {
        for (ver in versions) {
            ## Collect subunit/predictor directories from all datasets
            su_tags_all <- character(0)

            for (ds_nm in names(dataset_dirs)) {
                probe_dir <- file.path(dataset_dirs[[ds_nm]], st, "DPS", ver)
                if (!dir.exists(probe_dir)) next

                this_su_dirs <- list.dirs(probe_dir, full.names = FALSE, recursive = FALSE)
                this_su_tags <- this_su_dirs[grepl("^SU__", this_su_dirs)]
                if (length(this_su_tags)) {
                    su_tags_all <- union(su_tags_all, this_su_tags)
                }
            }

            if (!length(su_tags_all)) next

            ## For each predictor/subunit
            for (su_tag in su_tags_all) {
                su_name <- sub("^SU__", "", su_tag)

                ## Collect DPS_results.csv from all datasets
                per_ds <- lapply(names(dataset_dirs), function(ds_nm) {
                    fp <- file.path(
                        dataset_dirs[[ds_nm]],
                        st, "DPS", ver,
                        su_tag,
                        "DPS_results.csv"
                    )
                    .safe_read_dps(fp)
                })
                names(per_ds) <- names(dataset_dirs)

                ## Perform Stouffer combination
                comb <- .combine_one_predictor(per_ds)
                if (is.null(comb) || !nrow(comb)) next

                ## Write out meta results
                out_dir <- file.path(out_root, st, ver, su_name)
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                out_fp <- file.path(out_dir, "DPS_meta_stouffer.csv")
                data.table::fwrite(comb, out_fp)

                if (exists("log_msg", mode = "function")) {
                    try(
                        log_msg("  [DPS meta] stratum={st} ver={ver} su={su_name} -> {out_fp}"),
                        silent = TRUE
                    )
                }
            }
        }
    }

    invisible(TRUE)
}
