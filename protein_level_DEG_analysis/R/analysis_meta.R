## ============================================================================
## Meta-Analysis Functions for DEG Results
## Implements Stouffer's method for combining p-values across datasets
## ============================================================================

#' Stouffer's Z-score Method for Combining P-values
#'
#' Combines p-values from independent tests using Stouffer's Z-score method.
#'
#' @param pvals Numeric vector of p-values (0 < p <= 1)
#' @param weights Optional numeric vector of weights (e.g., sample sizes)
#' @param two_tailed Logical, whether p-values are two-tailed (default TRUE)
#'
#' @return List with z_score, p_value, n_studies

#' Helper: Convert P-values to Directional Z-scores
.meta_from_p_and_sign <- function(p, signv) {
    p <- pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - .Machine$double.eps)
    s <- sign(as.numeric(signv))
    z <- s * stats::qnorm(1 - p / 2)
    z[!is.finite(s)] <- NA_real_
    z
}

#' Meta-FDR Analysis for DEG
meta_fdr_stouffer_DEG <- function(dataset_dirs,
                                  strata = c("ALL", "TP53_mutant", "TP53_wild_type"),
                                  versions = c("BatchAdj"),
                                  output_prefix = "proteomic_DEG") {
    if (is.null(dataset_dirs) || !length(dataset_dirs)) {
        stop("dataset_dirs cannot be empty")
    }

    ds_ids <- names(dataset_dirs)

    deg_root_path <- function(ds, ver, st) {
        file.path(output_prefix, ds, "csn_gsea_results_TP53", st, "DEG", ver)
    }

    for (st in strata) {
        for (ver in versions) {
            log_msg("\n[Meta-DEG] Processing %s | %s", st, ver)

            preds_all <- unique(unlist(lapply(ds_ids, function(ds) {
                base <- deg_root_path(ds, ver, st)
                if (!dir.exists(base)) {
                    return(character(0))
                }
                subdirs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
                subdirs[nzchar(subdirs) & subdirs != "."]
            })))

            if (!length(preds_all)) {
                log_msg("  [Skip] No predictors found for %s | %s", st, ver)
                next
            }

            log_msg("  Found %d predictors", length(preds_all))

            for (su in preds_all) {
                lst <- lapply(ds_ids, function(ds) {
                    fp <- file.path(deg_root_path(ds, ver, st), su, "DEG_limma_cont.csv")
                    if (!file.exists(fp)) {
                        return(NULL)
                    }
                    dt <- tryCatch(data.table::fread(fp, na.strings = c("NA", "NaN", "")), error = function(e) NULL)
                    if (is.null(dt) || !nrow(dt)) {
                        return(NULL)
                    }

                    nm <- tolower(names(dt))
                    names(dt) <- nm

                    if (!("gene" %in% nm)) {
                        if ("pathway" %in% nm) names(dt)[match("pathway", nm)] <- "gene" else return(NULL)
                    }

                    need <- c("gene", "logfc", "p.value")
                    if (!all(need %in% names(dt))) {
                        return(NULL)
                    }

                    data.frame(
                        gene = as.character(dt$gene),
                        logfc = suppressWarnings(as.numeric(dt$logfc)),
                        pval = suppressWarnings(as.numeric(dt$p.value)),
                        stringsAsFactors = FALSE
                    )
                })

                keep <- lst[!vapply(lst, is.null, logical(1))]
                if (!length(keep)) next

                ds_ids_keep <- ds_ids[!vapply(lst, is.null, logical(1))]
                all <- do.call(rbind, Map(function(df, id) {
                    df$dataset <- id
                    as.data.frame(df, stringsAsFactors = FALSE)
                }, keep, ds_ids_keep))

                if (is.null(all) || !nrow(all)) next

                all <- all[is.finite(all$pval) & is.finite(all$logfc), , drop = FALSE]
                if (!nrow(all)) next

                all$Z <- .meta_from_p_and_sign(all$pval, all$logfc)

                z_tab <- split(all, all$gene)
                z_sum <- lapply(names(z_tab), function(g) {
                    z <- z_tab[[g]]$Z
                    k <- sum(is.finite(z))
                    Zmeta <- if (k) sum(z, na.rm = TRUE) / sqrt(k) else NA_real_
                    pmeta <- if (is.finite(Zmeta)) 2 * pnorm(-abs(Zmeta)) else NA_real_
                    c(gene = g, Z = Zmeta, p_meta = pmeta, n_ds = k)
                })

                Zdf <- data.frame(do.call(rbind, z_sum), stringsAsFactors = FALSE)
                Zdf$Z <- suppressWarnings(as.numeric(Zdf$Z))
                Zdf$p_meta <- suppressWarnings(as.numeric(Zdf$p_meta))
                Zdf$padj_meta <- p.adjust(Zdf$p_meta, method = "BH")

                Zdf$pathway <- Zdf$gene
                Zdf <- Zdf[, c("pathway", "Z", "p_meta", "padj_meta", "n_ds")]

                out_dir <- file.path(output_prefix, "csn_deg_pan_summary_TP53", "meta_fdr", st, ver, su)
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

                out_file <- file.path(out_dir, "DEG_meta_fdr.csv")
                data.table::fwrite(Zdf, out_file)

                log_msg("    [%s] %d genes → %s", su, nrow(Zdf), out_file)
            }
        }
    }

    invisible(TRUE)
}
