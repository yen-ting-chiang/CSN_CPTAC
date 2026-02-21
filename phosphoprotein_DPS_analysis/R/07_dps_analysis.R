## =========================================================
## 07_dps_analysis.R
## Core DPS (Differential Phosphosite) Analysis Functions
## =========================================================


.build_design_for_version <- function(version, pred_centered, sample_order,
                                      batch_all, purity_all, sa_all,
                                      USE_AGE_MISSING_INDICATOR = FALSE,
                                      ds_id = NULL, mat_for_sv = NULL) {
    # Consistent with run_predictor_analyses design: base covariates
    DF <- data.frame(
        pred_centered = as.numeric(pred_centered[sample_order]),
        row.names = sample_order, check.names = FALSE
    )
    if (!is.null(sa_all)) {
        if ("sex" %in% colnames(sa_all)) DF$sex <- factor(sa_all[sample_order, "sex"])
        if ("age" %in% colnames(sa_all)) DF$age <- suppressWarnings(as.numeric(sa_all[sample_order, "age"]))
        if (USE_AGE_MISSING_INDICATOR) {
            if ("age_missing" %in% colnames(sa_all)) DF$age_missing <- suppressWarnings(as.numeric(sa_all[sample_order, "age_missing"]))
            if ("age_z_imputed" %in% colnames(sa_all)) DF$age_z_imputed <- suppressWarnings(as.numeric(sa_all[sample_order, "age_z_imputed"]))
        }
    }
    if (identical(version, "BatchAdj")) {
        if (!is.null(batch_all)) DF$batch <- droplevels(batch_all[sample_order])
    }
    # Clean singular or empty columns
    DF <- coerce_covariates_safely(DF)
    DF
}

## ---- Small utility: align any vector/data.frame to specified sample order ----
.align_to_samples <- function(x, sam, what = "covariate") {
    if (is.null(x)) {
        return(NULL)
    }
    if (is.vector(x) || is.factor(x)) {
        # Allow length exactly equal to sam without names; else prefer align by name
        if (is.null(names(x))) {
            if (length(x) != length(sam)) stop(sprintf("[%s] Length mismatch: %d vs %d", what, length(x), length(sam)))
            names(x) <- sam
        }
        if (!all(sam %in% names(x))) stop(sprintf("[%s] Missing samples: %s", what, paste(setdiff(sam, names(x)), collapse = ", ")))
        out <- x[sam]
        return(out)
    } else {
        x <- as.data.frame(x)
        # Expect rownames(x) to be samples
        if (is.null(rownames(x))) {
            if (nrow(x) != length(sam)) stop(sprintf("[%s] Row count mismatch: %d vs %d and no rownames for alignment", what, nrow(x), length(sam)))
            rownames(x) <- sam
        }
        if (!all(sam %in% rownames(x))) stop(sprintf("[%s] Missing samples: %s", what, paste(setdiff(sam, rownames(x)), collapse = ", ")))
        out <- x[sam, , drop = FALSE]
        return(out)
    }
}

# Any object -> matrix; force set rownames; return NULL if 0 columns or row count doesn't match specified rows
.mk_mat_or_null <- function(x, rows) {
    if (is.null(x)) {
        return(NULL)
    }
    m <- as.matrix(x)
    if (is.null(rownames(m))) rownames(m) <- rows
    if (nrow(m) != length(rows) || ncol(m) == 0) {
        return(NULL)
    }
    m
}

## Safely convert vector to factor, and make NA explicit as "NA" (then add prefix later)
.factorize_with_explicit_NA <- function(x) {
    x_chr <- as.character(x)
    x_chr[!nzchar(x_chr) | is.na(x_chr)] <- "NA" # Empty string/NA both treated as "NA"
    factor(x_chr)
}

.run_dps_fit <- function(M, pred, version, ds_id, out_dir,
                         batch_all, purity_all, sa_all, mat_for_sv = NULL,
                         extra_covars = NULL) {
    ## 1) Align samples and predictors
    sample_order <- intersect(colnames(M), names(pred))
    pred <- pred[sample_order]
    M <- as.matrix(M[, sample_order, drop = FALSE])

    ## 2) Assemble design table from template (using existing .build_design_for_version)
    DF <- .build_design_for_version(
        version = version,
        pred_centered = pred, # Template already centered: use pred_centered to indicate
        sample_order = sample_order,
        batch_all = batch_all,
        purity_all = purity_all,
        sa_all = sa_all,
        USE_AGE_MISSING_INDICATOR = isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE)),
        ds_id = ds_id,
        mat_for_sv = mat_for_sv
    )

    ##  Merge extra_covars if provided 
    if (!is.null(extra_covars)) {
        extra_df <- .align_to_samples(extra_covars, sample_order, what = "extra_covars")
        extra_df <- as.data.frame(extra_df, check.names = FALSE)
        # Combine (be careful with name collisions, though usually distinctive)
        DF <- cbind(DF, extra_df)
    }

    ## 3) Select only template-allowed RHS variables; rename pred_centered to predictor as main coefficient
    # Include columns from extra_covars in RHS
    rhs_candidates <- c("sex", "age", "age_missing", "age_z_imputed", "batch", paste0("SV", 1:5))
    if (!is.null(extra_covars)) rhs_candidates <- c(rhs_candidates, colnames(extra_df))

    rhs <- intersect(rhs_candidates, colnames(DF))

    DF2 <- data.frame(
        predictor = DF$pred_centered,
        DF[, rhs, drop = FALSE],
        row.names = rownames(DF),
        check.names = FALSE
    )

    ## 4) Drop constant/single-level/zero-variance columns (always keep predictor); also drop NA rows and align M
    keep_cols <- rep(TRUE, ncol(DF2))
    names(keep_cols) <- colnames(DF2)
    if (ncol(DF2) > 1L) {
        keep_cols <- vapply(DF2, function(z) {
            if (is.factor(z)) {
                nz <- droplevels(z)
                (nlevels(nz) >= 2) && (sum(!is.na(nz)) >= 2)
            } else {
                vz <- suppressWarnings(stats::var(as.numeric(z), na.rm = TRUE))
                is.finite(vz) && vz > 0 && sum(is.finite(z)) >= 2
            }
        }, logical(1))
        if ("predictor" %in% names(keep_cols)) keep_cols["predictor"] <- TRUE
        
    }
    DF2 <- DF2[, keep_cols, drop = FALSE]

    row_keep <- stats::complete.cases(DF2)
    if (sum(row_keep) < 3L) {
        log_msg("[DPS|%s|%s] Design has insufficient usable samples (<3) -> skip", ds_id, version)
        return(invisible(NULL))
    }
    DF2 <- DF2[row_keep, , drop = FALSE]
    M <- M[, rownames(DF2), drop = FALSE]
    for (cn in colnames(DF2)) if (is.factor(DF2[[cn]])) DF2[[cn]] <- droplevels(DF2[[cn]])

    stopifnot(identical(colnames(M), rownames(DF2)))

    ## 5) Build design matrix (with intercept): ~ predictor + RHS
    make_X <- function(df, rhs_vars) {
        if (length(rhs_vars)) {
            stats::model.matrix(stats::as.formula(paste("~ predictor +", paste(rhs_vars, collapse = " + "))), data = df)
        } else {
            stats::model.matrix(~predictor, data = df)
        }
    }
    rhs_use <- setdiff(colnames(DF2), "predictor")
    des <- make_X(DF2, rhs_use)

    ## 6) Saturation protection: if residual df <= 0, sequentially remove covariates until df_res > 0
    rank_d <- qr(des)$rank
    res_df <- nrow(des) - rank_d
    if (res_df <= 0L) {
        #  Drop order: SVs first, then batch/clinical. Keep extra_covars (biology) as long as possible.
        drop_candidates <- intersect(
            c(paste0("SV", 5:1), "batch", "sex", "age_missing", "age_z_imputed", "age"),
            rhs_use
        )
        # If we need to drop extra_covars (e.g. CSN_SCORE) to run at all, that might be needed too, but put them last.
        if (!is.null(extra_covars)) drop_candidates <- c(drop_candidates, colnames(extra_df))

        for (cv in drop_candidates) {
            rhs_use <- setdiff(rhs_use, cv)
            des <- make_X(DF2, rhs_use)
            rank_d <- qr(des)$rank
            res_df <- nrow(des) - rank_d
            log_msg("[DPS|%s|%s] Saturation protection: remove covariate=%s -> df_res=%d", ds_id, version, cv, res_df)
            if (res_df > 0L) break
        }
    }
    if (res_df <= 0L) {
        log_msg("[DPS|%s|%s] No residual df after removing all optional covariates -> skip", ds_id, version)
        return(invisible(NULL))
    }
    if (nrow(des) != ncol(M)) {
        stop(sprintf(
            "[DPS|%s|%s] Design rows (%d) != matrix columns (%d): sample alignment failed",
            ds_id, version, nrow(des), ncol(M)
        ))
    }

    ## 7) Fit and output (using predictor coefficient as DPS statistic)
    fit <- limma::eBayes(limma::lmFit(M, des))
    j <- match("predictor", colnames(des))
    if (is.na(j)) stop("[DPS] Design matrix coefficient not found: predictor")
    tt <- limma::topTable(fit, coef = j, number = Inf, sort.by = "none")

    # Standardize output: site ID and observation count
    if ("ID" %in% colnames(tt)) {
        # If limma has the actual site ID (e.g. "AHNAK_S18") in tt$ID
        tt$site_id <- tt$ID
    } else {
        # Otherwise use rownames as site ID
        tt$site_id <- rownames(tt)
    }

    tt$n_obs <- rowSums(is.finite(M))

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(tt, file.path(out_dir, "DPS_results.csv"))
    invisible(tt)
}

run_dps_stratum <- function(ds_id, ds_dir, M_site_full, sample_keep, out_root, prot0_full) {
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

    sam <- intersect(colnames(M_site_full), sample_keep)
    if (length(sam) < (2L * min_per_group)) {
        log_msg("  [DPS] %s insufficient samples -> skip", basename(out_root))
        return(invisible(NULL))
    }
    M0 <- M_site_full[, sam, drop = FALSE]
    prot0 <- prot0_full[, sam, drop = FALSE]

    present_sub <- intersect(csn_subunits, rownames(prot0))
    if (!length(present_sub)) {
        log_msg("  [DPS] No CSN subunit -> skip")
        return(invisible(NULL))
    }

    csn_score <- build_csn_score_safe(
        prot0,
        subunits = present_sub, combine_7AB = TRUE,
        min_members = 5L, pca_min_samples = 10L
    )

    # Covariates / batch consistent with template
    batch_all <- get_batch_factor_phospho(ds_dir, sam)
    batch_all <- if (!is.null(batch_all)) droplevels(batch_all$fac[sam]) else NULL
    purity_all <- get_purity_covariate(ds_id, ds_dir, sam)
    sa_all <- get_sex_age_covariates(ds_dir, sam)

    run_passes <- getOption("csn.run_passes", c("BatchAdj"))
    for (version in run_passes) {
        out_ver <- file.path(out_root, version)
        dir.create(out_ver, recursive = TRUE, showWarnings = FALSE)

        ## (1) subunits
        for (su in present_sub) {
            log_msg("  [DPS|%s|%s] predictor=%s", version, basename(out_root), su)
            .run_dps_fit(
                M = M0, pred = prot0[su, ], version = version, ds_id = ds_id,
                out_dir = file.path(out_ver, paste0("SU__", safe_fs_name(su))),
                batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
                mat_for_sv = M0
            )
        }

        ## (2) CSN_SCORE
        if (sum(is.finite(csn_score)) >= (2L * min_per_group)) {
            log_msg("  [DPS|%s|%s] predictor=CSN_SCORE", version, basename(out_root))
            .run_dps_fit(
                M = M0, pred = csn_score, version = version, ds_id = ds_id,
                out_dir = file.path(out_ver, "SU__CSN_SCORE"),
                batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
                mat_for_sv = M0
            )
        }

        ## (3) RESIDUAL_<SU> 
        for (su in present_sub) {
            if (!sum(is.finite(csn_score))) next

            
            log_msg("  [DPS|%s|%s] predictor=RESIDUAL_%s (One-Step: ~ %s + CSN_SCORE)", version, basename(out_root), su, su)

            .run_dps_fit(
                M = M0, pred = prot0[su, ], version = version, ds_id = ds_id,
                out_dir = file.path(out_ver, paste0("SU__RESIDUAL_", safe_fs_name(su))),
                batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
                mat_for_sv = M0,
                extra_covars = data.frame(row.names = colnames(prot0), CSN_SCORE = csn_score)
            )
        }

        ## (4) Summarize all predictor results within this stratum x version into a wide table
        {
            
            pred_dirs <- list.dirs(out_ver, full.names = TRUE, recursive = FALSE)
            pred_dirs <- pred_dirs[grepl("^SU__", basename(pred_dirs))]

            wide_list <- list()

            for (pd in pred_dirs) {
                pred_label <- sub("^SU__", "", basename(pd)) # e.g. "GPS1", "CSN_SCORE", "RESIDUAL_GPS1"
                f_csv <- file.path(pd, "DPS_results.csv")
                if (!file.exists(f_csv)) next

                tt_local <- try(data.table::fread(f_csv), silent = TRUE)
                if (inherits(tt_local, "try-error") || !nrow(tt_local)) next

                
                if ("site_id" %in% colnames(tt_local)) {
                    site_vec <- tt_local$site_id
                } else if ("ID" %in% colnames(tt_local)) {
                    site_vec <- tt_local$ID
                } else {
                    site_vec <- rownames(tt_local)
                }

                tmp <- data.frame(
                    site_id = site_vec,
                    logFC = tt_local$logFC,
                    adj.P.Val = tt_local$adj.P.Val,
                    stringsAsFactors = FALSE
                )

                ## If same site_id appears multiple times for this predictor (e.g. isoform/duplicate rownames), keep only first row to avoid Cartesian explosion during merge
                tmp <- tmp[!duplicated(tmp$site_id), , drop = FALSE]

                ## Rename to <predictor>.logFC / <predictor>.adj.P.Val
                colnames(tmp)[colnames(tmp) == "logFC"] <- paste0(pred_label, ".logFC")
                colnames(tmp)[colnames(tmp) == "adj.P.Val"] <- paste0(pred_label, ".adj.P.Val")

                wide_list[[pred_label]] <- tmp
            }

            if (length(wide_list)) {
                summary_wide <- Reduce(
                    function(x, y) merge(x, y, by = "site_id", all = TRUE),
                    wide_list
                )
                data.table::fwrite(summary_wide, file.path(out_ver, "DPS_summary_wide.csv"))
            }
        }
    }
    invisible(NULL)
}
