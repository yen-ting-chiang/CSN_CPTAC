## ============================================================================
## Differential Expression Gene (DEG) Analysis Core Functions
## Limma-based continuous model analysis for proteomic data
## ============================================================================

# Internal helper: Build covariates data frame from sample annotation
.build_covars_df <- function(so, purity_all, sa_all) {
    USE_AGE_MISSING_INDICATOR <- isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))
    df <- data.frame(row.names = so, check.names = FALSE)
    if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
    if (!is.null(sa_all)) {
        if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
        if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
        if (USE_AGE_MISSING_INDICATOR) {
            if ("age_missing" %in% colnames(sa_all)) {
                df$age_missing <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
            }
            if ("age_z_imputed" %in% colnames(sa_all)) {
                df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
            }
        }
    }
    df
}

#' Run DEG Analysis for a Single Predictor
#'
#' Performs differential expression analysis using limma continuous model
#' with batch and covariate adjustment.
#'
#' @param predictor_name Name of the predictor (e.g., "Z_CSN_SCORE", "Z_GPS1")
#' @param predictor_vec Named numeric vector of predictor values (samples as names)
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param mat0 Raw expression matrix (genes x samples), no imputation
#' @param out_root Output root directory
#' @param batch_all Batch factor vector aligned to matrix columns
#' @param purity_all Tumor purity values (named numeric vector)
#' @param sa_all Sample annotation data frame (sex, age, age_missing, age_z_imputed)
#' @param min_per_group Minimum samples required (default from global)
#'
#' @return Invisible NULL; writes DEG results to file
#'
#' @details
#' Analysis steps:
#' 1. Align samples based on predictor availability
#' 2. Build covariate data frame (purity, sex, age)
#' 3. Add batch to design if available (BatchAdj version)
#' 4. Center predictor (mean-centering for interpretability)
#' 5. Clean covariates (remove single-level factors, all-NA columns)
#' 6. Build model.matrix with predictor + batch + covariates
#' 7. Filter genes with insufficient degrees of freedom
#' 8. Run limma::lmFit and eBayes (trend=TRUE)
#' 9. Extract topTable for predictor coefficient
#' 10. Write results to {out_root}/DEG/BatchAdj/{predictor_name}/DEG_limma_cont.csv
#'
#' Output columns: gene, logFC, t, P.Value, adj.P.Val, B
#'
#' @examples
#' run_deg_for_predictor(
#'     "Z_CSN_SCORE", csn_score, "brca_cptac_2020",
#'     ds_dir, mat0, out_root, batch, purity, sa
#' )
run_deg_for_predictor <- function(
  predictor_name,
  predictor_vec,
  ds_id, ds_dir,
  mat0,
  out_root,
  batch_all = NULL,
  purity_all = NULL,
  sa_all = NULL,
  min_per_group = NULL
) {
    # Backfill min_per_group from global settings
    if (is.null(min_per_group) || length(min_per_group) == 0L || !is.finite(min_per_group)) {
        if (exists("min_per_group", envir = .GlobalEnv, inherits = FALSE)) {
            min_per_group <- get("min_per_group", envir = .GlobalEnv, inherits = FALSE)
        } else {
            min_per_group <- 8L
        }
    }
    min_per_group <- as.integer(min_per_group)[1L]
    if (!is.finite(min_per_group) || min_per_group < 1L) min_per_group <- 8L

    ## 1) Align samples based on predictor
    if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat0)
    keep <- intersect(colnames(mat0), names(predictor_vec))
    pred_all <- suppressWarnings(as.numeric(predictor_vec[keep]))
    names(pred_all) <- keep
    fin <- is.finite(pred_all)

    if (sum(fin) < (2L * min_per_group)) {
        log_msg(
            "  [DEG-%s] predictor non-NA samples insufficient (%d < %d), skip",
            predictor_name, sum(fin), 2L * min_per_group
        )
        return(invisible(NULL))
    }

    sample_order <- keep[fin]
    pred <- pred_all[fin]
    M0 <- as.matrix(mat0[, sample_order, drop = FALSE])

    ## 2) Build covariates data frame
    df_covars0 <- .build_covars_df(sample_order, purity_all, sa_all)

    ## 3) BatchAdj version: Add batch/covars to design formula
    DF_ba <- df_covars0
    if (!is.null(batch_all)) DF_ba$batch <- factor(batch_all[sample_order])

    ok <- intersect(colnames(M0), rownames(DF_ba))
    DF_ba <- DF_ba[ok, , drop = FALSE]

    # Clean covariates: Normalize types + Remove single-level factors
    DF_ba <- coerce_covariates_safely(DF_ba)
    pred_ba <- pred[ok]
    names(pred_ba) <- ok
    M <- M0[, ok, drop = FALSE]

    # Center predictor (mean-centering)
    DF_ba$predictor <- pred_ba - mean(pred_ba, na.rm = TRUE)

    # Remove all-NA covariate columns
    if (ncol(DF_ba)) {
        all_na <- vapply(DF_ba, function(z) all(is.na(z)), logical(1))
        if (any(all_na)) DF_ba <- DF_ba[, !all_na, drop = FALSE]
    }

    # Clean again after removing all-NA columns
    if (ncol(DF_ba)) {
        DF_ba <- coerce_covariates_safely(DF_ba)
    }

    # Ensure predictor is still there
    if (!("predictor" %in% colnames(DF_ba))) {
        DF_ba$predictor <- pred_ba - mean(pred_ba, na.rm = TRUE)
    }

    # Build design matrix
    des <- stats::model.matrix(~ 1 + ., data = DF_ba)

    # Align expression matrix to samples kept in design matrix
    use <- rownames(des)
    M <- M[, use, drop = FALSE]
    pass_label <- "BatchAdj"

    ## DEG safety: Filter genes with insufficient degrees of freedom
    rnk <- qr(des)$rank
    need <- rnk + 1L # At least 1 residual degree of freedom
    nobs <- rowSums(is.finite(M))
    keep_rows <- nobs >= need

    if (!all(keep_rows)) {
        log_msg(
            "  [DEG|%s] Remove %d/%d genes (Valid samples < rank(des)+1 = %d; avoid NA df/weights)",
            pass_label, sum(!keep_rows), nrow(M), need
        )
        M <- M[keep_rows, , drop = FALSE]
    }

    if (nrow(M) == 0L) {
        log_msg("  [DEG|%s] No available genes after filtering, skip", pass_label)
        return(invisible(NULL))
    }

    ## Run limma
    fit <- limma::lmFit(M, design = des)
    eb <- limma::eBayes(fit, trend = TRUE)
    coef_name <- if ("predictor" %in% colnames(des)) "predictor" else tail(colnames(des), 1)
    tbl <- limma::topTable(eb, coef = coef_name, number = Inf, sort.by = "P")
    tbl$gene <- rownames(tbl)
    tbl <- tbl[, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "B")]

    ## Write output
    out_dir <- file.path(out_root, "DEG", pass_label, predictor_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(tbl, file.path(out_dir, "DEG_limma_cont.csv"))

    log_msg(
        "[DEG] %s | %s | %s: Wrote %s",
        ds_id, basename(out_root), predictor_name,
        file.path(out_dir, "DEG_limma_cont.csv")
    )

    invisible(NULL)
}

#' Run DEG Interaction Model (One-Step Approach)
#'
#' Performs DEG analysis for individual subunit while controlling for CSN
#' complex activity using a statistically correct one-step approach.
#'
#' @param subunit_name Subunit gene name (e.g., "GPS1", "COPS2")
#' @param subunit_vec Named numeric vector of subunit expression
#' @param csn_score_vec Named numeric vector of CSN complex scores
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param mat0 Raw expression matrix (genes x samples)
#' @param out_root Output root directory
#' @param batch_all Batch factor vector
#' @param purity_all Tumor purity values
#' @param sa_all Sample annotation data frame
#' @param min_per_group Minimum samples required
#'
#' @return Invisible NULL; writes results to file
#'
#' @details
#' This function implements the statistically correct one-step approach instead
#' of the problematic two-step residualization. It includes both subunit and
#' CSN_SCORE in the same limma model, allowing us to test the subunit's partial
#' effect while controlling for complex activity.
#'
#' Biological interpretation:
#' The coefficient for 'subunit' represents the effect of individual subunit
#' variation that is independent of CSN complex activity (same conceptual goal
#' as residualization but statistically correct).
#'
#' Model: gene ~ subunit + CSN_SCORE + batch + covariates
#'
#' Output saved to: {out_root}/DEG/BatchAdj/{subunit_name}_adj_CSN/DEG_limma_cont.csv
#'
#' @examples
#' run_deg_interaction_model(
#'     "GPS1", gps1_expr, csn_score, "brca_cptac_2020",
#'     ds_dir, mat0, out_root, batch, purity, sa
#' )
run_deg_interaction_model <- function(
  subunit_name,
  subunit_vec,
  csn_score_vec,
  ds_id, ds_dir,
  mat0,
  out_root,
  batch_all = NULL,
  purity_all = NULL,
  sa_all = NULL,
  min_per_group = NULL
) {
    # Backfill min_per_group
    if (is.null(min_per_group) || length(min_per_group) == 0L || !is.finite(min_per_group)) {
        if (exists("min_per_group", envir = .GlobalEnv, inherits = FALSE)) {
            min_per_group <- get("min_per_group", envir = .GlobalEnv, inherits = FALSE)
        } else {
            min_per_group <- 8L
        }
    }
    min_per_group <- as.integer(min_per_group)[1L]
    if (!is.finite(min_per_group) || min_per_group < 1L) min_per_group <- 8L

    ## 1) Align samples (must have both subunit and CSN_SCORE)
    if (is.null(names(subunit_vec))) names(subunit_vec) <- colnames(mat0)
    if (is.null(names(csn_score_vec))) names(csn_score_vec) <- colnames(mat0)

    common <- intersect(colnames(mat0), names(subunit_vec))
    common <- intersect(common, names(csn_score_vec))

    subunit_all <- suppressWarnings(as.numeric(subunit_vec[common]))
    csn_all <- suppressWarnings(as.numeric(csn_score_vec[common]))
    names(subunit_all) <- common
    names(csn_all) <- common

    # Both must be finite
    fin <- is.finite(subunit_all) & is.finite(csn_all)
    if (sum(fin) < (2L * min_per_group)) {
        log_msg(
            "  [DEG-INTERACTION-%s] Insufficient samples with both subunit and CSN_SCORE (%d < %d), skip",
            subunit_name, sum(fin), 2L * min_per_group
        )
        return(invisible(NULL))
    }

    sample_order <- common[fin]
    subunit <- subunit_all[fin]
    csn_score <- csn_all[fin]
    M0 <- as.matrix(mat0[, sample_order, drop = FALSE])

    ## 2) Build covariates
    df_covars0 <- .build_covars_df(sample_order, purity_all, sa_all)

    ## 3) ONE-STEP MODEL: Include both subunit and CSN_SCORE
    DF_ba <- df_covars0
    if (!is.null(batch_all)) DF_ba$batch <- factor(batch_all[sample_order])

    ok <- intersect(colnames(M0), rownames(DF_ba))
    DF_ba <- DF_ba[ok, , drop = FALSE]
    DF_ba <- coerce_covariates_safely(DF_ba)

    subunit_use <- subunit[ok]
    csn_use <- csn_score[ok]
    M <- M0[, ok, drop = FALSE]

    # Add both predictors (centered)
    DF_ba$subunit <- subunit_use - mean(subunit_use, na.rm = TRUE)
    DF_ba$CSN_SCORE <- csn_use - mean(csn_use, na.rm = TRUE)

    # Remove all-NA and single-level columns
    if (ncol(DF_ba)) {
        all_na <- vapply(DF_ba, function(z) all(is.na(z)), logical(1))
        if (any(all_na)) DF_ba <- DF_ba[, !all_na, drop = FALSE]
    }
    if (ncol(DF_ba)) {
        DF_ba <- coerce_covariates_safely(DF_ba)
    }

    # Ensure both predictors are still there
    if (!("subunit" %in% colnames(DF_ba))) {
        DF_ba$subunit <- subunit_use - mean(subunit_use, na.rm = TRUE)
    }
    if (!("CSN_SCORE" %in% colnames(DF_ba))) {
        DF_ba$CSN_SCORE <- csn_use - mean(csn_use, na.rm = TRUE)
    }

    # Build design matrix with BOTH subunit and CSN_SCORE
    des <- stats::model.matrix(~ 1 + ., data = DF_ba)
    use <- rownames(des)
    M <- M[, use, drop = FALSE]
    pass_label <- "BatchAdj"

    ## DEG safety: Filter genes with insufficient df
    rnk <- qr(des)$rank
    need <- rnk + 1L
    nobs <- rowSums(is.finite(M))
    keep_rows <- nobs >= need

    if (!all(keep_rows)) {
        log_msg(
            "  [DEG-INTERACTION|%s] Remove %d/%d genes (df insufficient)",
            pass_label, sum(!keep_rows), nrow(M)
        )
        M <- M[keep_rows, , drop = FALSE]
    }

    if (nrow(M) == 0L) {
        log_msg("  [DEG-INTERACTION|%s] No available genes, skip", pass_label)
        return(invisible(NULL))
    }

    ## Run limma
    fit <- limma::lmFit(M, design = des)
    eb <- limma::eBayes(fit, trend = TRUE)

    # Test the subunit coefficient (controlling for CSN_SCORE)
    coef_name <- "subunit"
    tbl <- limma::topTable(eb, coef = coef_name, number = Inf, sort.by = "P")
    tbl$gene <- rownames(tbl)
    tbl <- tbl[, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "B")]

    # Save with suffix indicating CSN_SCORE adjusted
    out_dir <- file.path(out_root, "DEG", pass_label, paste0(subunit_name, "_adj_CSN"))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(tbl, file.path(out_dir, "DEG_limma_cont.csv"))

    log_msg(
        "[DEG-INTERACTION] %s | %s | %s (controlling for CSN_SCORE): Wrote %s",
        ds_id, basename(out_root), subunit_name,
        file.path(out_dir, "DEG_limma_cont.csv")
    )

    invisible(NULL)
}

#' Summarize DEG Results Across Predictors
#'
#' Aggregates DEG results from multiple predictors into a wide-format summary
#' table with minimum adjusted p-value across predictors.
#'
#' @param out_root Output root directory containing DEG results
#'
#' @return Invisible NULL; writes summary files
#'
#' @details
#' Searches for DEG results in: {out_root}/DEG/BatchAdj/{predictor}/DEG_limma_cont.csv
#'
#' Creates wide-format table with columns:
#' - gene: Gene symbol
#' - min_padj: Minimum adjusted p-value across all predictors
#' - logFC_{predictor}: Log fold-change for each predictor
#' - padj_{predictor}: Adjusted p-value for each predictor
#'
#' Outputs:
#' - Summary_DEG_wide.csv
#' - Summary_DEG_wide.xlsx
#'
#' @examples
#' summarize_deg_across_predictors(out_root)
summarize_deg_across_predictors <- function(out_root) {
    log_msg("[Summary-DEG] Root directory: %s", out_root)
    sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

    ver <- "BatchAdj"
    base_dir <- file.path(out_root, "DEG", ver)

    if (!dir.exists(base_dir)) {
        log_msg("  [DEG-summary] %s | Directory not found (%s), skip", ver, base_dir)
        return(invisible(NULL))
    }

    preds <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
    if (!length(preds)) {
        log_msg("  [DEG-summary] %s | No predictor subdirectories, skip", ver)
        return(invisible(NULL))
    }

    tabs <- list()
    for (p in preds) {
        fp <- file.path(base_dir, p, "DEG_limma_cont.csv")
        if (!file.exists(fp)) next

        tb <- tryCatch(data.table::fread(fp, showProgress = FALSE), error = function(e) NULL)
        if (is.null(tb) || !nrow(tb)) next

        need <- c("gene", "logFC", "adj.P.Val")
        if (!all(need %in% names(tb))) next

        tb_slim <- tb[, ..need]
        data.table::setnames(tb_slim,
            old = c("gene", "logFC", "adj.P.Val"),
            new = c(
                "gene",
                paste0("logFC_", sfn(p)),
                paste0("padj_", sfn(p))
            )
        )
        tabs[[p]] <- tb_slim
    }

    if (!length(tabs)) {
        log_msg("  [DEG-summary] %s | Cannot find any available DEG_limma_cont.csv, skip", ver)
        return(invisible(NULL))
    }

    # Merge all tables
    df <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), tabs)

    # Calculate minimum padj across predictors
    padj_cols <- grep("^padj_", names(df), value = TRUE)
    if (length(padj_cols)) {
        suppressWarnings({
            df$min_padj <- apply(as.data.frame(df[, ..padj_cols]), 1, function(z) {
                z <- suppressWarnings(as.numeric(z))
                if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)
            })
        })
        data.table::setcolorder(df, c(
            "gene", "min_padj",
            setdiff(names(df), c("gene", "min_padj"))
        ))
        df <- df[order(df$min_padj), ]
    }

    # Write outputs
    out_dir <- file.path(out_root, "DEG", ver)
    base <- file.path(out_dir, "Summary_DEG_wide")

    data.table::fwrite(df, paste0(base, ".csv"))

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "DEG_wide")
    openxlsx::writeData(wb, "DEG_wide", df)
    openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)

    log_msg(
        "[DEG-summary] %s | predictors=%d | genes=%d: Wrote %s",
        ver, length(tabs), nrow(df), paste0(base, ".csv")
    )

    invisible(NULL)
}
