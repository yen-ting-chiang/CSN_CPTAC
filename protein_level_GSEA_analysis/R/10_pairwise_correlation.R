## ====================================================================
## 10_pairwise_correlation.R
##
## Purpose: CSN subunits pairwise correlation analysis
## Contains: Correlation calculation, residualization, heatmap visualization
## ====================================================================

#' Residualize a vector to covariates and batch effects
#'
#' @param y Named numeric vector
#' @param batch Batch factor (optional)
#' @param covars Data frame of covariates (optional)
#' @param min_n Minimum number of samples required
#' @return Named numeric vector of residuals
#' @export
residualize_to_covars <- function(y, batch = NULL, covars = NULL, min_n = 8L) {
    stopifnot(!is.null(names(y)))
    sam <- names(y)
    DF <- data.frame(row.names = sam, check.names = FALSE)
    if (!is.null(batch)) DF[["batch"]] <- droplevels(batch[sam])
    if (!is.null(covars)) {
        C <- as.data.frame(covars, check.names = FALSE)
        rn <- rownames(C)
        if (is.null(rn)) stop("[residualize_to_covars] covars need rownames=sample")
        C <- C[sam, , drop = FALSE]

        if (exists("coerce_covariates_safely", mode = "function")) C <- coerce_covariates_safely(C)
        for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
    }
    yv <- suppressWarnings(as.numeric(y[sam]))
    # Clear all NA/Constant/Single Level columns
    if (ncol(DF) > 0) {
        all_na <- vapply(DF, function(v) all(is.na(v)), logical(1))
        if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]
        is_const <- vapply(DF, function(v) {
            vv <- v[!is.na(v)]
            if (!length(vv)) TRUE else if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
        }, logical(1))
        if (any(is_const)) DF <- DF[, !is_const, drop = FALSE]
    }
    ok <- is.finite(yv) & (if (ncol(DF) == 0) TRUE else stats::complete.cases(DF))
    if (sum(ok) < min_n) {
        return(setNames(rep(NA_real_, length(y)), names(y)))
    }
    des <- if (ncol(DF) == 0) model.matrix(~1) else model.matrix(~ 1 + ., data = DF[ok, , drop = FALSE])
    fit <- lm.fit(x = des, y = yv[ok])
    out <- setNames(rep(NA_real_, length(y)), names(y))
    out[ok] <- fit$residuals
    out
}


#' Calculate pairwise correlations for CSN subunits in one dataset
#'
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory
#' @param out_root Output root directory
#' @param subunits CSN subunits to analyze
#' @param min_pairs Minimum number of pairs required
#' @param dpi DPI for output plots
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param point_size Point size for scatter plots
#' @return Data frame of correlation results (invisible)
#' @export
csn_pairwise_correlation_one_ds <- function(
  ds_id, ds_dir,
  out_root = "CSN_subunits_correlation_coefficient",
  subunits = get0("csn_subunits", ifnotfound = c("GPS1", "COPS2", "COPS3", "COPS4", "COPS5", "COPS6", "COPS7A", "COPS7B", "COPS8", "COPS9")),
  min_pairs = 10L,
  dpi = 600,
  width = 4.5, height = 4.5,
  point_size = 1.6
) {
    dir.create(out_root, showWarnings = FALSE)
    ds_out <- file.path(out_root, ds_id)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

    # Read matrix
    mat0 <- load_matrix_from_dataset_dir(ds_dir)
    present <- intersect(subunits, rownames(mat0))
    if (length(present) < 2) {
        log_msg("[pairwise] %s: available CSN subunits < 2, skip", ds_id)
        return(invisible(NULL))
    }
    sam_all <- colnames(mat0)

    # Covariate candidates (purity/sex/age)
    build_covars_df <- function(ds_id, ds_dir, sample_ids) {
        pur <- get_purity_covariate(ds_id, ds_dir, sample_ids)
        sa <- get_sex_age_covariates(ds_dir, sample_ids)
        df <- data.frame(
            purity = as.numeric(pur),
            sex = as.numeric(sa[, "sex"]),
            age = as.numeric(sa[, "age"]),
            row.names = sample_ids,
            check.names = FALSE
        )
        if (exists("coerce_covariates_safely", mode = "function")) {
            df <- coerce_covariates_safely(df)
        }
        df
    }
    cov0 <- build_covars_df(ds_id, ds_dir, sam_all)
    cov_raw <- cov0

    # Batch factor
    bi <- get_batch_factor(ds_dir, sam_all)
    batch <- NULL
    if (!is.null(bi) && !is.null(bi$fac)) {
        b <- bi$fac
        if (is.null(names(b))) {
            if (length(b) == length(sam_all)) {
                b <- stats::setNames(b, sam_all)
            } else {
                b <- NULL
            }
        }
        if (!is.null(b)) batch <- droplevels(b[sam_all])
    }

    # Imputation
    min_frac_complete <- get0("min_frac_complete", ifnotfound = 0.75)
    keep_rows <- rowMeans(is.finite(mat0)) >= min_frac_complete
    Mimp <- mat0[keep_rows, , drop = FALSE]
    set.seed(1234)
    Mimp <- imputeLCMD::impute.MinProb(Mimp, q = 0.01)

    # Generate all pairs
    pairs <- utils::combn(present, 2, simplify = FALSE)

    # Result collector
    all_rows <- list()

    # Helper for saving scatter plots (optional)
    .save_scatter <- function(df_xy, title, out_base) {
        gp <- ggplot2::ggplot(df_xy, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_point(size = point_size, alpha = 0.8) +
            ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
            ggplot2::labs(title = title, x = df_xy$gx[1], y = df_xy$gy[1]) +
            ggplot2::theme_classic(base_size = 11) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5),
                panel.background = ggplot2::element_rect(fill = "white", colour = NA),
                plot.background = ggplot2::element_rect(fill = "white", colour = NA)
            )
        MAKE_PLOTS <- get0("MAKE_PLOTS", ifnotfound = FALSE)
        if (isTRUE(MAKE_PLOTS)) {
            ggplot2::ggsave(paste0(out_base, ".tiff"), gp, width = width, height = height, dpi = dpi, compression = "lzw")
            ggplot2::ggsave(paste0(out_base, ".png"), gp, width = width, height = height, dpi = dpi)
            ggplot2::ggsave(paste0(out_base, ".jpg"), gp, width = width, height = height, dpi = dpi)
        }
    }

    # Main loop: Each pair of subunits
    for (p in pairs) {
        gx <- p[1]
        gy <- p[2]
        x <- as.numeric(mat0[gx, sam_all])
        names(x) <- sam_all
        y <- as.numeric(mat0[gy, sam_all])
        names(y) <- sam_all
        ok <- stats::complete.cases(x, y)
        if (sum(ok) < min_pairs) next

        ## Version 1: NoCovariate
        rr <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "pearson"))
        fit <- stats::lm(y ~ x, data = data.frame(x = x[ok], y = y[ok]))
        all_rows[[length(all_rows) + 1L]] <- data.frame(
            dataset = ds_id, version = "NoCovariate",
            gene_x = gx, gene_y = gy, n = sum(ok),
            pearson_r = as.numeric(rr$estimate), pearson_p = rr$p.value,
            R2 = summary(fit)$r.squared,
            slope = unname(stats::coef(fit)[["x"]]),
            intercept = unname(stats::coef(fit)[["(Intercept)"]]),
            stringsAsFactors = FALSE, check.names = FALSE
        )

        .save_scatter(
            df_xy = data.frame(x = x[ok], y = y[ok], gx = gx, gy = gy),
            title = sprintf("%s | %s vs %s (NoCovariate)", ds_id, gx, gy),
            out_base = file.path(ds_out, sprintf("%s_vs_%s_NoCovariate", safe_fs_name(gx), safe_fs_name(gy)))
        )

        ## Version 2: BatchAdj_covars
        cov_ba <- cov_raw
        batch_for_res <- NULL
        if (!is.null(batch)) {
            batch_for_res <- batch
        }

        xb <- residualize_to_covars(x, batch = batch_for_res, covars = cov_ba)
        yb <- residualize_to_covars(y, batch = batch_for_res, covars = cov_ba)
        ok3 <- stats::complete.cases(xb, yb)

        if (sum(ok3) >= min_pairs) {
            rr3 <- suppressWarnings(stats::cor.test(xb[ok3], yb[ok3], method = "pearson"))
            fit3 <- stats::lm(yb ~ xb, data = data.frame(xb = xb[ok3], yb = yb[ok3]))
            all_rows[[length(all_rows) + 1L]] <- data.frame(
                dataset = ds_id, version = "BatchAdj_covars",
                gene_x = gx, gene_y = gy, n = sum(ok3),
                pearson_r = as.numeric(rr3$estimate), pearson_p = rr3$p.value,
                R2 = summary(fit3)$r.squared,
                slope = unname(stats::coef(fit3)[["xb"]]),
                intercept = unname(stats::coef(fit3)[["(Intercept)"]]),
                stringsAsFactors = FALSE, check.names = FALSE
            )
            .save_scatter(
                df_xy = data.frame(x = xb[ok3], y = yb[ok3], gx = gx, gy = gy),
                title = sprintf("%s | %s vs %s (BatchAdj_covars)", ds_id, gx, gy),
                out_base = file.path(ds_out, sprintf("%s_vs_%s_BatchAdj_covars", safe_fs_name(gx), safe_fs_name(gy)))
            )
        }
    }

    # Summary and output
    if (!length(all_rows)) {
        log_msg("[pairwise] %s: No available pairs (insufficient samples or too many NAs)", ds_id)
        return(invisible(NULL))
    }
    RES <- do.call(rbind, all_rows)

    # BH correction by dataset + version
    RES$pearson_padj <- ave(RES$pearson_p, interaction(RES$dataset, RES$version, drop = TRUE),
        FUN = function(p) stats::p.adjust(p, method = "BH")
    )

    # Write CSV
    out_csv <- file.path(ds_out, "pairwise_correlations_all_versions.csv")
    data.table::fwrite(RES, out_csv)

    # Write Excel
    out_xlsx <- file.path(ds_out, "pairwise_correlations_all_versions.xlsx")
    wb <- openxlsx::createWorkbook()
    for (ver in unique(RES$version)) {
        openxlsx::addWorksheet(wb, ver)
        openxlsx::writeData(wb, ver, RES[RES$version == ver, ], withFilter = TRUE)
    }
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    log_msg("[pairwise] %s: completed. CSV: %s; XLSX: %s", ds_id, basename(out_csv), basename(out_xlsx))
    invisible(RES)
}


#' Run pairwise correlations for all datasets
#'
#' @param dataset_dirs_map Named vector of dataset directories
#' @param out_root Output root directory
#' @return Invisible TRUE
#' @export
run_csn_subunits_pairwise_correlations <- function(
  dataset_dirs_map = NULL,
  out_root = "CSN_subunits_correlation_coefficient"
) {
    if (is.null(dataset_dirs_map)) {
        if (exists("dataset_dirs_run")) {
            dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
        } else if (exists("dataset_dirs")) {
            dd <- get("dataset_dirs", inherits = TRUE)
            dataset_dirs_map <- dd[dir.exists(dd) & file.exists(file.path(dd, "data_protein_quantification.txt"))]
        } else {
            stop("cannot find dataset_dirs_run or dataset_dirs")
        }
    }
    dir.create(out_root, showWarnings = FALSE)
    for (ds in names(dataset_dirs_map)) {
        try(csn_pairwise_correlation_one_ds(ds_id = ds, ds_dir = dataset_dirs_map[[ds]], out_root = out_root), silent = FALSE)
    }
    invisible(TRUE)
}


# ---- End of 10_pairwise_correlation.R ----
