# =============================================================================
# BRCA_immune_main.R
# CSN Subunits -> BRCA Immune Score Correlation Analysis
# =============================================================================
#
# This script performs correlation analysis between CSN subunits (protein
# abundance) and immune/stromal scores in BRCA CPTAC dataset.
#
# Predictors:
#   - Individual CSN subunits (GPS1, COPS2-9)
#   - CSN_SCORE (PC1 of subunits)
#   - RESIDUAL_<subunit> (subunit expression regressed on CSN_SCORE)
#
# Outcomes (from data_clinical_patient.txt):
#   - ESTIMATE_IMMUNE_SCORE
#   - ESTIMATE_STROMAL_SCORE
#   - XCELL_IMMUNE_SCORE
#   - XCELL_STROMAL_SCORE
#   - CIBERSORT_ABSOLUTE_SCORE
#
# Output: BRCA_immune_analysis/<stratum>/<version>/immune_vs_predictor_*.csv
#
# =============================================================================

# =============================================================================
# WORKING DIRECTORY SETUP
# =============================================================================
# NOTE: This script should be run from the CSN_CPTAC project root directory.
# Please refer to README.md for data download instructions and setup guide.
#
# Example:
#   setwd("/path/to/CSN_CPTAC")
#   source("BRCA_CSN_immune_analysis/BRCA_immune_main.R")

message("====================================================================")
message("BRCA CSN Immune Score Correlation Analysis")
message("====================================================================\n")

# Verify we are in the correct directory
if (!file.exists("BRCA_CSN_immune_analysis/BRCA_immune_main.R")) {
    stop(
        "ERROR: Script must be run from CSN_CPTAC project root directory.\n",
        "Current directory: ", getwd(), "\n",
        "Please setwd() to the CSN_CPTAC root directory."
    )
}

message("Working directory: ", getwd(), "\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(janitor)
    library(glue)
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
    library(S4Vectors)
    library(limma)
    library(imputeLCMD)
    library(msigdbr)
    library(openxlsx)
    library(ComplexHeatmap)
    library(cowplot)
    library(matrixStats)
    library(yaml)
})

# =============================================================================
# SOURCE UTILITY MODULES
# =============================================================================
# Source all utilities from BRCA_CSN_immune_analysis/R/
# These must be sourced in order due to dependencies

script_dir <- "BRCA_CSN_immune_analysis/R"

source(file.path(script_dir, "01_config.R"))
source(file.path(script_dir, "02_utils_logging.R"))
source(file.path(script_dir, "03_utils_io.R"))
source(file.path(script_dir, "04_utils_csn_score.R"))
source(file.path(script_dir, "05_utils_batch.R"))
source(file.path(script_dir, "06_utils_tp53.R"))
source(file.path(script_dir, "07_utils_covariates.R"))
source(file.path(script_dir, "08_utils_correlation.R"))
source(file.path(script_dir, "09_utils_plotting.R"))
source(file.path(script_dir, "10_utils_audit.R"))

# =============================================================================
# INITIALIZATION
# =============================================================================

# Force serial execution for reproducibility
.force_serial_execution()

# Create run_info directory
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

datasets_root <- getwd()

dataset_ids <- c("brca_cptac_2020")

dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL", "TP53_mutant", "TP53_wild_type")

message("datasets_root = ", datasets_root)

# Filter to only existing datasets with protein data
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
    log_msg("Detected %d missing directories, will skip: %s", length(missing_dirs), paste(missing_dirs, collapse = ", "))
}

dataset_dirs_run <- dataset_dirs[
    dir.exists(dataset_dirs) &
        file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
]

if (!length(dataset_dirs_run)) {
    stop("dataset_dirs_run is empty, please check if directories and data_protein_quantification.txt exist")
}

log_msg("Available datasets this round (%d): %s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse = ", "))

# =============================================================================
# PLOT CONFIGURATION (Optional)
# =============================================================================

# Font sizes for plots
BRCA_TITLE_SIZE <- 8
BRCA_SUBTITLE_SIZE <- 5
BRCA_AXIS_TITLE_SIZE <- 8
BRCA_AXIS_TEXT_SIZE <- 8

# -----------------------------------------------------------------------------
# Scatter Plot Color Palette Selection
# Available palettes defined in R/09_utils_plotting.R:
# -----------------------------------------------------------------------------
BRCA_PLOT_COLORSET <- "CELL_GOLD" # <-- Change this to select scatter plot color

# Optional: Override individual point settings
# BRCA_POINT_COLOR <- "#2D2D2D"    # Custom point color
# BRCA_POINT_SIZE <- 1.6           # Custom point size

# Plot filters (can be customized)
BRCA_PLOT_FILTERS <- list(
    stratum   = c("ALL"),
    version   = c("RAW", "BatchAdj"),
    score     = c("XCELL_IMMUNE_SCORE"),
    predictor = c("CSN_SCORE", "COPS7A", "COPS7B"),
    method    = c("pearson")
)

# -----------------------------------------------------------------------------
# Heatmap Color Palette Selection
# Available palettes defined in R/09_utils_plotting.R:
# -----------------------------------------------------------------------------
BRCA_HM_COLORSET <- "BLUE_RED_CRIMSON" # <-- Change this to select heatmap color

# =============================================================================
# CORE ANALYSIS FUNCTION
# =============================================================================

#' Run BRCA immune correlation analysis
#'
#' Computes Pearson and Spearman partial correlations between CSN predictors
#' and immune/stromal scores, controlling for covariates.
#'
#' @param dataset_dirs_map Named vector of dataset directories
#' @param make_plots Whether to generate scatter plots
#' @param min_pairs Minimum sample pairs required for correlation
#' @return NULL (invisibly), outputs CSV files
.run_brca_immune_cor <- function(dataset_dirs_map = NULL, make_plots = FALSE, min_pairs = 8L) {
    ds_id <- "brca_cptac_2020"

    # Parse ds_dir
    ds_dir <- NULL
    if (!is.null(dataset_dirs_map) && ds_id %in% names(dataset_dirs_map)) {
        ds_dir <- dataset_dirs_map[[ds_id]]
    }
    if (is.null(ds_dir) && exists("dataset_dirs_run")) {
        dsm <- get("dataset_dirs_run", inherits = TRUE)
        if (ds_id %in% names(dsm)) ds_dir <- dsm[[ds_id]]
    }
    if (is.null(ds_dir) && exists("dataset_dirs")) {
        dsm <- get("dataset_dirs", inherits = TRUE)
        if (ds_id %in% names(dsm)) ds_dir <- dsm[[ds_id]]
    }
    if (is.null(ds_dir) || !dir.exists(ds_dir)) {
        message("[BRCA-immune] Cannot find brca_cptac_2020 directory, skipping")
        return(invisible(NULL))
    }

    out_root <- "BRCA_immune_analysis"
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

    # Plot configuration
    plot_filters <- get0("BRCA_PLOT_FILTERS", ifnotfound = list(
        stratum = NULL, version = NULL, score = NULL, predictor = NULL, method = NULL
    ))
    plot_formats <- toupper(get0("BRCA_PLOT_FORMATS", ifnotfound = c("TIFF")))
    plot_dpi <- get0("BRCA_PLOT_DPI", ifnotfound = 600L)
    plot_w <- get0("BRCA_PLOT_WIDTH", ifnotfound = 4.0)
    plot_h <- get0("BRCA_PLOT_HEIGHT", ifnotfound = 3.0)

    # Color palette
    plot_colorset <- toupper(get0("BRCA_PLOT_COLORSET", ifnotfound = "CELL_BLUE"))
    if (!plot_colorset %in% names(CELL_PALETTES)) plot_colorset <- "CELL_BLUE"

    # Read protein matrix & TP53
    mat0_full <- load_matrix_from_dataset_dir(ds_dir)
    tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))

    # Verify CSN subunits exist
    present_sub <- intersect(csn_subunits, rownames(mat0_full))
    if (!length(present_sub)) {
        message("[BRCA-immune] This dataset lacks CSN subunits, skipping")
        return(invisible(NULL))
    }

    # Read immune scores
    scores_all <- .read_brca_scores_from_patient(ds_dir, colnames(mat0_full))

    # Define strata
    strata_list <- list(
        ALL = colnames(mat0_full),
        TP53_mutant = names(tp53_status)[tp53_status == "TP53_mutant"],
        TP53_wild_type = names(tp53_status)[tp53_status == "TP53_wild_type"]
    )
    versions <- c("RAW", "BatchAdj")

    # ---------------------------------------------------------------------------
    # Main analysis loop: stratum x version
    # ---------------------------------------------------------------------------
    for (st_name in names(strata_list)) {
        keep <- intersect(colnames(mat0_full), strata_list[[st_name]])
        if (length(keep) < 4L) {
            message(sprintf("[BRCA-immune][%s] Too few samples (%d), skipping", st_name, length(keep)))
            next
        }
        mat0 <- mat0_full[, keep, drop = FALSE]

        # Build CSN_SCORE for this stratum
        csn <- build_csn_score_safe(mat0, present_sub, combine_7AB = TRUE)
        csn <- setNames(as.numeric(csn[colnames(mat0)]), colnames(mat0))

        # Outcome scores
        Y <- scores_all[keep, , drop = FALSE]

        # Covariates
        purity <- get_purity_covariate(ds_id, ds_dir, keep)
        sa <- get_sex_age_covariates(ds_dir, keep)
        cov0 <- cbind(sex = sa[, "sex"], age = sa[, "age"], purity = purity)
        rownames(cov0) <- keep
        bi <- get_batch_factor(ds_dir, keep)
        batch0 <- if (!is.null(bi)) droplevels(bi$fac[keep]) else NULL
        if (!is.null(batch0)) names(batch0) <- keep

        for (ver in versions) {
            if (ver == "RAW") {
                batch_use <- NULL
                cov_use <- cov0
            } else {
                batch_use <- batch0
                cov_use <- cov0
            }

            # Build predictors
            predictors <- list()
            for (su in present_sub) {
                v <- as.numeric(mat0[su, ])
                names(v) <- keep
                predictors[[su]] <- v
            }
            predictors[["CSN_SCORE"]] <- csn
            for (su in present_sub) {
                vv <- as.numeric(mat0[su, ])
                names(vv) <- keep
                predictors[[paste0("RESIDUAL_", su)]] <- residualize_vector(
                    y = vv, csn_score = csn, batch = batch_use, covars = cov_use
                )
            }

            rows <- list()
            plot_dir <- file.path(out_root, st_name, ver, "plots")
            if (isTRUE(make_plots)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
            plot_jobs <- list()

            # -----------------------------------------------------------------------
            # Correlation analysis for each score x predictor
            # -----------------------------------------------------------------------
            for (sc in colnames(Y)) {
                y0 <- suppressWarnings(as.numeric(Y[[sc]]))
                names(y0) <- keep

                for (pn in names(predictors)) {
                    x0 <- predictors[[pn]]
                    ok <- is.finite(x0) & is.finite(y0)
                    if (sum(ok) < min_pairs) next

                    # Build model data frame
                    df_model <- data.frame(y = y0[ok], x = x0[ok], check.names = FALSE)
                    if (!is.null(cov_use)) {
                        cov_sub <- as.data.frame(cov_use[ok, , drop = FALSE])
                        df_model <- cbind(df_model, cov_sub)
                    }
                    if (!is.null(batch_use)) {
                        df_model$batch <- droplevels(batch_use[ok])
                    }

                    # PEARSON: Linear model on raw values
                    form_p <- as.formula("y ~ .")
                    fit_p <- try(stats::lm(form_p, data = df_model), silent = TRUE)
                    stats_p <- .extract_partial_stats(fit_p, "x")

                    # For display: residualized plot
                    xy <- .residualize_pair(x0[ok], y0[ok],
                        batch = batch_use[ok],
                        covars = as.data.frame(cov_use[ok, , drop = FALSE])
                    )
                    lm_resid_p <- try(stats::lm(xy$y ~ xy$x), silent = TRUE)
                    intercept_p <- if (!inherits(lm_resid_p, "try-error")) unname(coef(lm_resid_p)[1]) else NA_real_
                    slope_p <- if (!inherits(lm_resid_p, "try-error")) unname(coef(lm_resid_p)[2]) else NA_real_
                    r2_p <- if (is.finite(stats_p$r)) stats_p$r^2 else NA_real_
                    eq_p <- sprintf("y = %.6f + %.6f * x", intercept_p, slope_p)

                    # SPEARMAN: Linear model on ranks
                    df_rank <- df_model
                    df_rank$y <- rank(df_model$y)
                    df_rank$x <- rank(df_model$x)
                    fit_s <- try(stats::lm(form_p, data = df_rank), silent = TRUE)
                    stats_s <- .extract_partial_stats(fit_s, "x")

                    xy_rank <- .residualize_pair(rank(x0[ok]), rank(y0[ok]),
                        batch = batch_use[ok],
                        covars = as.data.frame(cov_use[ok, , drop = FALSE])
                    )
                    lm_resid_s <- try(stats::lm(xy_rank$y ~ xy_rank$x), silent = TRUE)
                    intercept_s <- if (!inherits(lm_resid_s, "try-error")) unname(coef(lm_resid_s)[1]) else NA_real_
                    slope_s <- if (!inherits(lm_resid_s, "try-error")) unname(coef(lm_resid_s)[2]) else NA_real_
                    r2_s <- if (is.finite(stats_s$r)) stats_s$r^2 else NA_real_
                    eq_s <- sprintf("rank(y) = %.6f + %.6f * rank(x)", intercept_s, slope_s)

                    # Store results
                    rows[[length(rows) + 1L]] <- data.frame(
                        dataset = ds_id, stratum = st_name, version = ver,
                        score = sc, predictor = pn,
                        method = "pearson", n = sum(ok),
                        r = as.numeric(stats_p$r), p = as.numeric(stats_p$p),
                        padj = NA_real_, R2 = r2_p, equation = eq_p,
                        stringsAsFactors = FALSE, check.names = FALSE
                    )
                    rows[[length(rows) + 1L]] <- data.frame(
                        dataset = ds_id, stratum = st_name, version = ver,
                        score = sc, predictor = pn,
                        method = "spearman", n = sum(ok),
                        r = as.numeric(stats_s$r), p = as.numeric(stats_s$p),
                        padj = NA_real_, R2 = r2_s, equation = eq_s,
                        stringsAsFactors = FALSE, check.names = FALSE
                    )

                    # Collect plot jobs if needed
                    if (isTRUE(make_plots)) {
                        want_stratum <- is.null(plot_filters$stratum) || st_name %in% plot_filters$stratum
                        want_version <- is.null(plot_filters$version) || ver %in% plot_filters$version
                        want_score <- is.null(plot_filters$score) || sc %in% plot_filters$score
                        want_predictor <- is.null(plot_filters$predictor) || pn %in% plot_filters$predictor
                        if (want_stratum && want_version && want_score && want_predictor) {
                            dfp <- data.frame(x = xy$x, y = xy$y)
                            methods_want <- c("pearson", "spearman")
                            if (!is.null(plot_filters$method)) {
                                methods_want <- intersect(methods_want, tolower(plot_filters$method))
                            }
                            for (mm in methods_want) {
                                rr <- if (mm == "pearson") as.numeric(stats_p$r) else as.numeric(stats_s$r)
                                pp <- if (mm == "pearson") as.numeric(stats_p$p) else as.numeric(stats_s$p)
                                r2_val <- if (mm == "pearson") r2_p else r2_s
                                eq_val <- if (mm == "pearson") eq_p else eq_s
                                plot_jobs[[length(plot_jobs) + 1L]] <- list(
                                    dataset = ds_id, stratum = st_name, version = ver,
                                    score = sc, predictor = pn, method = mm,
                                    r = rr, p = pp, R2 = r2_val, equation = eq_val,
                                    df = dfp, n = sum(ok)
                                )
                            }
                        }
                    }
                }
            }

            # Output results
            out_dir <- file.path(out_root, st_name, ver)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

            if (length(rows)) {
                out_tbl <- do.call(rbind, rows)
                # BH adjustment by group
                out_tbl <- split(out_tbl, list(
                    out_tbl$dataset, out_tbl$stratum, out_tbl$version,
                    out_tbl$method, out_tbl$score
                ), drop = TRUE)
                out_tbl <- lapply(out_tbl, function(d) {
                    d$padj <- stats::p.adjust(d$p, method = "BH")
                    d
                })

                # Generate plots after BH adjustment
                if (isTRUE(make_plots) && length(plot_jobs)) {
                    ot <- do.call(rbind, out_tbl)
                    pal <- CELL_PALETTES[[plot_colorset]]
                    ln_col <- pal$ln
                    pt_col <- get0("BRCA_POINT_COLOR", ifnotfound = pal$pt)
                    pt_size <- as.numeric(get0("BRCA_POINT_SIZE", ifnotfound = 1.6))
                    ci_col <- .ci_fill(ln_col)

                    for (job in plot_jobs) {
                        if (!is.null(plot_filters$method) &&
                            !(tolower(job$method) %in% tolower(plot_filters$method))) {
                            next
                        }

                        idx <- which(ot$dataset == job$dataset & ot$stratum == job$stratum &
                            ot$version == job$version & ot$method == job$method &
                            ot$score == job$score & ot$predictor == job$predictor)
                        padj_txt <- "NA"
                        if (length(idx)) {
                            padj_val <- ot$padj[idx[1]]
                            if (is.finite(padj_val)) padj_txt <- format(padj_val, digits = 3, scientific = TRUE)
                        }

                        title_txt <- sprintf("%s vs %s", job$score, job$predictor)
                        sub_txt <- sprintf(
                            "%s | %s / %s / %s | r=%.3f, padj=%s, n=%d",
                            toupper(job$method), job$dataset, job$stratum, job$version,
                            job$r, padj_txt, job$n
                        )

                        p1 <- ggplot2::ggplot(job$df, ggplot2::aes(x = x, y = y)) +
                            ggplot2::geom_point(size = pt_size, alpha = 0.85, colour = pt_col) +
                            ggplot2::geom_smooth(
                                method = "lm", se = TRUE, linewidth = 0.8,
                                colour = ln_col, fill = ci_col
                            ) +
                            ggplot2::labs(
                                title = title_txt, subtitle = sub_txt,
                                x = sprintf("%s (residualized)", job$predictor),
                                y = sprintf("%s (residualized)", job$score)
                            ) +
                            .cell_theme(base_size = 10)

                        for (fmt in plot_formats) {
                            ext <- switch(toupper(fmt),
                                "TIFF" = "tiff",
                                "TIF" = "tiff",
                                "PNG" = "png",
                                "JPG" = "jpg",
                                "JPEG" = "jpg",
                                "PDF" = "pdf",
                                tolower(fmt)
                            )
                            fn <- sprintf(
                                "%s_vs_%s_%s_%s_%s.%s",
                                job$score, job$predictor, job$stratum, job$version, job$method, ext
                            )
                            ggplot2::ggsave(
                                filename = file.path(plot_dir, fn), plot = p1,
                                width = plot_w, height = plot_h, dpi = plot_dpi,
                                bg = "white", device = ext, limitsize = FALSE
                            )
                        }
                    }
                }

                out_tbl <- do.call(rbind, out_tbl)
                data.table::fwrite(out_tbl, file.path(out_dir, "immune_vs_predictor_correlations.csv"))
                message(sprintf(
                    "[BRCA-immune][%s/%s] Output %d results -> %s",
                    st_name, ver, nrow(out_tbl),
                    file.path(out_dir, "immune_vs_predictor_correlations.csv")
                ))
            } else {
                header <- data.frame(
                    dataset = character(0), stratum = character(0), version = character(0),
                    score = character(0), predictor = character(0), method = character(0),
                    correlation = numeric(0), p = numeric(0), padj = numeric(0),
                    R2 = numeric(0), equation = character(0),
                    stringsAsFactors = FALSE, check.names = FALSE
                )
                data.table::fwrite(header, file.path(out_dir, "immune_vs_predictor_correlations.csv"))
                message(sprintf("[BRCA-immune][%s/%s] No qualified pairs (rows=0); empty file output.", st_name, ver))
            }
        }

        # -------------------------------------------------------------------------
        # TP53 interaction analysis (ALL stratum only)
        # -------------------------------------------------------------------------
        if (identical(st_name, "ALL")) {
            tp53_num <- as.numeric(tp53_status[keep] == "TP53_mutant")
            names(tp53_num) <- keep

            for (ver in versions) {
                if (ver == "RAW") {
                    batch_use <- NULL
                    cov_use <- cov0
                } else {
                    batch_use <- batch0
                    cov_use <- cov0
                }
                rows_i <- list()

                for (sc in colnames(Y)) {
                    y0 <- suppressWarnings(as.numeric(Y[[sc]]))
                    names(y0) <- keep

                    for (pn in names(predictors)) {
                        x0 <- predictors[[pn]]
                        names(x0) <- keep
                        ok <- is.finite(x0) & is.finite(y0) & is.finite(tp53_num)
                        if (sum(ok) < min_pairs) next

                        dfm <- data.frame(
                            y = y0[ok], x = x0[ok], tp53 = tp53_num[ok],
                            cov_use[ok, , drop = FALSE]
                        )
                        if (!is.null(batch_use)) dfm$batch <- droplevels(batch_use[ok])

                        form <- as.formula(paste0(
                            "y ~ x * tp53 + ", paste(colnames(cov_use), collapse = " + "),
                            if (!is.null(batch_use)) " + batch" else ""
                        ))
                        fit <- try(stats::lm(form, data = dfm), silent = TRUE)
                        if (inherits(fit, "try-error")) next

                        sm <- summary(fit)
                        co <- coef(sm)
                        beta_int <- if ("x:tp53" %in% rownames(co)) co["x:tp53", "Estimate"] else NA_real_
                        p_int <- if ("x:tp53" %in% rownames(co)) co["x:tp53", "Pr(>|t|)"] else NA_real_
                        r2m <- as.numeric(sm$r.squared)
                        eq_line <- paste0(
                            "y = ", sprintf("%.6f", unname(coef(fit)[1])),
                            " + ", sprintf("%.6f", unname(coef(fit)["x"])), " * x",
                            " + ", sprintf("%.6f", unname(coef(fit)["tp53"])), " * tp53",
                            " + ", sprintf("%.6f", beta_int), " * x:tp53 + ..."
                        )
                        rows_i[[length(rows_i) + 1L]] <- data.frame(
                            dataset = ds_id, stratum = "ALL", version = ver,
                            score = sc, predictor = pn,
                            method = "lm_interaction",
                            beta_interaction = beta_int, p = p_int, padj = NA_real_,
                            R2 = r2m, equation = eq_line,
                            stringsAsFactors = FALSE, check.names = FALSE
                        )
                    }
                }

                out_dir <- file.path(out_root, "ALL", ver)
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

                if (length(rows_i)) {
                    out_i <- do.call(rbind, rows_i)
                    out_i$padj <- ave(out_i$p, out_i$dataset, out_i$stratum, out_i$version, out_i$score,
                        FUN = function(v) stats::p.adjust(v, method = "BH")
                    )
                    data.table::fwrite(out_i, file.path(out_dir, "immune_vs_predictor_interaction.csv"))
                    message(sprintf(
                        "[BRCA-immune][ALL/%s] Interaction: output %d rows -> %s",
                        ver, nrow(out_i), file.path(out_dir, "immune_vs_predictor_interaction.csv")
                    ))
                } else {
                    header_i <- data.frame(
                        dataset = character(0), stratum = character(0), version = character(0),
                        score = character(0), predictor = character(0), method = character(0),
                        beta_interaction = numeric(0), p = numeric(0), padj = numeric(0),
                        R2 = numeric(0), equation = character(0),
                        stringsAsFactors = FALSE, check.names = FALSE
                    )
                    data.table::fwrite(header_i, file.path(out_dir, "immune_vs_predictor_interaction.csv"))
                    message(sprintf("[BRCA-immune][ALL/%s] Interaction: no qualified results, empty file output.", ver))
                }
            }
        }
    }

    message("[BRCA-immune] Complete: output -> ", out_root)
}

# =============================================================================
# EXECUTE ANALYSIS
# =============================================================================

# Run correlation analysis (set make_plots = TRUE to generate scatter plots)
tryCatch(
    {
        .run_brca_immune_cor(dataset_dirs_map = dataset_dirs_run, make_plots = TRUE)
    },
    error = function(e) {
        calls <- sys.calls()
        message("[BRCA-immune] Call stack: ", paste(utils::capture.output(print(calls)), collapse = " | "))
        stop(e)
    }
)

# =============================================================================
# GENERATE HEATMAP
# =============================================================================

# Generate immune correlation heatmap
.plot_brca_immune_heatmap(
    stratum = "ALL",
    version = "BatchAdj",
    method = "pearson",
    hm_colorset = "BLUE_RED_CRIMSON",
    out_formats = c("TIFF")
)
