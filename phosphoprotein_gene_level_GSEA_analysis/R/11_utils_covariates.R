# =============================================================================
# 11_utils_covariates.R - Covariate Handling (Complete Implementation)
# =============================================================================
# This file contains functions for handling covariates (purity, sex, age)
# with dataset-specific logic for CPTAC datasets.
# =============================================================================

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Safe z-score: estimate mean/sd only from finite values, keep NA, no mean imputation
.z_no_impute <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    fin <- is.finite(v)
    mu <- mean(v[fin], na.rm = TRUE)
    sdv <- stats::sd(v[fin], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out <- (v - mu) / sdv
    out[!fin] <- NA_real_
    out
}

# Column name normalization
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))

# Z-score with mean imputation
.zscore <- function(v) {
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    (v - mu) / sdv
}

# Convert to 0-1 range (auto-detect percentage vs fraction)
.to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
    pmin(pmax(v, 0), 1)
}

# For PAAD: get median from semicolon-separated values
.median_from_semicolon <- function(x_chr) {
    vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
    vv <- vv[is.finite(vv)]
    if (!length(vv)) {
        return(NA_real_)
    }
    stats::median(vv)
}

# Format batch level sample counts as string
.format_batch_sizes <- function(fac) {
    if (is.null(fac)) {
        return(NA_character_)
    }
    tb <- sort(table(fac), decreasing = TRUE)
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}

# -----------------------------------------------------------------------------
# Get Purity by Dataset (Named Numeric, 0~1)
# Dataset-specific logic for CPTAC datasets
# -----------------------------------------------------------------------------
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    samp <- if (file.exists(samp_fp)) suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) else NULL
    if (!is.null(samp)) {
        samp <- as.data.frame(samp)
        names(samp) <- .norm_names(names(samp))
    }
    purity <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

    if (ds_id == "brca_cptac_2020") {
        pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
        if (file.exists(pat_fp) && !is.null(samp)) {
            pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
            names(pat) <- .norm_names(names(pat))
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            pid_in_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
            pid_in_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
            if (!is.na(sid) && !is.na(pid_in_samp) && !is.na(pid_in_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat)) {
                map_pt <- setNames(as.character(samp[[pid_in_samp]]), samp[[sid]])
                pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_in_pat]]))
                purity[] <- unname(pt_pur[map_pt[sample_ids]])
            }
        }
    } else if (ds_id == "luad_cptac_2020") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp)) {
                purity[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "lusc_cptac_2021") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp)) {
                purity[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "paad_cptac_2021") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
                medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
                purity[] <- .to01(medv[match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "ucec_cptac_2020") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid)) {
                pc <- suppressWarnings(as.numeric(samp[["PURITY_CANCER"]]))
                pi <- suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]]))
                ps <- suppressWarnings(as.numeric(samp[["PURITY_STROMA"]]))
                idx <- match(sample_ids, samp[[sid]])
                purity_calc <- ifelse(is.finite(pc[idx]), pc[idx], 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0))
                purity[] <- pmin(pmax(purity_calc, 0), 1)
            }
        }
    }
    purity
}

# -----------------------------------------------------------------------------
# Get Sex/Age Covariates (Sample Aligned)
# Sex: 0/1 (Female/Male); Age: z-score
# -----------------------------------------------------------------------------
get_sex_age_covariates <- function(ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")

    if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
        log_msg("[covars] Missing sample/patient file, using NA for sex/age")
        out <- cbind(
            sex = rep(NA_real_, length(sample_ids)),
            age = rep(NA_real_, length(sample_ids))
        )
        rownames(out) <- sample_ids
        return(out)
    }

    .NN <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))

    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .NN(names(samp))
    names(pat) <- .NN(names(pat))

    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
    if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) {
        stop("[get_sex_age_covariates] Cannot establish sample-patient mapping (missing SAMPLE_ID/PATIENT_ID)")
    }
    map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

    ## SEX: Male=1, Female=0
    sex_col <- intersect(c("SEX", "GENDER"), names(pat))[1]
    sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    if (!is.na(sex_col)) {
        raw <- toupper(as.character(pat[[sex_col]]))
        val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
        names(val) <- as.character(pat[[pid_pat]])
        sex[] <- unname(val[map_pt[sample_ids]])
    }

    ## AGE: main analysis does not impute; sensitivity can use missing-indicator + z after mean imputation
    age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
    age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    age_missing <- age_z_imputed <- NULL
    if (!is.na(age_col)) {
        v <- suppressWarnings(as.numeric(pat[[age_col]]))
        names(v) <- as.character(pat[[pid_pat]])
        age_raw <- unname(v[map_pt[sample_ids]])
        age[] <- .z_no_impute(age_raw)
        if (!exists("USE_AGE_MISSING_INDICATOR", inherits = FALSE)) USE_AGE_MISSING_INDICATOR <- FALSE
        if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
            age_missing <- as.numeric(is.na(age_raw))
            age_z_imputed <- .zscore(age_raw)
        }
    }

    cov_sex <- mean(is.finite(sex)) * 100
    cov_age <- mean(is.finite(age)) * 100
    if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
        cov_age_imp <- mean(is.finite(age_z_imputed)) * 100
        cov_age_mis <- mean(is.finite(age_missing)) * 100
        log_msg(
            "    covariates coverage: sex %.1f%%, age(NA-as-NA) %.1f%%, age_z_imputed %.1f%%, age_missing %.1f%%",
            cov_sex, cov_age, cov_age_imp, cov_age_mis
        )
    } else {
        log_msg("    covariates coverage: sex %.1f%%, age(NA-as-NA) %.1f%%", cov_sex, cov_age)
    }

    out <- data.frame(
        sex = as.numeric(sex), age = as.numeric(age),
        row.names = sample_ids, check.names = FALSE
    )
    out$age_z <- out$age
    if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
        out$age_missing <- as.numeric(age_missing)
        out$age_z_imputed <- as.numeric(age_z_imputed)
    }
    out
}

# -----------------------------------------------------------------------------
# Select Covariates Safely (With Coverage Checks, Correlation Filtering, TP53 Exclusion)
# -----------------------------------------------------------------------------
select_covars_safely <- function(
  df, # candidate covariates; rownames=sample_order
  sample_order,
  label = "covars",
  y = NULL, # continuous predictor; can be NULL
  min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
  max_abs_cor = 0.30,
  min_pairs = 20L
) {
    logf <- function(...) if (exists("log_msg", mode = "function")) try(log_msg(...), silent = TRUE)
    ## column-by-column drop reason collector
    drop_log <- list()
    .append_drop <- function(step, col, reason, metric = NA_real_, threshold = NA_real_, extra = NA_character_) {
        drop_log[[length(drop_log) + 1L]] <<- data.frame(
            step = step,
            covariate = col,
            reason = reason,
            metric = metric,
            threshold = threshold,
            extra = extra,
            stringsAsFactors = FALSE
        )
    }
    if (is.null(df) || !nrow(df)) {
        return(NULL)
    }


    # 1) Align sample order
    if (is.null(rownames(df))) {
        logf("  [covars-%s] df has no rownames, skip", label)
        return(NULL)
    }
    so <- as.character(sample_order)
    logf("  [covars-%s] before align: C_dim=%d x %d; has_rownames=%s", label, NROW(df), NCOL(df), !is.null(rownames(df)))
    df <- df[so, , drop = FALSE]

    ## [POLICY] Never include TP53 as a covariate (ALL / MT / WT strata)
    tp_cols_idx <- grep("^TP53($|_|)", colnames(df), ignore.case = FALSE)
    if (length(tp_cols_idx) > 0L) {
        tp_cols <- colnames(df)[tp_cols_idx]
        for (cc in tp_cols) {
            .append_drop(
                step = "policy", col = cc, reason = "exclude_TP53_as_covariate",
                metric = NA_real_, threshold = NA_real_, extra = "global policy"
            )
        }
        df <- df[, setdiff(colnames(df), tp_cols), drop = FALSE]
        logf(sprintf("  [covars-%s] policy: drop TP53 columns from covariates -> %s", label, paste(tp_cols, collapse = ",")))
    }

    logf("  [covars-%s] after  align: C_dim=%d x %d", label, NROW(df), NCOL(df))


    # 2) Always use data.frame; preserve factors
    df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)

    # 3) coverage: use global minimum threshold for unnamed columns (safe indexing, do not use [[cn]])
    global_min <- if (length(min_cov_named)) min(min_cov_named, na.rm = TRUE) else 0
    keep_cov <- rep(TRUE, ncol(df))
    names(keep_cov) <- colnames(df)

    for (cn in colnames(df)) {
        v <- df[[cn]]
        # Unified coverage definition: use is.finite for numeric; use !is.na for non-numeric
        cover <- if (is.numeric(v)) mean(is.finite(v)) else mean(!is.na(v))

        thr <- if (!is.null(names(min_cov_named)) && (cn %in% names(min_cov_named))) {
            min_cov_named[[cn]] # Safe: only use [[cn]] when name exists
        } else {
            global_min
        }

        if (is.na(cover) || cover < thr) {
            logf(
                "  [covars-%s] drop %s (coverage=%s < %s)", label, cn,
                sprintf("%.0f%%", 100 * ifelse(is.na(cover), 0, cover)),
                sprintf("%.0f%%", 100 * thr)
            )
            keep_cov[cn] <- FALSE
            .append_drop("coverage", cn, "low_coverage",
                metric = ifelse(is.na(cover), 0, cover), threshold = thr
            )
        }
    }
    df <- df[, keep_cov, drop = FALSE]
    if (!ncol(df)) {
        return(NULL)
    }

    
    keep_cov <- rep(TRUE, ncol(df))
    names(keep_cov) <- colnames(df)
    

    # 4) Correlation filtering with biology (only for numeric columns)
    if (!is.null(y)) {
        y <- suppressWarnings(as.numeric(y))
        ## Do not apply rho gate to these covariates
        
        skip_rho_gate <- c("purity", "age", "sex", "CSN_SCORE")

        for (cn in colnames(df)) {
            v <- df[[cn]]
            ## Only apply |rho| filtering when numeric and not in skip list
            if (is.numeric(v) && !(cn %in% skip_rho_gate)) {
                fin <- is.finite(v) & is.finite(y)
                if (sum(fin) >= min_pairs) {
                    r <- suppressWarnings(stats::cor(v[fin], y[fin], method = "spearman"))
                    if (is.finite(r) && abs(r) >= max_abs_cor) {
                        logf("  [covars-%s] drop %s (|rho|=%.3f >= %.2f vs biology)", label, cn, abs(r), max_abs_cor)
                        keep_cov[cn] <- FALSE
                    }
                }
            }
        }
        df <- df[, keep_cov, drop = FALSE]
        if (!ncol(df)) {
            return(NULL)
        }
    }


    # 5) Return: row order = sample_order; for model.matrix expansion
    stopifnot(nrow(df) == length(so), identical(rownames(df), so))
    attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
    df
}

# -----------------------------------------------------------------------------
# Audit Covariates Coverage (Log Coverage Stats)
# -----------------------------------------------------------------------------
audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL) {
    dir.create(file.path("run_info", "covars_audit"), recursive = TRUE, showWarnings = FALSE)

    cov_df <- if (is.null(covars)) NULL else as.data.frame(covars, check.names = FALSE)
    cov_nms <- if (!is.null(cov_df)) colnames(cov_df) else character(0)

    get_cov <- function(k) {
        if (is.null(cov_df) || !(k %in% names(cov_df))) {
            return(NA_real_)
        }
        v <- cov_df[[k]]
        if (is.numeric(v)) {
            mean(is.finite(v)) * 100
        } else {
            mean(!is.na(v)) * 100 # factor/string -> calculate non-NA ratio
        }
    }
    cov_purity <- get_cov("purity")
    cov_sex <- get_cov("sex")
    cov_age <- get_cov("age")

    line <- sprintf(
        "  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s",
        tag,
        if (length(cov_nms)) paste(cov_nms, collapse = ",") else "NULL",
        cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
        if (!is.null(batch)) nlevels(batch) else 0L
    )
    log_msg(line)

    # Append one row to summary table
    df <- data.frame(
        dataset = ds_id, stratum = stratum, subunit = su, tag = tag,
        batch_levels = if (!is.null(batch)) nlevels(batch) else 0L,
        batch_sizes = if (!is.null(batch)) {
            paste(sprintf(
                "%s=%d", names(sort(table(batch), decreasing = TRUE)),
                as.integer(sort(table(batch), decreasing = TRUE))
            ), collapse = "; ")
        } else {
            NA_character_
        },
        covars_cols = if (length(cov_nms)) paste(cov_nms, collapse = ";") else "NULL",
        purity_cov = cov_purity, sex_cov = cov_sex, age_cov = cov_age,
        stringsAsFactors = FALSE
    )
    fp <- file.path("run_info", "covars_audit", "audit_rows.csv")
    data.table::fwrite(df, fp, append = file.exists(fp))
}
