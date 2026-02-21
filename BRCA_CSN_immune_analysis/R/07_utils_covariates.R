# =============================================================================
# 07_utils_covariates.R
# Covariate Extraction Utilities (Purity, Sex, Age)
# =============================================================================

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Normalize column names
#'
#' Converts to uppercase and replaces non-alphanumeric characters with underscore.
#'
#' @param x Character vector of names
#' @return Normalized names
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))

#' Z-score normalization with mean imputation
#'
#' Computes z-scores, imputing non-finite values with the mean.
#'
#' @param v Numeric vector
#' @return Z-scored vector
.zscore <- function(v) {
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    (v - mu) / sdv
}

#' Z-score without imputation
#'
#' Computes z-scores but keeps NA values as NA.
#'
#' @param v Numeric vector
#' @return Z-scored vector with NA preserved
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

#' Convert values to 0-1 range
#'
#' Handles percentage values (>100) by dividing by 100.
#'
#' @param v Numeric vector
#' @return Values clamped to 0-1
.to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
    pmin(pmax(v, 0), 1)
}

#' Extract median from semicolon-separated values
#'
#' Used for PAAD dataset where cellularity is semicolon-separated.
#'
#' @param x_chr Character string with semicolon-separated values
#' @return Median value or NA
.median_from_semicolon <- function(x_chr) {
    vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
    vv <- vv[is.finite(vv)]
    if (!length(vv)) {
        return(NA_real_)
    }
    stats::median(vv)
}

# -----------------------------------------------------------------------------
# Purity Covariate
# -----------------------------------------------------------------------------

#' Get tumor purity covariate by dataset
#'
#' Extracts purity values from clinical data, handling dataset-specific
#' column names and formats.
#'
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param sample_ids Character vector of sample IDs
#' @return Named numeric vector of purity values (0-1)
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
# Sex and Age Covariates
# -----------------------------------------------------------------------------

#' Get sex and age covariates
#'
#' Extracts sex (0=Female, 1=Male) and age (z-scored) from patient data.
#'
#' @param ds_dir Dataset directory path
#' @param sample_ids Character vector of sample IDs
#' @return Data frame with sex, age, and optional age_missing/age_z_imputed
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
        stop("[get_sex_age_covariates] Cannot build sample<->patient mapping (missing SAMPLE_ID/PATIENT_ID)")
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

    ## AGE: main analysis no imputation
    age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX", "AGE_YEARS"), names(pat))[1]
    age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    age_missing <- age_z_imputed <- NULL
    if (!is.na(age_col)) {
        v <- suppressWarnings(as.numeric(pat[[age_col]]))
        names(v) <- as.character(pat[[pid_pat]])
        age_raw <- unname(v[map_pt[sample_ids]])
        age[] <- .z_no_impute(age_raw) # NA preserved (main analysis)
        if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
            age_missing <- as.numeric(is.na(age_raw))
            age_z_imputed <- .zscore(age_raw) # mean-imputed z (sensitivity)
        }
    }

    cov_sex <- mean(is.finite(sex)) * 100
    cov_age <- mean(is.finite(age)) * 100
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
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
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
        out$age_missing <- as.numeric(age_missing)
        out$age_z_imputed <- as.numeric(age_z_imputed)
    }
    out
}
