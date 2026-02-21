## ============================================================================
## Stratum-Level Analysis Workflow (Modularized)
## Manages complete analysis for one stratum (ALL/TP53_mutant/TP53_wild_type)
## ============================================================================

#' Prepare Stratum Data
#'
#' Initial data preparation: sample subsetting, log transformation, filtering.
#'
#' @param mat0_full Full expression matrix (all samples)
#' @param sample_keep Sample IDs to keep for this stratum
#' @param min_frac_complete Minimum fraction of non-NA values per gene
#'
#' @return List with $mat0 (raw matrix), $mat (imputed/filtered), $present_sub (CSN subunits)
#'
#' @keywords internal
.prepare_stratum_data <- function(mat0_full, sample_keep, min_frac_complete) {
    # Subset samples
    keep <- intersect(colnames(mat0_full), sample_keep)

    if (length(keep) < 4) {
        log_msg("  [Skip] Samples too few: %d", length(keep))
        return(NULL)
    }

    mat0 <- mat0_full[, keep, drop = FALSE]

    # Log transformation if needed
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) {
        mat0 <- log2(mat0 + 1)
    }

    # Imputation and filtering for limma
    mat <- impute_and_filter(mat0, min_frac = min_frac_complete)

    # Check CSN subunits
    present_sub <- intersect(csn_subunits, rownames(mat0))

    if (!length(present_sub)) {
        log_msg("  [Skip] No CSN subunit in this stratum")
        return(NULL)
    }

    list(mat0 = mat0, mat = mat, present_sub = present_sub)
}

#' Calculate CSN Score and Prepare Covariates
#'
#' Computes CSN complex score and gathers all covariate information.
#'
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory
#' @param mat0 Raw expression matrix for this stratum
#' @param present_sub CSN subunits present in data
#' @param stratum Stratum name
#'
#' @return List with $csn_score, $batch_all, $purity_all, $sa_all, $tp53_num_all, $is_ALL
#'
#' @keywords internal
.calculate_predictors_and_covariates <- function(ds_id, ds_dir, mat0, present_sub, stratum) {
    # Audit CSN score feasibility
    audit_csn_score_feasibility_safe(
        ds_id = ds_id,
        stratum = stratum,
        mat0 = mat0,
        present_sub = present_sub,
        min_members = 5L,
        pca_min_samples = 10L,
        min_per_group = min_per_group,
        out_dir = file.path("run_info", "csn_score_audit")
    )

    # Calculate CSN score
    csn_score <- build_csn_score_safe(
        mat0,
        subunits = present_sub,
        combine_7AB = TRUE,
        min_members = 5L,
        pca_min_samples = 10L
    )

    # Get batch information
    bi_all <- get_batch_factor(ds_dir, colnames(mat0))
    batch_all <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL

    # Get covariates
    purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
    sa_all <- get_sex_age_covariates(ds_dir, colnames(mat0))

    # TP53 status (only for ALL stratum)
    is_ALL <- identical(stratum, "ALL")
    tp53_num_all <- NULL

    if (is_ALL) {
        tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
        tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
        names(tp53_num_all) <- colnames(mat0)
    }

    list(
        csn_score = csn_score,
        batch_all = batch_all,
        purity_all = purity_all,
        sa_all = sa_all,
        tp53_num_all = tp53_num_all,
        is_ALL = is_ALL
    )
}

#' Get Purity Covariate
#'
#' Extracts tumor purity values from clinical data.
#'
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory
#' @param sample_ids Sample IDs to extract
#'
#' @return Named numeric vector of purity values, or NULL
#'
#' @keywords internal
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
    # Try to read clinical patient file for purity
    fp <- file.path(ds_dir, "data_clinical_patient.txt")

    if (!file.exists(fp)) {
        return(NULL)
    }

    clin <- tryCatch(
        readr::read_tsv(fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )

    if (is.null(clin)) {
        return(NULL)
    }

    # Look for purity column
    purity_cols <- c(
        "TUMOR_PURITY", "Tumor_Purity", "tumor_purity",
        "PURITY", "Purity", "purity"
    )
    hit <- intersect(purity_cols, colnames(clin))

    if (!length(hit)) {
        return(NULL)
    }

    # Match to samples (patient ID might differ from sample ID)
    patient_id_cols <- c("PATIENT_ID", "Patient_ID", "patient_id", "PATIENT", "Patient")
    pid_col <- intersect(patient_id_cols, colnames(clin))

    if (!length(pid_col)) {
        return(NULL)
    }

    clin_df <- data.frame(
        patient_id = as.character(clin[[pid_col[1]]]),
        purity = suppressWarnings(as.numeric(clin[[hit[1]]])),
        stringsAsFactors = FALSE
    )

    # Try to match patient ID to sample ID (often sample ID contains patient ID)
    purity_vec <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

    for (i in seq_along(sample_ids)) {
        sid <- sample_ids[i]
        # Simple matching: exact match or patient ID is substring of sample ID
        matches <- clin_df$patient_id[clin_df$patient_id == sid | grepl(clin_df$patient_id, sid, fixed = TRUE)]
        if (length(matches) > 0) {
            purity_vec[sid] <- clin_df$purity[clin_df$patient_id == matches[1]]
        }
    }

    purity_vec
}

#' Get Sex and Age Covariates
#'
#' Extracts sex and age information from clinical data, with optional
#' age missing indicator.
#'
#' @param ds_dir Dataset directory
#' @param sample_ids Sample IDs to extract
#'
#' @return Data frame with columns: sex, age, [age_missing, age_z_imputed]
#'
#' @keywords internal
get_sex_age_covariates <- function(ds_dir, sample_ids) {
    fp <- file.path(ds_dir, "data_clinical_sample.txt")

    if (!file.exists(fp)) {
        return(data.frame(
            sex = rep(NA_character_, length(sample_ids)),
            age = rep(NA_real_, length(sample_ids)),
            row.names = sample_ids,
            stringsAsFactors = FALSE
        ))
    }

    clin <- tryCatch(
        readr::read_tsv(fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )

    if (is.null(clin)) {
        return(data.frame(
            sex = rep(NA_character_, length(sample_ids)),
            age = rep(NA_real_, length(sample_ids)),
            row.names = sample_ids,
            stringsAsFactors = FALSE
        ))
    }

    # Find sample ID column
    id_cols <- c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample")
    id_hit <- intersect(id_cols, colnames(clin))

    if (!length(id_hit)) {
        return(data.frame(
            sex = rep(NA_character_, length(sample_ids)),
            age = rep(NA_real_, length(sample_ids)),
            row.names = sample_ids,
            stringsAsFactors = FALSE
        ))
    }

    clin$SAMPLE_ID <- as.character(clin[[id_hit[1]]])

    # Find sex column
    sex_cols <- c("SEX", "Sex", "sex", "GENDER", "Gender", "gender")
    sex_hit <- intersect(sex_cols, colnames(clin))
    sex_values <- if (length(sex_hit)) as.character(clin[[sex_hit[1]]]) else rep(NA_character_, nrow(clin))

    # Find age column
    age_cols <- c("AGE", "Age", "age", "AGE_AT_DIAGNOSIS", "Age_at_Diagnosis")
    age_hit <- intersect(age_cols, colnames(clin))
    age_values <- if (length(age_hit)) suppressWarnings(as.numeric(clin[[age_hit[1]]])) else rep(NA_real_, nrow(clin))

    # Align to sample_ids
    clin_aligned <- data.frame(
        SAMPLE_ID = clin$SAMPLE_ID,
        sex = sex_values,
        age = age_values,
        stringsAsFactors = FALSE
    )

    clin_aligned <- clin_aligned[match(sample_ids, clin_aligned$SAMPLE_ID), ]
    rownames(clin_aligned) <- sample_ids

    # Age missing indicator (optional)
    USE_AGE_MISSING_INDICATOR <- isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))

    if (USE_AGE_MISSING_INDICATOR) {
        clin_aligned$age_missing <- as.numeric(is.na(clin_aligned$age))

        # Z-score imputed age (use mean age where NA)
        age_mean <- mean(clin_aligned$age, na.rm = TRUE)
        age_sd <- sd(clin_aligned$age, na.rm = TRUE)
        if (!is.finite(age_mean)) age_mean <- 0
        if (!is.finite(age_sd) || age_sd == 0) age_sd <- 1

        age_imputed <- clin_aligned$age
        age_imputed[is.na(age_imputed)] <- age_mean
        clin_aligned$age_z_imputed <- (age_imputed - age_mean) / age_sd
    }

    clin_aligned[, setdiff(names(clin_aligned), "SAMPLE_ID"), drop = FALSE]
}

#' Run All DEG Analyses for Stratum
#'
#' Executes DEG analysis for all predictors: individual subunits, CSN_SCORE,
#' and interaction models.
#'
#' @param ds_id Dataset ID
#' @param ds_dir Dataset directory
#' @param mat0 Raw expression matrix
#' @param out_root Output root directory
#' @param present_sub CSN subunits present
#' @param csn_score CSN complex score vector
#' @param batch_all Batch factor
#' @param purity_all Purity vector
#' @param sa_all Sample annotation data frame
#'
#' @return Invisible NULL
#'
#' @keywords internal
.run_all_deg_analyses <- function(ds_id, ds_dir, mat0, out_root, present_sub,
                                  csn_score, batch_all, purity_all, sa_all) {
    log_msg("[DEG] Starting all predictor analyses")

    # 1) Individual subunits
    for (su in present_sub) {
        run_deg_for_predictor(
            predictor_name = su,
            predictor_vec = mat0[su, ],
            ds_id = ds_id,
            ds_dir = ds_dir,
            mat0 = mat0,
            out_root = out_root,
            batch_all = batch_all,
            purity_all = purity_all,
            sa_all = sa_all
        )
    }

    # 2) CSN complex score
    if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
        run_deg_for_predictor(
            predictor_name = "CSN_SCORE",
            predictor_vec = csn_score,
            ds_id = ds_id,
            ds_dir = ds_dir,
            mat0 = mat0,
            out_root = out_root,
            batch_all = batch_all,
            purity_all = purity_all,
            sa_all = sa_all
        )
    } else {
        log_msg("  [DEG-CSN_SCORE] Insufficient non-NA samples, skip")
    }

    # 3) Interaction models (subunit controlling for CSN)
    min_n_resid <- opt("min_n_resid", 10L)
    csn_nonNA <- sum(is.finite(csn_score))

    if (csn_nonNA < min_n_resid) {
        log_msg(
            "  [DEG-INTERACTION] CSN score non-NA samples insufficient (%d < %d), skip all interaction models",
            csn_nonNA, min_n_resid
        )
    } else {
        log_msg("  [DEG-INTERACTION] Running one-step models (subunit + CSN_SCORE + covariates)")

        for (su in present_sub) {
            su_nonNA <- sum(is.finite(mat0[su, ]))

            if (su_nonNA < (2 * min_per_group)) {
                log_msg("  [%s_adj_CSN] Insufficient non-NA samples, skip", su)
                next
            }

            run_deg_interaction_model(
                subunit_name = su,
                subunit_vec = mat0[su, ],
                csn_score_vec = csn_score,
                ds_id = ds_id,
                ds_dir = ds_dir,
                mat0 = mat0,
                out_root = out_root,
                batch_all = batch_all,
                purity_all = purity_all,
                sa_all = sa_all,
                min_per_group = min_per_group
            )
        }
    }

    invisible(NULL)
}

#' Run Complete Analysis for One Stratum
#'
#' Main orchestration function that manages the complete analysis workflow
#' for a single stratum (ALL, TP53_mutant, or TP53_wild_type).
#'
#' @param ds_id Dataset identifier
#' @param ds_dir Dataset directory path
#' @param mat0_full Full expression matrix (all samples)
#' @param sample_keep Character vector of sample IDs to analyze in this stratum
#' @param out_root Output directory for this stratum
#' @param genesets_by_group Named list of gene sets (for future GSEA, currently unused in DEG mode)
#'
#' @return Invisible NULL
#'
#' @details
#' Analysis workflow:
#' 1. Prepare data: subset samples, log transform, filter/impute
#' 2. Calculate CSN score and gather covariates (batch, purity, sex, age)
#' 3. Run DEG analyses:
#'    - Individual subunits (e.g., GPS1, COPS2, ...)
#'    - CSN complex score (PC1)
#'    - Interaction models (subunit effect controlling for CSN complex)
#' 4. All results saved to: {out_root}/DEG/BatchAdj/{predictor}/DEG_limma_cont.csv
#'
#' @examples
#' run_one_stratum(
#'     "brca_cptac_2020", ds_dir, mat0_full,
#'     all_samples, out_root, genesets_by_group
#' )
run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group) {
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
    log_msg("  -- stratum: %s | N(sample_keep)=%d", basename(out_root), length(sample_keep))

    ## Step 1: Prepare data
    data_prep <- .prepare_stratum_data(mat0_full, sample_keep, min_frac_complete)
    if (is.null(data_prep)) {
        return(invisible(NULL))
    }

    mat0 <- data_prep$mat0
    mat <- data_prep$mat
    present_sub <- data_prep$present_sub

    ## Step 2: Calculate predictors and covariates
    predictors <- .calculate_predictors_and_covariates(
        ds_id, ds_dir, mat0, present_sub, basename(out_root)
    )

    ## Step 3: Run all DEG analyses
    .run_all_deg_analyses(
        ds_id = ds_id,
        ds_dir = ds_dir,
        mat0 = mat0,
        out_root = out_root,
        present_sub = present_sub,
        csn_score = predictors$csn_score,
        batch_all = predictors$batch_all,
        purity_all = predictors$purity_all,
        sa_all = predictors$sa_all
    )

    log_msg("[Complete] Stratum analysis finished: %s", basename(out_root))
    invisible(NULL)
}
