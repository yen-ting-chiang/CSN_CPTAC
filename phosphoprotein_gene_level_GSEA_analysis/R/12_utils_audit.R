# =============================================================================
# 12_utils_audit.R - Dataset Audit Utilities
# =============================================================================
# This file contains functions for auditing dataset covariates and batch info.
# =============================================================================

# -----------------------------------------------------------------------------
# Inspect PLEX Column 
# -----------------------------------------------------------------------------
inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
    # Read samples actually used (based on protein matrix)
    mat0 <- load_phospho_matrix_from_dataset_dir(ds_dir)
    sample_ids <- colnames(mat0)

    # Read clinical table and align samples
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    stopifnot(file.exists(meta_fp))
    meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"), names(meta))
    stopifnot(length(id_cols) > 0)
    meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
    meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]

    # Get raw and cleaned PLEX
    raw <- as.character(meta[[col]])
    clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)

    # Summary
    cat("\n==== ", basename(ds_dir), " | Column: ", col, " ====\n", sep = "")
    cat("[Raw PLEX levels]:\n")
    print(sort(table(raw), decreasing = TRUE))
    cat("\n[Raw values containing '|' (first few)]:\n")
    print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))

    cat("\n[Cleaned PLEX levels] (pipe_policy = ", pipe_policy,
        ", min_per_level = ", min_per_level, "):\n",
        sep = ""
    )
    print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))

    # Output mapping table (first 10 rows as demo)
    df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
    cat("\n[Sample mapping table (first 10 rows)]\n")
    print(utils::head(df_map, 10))

    invisible(list(
        raw_counts = sort(table(raw), decreasing = TRUE),
        clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
        map = df_map
    ))
}

# -----------------------------------------------------------------------------
# Audit Purity Source and Coverage by Dataset Rules
# -----------------------------------------------------------------------------
.audit_purity_for_dataset <- function(ds_id, ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    has_samp <- file.exists(samp_fp)
    has_pat <- file.exists(pat_fp)

    purity_col <- "NONE"
    purity_vec <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

    samp <- pat <- NULL
    sid <- pid_samp <- pid_pat <- NA

    if (has_samp) {
        samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
        names(samp) <- .norm_names(names(samp))
        sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
        pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
    }
    if (has_pat) {
        pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
        names(pat) <- .norm_names(names(pat))
        pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
    }

    # Build sample->patient mapping (when needed)
    map_pt <- NULL
    if (!is.na(sid) && !is.na(pid_samp)) map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

    ds <- ds_id

    if (ds == "brca_cptac_2020" && has_pat && !is.na(pid_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat) && !is.null(map_pt)) {
        purity_col <- "PATIENT:ESTIMATE_TUMORPURITY"
        pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_pat]]))
        purity_vec[] <- unname(pt_pur[map_pt[sample_ids]])
    } else if (ds == "luad_cptac_2020" && has_samp && !is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp)) {
        purity_col <- "SAMPLE:TUMOR_PURITY_BYESTIMATE_RNASEQ"
        purity_vec[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][match(sample_ids, samp[[sid]])])
    } else if (ds == "lusc_cptac_2021" && has_samp && !is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp)) {
        purity_col <- "SAMPLE:ESTIMATE_TUMORPURITY"
        purity_vec[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][match(sample_ids, samp[[sid]])])
    } else if (ds == "paad_cptac_2021" && has_samp && !is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
        purity_col <- "SAMPLE:NEOPLASTIC_CELLULARITY(median;)"
        medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
        purity_vec[] <- .to01(medv[match(sample_ids, samp[[sid]])])
    } else if (ds == "ucec_cptac_2020" && has_samp && !is.na(sid)) {
        if ("PURITY_CANCER" %in% names(samp)) {
            purity_col <- "SAMPLE:PURITY_CANCER"
            purity_vec[] <- .to01(samp[["PURITY_CANCER"]][match(sample_ids, samp[[sid]])])
        } else {
            # Try to infer from immune/stroma (if columns exist)
            pi <- if ("PURITY_IMMUNE" %in% names(samp)) suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]])) else NULL
            ps <- if ("PURITY_STROMA" %in% names(samp)) suppressWarnings(as.numeric(samp[["PURITY_STROMA"]])) else NULL
            if (!is.null(pi) || !is.null(ps)) {
                idx <- match(sample_ids, samp[[sid]])
                purity_est <- 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0)
                purity_vec[] <- .to01(purity_est)
                purity_col <- "SAMPLE:1-(PURITY_IMMUNE+PURITY_STROMA)"
            }
        }
    } else {
        purity_col <- "NONE"
    }

    nonNA <- mean(is.finite(purity_vec)) * 100
    pur_min <- if (any(is.finite(purity_vec))) min(purity_vec, na.rm = TRUE) else NA_real_
    pur_med <- if (any(is.finite(purity_vec))) stats::median(purity_vec, na.rm = TRUE) else NA_real_
    pur_max <- if (any(is.finite(purity_vec))) max(purity_vec, na.rm = TRUE) else NA_real_

    list(
        column = purity_col,
        nonNA_pct = nonNA,
        min = pur_min, median = pur_med, max = pur_max,
        vec = purity_vec
    )
}

# -----------------------------------------------------------------------------
# Audit One Dataset: Sex/Age + Purity + Batch
# -----------------------------------------------------------------------------
audit_one_dataset_sa_batch <- function(ds_dir, pipe_policy = "NA", min_per_level = 2) {
    ds_id <- basename(ds_dir)
    # 1) Get protein matrix samples
    m <- load_phospho_matrix_from_dataset_dir(ds_dir)
    sample_ids <- colnames(m)

    # 2) Only get SEX/AGE from patient file
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    has_samp <- file.exists(samp_fp)
    has_pat <- file.exists(pat_fp)

    sex_nonNA <- age_nonNA <- NA_real_
    sex_M <- sex_F <- NA_integer_
    age_min <- age_med <- age_max <- NA_real_
    sex_col <- age_col <- "MISSING"
    sex_src <- age_src <- if (has_pat) "patient" else "NONE"

    if (has_samp && has_pat) {
        samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
        pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
        names(samp) <- .norm_names(names(samp))
        names(pat) <- .norm_names(names(pat))

        sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
        pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
        pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]

        if (!is.na(sid) && !is.na(pid_samp) && !is.na(pid_pat)) {
            map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

            # SEX
            if ("SEX" %in% names(pat)) {
                sex_col <- "SEX"
                raw <- toupper(as.character(pat$SEX))
                val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
                names(val) <- as.character(pat[[pid_pat]])
                sex_vec <- unname(val[map_pt[sample_ids]])
                sex_nonNA <- mean(is.finite(sex_vec)) * 100
                sex_M <- sum(sex_vec == 1, na.rm = TRUE)
                sex_F <- sum(sex_vec == 0, na.rm = TRUE)
            }

            # AGE (raw value statistics; z-score done in main flow)
            if ("AGE" %in% names(pat)) {
                age_col <- "AGE"
                v <- suppressWarnings(as.numeric(pat$AGE))
                names(v) <- as.character(pat[[pid_pat]])
                age_vec <- unname(v[map_pt[sample_ids]])
                age_nonNA <- mean(is.finite(age_vec)) * 100
                if (any(is.finite(age_vec))) {
                    age_min <- min(age_vec, na.rm = TRUE)
                    age_med <- stats::median(age_vec, na.rm = TRUE)
                    age_max <- max(age_vec, na.rm = TRUE)
                }
            }
        }
    }

    # 3) purity audit
    pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)

    # 4) batch check (using existing get_batch_factor_phospho)
    bi <- get_batch_factor_phospho(ds_dir, sample_ids,
        pipe_policy   = pipe_policy,
        min_per_level = min_per_level
    )
    batch_col <- if (!is.null(bi)) bi$name else "NONE"
    batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
    batch_nonNA <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
    batch_sizes <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_

    # Print one-line summary (added policy/batch_min/limma_min)
    line <- sprintf(
        "[Audit] %-24s | SEX:%-7s nonNA=%5.1f%% (M=%d,F=%d) | AGE:%-7s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | PURITY:%-35s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | BATCH:%-18s levels=%-2s nonNA=%5.1f%% [%s] | policy=%s | batch_min=%d | limma_min=%d",
        ds_id,
        sex_col, sex_nonNA %||% NaN, sex_M %||% NA_integer_, sex_F %||% NA_integer_,
        age_col, age_nonNA %||% NaN,
        ifelse(is.finite(age_min), round(age_min, 1), "NA"),
        ifelse(is.finite(age_med), round(age_med, 1), "NA"),
        ifelse(is.finite(age_max), round(age_max, 1), "NA"),
        pur$column, pur$nonNA_pct %||% NaN,
        ifelse(is.finite(pur$min), round(pur$min, 3), "NA"),
        ifelse(is.finite(pur$median), round(pur$median, 3), "NA"),
        ifelse(is.finite(pur$max), round(pur$max, 3), "NA"),
        batch_col, as.integer(batch_levels), batch_nonNA %||% NaN, batch_sizes %||% "NA",
        pipe_policy, min_per_level, min_per_group
    )
    if (exists("log_msg", mode = "function")) log_msg(line) else cat(line, "\n")

    data.frame(
        dataset = ds_id,
        sex_column = sex_col, sex_nonNA_pct = round(sex_nonNA, 1), sex_M = sex_M, sex_F = sex_F,
        age_column = age_col, age_nonNA_pct = round(age_nonNA, 1),
        age_min = ifelse(is.finite(age_min), round(age_min, 1), NA),
        age_median = ifelse(is.finite(age_med), round(age_med, 1), NA),
        age_max = ifelse(is.finite(age_max), round(age_max, 1), NA),
        purity_column = pur$column, purity_nonNA_pct = round(pur$nonNA_pct, 1),
        purity_min = ifelse(is.finite(pur$min), round(pur$min, 3), NA),
        purity_median = ifelse(is.finite(pur$median), round(pur$median, 3), NA),
        purity_max = ifelse(is.finite(pur$max), round(pur$max, 3), NA),
        batch_column = batch_col, batch_levels = as.integer(batch_levels),
        batch_nonNA_pct = round(batch_nonNA, 1), batch_sizes = batch_sizes,
        batch_policy = pipe_policy,
        batch_min_per_level = as.integer(min_per_level),
        limma_min_per_group = as.integer(min_per_group),
        stringsAsFactors = FALSE
    )
}

# -----------------------------------------------------------------------------
# Audit All Datasets: Sex/Age + Purity + Batch
# -----------------------------------------------------------------------------
audit_all_datasets_sa_batch <- function(dataset_dirs, pipe_policy = BATCH_PIPE_POLICY, min_per_level = 2) {
    out <- lapply(dataset_dirs, function(p) {
        tryCatch(audit_one_dataset_sa_batch(p, pipe_policy, min_per_level),
            error = function(e) {
                if (exists("log_msg", mode = "function")) log_msg("[Audit] %s ERROR: %s", basename(p), e$message)
                data.frame(
                    dataset = basename(p),
                    sex_column = "ERROR", sex_nonNA_pct = NA, sex_M = NA, sex_F = NA,
                    age_column = "ERROR", age_nonNA_pct = NA, age_min = NA, age_median = NA, age_max = NA,
                    purity_column = "ERROR", purity_nonNA_pct = NA, purity_min = NA, purity_median = NA, purity_max = NA,
                    batch_column = "ERROR", batch_levels = NA, batch_nonNA_pct = NA, batch_sizes = NA,
                    stringsAsFactors = FALSE
                )
            }
        )
    })
    dplyr::bind_rows(out)
}

# -----------------------------------------------------------------------------
# Utility: Check if Stratum Should Run
# -----------------------------------------------------------------------------
.should_run <- function(tag) {
    if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) {
        return(TRUE)
    }
    tag %in% only_strata
}
