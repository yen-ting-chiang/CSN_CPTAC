## ====================================================================
## 02_data_loading.R
##
## Purpose: Data loading functions for CSN CPTAC GSEA analysis
## Contains: protein matrix loading, clinical data reading, case list parsing
## ====================================================================

#' Read case list from cBioPortal case_lists file
#'
#' @param path_file Path to case list file
#' @return Character vector of case IDs
#' @export
read_case_list <- function(path_file) {
    if (!file.exists(path_file)) {
        return(character(0))
    }
    x <- readLines(path_file, warn = FALSE, encoding = "UTF-8")
    line <- x[grepl("^case_list_ids:", x)]
    if (!length(line)) {
        return(character(0))
    }
    ids <- sub("^case_list_ids:\\s*", "", line[1])
    ids <- unlist(strsplit(ids, "[,\\s]+"))
    unique(ids[nchar(ids) > 0])
}


#' Load protein quantification matrix from dataset directory
#'
#' @param dir Dataset directory path
#' @return Numeric matrix with genes as rows, samples as columns
#' @export
load_matrix_from_dataset_dir <- function(dir) {
    fp <- file.path(dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(glue::glue("File not found: {fp}"))

    if (exists("log_msg", mode = "function")) {
        log_msg("Reading protein matrix: {basename(fp)}")
    }

    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))

    # Identify gene column
    gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol", "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol")
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]
    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene)

    # Identify sample columns (exclude metadata columns)
    not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.", "ENTREZ_GENE_ID", "Description", "Gene_Name", "GeneName", "Gene_Symbol")
    sample_cols_all <- setdiff(names(dat), not_sample)

    # Filter by case list if available
    case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)
    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10 && exists("log_msg", mode = "function")) {
            log_msg("hint: case_list intersection too small ({length(inter)}); using all sample columns")
        }
    } else {
        sample_cols <- sample_cols_all
    }

    # Extract matrix
    m <- dat %>%
        dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
        janitor::remove_empty("cols")
    rn <- m$Gene
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn

    if (ncol(m) == 0) stop("0 sample columns read")

    # Handle duplicate genes (average)
    if (anyDuplicated(rownames(m))) {
        if (exists("log_msg", mode = "function")) {
            log_msg("Duplicate genes detected; averaging duplicate rows")
        }
        m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
    }

    if (exists("log_msg", mode = "function")) {
        log_msg("Matrix dimension: {nrow(m)} genes × {ncol(m)} samples")
    }

    m
}


#' Impute missing values and filter low-completeness genes
#'
#' @param mat Numeric matrix
#' @param min_frac Minimum fraction of non-NA values per gene
#' @return Imputed and filtered matrix
#' @export
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(!is.na(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    if (any(is.na(m))) {
        set.seed(1234) # Fixed seed for reproducibility
        m <- imputeLCMD::impute.MinProb(m, q = 0.01)
    }
    m
}


#' Get TP53 mutation status for samples
#'
#' @param ds_dir Dataset directory
#' @param sample_ids Character vector of sample IDs
#' @return Named character vector: "TP53_wild_type" or "TP53_mutant"
#' @export
get_tp53_status <- function(ds_dir, sample_ids) {
    # Define protein-altering variant classes
    TP53_KEEP_CLASSES <- c(
        "MISSENSE_MUTATION", "NONSENSE_MUTATION",
        "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
        "IN_FRAME_DEL", "IN_FRAME_INS",
        "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
    )

    # Helper: normalize variant classification
    normalize_vc <- function(x) {
        x <- toupper(trimws(as.character(x)))
        gsub("[^A-Z0-9]+", "_", x)
    }

    status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

    mut_fp <- file.path(ds_dir, "data_mutations.txt")
    if (!file.exists(mut_fp)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("  [TP53] data_mutations.txt not found; all samples considered wild-type")
        }
        return(status)
    }

    mutation_df <- tryCatch(
        readr::read_tsv(mut_fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )
    if (is.null(mutation_df)) {
        return(status)
    }

    req <- c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")
    miss <- setdiff(req, colnames(mutation_df))
    if (length(miss)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("  [TP53] data_mutations.txt missing columns: %s → treated as wild-type", paste(miss, collapse = ", "))
        }
        return(status)
    }

    tp53_df <- subset(mutation_df, Hugo_Symbol == "TP53",
        select = c("Variant_Classification", "Tumor_Sample_Barcode")
    )
    if (!nrow(tp53_df)) {
        return(status)
    }

    vc_norm <- normalize_vc(tp53_df$Variant_Classification)

    # Keep only protein-altering mutations
    keep <- vc_norm %in% TP53_KEEP_CLASSES |
        grepl("^FRAME_SHIFT_", vc_norm) |
        grepl("^IN_FRAME_", vc_norm)

    if (!any(keep)) {
        return(status)
    }

    tp53_samples <- toupper(unique(tp53_df$Tumor_Sample_Barcode[keep]))
    sid_up <- toupper(sample_ids)
    status[sid_up %in% tp53_samples] <- "TP53_mutant"
    status
}


#' Get sex and age covariates from clinical data
#'
#' @param ds_dir Dataset directory
#' @param sample_ids Character vector of sample IDs
#' @return Data frame with sex, age, age_z columns
#' @export
get_sex_age_covariates <- function(ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")

    sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    age_missing <- NULL
    age_z_imputed <- NULL

    if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("    sex/age: missing clinical files → all NA")
        }
        return(data.frame(
            sex = as.numeric(sex), age = as.numeric(age), age_z = as.numeric(age),
            row.names = sample_ids, check.names = FALSE
        ))
    }

    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()

    names(samp) <- .norm_names(names(samp))
    names(pat) <- .norm_names(names(pat))

    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]

    if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("    sex/age: ID columns not found → all NA")
        }
        return(data.frame(
            sex = as.numeric(sex), age = as.numeric(age), age_z = as.numeric(age),
            row.names = sample_ids, check.names = FALSE
        ))
    }

    # Create sample → patient mapping
    map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

    # SEX: 1=Male, 0=Female
    if ("SEX" %in% names(pat)) {
        raw <- toupper(as.character(pat$SEX))
        val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
        names(val) <- as.character(pat[[pid_pat]])
        sex[] <- unname(val[map_pt[sample_ids]])
    }

    # AGE
    if ("AGE" %in% names(pat)) {
        v <- suppressWarnings(as.numeric(pat$AGE))
        names(v) <- as.character(pat[[pid_pat]])
        age_raw <- unname(v[map_pt[sample_ids]])
        age[] <- .z_no_impute(age_raw) # NA preserved

        # Optional: age missing indicator
        if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
            age_missing <- as.numeric(is.na(age_raw))
            age_z_imputed <- .zscore(age_raw) # NA imputed
        }
    }

    cov_sex <- mean(is.finite(sex)) * 100
    cov_age <- mean(is.finite(age)) * 100
    if (exists("log_msg", mode = "function")) {
        if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE)) && !is.null(age_z_imputed)) {
            cov_age_imp <- mean(is.finite(age_z_imputed)) * 100
            cov_age_mis <- mean(is.finite(age_missing)) * 100
            log_msg(
                "    covariates coverage: sex %.1f%%, age(NA-as-NA) %.1f%%, age_z_imputed %.1f%%, age_missing %.1f%%",
                cov_sex, cov_age, cov_age_imp, cov_age_mis
            )
        } else {
            log_msg("    covariates coverage: sex %.1f%%, age(NA-as-NA) %.1f%%", cov_sex, cov_age)
        }
    }

    out <- data.frame(
        sex = as.numeric(sex), age = as.numeric(age),
        row.names = sample_ids, check.names = FALSE
    )
    out$age_z <- out$age

    if (!is.null(age_missing)) {
        out$age_missing <- as.numeric(age_missing)
        out$age_z_imputed <- as.numeric(age_z_imputed)
    }

    out
}


# ---- End of 02_data_loading.R ----
