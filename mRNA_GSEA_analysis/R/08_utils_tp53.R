# =============================================================================
# 08_utils_tp53.R - TP53 Stratification Utilities
# =============================================================================
# This module provides functions for determining TP53 mutation status
# and stratifying samples for analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# TP53 Protein-Altering Variant Classes
# -----------------------------------------------------------------------------
#' Define protein-altering variant classes for TP53
TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION", "NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
    "IN_FRAME_DEL", "IN_FRAME_INS",
    "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

# -----------------------------------------------------------------------------
# Normalize Variant Classification
# -----------------------------------------------------------------------------
#' Normalize variant classification string for comparison
#' @param x Character vector of variant classifications
#' @return Normalized uppercase strings with underscores
normalize_vc <- function(x) {
    x <- toupper(trimws(as.character(x)))
    gsub("[^A-Z0-9]+", "_", x)
}

# -----------------------------------------------------------------------------
# Get TP53 Status
# -----------------------------------------------------------------------------
#' Determine TP53 mutation status for samples
#' @param ds_dir Dataset directory path
#' @param sample_ids Sample IDs to classify
#' @return Named character vector with "TP53_mutant" or "TP53_wild_type"
get_tp53_status <- function(ds_dir, sample_ids) {
    # Default all wild-type
    status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

    mut_fp <- file.path(ds_dir, "data_mutations.txt")
    if (!file.exists(mut_fp)) {
        log_msg("  [TP53] data_mutations.txt not found, treating all as wild-type/ALL available")
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
        log_msg("  [TP53] data_mutations.txt missing columns: %s -> treating as wild-type", paste(miss, collapse = ", "))
        return(status)
    }

    tp53_df <- subset(mutation_df, Hugo_Symbol == "TP53",
        select = c("Variant_Classification", "Tumor_Sample_Barcode")
    )
    if (!nrow(tp53_df)) {
        return(status)
    }

    vc_norm <- normalize_vc(tp53_df$Variant_Classification)

    # Strictly use protein-altering categories
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

# -----------------------------------------------------------------------------
# Get Protein-Altering Variant Classes
# -----------------------------------------------------------------------------
#' Return list of protein-altering variant classes
#' @return Character vector
.protein_altering_vc <- function() {
    c(
        "MISSENSE_MUTATION", "NONSENSE_MUTATION",
        "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
        "IN_FRAME_DEL", "IN_FRAME_INS",
        "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
    )
}

# -----------------------------------------------------------------------------
# Summarize TP53 Counts for Dataset
# -----------------------------------------------------------------------------
#' Generate TP53 sample count summary for a dataset
#' @param ds_dir Dataset directory path
#' @return List with binary and class-level counts
summarize_tp53_counts_for_dataset <- function(ds_dir) {
    ds_id <- basename(ds_dir)

    # Use protein matrix samples as the population
    M <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
    if (inherits(M, "try-error")) {
        return(list(binary = NULL, class_long = NULL))
    }
    sample_ids <- colnames(M)
    tp53_status <- get_tp53_status(ds_dir, sample_ids)

    tb_bin <- table(tp53_status)
    bin_row <- data.frame(
        dataset = ds_id,
        total_samples = length(sample_ids),
        WT_n = as.integer(tb_bin["TP53_wild_type"]),
        MUT_n = as.integer(tb_bin["TP53_mutant"]),
        stringsAsFactors = FALSE
    )

    # Read MAF (original Variant_Classification distribution)
    mut_fp <- file.path(ds_dir, "data_mutations.txt")
    class_df <- NULL
    if (file.exists(mut_fp)) {
        mutation_df <- try(readr::read_tsv(mut_fp, comment = "#", show_col_types = FALSE), silent = TRUE)
        if (!inherits(mutation_df, "try-error") &&
            all(c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode") %in% colnames(mutation_df))) {
            tp53 <- subset(mutation_df, Hugo_Symbol == "TP53",
                select = c("Variant_Classification", "Tumor_Sample_Barcode")
            )
            if (nrow(tp53) > 0) {
                tp53$Variant_Classification <- normalize_vc(tp53$Variant_Classification)
                tp53$Tumor_Sample_Barcode <- toupper(tp53$Tumor_Sample_Barcode)

                # Only count samples in the protein matrix
                tp53 <- tp53[tp53$Tumor_Sample_Barcode %in% toupper(sample_ids), ]
                class_counts <- data.frame(
                    Variant_Classification = names(table(tp53$Variant_Classification)),
                    sample_n = as.integer(table(tp53$Variant_Classification)),
                    stringsAsFactors = FALSE
                )

                any_samples <- unique(tp53$Tumor_Sample_Barcode)
                prot_alt <- .protein_altering_vc()
                prot_alt_samples <- unique(tp53$Tumor_Sample_Barcode[tp53$Variant_Classification %in% prot_alt])

                add_rows <- data.frame(
                    Variant_Classification = c("Any_TP53_mutation", "protein_altering"),
                    sample_n = c(length(any_samples), length(prot_alt_samples)),
                    stringsAsFactors = FALSE
                )

                class_df <- dplyr::bind_rows(class_counts, add_rows)
                class_df$dataset <- ds_id
                class_df <- class_df[, c("dataset", "Variant_Classification", "sample_n")]
            }
        }
    }

    # Write individual files
    dir.create(file.path("run_info", "tp53_status"), recursive = TRUE, showWarnings = FALSE)
    if (!is.null(class_df)) {
        data.table::fwrite(
            class_df,
            file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
        )
    } else {
        data.table::fwrite(
            data.frame(dataset = ds_id, Variant_Classification = NA_character_, sample_n = NA_integer_),
            file.path("run_info", "tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
        )
    }

    list(binary = bin_row, class_long = class_df)
}

# -----------------------------------------------------------------------------
# Summarize TP53 Counts All Datasets
# -----------------------------------------------------------------------------
#' Generate TP53 summary across all datasets
#' @param dataset_dirs Named list/vector of dataset directories
summarize_tp53_counts_all_datasets <- function(dataset_dirs) {
    all_bin <- list()
    all_class <- list()
    k <- 1L
    j <- 1L
    for (ds in names(dataset_dirs)) {
        ds_dir <- dataset_dirs[[ds]]
        res <- try(summarize_tp53_counts_for_dataset(ds_dir), silent = TRUE)
        if (!inherits(res, "try-error")) {
            if (!is.null(res$binary)) {
                all_bin[[k]] <- res$binary
                k <- k + 1L
            }
            if (!is.null(res$class_long)) {
                all_class[[j]] <- res$class_long
                j <- j + 1L
            }
        }
    }
    if (length(all_bin)) {
        df_bin <- dplyr::bind_rows(all_bin)
        data.table::fwrite(df_bin, file.path("run_info", "tp53_status", "tp53_binary_counts_by_dataset.csv"))
    }
    if (length(all_class)) {
        df_class <- dplyr::bind_rows(all_class)
        data.table::fwrite(df_class, file.path("run_info", "tp53_status", "tp53_class_sample_counts_long.csv"))
    }
}
