# =============================================================================
# 06_utils_tp53.R
# TP53 Mutation Status Utilities
# =============================================================================

# -----------------------------------------------------------------------------
# Variant Classification Normalization
# -----------------------------------------------------------------------------

#' Normalize variant classification string
#'
#' Converts to uppercase and replaces non-alphanumeric characters with underscore.
#'
#' @param x Character vector of variant classifications
#' @return Normalized character vector
normalize_vc <- function(x) {
    x <- toupper(trimws(as.character(x)))
    gsub("[^A-Z0-9]+", "_", x)
}

# -----------------------------------------------------------------------------
# TP53 Status Detection
# -----------------------------------------------------------------------------

#' Get TP53 mutation status for samples
#'
#' Reads mutation data and classifies samples as TP53_mutant or TP53_wild_type
#' based on protein-altering variants.
#'
#' @param ds_dir Dataset directory path
#' @param sample_ids Character vector of sample IDs
#' @return Named character vector with TP53 status for each sample
get_tp53_status <- function(ds_dir, sample_ids) {
    # Default all to wild-type
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

    # Strictly use protein-altering classes; also tolerate FRAME_SHIFT_* / IN_FRAME_* prefixes
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
