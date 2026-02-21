# =============================================================================
# 08_utils_tp53.R - TP53 Mutation Status Handling
# =============================================================================
# This file contains functions for parsing TP53 mutation status from cBioPortal data.
# =============================================================================

# -----------------------------------------------------------------------------
# TP53 Keep Classes (Protein-altering Mutations)
# -----------------------------------------------------------------------------
TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION", "NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
    "IN_FRAME_DEL", "IN_FRAME_INS",
    "SPLICE_SITE", "TRANSLATION_START_SITE",
    "NONSTOP_MUTATION"
)

# -----------------------------------------------------------------------------
# Normalize Variant Classification
# -----------------------------------------------------------------------------
normalize_vc <- function(x) {
    x <- toupper(trimws(as.character(x)))
    gsub("[^A-Z0-9]+", "_", x)
}

# -----------------------------------------------------------------------------
# Get TP53 Mutation Status for Each Sample
# -----------------------------------------------------------------------------
get_tp53_status <- function(ds_dir, sample_ids) {
    # Default all wild-type
    status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)

    mut_fp <- file.path(ds_dir, "data_mutations.txt")
    if (!file.exists(mut_fp)) {
        log_msg("  [TP53] Cannot find data_mutations.txt, treating entire batch as wild-type/ALL available")
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
