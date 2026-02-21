## =========================================================
## 05_tp53_status.R
## TP53 Mutation Status Functions for DPS Analysis
## =========================================================

## ===== TP53 protein-altering variant classes =====
TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION", "NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
    "IN_FRAME_DEL", "IN_FRAME_INS",
    "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

## ===== Normalize variant classification =====
normalize_vc <- function(x) {
    x <- toupper(trimws(as.character(x)))
    gsub("[^A-Z0-9]+", "_", x)
}

## ===== Get TP53 mutation status =====
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

    # Strictly use protein-altering classes
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

## ===== Summarize TP53 counts for one dataset =====
summarize_tp53_counts_for_dataset <- function(ds_dir) {
    ds_id <- basename(ds_dir)

    # Use phospho matrix samples as "population"
    M <- try(load_phospho_matrix_from_dataset_dir(ds_dir), silent = TRUE)
    if (inherits(M, "try-error")) {
        log_msg("[TP53-audit] %s: Cannot read matrix, skipping", ds_id)
        return(NULL)
    }
    sample_ids <- colnames(M)
    sid_up <- toupper(sample_ids)
    n_all <- length(sample_ids)

    # Get mutant/wild type
    status <- get_tp53_status(ds_dir, sample_ids)
    tb_bin <- table(factor(status, levels = c("TP53_wild_type", "TP53_mutant")))
    bin_row <- data.frame(
        dataset = ds_id,
        in_matrix_n = n_all,
        WT_n = as.integer(tb_bin["TP53_wild_type"]),
        MUT_n = as.integer(tb_bin["TP53_mutant"]),
        stringsAsFactors = FALSE
    )

    # Read MAF
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

                # Only count samples "in protein matrix"
                tp53 <- tp53[tp53$Tumor_Sample_Barcode %in% sid_up, , drop = FALSE]

                # Each classification -> how many "unique samples" hit
                class_counts <- tp53 |>
                    dplyr::group_by(Variant_Classification) |>
                    dplyr::summarise(sample_n = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
                    dplyr::arrange(dplyr::desc(sample_n), Variant_Classification)

                # Add summary rows
                any_samples <- unique(tp53$Tumor_Sample_Barcode)

                keep_flag <- tp53$Variant_Classification %in% TP53_KEEP_CLASSES |
                    grepl("^FRAME_SHIFT_", tp53$Variant_Classification) |
                    grepl("^IN_FRAME_", tp53$Variant_Classification)
                prot_alt_samples <- unique(tp53$Tumor_Sample_Barcode[keep_flag])

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

## ===== Summarize TP53 counts for all datasets =====
summarize_tp53_counts_all_datasets <- function(dataset_dirs) {
    all_bin <- list()
    all_class <- list()
    k <- 1L
    j <- 1L
    for (ds in names(dataset_dirs)) {
        ds_dir <- dataset_dirs[[ds]]
        if (!dir.exists(ds_dir)) next
        log_msg("[TP53-audit] Starting: %s", ds)
        res <- summarize_tp53_counts_for_dataset(ds_dir)
        if (is.null(res)) next
        all_bin[[k]] <- res$binary
        k <- k + 1L
        if (!is.null(res$class_long)) {
            all_class[[j]] <- res$class_long
            j <- j + 1L
        }
    }

    if (length(all_bin)) {
        bin_df <- dplyr::bind_rows(all_bin) |>
            dplyr::mutate(Any_TP53_mutation_n = in_matrix_n - WT_n)
        data.table::fwrite(bin_df, file.path("run_info", "tp53_status", "tp53_binary_counts_by_dataset.csv"))
        log_msg("[TP53-audit] Written: tp53_binary_counts_by_dataset.csv")
    }
    if (length(all_class)) {
        class_df <- dplyr::bind_rows(all_class)
        data.table::fwrite(class_df, file.path("run_info", "tp53_status", "tp53_class_sample_counts_long.csv"))
        log_msg("[TP53-audit] Written: tp53_class_sample_counts_long.csv")
    }
}
