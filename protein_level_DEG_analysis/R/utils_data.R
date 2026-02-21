## ============================================================================
## Data Loading and Filtering Utilities
## Functions for reading protein quantification data and case lists
## ============================================================================

#' Read Case List from cBioPortal Format File
#'
#' Parses a case list file (e.g., cases_protein_quantification.txt) and
#' extracts sample IDs.
#'
#' @param path_file Path to the case list file
#' @return Character vector of sample IDs, or empty vector if file not found
#'
#' @examples
#' ids <- read_case_list("brca_cptac_2020/case_lists/cases_protein_quantification.txt")
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

#' Load Protein Quantification Matrix from Dataset Directory
#'
#' Reads the data_protein_quantification.txt file from a cBioPortal dataset
#' directory and converts it to a numeric matrix with genes as rows and
#' samples as columns.
#'
#' @param dir Path to the dataset directory
#' @return Numeric matrix (genes x samples) with gene symbols as rownames
#'
#' @details
#' - Automatically detects gene symbol column (looks for Hugo_Symbol, Gene, etc.)
#' - Filters samples based on case list if available
#' - Averages duplicate gene entries
#' - Removes non-numeric columns
#'
#' @examples
#' mat <- load_matrix_from_dataset_dir("brca_cptac_2020")
load_matrix_from_dataset_dir <- function(dir) {
    fp <- file.path(dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(glue::glue("File not found: {fp}"))

    log_msg("Reading protein matrix: {basename(fp)}")
    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))

    # Identify gene column
    gene_cols <- c(
        "Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol",
        "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol"
    )
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]

    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene) # Remove Entrez ID if present

    # Identify sample columns (exclude metadata columns)
    not_sample <- c(
        "Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.",
        "ENTREZ_GENE_ID", "Description", "Gene_Name",
        "GeneName", "Gene_Symbol"
    )
    sample_cols_all <- setdiff(names(dat), not_sample)

    # Filter by case list if available
    case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)

    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10) {
            log_msg("Hint: case_list intersection too small ({length(inter)}), using all sample columns instead")
        }
    } else {
        sample_cols <- sample_cols_all
    }

    # Create matrix
    m <- dat %>%
        dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
        janitor::remove_empty("cols")

    rn <- m$Gene
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn

    if (ncol(m) == 0) stop("Read 0 sample columns")

    # Handle duplicate genes by averaging
    if (anyDuplicated(rownames(m))) {
        log_msg("Detected duplicate genes, averaging duplicate rows")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) /
            as.vector(table(rownames(m)))
    }

    log_msg("Matrix dimensions: {nrow(m)} genes x {ncol(m)} samples")
    m
}

#' Impute Missing Values and Filter Low-Coverage Genes
#'
#' Filters genes with too many missing values, then imputes remaining
#' NAs using MinProb method (left-censored missing data imputation).
#'
#' @param mat Numeric matrix (genes x samples)
#' @param min_frac Minimum fraction of non-NA values required per gene
#' @return Imputed matrix with low-coverage genes removed
#'
#' @details
#' Uses imputeLCMD::impute.MinProb with quantile = 0.01 and fixed seed
#' for reproducibility.
#'
#' @examples
#' mat_clean <- impute_and_filter(mat, min_frac = 0.75)
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(!is.na(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]

    if (any(is.na(m))) {
        # NOTE: Rely on global set.seed(1234) from config.R
        # Do NOT set seed here to match original script behavior
        m <- imputeLCMD::impute.MinProb(m, q = 0.01)
    }

    m
}

#' Get Dataset Directories Mapping
#'
#' Creates a named vector mapping dataset IDs to their full paths,
#' and filters to only include existing directories with protein data.
#'
#' @param dataset_ids Character vector of dataset IDs
#' @param datasets_root Root directory containing all datasets
#' @return Named character vector of paths to valid datasets
#'
#' @details
#' Only returns datasets that:
#' - Have an existing directory
#' - Contain a data_protein_quantification.txt file
#'
#' @examples
#' dirs <- get_dataset_dirs(
#'     c("brca_cptac_2020", "luad_cptac_2020"),
#'     getwd()
#' )
get_dataset_dirs <- function(dataset_ids, datasets_root = getwd()) {
    dataset_dirs <- setNames(
        file.path(datasets_root, dataset_ids),
        dataset_ids
    )

    # Check for missing directories
    missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
    if (length(missing_dirs)) {
        log_msg(
            "Detected %d missing folders, skipping: %s",
            length(missing_dirs), paste(missing_dirs, collapse = ", ")
        )
    }

    # Filter to existing directories with protein data
    dataset_dirs_run <- dataset_dirs[
        dir.exists(dataset_dirs) &
            file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
    ]

    if (!length(dataset_dirs_run)) {
        stop("No valid datasets found. Check if folders and data_protein_quantification.txt exist")
    }

    log_msg(
        "Available datasets this round (%d): %s",
        length(dataset_dirs_run),
        paste(names(dataset_dirs_run), collapse = ", ")
    )

    dataset_dirs_run
}
