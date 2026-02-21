# =============================================================================
# 03_utils_io.R
# File Reading and I/O Utilities
# =============================================================================

# -----------------------------------------------------------------------------
# Case List Reading
# -----------------------------------------------------------------------------

#' Read case list from cBioPortal-style file
#'
#' @param path_file Path to case list file
#' @return Character vector of sample IDs
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

# -----------------------------------------------------------------------------
# Protein Matrix Loading
# -----------------------------------------------------------------------------

#' Load protein quantification matrix from dataset directory
#'
#' Reads data_protein_quantification.txt, handles gene column naming,
#' filters samples, and averages duplicate genes.
#'
#' @param dir Dataset directory path
#' @return Numeric matrix with genes as rows, samples as columns
load_matrix_from_dataset_dir <- function(dir) {
    fp <- file.path(dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(glue::glue("File not found: {fp}"))
    log_msg("Reading protein matrix: {basename(fp)}")
    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))
    gene_cols <- c(
        "Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol",
        "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol"
    )
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]
    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene)
    not_sample <- c(
        "Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.", "ENTREZ_GENE_ID",
        "Description", "Gene_Name", "GeneName", "Gene_Symbol"
    )
    sample_cols_all <- setdiff(names(dat), not_sample)
    # If folder has case_list, try to match
    case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)
    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10) log_msg("Note: case_list intersection too small ({length(inter)}), using all sample columns")
    } else {
        sample_cols <- sample_cols_all
    }
    m <- dat %>%
        dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
        janitor::remove_empty("cols")
    rn <- m$Gene
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn
    if (ncol(m) == 0) stop("Read 0 sample columns")
    if (anyDuplicated(rownames(m))) {
        log_msg("Detected duplicate genes, averaging duplicate rows")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
    }
    log_msg("Matrix dimensions: {nrow(m)} genes x {ncol(m)} samples")
    m
}

# -----------------------------------------------------------------------------
# Safe CSV Reading
# -----------------------------------------------------------------------------

#' Safely read CSV file with multiple fallback methods
#'
#' Tries data.table::fread, then readr::read_csv, then base R read.csv.
#' Handles Windows path separators.
#'
#' @param path Path to CSV file
#' @return tibble containing CSV data
.read_csv_safe <- function(path) {
    p <- path

    # Windows: if original string doesn't exist, try replacing / with \\
    if (!file.exists(p) && .Platform$OS.type == "windows") {
        p2 <- gsub("/", "\\\\", p, fixed = TRUE)
        if (file.exists(p2)) p <- p2
    }

    if (!file.exists(p)) stop("File not found: ", path)

    # 1) Try data.table::fread
    out <- try(suppressMessages(data.table::fread(p, check.names = FALSE)), silent = TRUE)
    if (!inherits(out, "try-error")) {
        return(tibble::as_tibble(out))
    }

    # 2) Try readr::read_csv
    if (requireNamespace("readr", quietly = TRUE)) {
        out <- try(suppressMessages(readr::read_csv(
            p,
            show_col_types = FALSE, guess_max = 100000
        )), silent = TRUE)
        if (!inherits(out, "try-error")) {
            return(out)
        }
    }

    # 3) Last resort: base R
    out <- try(utils::read.csv(p, check.names = FALSE), silent = TRUE)
    if (inherits(out, "try-error")) stop(out)
    tibble::as_tibble(out)
}
