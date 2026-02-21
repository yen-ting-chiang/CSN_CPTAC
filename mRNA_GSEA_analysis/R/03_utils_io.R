# =============================================================================
# 03_utils_io.R - File I/O Utilities
# =============================================================================
# This module provides file reading and writing utilities for the mRNA GSEA
# analysis pipeline, including matrix loading and case list parsing.
# =============================================================================

# -----------------------------------------------------------------------------
# Read Case List
# -----------------------------------------------------------------------------
#' Parse cBioPortal case list file to extract sample IDs
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
# Load RNA Expression Matrix
# -----------------------------------------------------------------------------
#' Load mRNA expression matrix from dataset directory
#' @param dir Dataset directory path
#' @return Numeric matrix with genes as rows, samples as columns
load_matrix_from_dataset_dir <- function(dir) {
    # Try the following filenames in order (pick first hit)
    candidates <- c(
        "data_mrna_seq_fpkm.txt",
        "data_mrna_seq_v2_rsem.txt",
        "data_mrna_seq_rpkm.txt",
        "data_mrna_seq_rsem.txt"
    )
    hit_flags <- file.exists(file.path(dir, candidates))
    if (!any(hit_flags)) {
        stop(glue::glue(
            "RNA file not found: {dir}\nTried: {paste(candidates, collapse = ', ')}"
        ))
    }
    fp <- file.path(dir, candidates[which(hit_flags)[1]])
    log_msg("Reading RNA matrix: {basename(fp)}")

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

    # RNA-specific case list: prefer matching data filename, fallback to common names
    bn <- basename(fp)
    case_candidates <- c(
        file.path(dir, "case_lists", sub("^data_", "cases_", bn)),
        file.path(dir, "case_lists", "cases_mrna_seq_fpkm.txt"),
        file.path(dir, "case_lists", "cases_mrna_seq_v2_rsem.txt"),
        file.path(dir, "case_lists", "cases_mrna_seq_rpkm.txt"),
        file.path(dir, "case_lists", "cases_mrna_seq_rsem.txt")
    )
    case_hit <- case_candidates[file.exists(case_candidates)]
    keep_ids <- if (length(case_hit)) read_case_list(case_hit[1]) else character(0)

    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10) log_msg("Note: case_list intersection too small (%d), using all sample columns instead", length(inter))
    } else {
        sample_cols <- sample_cols_all
    }

    m <- dat %>%
        dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
        janitor::remove_empty("cols")
    rn <- as.character(m$Gene)
    ok <- !is.na(rn) & nzchar(trimws(rn))
    m <- m[ok, , drop = FALSE]
    rn <- rn[ok]
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn

    if (ncol(m) == 0) stop("Read 0 sample columns")
    if (anyDuplicated(rownames(m))) {
        log_msg("Detected duplicate genes, averaging duplicate rows")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
    }
    log_msg("RNA matrix dimensions: %d genes x %d samples", nrow(m), ncol(m))
    m
}

# -----------------------------------------------------------------------------
# Load Protein Matrix (for Predictors)
# -----------------------------------------------------------------------------
#' Load protein quantification matrix from dataset directory
#' @param dir Dataset directory path
#' @return Numeric matrix with genes as rows, samples as columns
load_matrix_from_dataset_dir_protein <- function(dir) {
    fp <- file.path(dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(glue::glue("File not found: {fp}"))
    log_msg("Reading protein matrix (for predictors): {basename(fp)}")

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

    case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)

    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10) log_msg("Note: case_list intersection too small (%d), using all sample columns instead", length(inter))
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
    log_msg("predictor-protein matrix dimensions: %d genes x %d samples", nrow(m), ncol(m))
    m
}

# -----------------------------------------------------------------------------
# Impute and Filter Matrix
# -----------------------------------------------------------------------------
#' Filter genes by non-NA proportion and impute remaining missing values
#' @param mat Numeric matrix
#' @param min_frac Minimum fraction of non-NA values per gene (default 0.75)
#' @return Filtered and imputed matrix
impute_and_filter <- function(mat, min_frac = 0.75) {
    # Filter by gene's non-NA proportion across samples
    keep <- rowMeans(is.finite(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]

    # If still some missing values, use gene-median imputation
    if (anyNA(m)) {
        meds <- apply(m, 1, function(v) median(v[is.finite(v)], na.rm = TRUE))
        for (i in seq_len(nrow(m))) {
            vi <- m[i, ]
            vi[!is.finite(vi)] <- meds[i]
            m[i, ] <- vi
        }
    }
    m
}

# -----------------------------------------------------------------------------
# Write Gene Set Manifest
# -----------------------------------------------------------------------------
#' Export gene set manifest to CSV
#' @param genesets_by_group Named list of gene sets
#' @param out_csv Output file path
#' @return msigdbr package version
write_geneset_manifest <- function(genesets_by_group,
                                   out_csv = file.path("run_info", "geneset_manifest.csv")) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    ver <- tryCatch(as.character(utils::packageVersion("msigdbr")), error = function(e) NA_character_)
    rows <- lapply(names(genesets_by_group), function(g) {
        gs <- genesets_by_group[[g]]
        if (is.null(gs) || !length(gs)) {
            return(NULL)
        }
        data.frame(
            group = g,
            pathway = names(gs),
            genes_n = vapply(gs, function(v) length(unique(v)), integer(1)),
            stringsAsFactors = FALSE
        )
    })
    df <- dplyr::bind_rows(rows)
    if (!is.null(df)) data.table::fwrite(df, out_csv)
    return(ver)
}

# -----------------------------------------------------------------------------
# Ensure Matrix or NULL
# -----------------------------------------------------------------------------
#' Validate and convert to matrix, returning NULL if invalid
#' @param x Input object
#' @return Matrix or NULL
.ensure_mat_or_null <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m)) {
        return(NULL)
    }
    if (is.null(dim(m))) {
        m <- matrix(m, ncol = 1)
    }
    if (ncol(m) == 0) {
        return(NULL)
    }
    m
}
