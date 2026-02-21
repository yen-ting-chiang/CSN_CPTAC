# =============================================================================
# Data Loading Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for loading data matrices.
# =============================================================================

# -----------------------------------------------------------------------------
# read_case_list()
# -----------------------------------------------------------------------------
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
# load_matrix_from_dataset_dir()
# -----------------------------------------------------------------------------
load_matrix_from_dataset_dir <- function(dir) {
    fp <- file.path(dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(glue::glue("Cannot find file: {fp}"))
    log_msg("Reading protein matrix: {basename(fp)}")
    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))
    gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol", "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol")
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]
    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene)
    not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.", "ENTREZ_GENE_ID", "Description", "Gene_Name", "GeneName", "Gene_Symbol")
    sample_cols_all <- setdiff(names(dat), not_sample)

    case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)
    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10) log_msg("Note: The intersection of case_list is too small ({length(inter)}), use all sample fields instead.")
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
    if (ncol(m) == 0) stop("0 sample columns read")
    if (anyDuplicated(rownames(m))) {
        log_msg("Duplicate genes detected; averaging duplicate rows.")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
    }
    log_msg("Matrix: {nrow(m)} genes x {ncol(m)} samples")
    m
}


# -----------------------------------------------------------------------------
# impute_and_filter()
# -----------------------------------------------------------------------------
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(!is.na(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    if (any(is.na(m))) {
        set.seed(1234)
        m <- imputeLCMD::impute.MinProb(m, q = 0.01)
    }
    m
}


# -----------------------------------------------------------------------------
# load_phosphosite_matrix_from_dataset_dir()
# -----------------------------------------------------------------------------
load_phosphosite_matrix_from_dataset_dir <- function(ds_dir,
                                                     protein_adjust = TRUE,
                                                     min_nonNA_row_frac = 0) {
    fp <- file.path(ds_dir, "data_phosphoprotein_quantification.txt")
    if (!file.exists(fp)) stop("Phosphoprotein file not found: ", fp)
    log_msg("Reading phospho (site-level) matrix: data_phosphoprotein_quantification.txt")

    df <- suppressMessages(vroom::vroom(fp, delim = "\t", col_types = vroom::cols(
        .default = "d",
        ENTITY_STABLE_ID = "c",
        NAME = "c",
        DESCRIPTION = "c",
        GENE_SYMBOL = "c",
        PHOSPHOSITES = "c", PHOSPHOSITE = "c"
    )))
    df <- as.data.frame(df, check.names = FALSE)


    anno_cols <- intersect(c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES", "PHOSPHOSITE"), colnames(df))
    samp_cols <- setdiff(colnames(df), anno_cols)
    stopifnot(length(samp_cols) > 0)


    site_id <- .infer_site_id(df)
    keep <- !is.na(site_id) & nzchar(site_id)
    df <- df[keep, , drop = FALSE]
    site_id <- site_id[keep]

    M_site <- as.matrix(df[, samp_cols, drop = FALSE])
    storage.mode(M_site) <- "numeric"
    rownames(M_site) <- site_id
    colnames(M_site) <- samp_cols


    if (min_nonNA_row_frac > 0) {
        ok <- rowMeans(is.finite(M_site)) >= min_nonNA_row_frac
        M_site <- M_site[ok, , drop = FALSE]
    }


    if (isTRUE(protein_adjust)) {
        prot_fp <- file.path(ds_dir, "data_protein_quantification.txt")
        if (!file.exists(prot_fp)) stop("[phospho-site] protein_adjust=TRUE but protein file not found.")
        log_msg("Site-level protein adjustment: reading protein matrix data_protein_quantification.txt")

        prot_df <- suppressMessages(vroom::vroom(prot_fp, delim = "\t"))
        prot_df <- as.data.frame(prot_df, check.names = FALSE)


        gcol <- which(tolower(colnames(prot_df)) %in% tolower(c("Composite.Element.REF", "Gene", "GENE", "GENE_SYMBOL", "GENE_NAME")))[1]
        if (is.na(gcol)) stop("[phospho-site] protein_adjust=TRUE, but the protein file is missing the gene field.")
        prot_genes <- as.character(prot_df[[gcol]])
        prot_mat <- as.matrix(prot_df[, setdiff(colnames(prot_df), colnames(prot_df)[gcol]), drop = FALSE])
        rownames(prot_mat) <- prot_genes
        storage.mode(prot_mat) <- "numeric"


        gene_of_site <- as.character(df$GENE_SYMBOL)


        common <- intersect(colnames(M_site), colnames(prot_mat))
        if (length(common) >= 2) {
            M_site <- M_site[, common, drop = FALSE]
            prot_mat <- prot_mat[, common, drop = FALSE]


            idx <- match(gene_of_site, rownames(prot_mat))
            hasp <- !is.na(idx)
            if (any(hasp)) {
                M_site[hasp, ] <- M_site[hasp, , drop = FALSE] - prot_mat[idx[hasp], , drop = FALSE]
            }
        } else {
            log_msg("(WARNING) The intersection of phospho and protein samples is < 2; protein correction is skipped.")
        }
    }


    M_site
}
