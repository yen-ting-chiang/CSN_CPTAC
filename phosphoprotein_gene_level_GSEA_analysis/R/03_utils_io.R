# =============================================================================
# 03_utils_io.R - File I/O and Data Loading Utilities
# =============================================================================
# This file contains functions for reading data files and writing outputs.
# =============================================================================

# -----------------------------------------------------------------------------
# Read Case List from cBioPortal Format
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
# Load Protein Quantification Matrix from Dataset Directory
# -----------------------------------------------------------------------------
load_matrix_from_dataset_dir <- function(dir) {
    fp <- file.path(dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(glue::glue("File not found: {fp}"))
    log_msg("Reading protein matrix: {basename(fp)}")
    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))
    gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol", "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol")
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]
    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene)
    not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.", "ENTREZ_GENE_ID", "Description", "Gene_Name", "GeneName", "Gene_Symbol")
    sample_cols_all <- setdiff(names(dat), not_sample)
    # If folder has case_list, try to match as much as possible
    case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)
    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
        if (length(inter) < 10) log_msg("Hint: case_list intersection too small ({length(inter)}), using all sample columns")
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
# Load Phosphoprotein Quantification Matrix from Dataset Directory
# -----------------------------------------------------------------------------
load_phospho_matrix_from_dataset_dir <- function(dir,
                                                 site_collapse = c("median", "max"),
                                                 protein_adjust = FALSE,
                                                 protein_fp = file.path(dir, "data_protein_quantification.txt")) {
    site_collapse <- match.arg(site_collapse)
    fp <- file.path(dir, "data_phosphoprotein_quantification.txt")
    if (!file.exists(fp)) stop("[phospho] File does not exist: ", fp)

    suppressMessages({
        df <- readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE)
    })

    # Compatible columns: ENTITY_STABLE_ID, NAME, DESCRIPTION, GENE_SYMBOL, PHOSPHOSITES, <samples...>
    meta_cols <- intersect(
        colnames(df),
        c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES")
    )
    sample_cols <- setdiff(colnames(df), meta_cols)
    if (!length(sample_cols)) stop("[phospho] Cannot find sample columns")

    # --- site x sample matrix (GENE_SYMBOL as group key), values are already log2 ratio (CPTAC standard)
    M_site <- as.matrix(df[, sample_cols, drop = FALSE])
    storage.mode(M_site) <- "numeric"
    rownames(M_site) <- df$GENE_SYMBOL

    # --- Summarize to gene x sample (consistent with protein version interface)
    split_idx <- split(seq_len(nrow(M_site)), rownames(M_site))
    agg_fun <- if (site_collapse == "median") {
        function(x) stats::median(x, na.rm = TRUE)
    } else {
        function(x) suppressWarnings(max(x, na.rm = TRUE))
    }
    M_gene <- do.call(rbind, lapply(split_idx, function(ii) {
        apply(M_site[ii, , drop = FALSE], 2, agg_fun)
    }))
    M_gene <- M_gene[order(rownames(M_gene)), , drop = FALSE]

    # --- Optional: adjust by protein abundance (remove same-gene protein level effect)
    if (isTRUE(protein_adjust) && file.exists(protein_fp)) {
        suppressMessages({
            prot_df <- readr::read_tsv(protein_fp, show_col_types = FALSE, progress = FALSE)
        })
        # Protein file sample columns: exclude common annotation columns
        prot_meta <- intersect(
            colnames(prot_df),
            c(
                "Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "Description",
                "GENE", "GENE_ID", "GENE.STABLE.ID", "Composite.Element.REF"
            )
        )
        prot_sample <- setdiff(colnames(prot_df), prot_meta)
        # Try to find gene symbol column (including Composite.Element.REF; case-insensitive)
        cand <- c("Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "GENE", "Composite.Element.REF")
        cn <- colnames(prot_df)
        lc <- tolower(cn)
        hit <- match(tolower(cand), lc, nomatch = NA_integer_)
        hit <- hit[!is.na(hit)]
        gcol <- if (length(hit)) cn[hit[1]] else NA_character_
        if (is.na(gcol)) stop("[phospho] protein_adjust=TRUE but protein file lacks gene column")
        # Ensure sample columns do not contain gene column
        prot_sample <- setdiff(prot_sample, gcol)
        M_prot <- as.matrix(prot_df[, prot_sample, drop = FALSE])
        storage.mode(M_prot) <- "numeric"
        rownames(M_prot) <- prot_df[[gcol]]
        rownames(M_prot) <- sub("\\|.*$", "", rownames(M_prot))

        # Intersection of genes and samples
        g_common <- intersect(rownames(M_gene), rownames(M_prot))
        s_common <- intersect(colnames(M_gene), colnames(M_prot))
        if (length(g_common) >= 1 && length(s_common) >= 3) {
            Y <- M_gene[g_common, s_common, drop = FALSE]
            X <- M_prot[g_common, s_common, drop = FALSE]
            Y_adj <- Y
            for (g in g_common) {
                y <- as.numeric(Y[g, ])
                x <- as.numeric(X[g, ])
                ok <- is.finite(y) & is.finite(x)
                if (sum(ok) >= 3) {
                    fit <- stats::lm(y[ok] ~ x[ok])
                    res <- rep(NA_real_, length(y))
                    res[ok] <- stats::residuals(fit)
                    Y_adj[g, ] <- res
                } else {
                    Y_adj[g, ] <- y
                }
            }
            M_gene[g_common, s_common] <- Y_adj
            attr(M_gene, "protein_adjusted") <- TRUE
        }
    }

    # Sample name minor cleanup (consistent with protein version)
    colnames(M_gene) <- gsub("\\s+", "", colnames(M_gene))
    if (exists("log_msg")) log_msg("[phospho] Matrix dimensions: %d genes x %d samples", nrow(M_gene), ncol(M_gene))
    return(M_gene)
}

# -----------------------------------------------------------------------------
# Write Gene Set Manifest to CSV
# -----------------------------------------------------------------------------
write_geneset_manifest <- function(genesets_by_group, out_csv = file.path("run_info", "geneset_manifest.csv")) {
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
# Clean Matrix for PCA (remove Inf->NA, impute NA with row median)
# -----------------------------------------------------------------------------
.clean_for_pca <- function(X, min_samples = 10L, min_genes = 5L) {
    X <- as.matrix(X)
    X[!is.finite(X)] <- NA
    keep_rows <- rowSums(is.finite(X)) >= min_samples
    if (!any(keep_rows)) {
        return(NULL)
    }
    X <- X[keep_rows, , drop = FALSE]
    if (nrow(X) < min_genes) {
        return(NULL)
    }
    if (anyNA(X)) {
        med <- apply(X, 1, function(r) median(r[is.finite(r)], na.rm = TRUE))
        for (i in seq_len(nrow(X))) {
            xi <- X[i, ]
            xi[!is.finite(xi)] <- med[i]
            X[i, ] <- xi
        }
    }
    X
}

# -----------------------------------------------------------------------------
# Impute and Filter Matrix
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
# Align Vector/Data Frame to Sample Order
# -----------------------------------------------------------------------------
.align_to_samples <- function(x, sam, what = "covariate") {
    if (is.null(x)) {
        return(NULL)
    }
    if (is.vector(x) || is.factor(x)) {
        # Allow length exactly equal to sam without names; otherwise use names to align
        if (is.null(names(x))) {
            if (length(x) != length(sam)) stop(sprintf("[%s] Length mismatch: %d vs %d", what, length(x), length(sam)))
            names(x) <- sam
        }
        if (!all(sam %in% names(x))) stop(sprintf("[%s] Missing samples: %s", what, paste(setdiff(sam, names(x)), collapse = ", ")))
        out <- x[sam]
        return(out)
    } else {
        x <- as.data.frame(x)
        # Expect rownames(x) to be samples
        if (is.null(rownames(x))) {
            if (nrow(x) != length(sam)) stop(sprintf("[%s] Row count mismatch: %d vs %d and no rownames for alignment", what, nrow(x), length(sam)))
            rownames(x) <- sam
        }
        if (!all(sam %in% rownames(x))) stop(sprintf("[%s] Missing samples: %s", what, paste(setdiff(sam, rownames(x)), collapse = ", ")))
        out <- x[sam, , drop = FALSE]
        return(out)
    }
}

# -----------------------------------------------------------------------------
# Convert Object to Matrix or NULL
# -----------------------------------------------------------------------------
.mk_mat_or_null <- function(x, rows) {
    if (is.null(x)) {
        return(NULL)
    }
    m <- as.matrix(x)
    if (is.null(rownames(m))) rownames(m) <- rows
    if (ncol(m) == 0 || nrow(m) != length(rows)) {
        return(NULL)
    }
    m
}
