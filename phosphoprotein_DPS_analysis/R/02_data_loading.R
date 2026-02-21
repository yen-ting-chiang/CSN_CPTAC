## =========================================================
## 02_data_loading.R
## Data Loading Functions for DPS Analysis
## =========================================================

## ===== Read case list from file =====
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

## ===== Load protein quantification matrix =====
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
    # If directory has case_list, try to match
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

## ===== Load phosphoprotein matrix  =====
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

    # Compatible columns
    meta_cols <- intersect(
        colnames(df),
        c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES")
    )
    sample_cols <- setdiff(colnames(df), meta_cols)
    if (!length(sample_cols)) stop("[phospho] Cannot find sample columns")

    # Site x sample matrix
    M_site <- as.matrix(df[, sample_cols, drop = FALSE])
    storage.mode(M_site) <- "numeric"
    rownames(M_site) <- df$GENE_SYMBOL

    # Aggregate to gene x sample
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

    # Optional: adjust by protein abundance
    if (isTRUE(protein_adjust) && file.exists(protein_fp)) {
        suppressMessages({
            prot_df <- readr::read_tsv(protein_fp, show_col_types = FALSE, progress = FALSE)
        })
        prot_meta <- intersect(
            colnames(prot_df),
            c(
                "Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "Description",
                "GENE", "GENE_ID", "GENE.STABLE.ID", "Composite.Element.REF"
            )
        )
        prot_sample <- setdiff(colnames(prot_df), prot_meta)
        cand <- c("Hugo_Symbol", "GENE_SYMBOL", "GeneSymbol", "GENE", "Composite.Element.REF")
        cn <- colnames(prot_df)
        lc <- tolower(cn)
        hit <- match(tolower(cand), lc, nomatch = NA_integer_)
        hit <- hit[!is.na(hit)]
        gcol <- if (length(hit)) cn[hit[1]] else NA_character_
        if (is.na(gcol)) stop("[phospho] protein_adjust=TRUE but protein file missing gene column")
        prot_sample <- setdiff(prot_sample, gcol)
        M_prot <- as.matrix(prot_df[, prot_sample, drop = FALSE])
        storage.mode(M_prot) <- "numeric"
        rownames(M_prot) <- prot_df[[gcol]]
        rownames(M_prot) <- sub("\\|.*$", "", rownames(M_prot))

        # Gene and sample intersection
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

    # Sample name cleanup
    colnames(M_gene) <- gsub("\\s+", "", colnames(M_gene))
    if (exists("log_msg")) log_msg("[phospho] Matrix dimensions: %d genes x %d samples", nrow(M_gene), ncol(M_gene))
    return(M_gene)
}

## ===== PTMsigDB helper functions =====
.std_site_id <- function(x) {
    # Standardize to GENE_S123 / GENE_T123 / GENE_Y123
    x <- toupper(as.character(x))
    x <- gsub("[:\\s]", "_", x)
    x <- gsub("_([STY])(\\d+)[A-Z]*$", "_\\1\\2", x)
    x
}

.make_site_id <- function(gene, name, phosphosites) {
    gene <- as.character(gene)
    name <- as.character(name)
    phosphosites <- as.character(phosphosites)

    row_fun <- function(g, n, p) {
        cand <- NA_character_
        if (!is.na(n) && nzchar(n) && grepl("_[STY]\\d+", n, ignore.case = TRUE)) cand <- n
        if (is.na(cand) || !nzchar(cand)) {
            ps <- NA_character_
            if (!is.na(p) && nzchar(p)) {
                m <- regmatches(p, regexpr("([STY])(\\d+)", p, ignore.case = TRUE))
                if (length(m) == 1 && nzchar(m)) ps <- m
            }
            if (!is.na(g) && nzchar(g) && !is.na(ps) && nzchar(ps)) cand <- paste0(g, "_", ps)
        }
        .std_site_id(cand)
    }

    vapply(seq_along(gene), function(i) row_fun(gene[i], name[i], phosphosites[i]), character(1))
}

.infer_site_id <- function(df) {
    g <- if ("GENE_SYMBOL" %in% names(df)) as.character(df$GENE_SYMBOL) else rep(NA_character_, nrow(df))
    nm <- if ("NAME" %in% names(df)) as.character(df$NAME) else rep(NA_character_, nrow(df))
    p1 <- if ("PHOSPHOSITES" %in% names(df)) as.character(df$PHOSPHOSITES) else rep(NA_character_, nrow(df))
    p2 <- if ("PHOSPHOSITE" %in% names(df)) as.character(df$PHOSPHOSITE) else rep(NA_character_, nrow(df))
    ent <- if ("ENTITY_STABLE_ID" %in% names(df)) as.character(df$ENTITY_STABLE_ID) else rep(NA_character_, nrow(df))

    # 1) from NAME
    id_from_name <- ifelse(!is.na(nm) & grepl("_[STY]\\d+", nm, ignore.case = TRUE),
        sub(":.*$", "", nm), NA_character_
    )

    # helper: extract first site token
    site_token <- function(x) {
        y <- NA_character_
        if (!is.na(x) && nzchar(x)) {
            m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
            if (length(m) == 1 && nzchar(m)) y <- toupper(m)
        }
        y
    }

    # 2) from GENE + PHOSPHOSITES/PHOSPHOSITE
    ps <- vapply(seq_along(p1), function(i) {
        s <- if (!is.na(p1[i]) && nzchar(p1[i])) p1[i] else p2[i]
        site_token(s)
    }, character(1))
    id_from_gp <- ifelse(!is.na(g) & nzchar(g) & !is.na(ps) & nzchar(ps),
        paste0(g, "_", ps), NA_character_
    )

    # 3) from ENTITY_STABLE_ID
    ent_base <- toupper(sub(":.*$", "", ent))
    ent_gene <- ifelse(!is.na(g) & nzchar(g), toupper(g), sub("_.*$", "", ent_base))
    ent_site <- vapply(ent_base, function(x) {
        m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
        if (length(m) == 1 && nzchar(m)) toupper(m) else NA_character_
    }, character(1))
    id_from_ent <- ifelse(!is.na(ent_gene) & nzchar(ent_gene) & !is.na(ent_site) & nzchar(ent_site),
        paste0(ent_gene, "_", ent_site), NA_character_
    )

    # Backfill in order
    cand <- id_from_name
    need <- is.na(cand) | !nzchar(cand)
    cand[need] <- id_from_gp[need]
    need <- is.na(cand) | !nzchar(cand)
    cand[need] <- id_from_ent[need]

    .std_site_id(cand)
}

## ===== Read PTMsigDB GMT =====
read_ptmsigdb_gmt <- function(fp) {
    if (missing(fp) || !nzchar(fp) || !file.exists(fp)) {
        stop("[PTMsigDB] File does not exist: ", fp %||% "<empty>")
    }

    std_id <- function(x) {
        x <- toupper(trimws(x))
        x <- sub(":.*$", "", x)
        x <- gsub("[^A-Z0-9_]", "", x)
        x
    }

    # If .xlsx
    if (grepl("\\.xlsx$", fp, ignore.case = TRUE)) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
            stop("[PTMsigDB] Package readxl required to read .xlsx")
        }
        df <- readxl::read_xlsx(fp)
        req <- c("signature", "site.annotation")
        if (!all(req %in% names(df))) {
            stop("[PTMsigDB] .xlsx missing required columns: ", paste(setdiff(req, names(df)), collapse = ", "))
        }
        df$gene_site <- std_id(df$site.annotation)
        by_sig <- split(df$gene_site, df$signature)
        out <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
        out[lengths(out) == 0] <- NULL
        return(out)
    }

    # If .gmt
    lines <- readr::read_lines(fp)
    res <- vector("list", length(lines))
    nm <- character(length(lines))

    for (i in seq_along(lines)) {
        fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
        if (length(fields) < 2) next
        sig <- fields[1]
        nm[i] <- sig

        toks <- if (length(fields) >= 3) fields[-(1:2)] else character(0)
        if (length(toks) == 0 && grepl("\\|", fields[2], fixed = TRUE)) {
            t2 <- strsplit(fields[2], "\\|")[[1]]
            toks <- t2[!grepl("^n=\\d+$", t2, ignore.case = TRUE)]
        }

        toks <- std_id(toks)
        keep <- grepl("^[A-Z0-9._-]+_[STY][0-9]+$", toks)
        toks <- unique(toks[keep & nzchar(toks)])
        res[[i]] <- toks
    }

    names(res) <- nm
    res[lengths(res) == 0] <- NULL
    if (!length(res)) {
        stop(
            "[PTMsigDB] Parsing result is empty. This GMT is likely UniProt version (e.g., O14974;T696).",
            " Please use source with gene-symbol (.xlsx or gene version GMT)."
        )
    }
    res
}

## ===== Load site-level phosphosite matrix =====
load_phosphosite_matrix_from_dataset_dir <- function(ds_dir,
                                                     protein_adjust = TRUE,
                                                     min_nonNA_row_frac = 0) {
    fp <- file.path(ds_dir, "data_phosphoprotein_quantification.txt")
    if (!file.exists(fp)) stop("Cannot find phosphoprotein file: ", fp)
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

    # Detect annotation columns + sample columns
    anno_cols <- intersect(c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES", "PHOSPHOSITE"), colnames(df))
    samp_cols <- setdiff(colnames(df), anno_cols)
    stopifnot(length(samp_cols) > 0)

    # Generate site_id and normalize
    site_id <- .infer_site_id(df)
    keep <- !is.na(site_id) & nzchar(site_id)
    df <- df[keep, , drop = FALSE]
    site_id <- site_id[keep]

    M_site <- as.matrix(df[, samp_cols, drop = FALSE])
    storage.mode(M_site) <- "numeric"
    rownames(M_site) <- site_id
    colnames(M_site) <- samp_cols

    # Filter by minimum non-missing ratio
    if (min_nonNA_row_frac > 0) {
        ok <- rowMeans(is.finite(M_site)) >= min_nonNA_row_frac
        M_site <- M_site[ok, , drop = FALSE]
    }

    # Protein adjustment
    if (isTRUE(protein_adjust)) {
        prot_fp <- file.path(ds_dir, "data_protein_quantification.txt")
        if (!file.exists(prot_fp)) stop("[phospho-site] protein_adjust=TRUE but protein file not found")
        log_msg("Site-level protein adjustment: reading protein matrix: data_protein_quantification.txt")

        prot_df <- suppressMessages(vroom::vroom(prot_fp, delim = "\t"))
        prot_df <- as.data.frame(prot_df, check.names = FALSE)

        gcol <- which(tolower(colnames(prot_df)) %in% tolower(c("Composite.Element.REF", "Gene", "GENE", "GENE_SYMBOL", "GENE_NAME")))[1]
        if (is.na(gcol)) stop("[phospho-site] protein_adjust=TRUE but protein file missing gene column")
        prot_genes <- as.character(prot_df[[gcol]])
        prot_mat <- as.matrix(prot_df[, setdiff(colnames(prot_df), colnames(prot_df)[gcol]), drop = FALSE])
        rownames(prot_mat) <- prot_genes
        storage.mode(prot_mat) <- "numeric"

        gene_of_site <- as.character(df$GENE_SYMBOL)
        gene_of_site <- gene_of_site[keep]

        # Subtract protein from phosphosite
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
            log_msg("(Warning) phospho and protein sample intersection < 2; skipping protein adjustment")
        }
    }

    M_site
}
