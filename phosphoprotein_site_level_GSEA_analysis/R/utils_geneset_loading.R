# =============================================================================
# Gene Set Loading Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains functions for loading PTMsigDB gene sets.
# =============================================================================


.std_site_id <- function(x) {
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

    # helper
    site_token <- function(x) {
        y <- NA_character_
        if (!is.na(x) && nzchar(x)) {
            m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
            if (length(m) == 1 && nzchar(m)) y <- toupper(m)
        }
        y
    }

    # 2) from GENE + (PHOSPHOSITES/PHOSPHOSITE)
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


    cand <- id_from_name
    need <- is.na(cand) | !nzchar(cand)
    cand[need] <- id_from_gp[need]
    need <- is.na(cand) | !nzchar(cand)
    cand[need] <- id_from_ent[need]

    .std_site_id(cand)
}


read_ptmsigdb_gmt <- function(fp) {
    if (missing(fp) || !nzchar(fp) || !file.exists(fp)) {
        stop("[PTMsigDB] File does not exist:", fp %||% "<empty>")
    }

    std_id <- function(x) {
        x <- toupper(trimws(x))
        x <- sub(":.*$", "", x)
        x <- gsub("[^A-Z0-9_]", "", x)
        x
    }


    if (grepl("\\.xlsx$", fp, ignore.case = TRUE)) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
            stop("[PTMsigDB] requires the readxl package to read .xlsx files.")
        }
        df <- readxl::read_xlsx(fp)
        req <- c("signature", "site.annotation")
        if (!all(req %in% names(df))) {
            stop("[PTMsigDB].xlsx is missing necessary fields:", paste(setdiff(req, names(df)), collapse = ", "))
        }

        df$gene_site <- std_id(df$site.annotation)
        by_sig <- split(df$gene_site, df$signature)
        out <- lapply(by_sig, function(v) unique(v[nzchar(v)]))
        out[lengths(out) == 0] <- NULL
        return(out)
    }


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
            "The [PTMsigDB] parsing result is empty. This GMT is likely a UniProt version (e.g., O14974; T696).",
            " Please use a source with gene-symbol instead (.xlsx or gene version of GMT)."
        )
    }
    res
}
