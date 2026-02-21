## ====================================================================
## 07_msigdb_utils.R
##
## Purpose: MSigDB gene set preparation and manifest generation
## Contains: Gene set loading, collection filtering, manifest writing
## ====================================================================

#' Write gene set manifest to CSV
#'
#' @param genesets_by_group List of gene sets organized by group
#' @param out_csv Output CSV file path
#' @return MSigDB version string or NA
#' @export
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


#' Get subcollection column name from MSigDB data frame
#'
#' @param df MSigDB data frame
#' @return Subcollection column name or NULL
#' @export
.get_subcat_col <- function(df) {
    cand <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
    hit <- cand[cand %in% names(df)]
    if (length(hit)) hit[1] else NULL
}


#' Fetch MSigDB data for a specific category
#'
#' @param cat Category code (e.g., "H", "C2")
#' @return MSigDB data frame or try-error
#' @export
.fetch_msig <- function(cat) {
    try(msigdbr::msigdbr(species = "Homo sapiens", category = toupper(cat)), silent = TRUE)
}


#' Parse gene set group token (e.g., "C2:CP:REACTOME")
#'
#' @param tok Token string
#' @return List with cat, sub, raw components
#' @export
.parse_group_token <- function(tok) {
    t <- toupper(trimws(tok))
    sp <- unlist(strsplit(t, ":", fixed = TRUE))
    list(cat = sp[1], sub = if (length(sp) >= 2) paste(sp[-1], collapse = ":") else NULL, raw = t)
}


#' Create group label from category and subcollection
#'
#' @param cat Category code
#' @param df MSigDB data frame
#' @param want_sub Desired subcollection
#' @return Group label string
#' @export
.make_group_label <- function(cat, df, want_sub) {
    subcol <- .get_subcat_col(df)
    if (!is.null(want_sub)) {
        if (!is.null(subcol)) {
            uniq <- unique(df[[subcol]])
            hit <- uniq[grepl(want_sub, uniq, ignore.case = TRUE)]
            if (length(hit) >= 1) {
                return(sprintf("%s__%s", cat, hit[1]))
            }
        }
        return(sprintf("%s__%s", cat, safe_fs_name(want_sub)))
    }
    return(cat)
}


#' Load and organize MSigDB gene sets
#'
#' @param geneset_groups_to_run Character vector of collection codes
#' @return List of gene sets organized by group
#' @export
load_msigdb_genesets <- function(geneset_groups_to_run = c("H")) {
    if (exists("log_msg", mode = "function")) {
        log_msg(
            "Prepare MSigDB gene sets; collections: %s",
            paste(geneset_groups_to_run, collapse = ", ")
        )
    }

    # Check MSigDB availability and write collection manifest
    df0 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)

    if (!is.null(df0) && NROW(df0) > 0) {
        df0 <- as.data.frame(df0)

        # Detect actual field names
        cat_candidates <- c("gs_collection", "gs_cat", "gs_category", "category", "collection")
        sub_candidates <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
        name_candidates <- c("gs_name", "geneset_name", "set_name")

        cat_hits <- intersect(cat_candidates, names(df0))
        sub_hits <- intersect(sub_candidates, names(df0))
        name_hits <- intersect(name_candidates, names(df0))

        if (length(cat_hits) == 0 || length(name_hits) == 0) {
            if (exists("log_msg", mode = "function")) {
                log_msg("(Omitted) Required field (collection or gs_name) not found")
            }
        } else {
            catcol <- cat_hits[1]
            subcol <- if (length(sub_hits) >= 1) sub_hits[1] else NULL
            namecol <- name_hits[1]

            # Count unique gene sets per collection
            tmp <- data.frame(
                gs_cat = as.character(df0[[catcol]]),
                gs_subcat = if (!is.null(subcol)) as.character(df0[[subcol]]) else NA_character_,
                gs_name = as.character(df0[[namecol]]),
                stringsAsFactors = FALSE
            )

            counts <- stats::aggregate(gs_name ~ gs_cat + gs_subcat,
                data = tmp,
                FUN = function(x) length(unique(x))
            )
            names(counts)[3] <- "n_genesets"

            # Sort
            ord <- order(counts$gs_cat, ifelse(is.na(counts$gs_subcat), "", counts$gs_subcat))
            counts <- counts[ord, , drop = FALSE]

            # Write manifest
            dir.create("run_info", showWarnings = FALSE, recursive = TRUE)
            readr::write_csv(counts, file.path("run_info", "msigdb_available_collections.csv"))
            if (exists("log_msg", mode = "function")) {
                log_msg("Available collections (with counts) → run_info/msigdb_available_collections.csv")
            }
        }
    } else {
        if (exists("log_msg", mode = "function")) {
            log_msg("(Omitted) msigdbr() returns null or call failed")
        }
    }

    # Load requested gene set groups
    genesets_by_group <- list()

    for (tok in geneset_groups_to_run) {
        pp <- .parse_group_token(tok)
        cat <- pp$cat
        df <- suppressMessages(.fetch_msig(cat))
        if (inherits(df, "try-error") || is.null(df) || !nrow(df)) {
            if (exists("log_msg", mode = "function")) {
                log_msg("  Note: msigdbr cannot provide category %s; skip %s", cat, pp$raw)
            }
            next
        }

        if (!is.null(pp$sub)) {
            subcol <- .get_subcat_col(df)
            subtxt <- if (!is.null(subcol)) as.character(df[[subcol]]) else ""
            # Match only at beginning of subclass string
            esc <- function(s) gsub("([\\^\\$\\.|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])", "\\\\\\1", s)
            pattern <- paste0("^", esc(pp$sub), "($|:)")
            keep <- grepl(pattern, subtxt, ignore.case = TRUE)
            df <- df[keep, , drop = FALSE]
            if (!nrow(df)) {
                if (exists("log_msg", mode = "function")) {
                    log_msg("  Note: %s does not find subset %s; skip", cat, pp$sub)
                }
                next
            }
            grp <- .make_group_label(cat, df, pp$sub)
        } else {
            grp <- cat
        }

        genesets_by_group[[grp]] <- lapply(split(df$gene_symbol, df$gs_name), unique)
        if (exists("log_msg", mode = "function")) {
            log_msg("  Loaded %s: %d sets", grp, length(genesets_by_group[[grp]]))
        }
    }

    if (!length(genesets_by_group)) {
        stop("No gene sets available. Check GENESET_GROUPS_TO_RUN setting.")
    }

    if (exists("log_msg", mode = "function")) {
        log_msg("Final gene set groups: %s", paste(names(genesets_by_group), collapse = ", "))
    }

    genesets_by_group
}


# ---- End of 07_msigdb_utils.R ----
