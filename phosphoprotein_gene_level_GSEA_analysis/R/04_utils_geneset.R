# =============================================================================
# 04_utils_geneset.R - Gene Set Loading and Handling
# =============================================================================
# This file contains functions for loading MSigDB gene sets.
# =============================================================================

# -----------------------------------------------------------------------------
# Get Subcategory Column Name from MSigDB Data Frame
# -----------------------------------------------------------------------------
.get_subcat_col <- function(df) {
    cand <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
    hit <- cand[cand %in% names(df)]
    if (length(hit)) hit[1] else NULL
}

# -----------------------------------------------------------------------------
# Fetch MSigDB Gene Sets by Category
# -----------------------------------------------------------------------------
.fetch_msig <- function(cat) {
    try(msigdbr::msigdbr(species = "Homo sapiens", category = toupper(cat)), silent = TRUE)
}

# -----------------------------------------------------------------------------
# Parse Gene Set Group Token 
# -----------------------------------------------------------------------------
.parse_group_token <- function(tok) {
    t <- toupper(trimws(tok))
    sp <- unlist(strsplit(t, ":", fixed = TRUE))
    list(cat = sp[1], sub = if (length(sp) >= 2) paste(sp[-1], collapse = ":") else NULL, raw = t)
}

# -----------------------------------------------------------------------------
# Make Group Label from Category and Subcategory
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Build Gene Sets by Group (Main Entry Point)
# -----------------------------------------------------------------------------
build_genesets_by_group <- function(groups_to_run) {
    log_msg(
        "Preparing MSigDB gene sets; will build according to GENESET_GROUPS_TO_RUN: %s",
        paste(groups_to_run, collapse = ", ")
    )

    genesets_by_group <- list()

    for (tok in groups_to_run) {
        pp <- .parse_group_token(tok)
        cat <- pp$cat
        df <- suppressMessages(.fetch_msig(cat))
        if (inherits(df, "try-error") || is.null(df) || !nrow(df)) {
            log_msg("  Note: msigdbr cannot provide category %s; skipping %s", cat, pp$raw)
            next
        }

        if (!is.null(pp$sub)) {
            subcol <- .get_subcat_col(df)
            subtxt <- if (!is.null(subcol)) as.character(df[[subcol]]) else ""
            # Only match at the "beginning" of subcategory string 
            # And allow followed by ":" or end
            esc <- function(s) gsub("([\\^\\$\\.|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])", "\\\\\\1", s)
            pattern <- paste0("^", esc(pp$sub), "($|:)")
            keep <- grepl(pattern, subtxt, ignore.case = TRUE)
            df <- df[keep, , drop = FALSE]
            if (!nrow(df)) {
                log_msg("  Note: %s cannot find any subcollection for %s; skipping", cat, pp$sub)
                next
            }
            grp <- .make_group_label(cat, df, pp$sub)
        } else {
            grp <- cat
        }

        genesets_by_group[[grp]] <- lapply(split(df$gene_symbol, df$gs_name), unique)
        log_msg("  Obtained %s: %d sets", grp, length(genesets_by_group[[grp]]))
    }

    if (!length(genesets_by_group)) {
        stop("No available gene-set. Please check GENESET_GROUPS_TO_RUN settings.")
    }
    log_msg("Final gene-set groups: %s", paste(names(genesets_by_group), collapse = ", "))

    return(genesets_by_group)
}

# -----------------------------------------------------------------------------
# Output MSigDB Available Collections (Optional Diagnostic)
# -----------------------------------------------------------------------------
output_msigdb_collections <- function() {
    df0 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)

    if (!is.null(df0) && NROW(df0) > 0) {
        df0 <- as.data.frame(df0) # avoid tibble strict subsetting behavior

        # Detect actual column names (2025.1 uses gs_collection / gs_subcollection)
        cat_candidates <- c("gs_collection", "gs_cat", "gs_category", "category", "collection")
        sub_candidates <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
        name_candidates <- c("gs_name", "geneset_name", "set_name")

        cat_hits <- intersect(cat_candidates, names(df0))
        sub_hits <- intersect(sub_candidates, names(df0))
        name_hits <- intersect(name_candidates, names(df0))

        if (length(cat_hits) == 0 || length(name_hits) == 0) {
            log_msg("(Skipped) Required columns not found (collection or gs_name); collections list not output.")
        } else {
            catcol <- cat_hits[1]
            subcol <- if (length(sub_hits) >= 1) sub_hits[1] else NULL
            namecol <- name_hits[1]

            # Build data frame and count unique geneset (gs_name) counts
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

            # Sort (treat NA subcategory as empty string to avoid order encountering NA)
            ord <- order(counts$gs_cat, ifelse(is.na(counts$gs_subcat), "", counts$gs_subcat))
            counts <- counts[ord, , drop = FALSE]

            # Output list to run_info
            dir.create("run_info", showWarnings = FALSE, recursive = TRUE)
            readr::write_csv(counts, file.path("run_info", "msigdb_available_collections.csv"))
            log_msg("Output msigdb available collections (with geneset counts) to run_info/msigdb_available_collections.csv")
        }
    } else {
        log_msg("(Skipped) msigdbr() returned empty or call failed; collections list not output.")
    }
}
