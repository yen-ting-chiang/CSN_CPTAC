## ====================================================================
## 01_utility_functions.R
##
## Purpose: Utility and helper functions for CSN CPTAC GSEA analysis
## Contains: logging, file operations, safe naming, operators, etc.
## ====================================================================

# ---- Logging Functions ----

#' Log message with timestamp
#'
#' @param text Message text (supports sprintf and glue formats)
#' @param ... Additional arguments for sprintf/glue
#' @param .envir Environment for glue evaluation
#' @return NULL (prints to console)
#' @export
log_msg <- function(text, ..., .envir = parent.frame()) {
    ts <- format(Sys.time(), "%H:%M:%S")
    msg <- tryCatch(
        {
            if (grepl("\\{[^}]+\\}", text)) {
                glue::glue(text, ..., .envir = .envir)
            } else if (grepl("%", text)) {
                do.call(sprintf, c(list(fmt = text), list(...)))
            } else {
                if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
            }
        },
        error = function(e) text
    )
    cat(sprintf("[%s] %s\n", ts, msg))
}


# ---- Operators ----

#' Fallback operator (returns b if a is NULL)
#'
#' @param a First value
#' @param b Fallback value
#' @return a if not NULL, otherwise b
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Safe option getter
#'
#' @param nm Variable name to get
#' @param default Default value if variable doesn't exist
#' @return Value of variable or default
#' @export
opt <- function(nm, default) {
    if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
}


# ---- File and String Operations ----

#' Create safe filesystem name
#'
#' @param s String to convert
#' @return Filesystem-safe string
#' @export
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}


#' Normalize column names (uppercase, remove special chars)
#'
#' @param x Character vector
#' @return Normalized names
#' @export
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))


# ---- Data Conversion Functions ----

#' Convert values to 0-1 range (handles percentage inputs)
#'
#' @param v Numeric vector
#' @return Values clipped to [0,1]
#' @export
.to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    # If most values > 1, assume percentage (divide by 100)
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
    pmin(pmax(v, 0), 1)
}


#' Extract median from semicolon-separated string
#'
#' @param x_chr Character string with values separated by ";"
#' @return Median value or NA
#' @export
.median_from_semicolon <- function(x_chr) {
    vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
    vv <- vv[is.finite(vv)]
    if (!length(vv)) {
        return(NA_real_)
    }
    stats::median(vv)
}


# ---- Statistical Helper Functions ----

#' Z-score transformation (NA preserved)
#'
#' @param x Numeric vector
#' @return Z-scores with NA preserved
#' @export
.z_no_impute <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    fin <- is.finite(x)
    if (sum(fin) < 2) {
        return(rep(NA_real_, length(x)))
    }

    mu <- mean(x[fin], na.rm = TRUE)
    sdv <- stats::sd(x[fin], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1

    out <- rep(NA_real_, length(x))
    out[fin] <- (x[fin] - mu) / sdv
    out
}


#' Z-score transformation (NA imputed with mean)
#'
#' @param x Numeric vector
#' @return Z-scores with NA imputed to mean (z=0)
#' @export
.zscore <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    fin <- is.finite(x)
    if (sum(fin) < 2) {
        return(rep(0, length(x)))
    }

    mu <- mean(x[fin], na.rm = TRUE)
    sdv <- stats::sd(x[fin], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1

    x[!fin] <- mu # Impute NA with mean
    (x - mu) / sdv
}


# ---- Batch Operations ----

#' Format batch factor sizes as string
#'
#' @param fac Factor variable
#' @return String with level counts (e.g., "level1=10; level2=8")
#' @export
.format_batch_sizes <- function(fac) {
    if (is.null(fac)) {
        return(NA_character_)
    }
    tb <- sort(table(fac), decreasing = TRUE)
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}


# ---- Covariate Safety Functions ----

#' Coerce covariates to appropriate types and drop single-level factors
#'
#' @param df Data frame with covariates
#' @return Cleaned data frame
#' @export
coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)

    for (cn in colnames(df)) {
        v <- df[[cn]]
        if (is.factor(v) || is.character(v) || is.logical(v)) {
            v <- factor(v)
            lv <- levels(droplevels(v[!is.na(v)]))
            if (length(lv) <= 1) {
                keep[cn] <- FALSE
                if (exists("log_msg", mode = "function")) {
                    try(log_msg("  [covars] drop single-level factor: %s", cn), silent = TRUE)
                }
            } else {
                df[[cn]] <- v
            }
        } else {
            df[[cn]] <- suppressWarnings(as.numeric(v))
        }
    }
    df <- df[, keep, drop = FALSE]
    df
}


# ---- GSEA Helper Functions ----

#' Ensure stats vector has correct gene names
#'
#' @param stats Numeric vector of statistics
#' @param gene_names Character vector of gene names
#' @param label Optional label for error messages
#' @return Named numeric vector
#' @export
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
    if (is.null(stats)) {
        return(NULL)
    }
    v <- suppressWarnings(as.numeric(stats))
    if (length(v) != length(gene_names)) {
        stop(sprintf(
            "[ensure-names%s] stats(%d) and gene_names(%d) length mismatch",
            if (!is.null(label)) paste0("-", label) else "",
            length(v), length(gene_names)
        ))
    }
    names(v) <- as.character(gene_names)
    v
}


#' Clean and sort GSEA stats (finite values, remove duplicates)
#'
#' @param stats Named numeric vector
#' @param gene_names Optional gene names
#' @param label Label for logging
#' @param min_n Minimum number of genes required
#' @return Sorted named numeric vector or NULL
#' @export
._finite_rank_stats <- function(stats, gene_names = NULL, label = "", min_n = 20L) {
    nm <- names(stats)
    stats <- suppressWarnings(as.numeric(stats))
    if (is.null(nm) && !is.null(gene_names) && length(stats) == length(gene_names)) {
        nm <- gene_names
    }
    names(stats) <- nm

    ok <- is.finite(stats) & !is.na(stats)
    if (sum(!ok) > 0L) {
        if (!is.null(label) && nzchar(label)) {
            if (exists("log_msg", mode = "function")) {
                log_msg("[gsea-%s] drop non-finite stats: %d", label, sum(!ok))
            }
        }
    }
    stats <- stats[ok]
    nm <- nm[ok]
    names(stats) <- nm

    # Names must exist
    if (is.null(nm)) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[gsea-%s] stats has no names after filtering → skip", label)
        }
        return(NULL)
    }

    # Remove duplicate names (keep first)
    dup <- duplicated(nm)
    if (any(dup)) {
        if (!is.null(label) && nzchar(label)) {
            if (exists("log_msg", mode = "function")) {
                log_msg("[gsea-%s] drop duplicated gene names in stats: %d", label, sum(dup))
            }
        }
        stats <- stats[!dup]
        nm <- nm[!dup]
    }
    names(stats) <- nm

    if (length(stats) < min_n) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[gsea-%s] too few finite stats after filtering: %d < %d → skip", label, length(stats), min_n)
        }
        return(NULL)
    }
    stats[order(stats, decreasing = TRUE)]
}


#' Intersect pathways with universe and filter by size
#'
#' @param pathways List of gene sets
#' @param universe Character vector of genes
#' @param minSize Minimum pathway size
#' @param maxSize Maximum pathway size
#' @param label Label for logging
#' @return Filtered list of pathways or NULL
#' @export
._intersect_and_filter_pathways <- function(pathways, universe, minSize, maxSize, label = "") {
    pw <- lapply(pathways, function(gs) unique(intersect(gs, universe)))
    lens <- vapply(pw, length, integer(1))
    keep <- which(lens >= minSize & lens <= maxSize)
    if (length(keep) == 0L) {
        if (exists("log_msg", mode = "function")) {
            log_msg("[gsea-%s] no pathways within size bounds after intersect (min=%d, max=%d)", label, minSize, maxSize)
        }
        return(NULL)
    }
    pw[keep]
}


# ---- Serialization Control ----

#' Force single-threaded execution across all packages
#'
#' @return NULL (side effects: sets options)
#' @export
.force_serial_execution <- function() {
    # base R / data.table
    options(mc.cores = 1L)
    if ("package:data.table" %in% search() &&
        exists("setDTthreads", asNamespace("data.table"))) {
        data.table::setDTthreads(1L)
    }

    # BiocParallel
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::register(BiocParallel::SerialParam())
    }

    # foreach / doParallel
    if (requireNamespace("foreach", quietly = TRUE)) {
        foreach::registerDoSEQ()
    }

    # future
    if (requireNamespace("future", quietly = TRUE)) {
        future::plan(future::sequential)
    }
    Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")

    options(.fgsea_nproc = 1L)

    invisible(NULL)
}


# ---- Design Alignment Audit ----

#' Audit alignment between matrix samples and design matrices
#'
#' @param tag Label for this audit
#' @param samples Character vector of sample IDs
#' @param mod_interest Design matrix for interest variables
#' @param mod_nuisance Design matrix for nuisance variables
#' @param out_dir Optional output directory for CSV
#' @return Invisible data frame with alignment results
#' @export
audit_design_alignment <- function(tag, samples, mod_interest, mod_nuisance, out_dir = NULL) {
    safe_rn <- function(X) if (is.null(X)) rep(NA_character_, length(samples)) else rownames(as.matrix(X))
    df <- data.frame(
        pos = seq_along(samples),
        M_sample = samples,
        interest_rn = safe_rn(mod_interest),
        nuisance_rn = safe_rn(mod_nuisance),
        interest_ok = if (!is.null(mod_interest)) safe_rn(mod_interest) == samples else NA,
        nuisance_ok = if (!is.null(mod_nuisance)) safe_rn(mod_nuisance) == samples else NA,
        stringsAsFactors = FALSE
    )
    n_i <- sum(df$interest_ok, na.rm = TRUE)
    n_n <- sum(df$nuisance_ok, na.rm = TRUE)

    if (exists("log_msg", mode = "function")) {
        log_msg(
            "[align:%s] interest match = %d/%d, nuisance match = %d/%d, both = %d",
            tag, n_i, length(samples), n_n, length(samples), sum(df$interest_ok & df$nuisance_ok, na.rm = TRUE)
        )
    }

    if (!is.null(out_dir)) {
        fn <- file.path(out_dir, sprintf("ALIGN_%s.csv", gsub("[^A-Za-z0-9]+", "_", tag)))
        utils::write.csv(df, fn, row.names = FALSE)
        if (exists("log_msg", mode = "function")) {
            log_msg("[align:%s] wrote %s", tag, fn)
        }
    }

    head_show <- utils::head(df[, c("pos", "M_sample", "interest_rn", "nuisance_rn", "interest_ok", "nuisance_ok")], 6)
    if (exists("log_msg", mode = "function")) {
        log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse = "\n"))
    }

    invisible(df)
}


# ---- End of 01_utility_functions.R ----
