# =============================================================================
# Utility Helper Functions for CSN CPTAC Phosphosite GSEA Analysis
# =============================================================================
# This file contains general utility functions extracted from the main script.
# =============================================================================

# -----------------------------------------------------------------------------
# opt()
# -----------------------------------------------------------------------------
if (!exists("opt", mode = "function")) {
    opt <- function(nm, default) {
        if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
    }
}

# -----------------------------------------------------------------------------
# coerce_covariates_safely()
# -----------------------------------------------------------------------------
if (!exists("coerce_covariates_safely", mode = "function")) {
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
                    if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
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
}

# -----------------------------------------------------------------------------
# .force_serial_execution()
# -----------------------------------------------------------------------------
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
}

# -----------------------------------------------------------------------------
# log_msg() 
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# %||% operator 
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# safe_fs_name() 
# -----------------------------------------------------------------------------
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}

# -----------------------------------------------------------------------------
# .z_no_impute() 
# Z-score transformation without imputation (keeps NA for non-finite values)
# -----------------------------------------------------------------------------
.z_no_impute <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    fin <- is.finite(v)
    mu <- mean(v[fin], na.rm = TRUE)
    sdv <- stats::sd(v[fin], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out <- (v - mu) / sdv
    out[!fin] <- NA_real_
    out
}

# -----------------------------------------------------------------------------
# .norm_names() 
# Normalize names to uppercase with underscores
# -----------------------------------------------------------------------------
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))

# -----------------------------------------------------------------------------
# .zscore() 
# Z-score transformation with mean imputation for non-finite values
# -----------------------------------------------------------------------------
.zscore <- function(v) {
    v <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    (v - mu) / sdv
}

# -----------------------------------------------------------------------------
# .to01() 
# Convert values to 0-1 range (handles percentage vs fraction)
# -----------------------------------------------------------------------------
.to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
    pmin(pmax(v, 0), 1)
}

# -----------------------------------------------------------------------------
# .median_from_semicolon() 
# Extract median from semicolon-separated values (for PAAD dataset)
# -----------------------------------------------------------------------------
.median_from_semicolon <- function(x_chr) {
    vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
    vv <- vv[is.finite(vv)]
    if (!length(vv)) {
        return(NA_real_)
    }
    stats::median(vv)
}
