# =============================================================================
# 02_utils_logging.R - Logging and Helper Utilities
# =============================================================================
# This file contains logging functions and basic helper utilities.
# =============================================================================

# -----------------------------------------------------------------------------
# Global Helper: opt() - Get option with default value
# -----------------------------------------------------------------------------
if (!exists("opt", mode = "function")) {
    opt <- function(nm, default) {
        if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
    }
}

# -----------------------------------------------------------------------------
# Null Coalescing Operator
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# Logging Function
# -----------------------------------------------------------------------------
log_msg <- function(text, ..., .envir = parent.frame()) {
    ts <- format(Sys.time(), "%H:%M:%S")
    msg <- tryCatch(
        {
            if (grepl("\\{[^}]+\\}", text)) {
                # has { } -> use glue
                glue::glue(text, ..., .envir = .envir)
            } else if (grepl("%", text)) {
                # has % -> use sprintf
                do.call(sprintf, c(list(fmt = text), list(...)))
            } else {
                # plain text; if there are extra arguments, also use sprintf format
                if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
            }
        },
        error = function(e) text
    )
    cat(sprintf("[%s] %s\n", ts, msg))
}

# -----------------------------------------------------------------------------
# Safe Filesystem Name Helper
# -----------------------------------------------------------------------------
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}

# -----------------------------------------------------------------------------
# Coerce Covariates Safely
# -----------------------------------------------------------------------------
# Converts columns to appropriate types and drops single-level factors
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
                if (length(lv) <= 1) { # single-level factor -> drop to avoid contrasts error
                    keep[cn] <- FALSE
                    if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
                } else {
                    df[[cn]] <- v
                }
            } else {
                df[[cn]] <- suppressWarnings(as.numeric(v)) # keep numeric as numeric
            }
        }
        df <- df[, keep, drop = FALSE]
        df
    }
}
