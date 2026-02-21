# =============================================================================
# 02_utils_logging.R - Logging and Helper Functions
# =============================================================================
# This module provides logging utilities and common helper functions
# used throughout the mRNA GSEA analysis pipeline.
# =============================================================================

# -----------------------------------------------------------------------------
# Option Retrieval Helper
# -----------------------------------------------------------------------------
#' Retrieve a value by name, with default fallback
#' @param nm Name of the variable to retrieve
#' @param default Default value if variable not found
#' @return The value of the variable or the default
if (!exists("opt", mode = "function")) {
    opt <- function(nm, default) {
        if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
    }
}

# -----------------------------------------------------------------------------
# Null Coalescing Operator
# -----------------------------------------------------------------------------
#' Return first argument if not NULL, otherwise second
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# Safe Filesystem Naming
# -----------------------------------------------------------------------------
#' Convert string to safe filesystem name
#' @param s Input string
#' @return String with unsafe characters replaced by underscores
safe_fs_name <- function(s) {
    s <- gsub('[<>:"/\\\\|?*]', "_", s)
    s <- gsub("\\s+", "_", s)
    s
}

# -----------------------------------------------------------------------------
# Timestamped Logging
# -----------------------------------------------------------------------------
#' Log a message with timestamp
#' @param text Message text (supports glue {} and sprintf % formatting)
#' @param ... Additional arguments for formatting
#' @param .envir Environment for glue evaluation
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
                # plain text; if extra arguments exist, treat as sprintf format
                if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
            }
        },
        error = function(e) text
    )
    cat(sprintf("[%s] %s\n", ts, msg))
}

# -----------------------------------------------------------------------------
# Force Serial Execution
# -----------------------------------------------------------------------------
#' Disable all parallel backends to ensure reproducibility
#' Affects: base R, data.table, BiocParallel, foreach, doParallel, future
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

    # fgsea: set option for wrapper to read
    options(.fgsea_nproc = 1L)
}

# Execute serial mode by default
.force_serial_execution()

# -----------------------------------------------------------------------------
# Directory Setup
# -----------------------------------------------------------------------------
#' Create run_info directory for logging and audit files
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Multiple Testing Policy Documentation
# -----------------------------------------------------------------------------
#' Write analysis notes about multiple testing policy
.write_analysis_notes <- function() {
    cat(
        paste(
            "Multiple-testing policy:",
            " - Within each gene-set group and per statistic, we control FDR using fgsea padj (per-dataset/stratum).",
            " - Pan-cancer summaries aggregate directions and counts across datasets/strata (descriptive).",
            " - Additionally, we provide a simple meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
            "   followed by Benjamini-Hochberg to report meta-level FDR (see csn_gsea_pan_summary_TP53/meta_fdr).",
            sep = "\n"
        ),
        file = file.path("run_info", "analysis_notes.txt")
    )
}
