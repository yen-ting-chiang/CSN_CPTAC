# =============================================================================
# 02_utils_logging.R
# Logging and Execution Control Utilities
# =============================================================================

# -----------------------------------------------------------------------------
# Logging Function
# -----------------------------------------------------------------------------

#' Flexible logging function with timestamp
#'
#' Supports glue-style {variable} interpolation, sprintf-style %s formatting,
#' or plain text output. Automatically detects the format based on the text.
#'
#' @param text Message text (may contain {var} or %s placeholders)
#' @param ... Additional arguments for sprintf or glue
#' @param .envir Environment for glue interpolation
#' @return NULL (invisibly), prints message to console
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
                # plain text; if extra args, treat as sprintf format
                if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
            }
        },
        error = function(e) text
    )
    cat(sprintf("[%s] %s\n", ts, msg))
}

# -----------------------------------------------------------------------------
# Serial Execution Control
# -----------------------------------------------------------------------------

#' Force single-threaded execution across common parallel frameworks
#'
#' Disables parallelism in base R, data.table, BiocParallel, foreach,
#' doParallel, and future packages to ensure reproducibility.
#'
#' @return NULL (invisibly)
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
    if (requireNamespace("doParallel", quietly = TRUE)) {}

    # future
    if (requireNamespace("future", quietly = TRUE)) {
        future::plan(future::sequential)
    }
    Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")
}
