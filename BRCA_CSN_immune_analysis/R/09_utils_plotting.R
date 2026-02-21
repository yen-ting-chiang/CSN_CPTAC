# =============================================================================
# 09_utils_plotting.R
# Plotting Utilities for BRCA Immune Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Color Palettes
# -----------------------------------------------------------------------------

#' Cell-style scatter plot color palettes
CELL_PALETTES <- list(
    # Original four sets
    CELL_BLUE        = list(pt = "#2D2D2D", ln = "#1F77B4"), # dark gray point / blue line
    CELL_TEAL        = list(pt = "#2D2D2D", ln = "#1B9E77"), # dark gray point / teal line
    CELL_ORANGE      = list(pt = "#2D2D2D", ln = "#D55E00"), # dark gray point / orange line
    CELL_PURPLE      = list(pt = "#2D2D2D", ln = "#6A3D9A"), # dark gray point / purple line
    # Brighter series
    CELL_BRIGHT_BLUE = list(pt = "#2D2D2D", ln = "#0077FF"),
    CELL_AQUA        = list(pt = "#2D2D2D", ln = "#00D5C6"),
    CELL_MINT        = list(pt = "#2D2D2D", ln = "#00C853"),
    CELL_MAGENTA     = list(pt = "#2D2D2D", ln = "#E91E63"),
    CELL_SCARLET     = list(pt = "#2D2D2D", ln = "#FF3B30"),
    CELL_GOLD        = list(pt = "#2D2D2D", ln = "#F6B400")
)

#' Heatmap color palettes
.BRCA_HM_PALETTES <- list(
    GSEA_DEFAULT      = c(neg = "#053061", mid = "#FFFFFF", pos = "#67001F"),
    BLUE_RED          = c(neg = "#2166AC", mid = "#F7F7F7", pos = "#B2182B"),
    BLUE_ORANGE       = c(neg = "#2B8CBE", mid = "#F7F7F7", pos = "#E34A33"),
    GREEN_MAGENTA     = c(neg = "#1B7837", mid = "#F7F7F7", pos = "#762A83"),
    BLUE_RED_DEEP     = c(neg = "#08306B", mid = "#F7F7F7", pos = "#7F0000"),
    BLUE_RED_BRIGHT   = c(neg = "#1F78B4", mid = "#FFFFFF", pos = "#E31A1C"),
    BLUE_RED_LIGHT    = c(neg = "#6BAED6", mid = "#F7F7F7", pos = "#FB6A4A"),
    BLUE_RED_TWILIGHT = c(neg = "#2C7FB8", mid = "#EEEEEE", pos = "#D7301F"),
    BLUE_RED_CRIMSON  = c(neg = "#0B3C5D", mid = "#FAFAFA", pos = "#B80C09")
)

#' Predictor ordering for heatmap x-axis
.pred_order_all <- c(
    "NES_CSN_SCORE", "NES_GPS1",
    "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
    "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
    "NES_RESIDUAL_GPS1",
    "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
    "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
)

# -----------------------------------------------------------------------------
# Theme and Helper Functions
# -----------------------------------------------------------------------------

#' Lighter CI fill color
#' @param col Color to lighten
#' @return Color with alpha = 0.25
.ci_fill <- function(col) tryCatch(grDevices::adjustcolor(col, alpha.f = 0.25), error = function(e) col)

#' Cell-style ggplot2 theme
#'
#' Clean theme with no outer frame, thin black axis lines, small text.
#'
#' @param base_size Base font size
#' @return ggplot2 theme object
.cell_theme <- function(base_size = 10) {
    title_sz <- get0("BRCA_TITLE_SIZE", ifnotfound = base_size * 1.0)
    subtitle_sz <- get0("BRCA_SUBTITLE_SIZE", ifnotfound = base_size * 0.9)
    axis_title_sz <- get0("BRCA_AXIS_TITLE_SIZE", ifnotfound = base_size * 0.9)
    axis_text_sz <- get0("BRCA_AXIS_TEXT_SIZE", ifnotfound = base_size * 0.85)

    ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            panel.border = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(linewidth = 0.6, colour = "black"),
            axis.ticks = ggplot2::element_line(linewidth = 0.6, colour = "black"),
            plot.title = ggplot2::element_text(
                face = "bold", size = title_sz,
                margin = ggplot2::margin(b = 2)
            ),
            plot.subtitle = ggplot2::element_text(
                size = subtitle_sz,
                margin = ggplot2::margin(b = 4)
            ),
            axis.title = ggplot2::element_text(size = axis_title_sz),
            axis.text = ggplot2::element_text(size = axis_text_sz)
        )
}

# -----------------------------------------------------------------------------
# Heatmap Plotting
# -----------------------------------------------------------------------------

#' Plot BRCA immune correlation heatmap
#'
#' Creates a heatmap of immune scores vs predictors correlation.
#'
#' @param dataset_id Dataset identifier
#' @param dataset_dirs_map Named vector of dataset directories
#' @param stratum Stratification (ALL, TP53_mutant, TP53_wild_type)
#' @param version Covariate version (RAW, BatchAdj)
#' @param method Correlation method (pearson, spearman)
#' @param hm_colorset Color palette name
#' @param out_formats Output file formats
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
#' @return List with plot and file_base
.plot_brca_immune_heatmap <- function(
  dataset_id = "brca_cptac_2020",
  dataset_dirs_map = get0("dataset_dirs_run", ifnotfound = NULL),
  stratum = c("ALL", "TP53_mutant", "TP53_wild_type"),
  version = c("RAW", "BatchAdj"),
  method = c("pearson", "spearman"),
  hm_colorset = get0("BRCA_HM_COLORSET", ifnotfound = "GSEA_DEFAULT"),
  out_formats = get0("BRCA_PLOT_FORMATS", ifnotfound = c("TIFF")),
  width = 7.0, height = 3.2, dpi = get0("BRCA_PLOT_DPI", ifnotfound = 600)
) {
    stratum <- match.arg(stratum)
    version <- match.arg(version)
    method <- match.arg(method)

    # 1) Find immune_vs_predictor_correlations.csv
    ds_root <- if (!is.null(dataset_dirs_map) && dataset_id %in% names(dataset_dirs_map)) {
        dataset_dirs_map[[dataset_id]]
    } else {
        NULL
    }
    if (is.null(ds_root) || !dir.exists(ds_root)) {
        stop("Cannot find dataset directory: ", dataset_id, " (dataset_dirs_map not set or path does not exist)")
    }

    # Try several candidate paths based on actual output location
    candidates <- c(
        file.path("BRCA_immune_analysis", stratum, version, "immune_vs_predictor_correlations.csv"),
        file.path(ds_root, "BRCA_immune_analysis", stratum, version, "immune_vs_predictor_correlations.csv")
    )
    csv_file <- NA_character_
    for (p in candidates) {
        if (file.exists(p)) {
            csv_file <- p
            break
        }
    }
    if (!is.character(csv_file) || is.na(csv_file) || !nzchar(csv_file)) {
        stop("Cannot find immune correlation CSV: ", paste(candidates, collapse = " | "))
    }

    dt <- data.table::fread(csv_file, na.strings = c("NA", "NaN", ""))
    if (!nrow(dt)) stop("File is empty: ", csv_file)

    # 2) Filter by method (PEARSON / SPEARMAN)
    keep <- grepl(method, dt$method, ignore.case = TRUE)
    dt <- dt[keep, ]
    if (!nrow(dt)) stop("No rows matching method: ", method, " (check CSV's method column)")

    # 3) padj column selection
    padj_col <- if (tolower(method) == "pearson" && "padj_pearson" %in% names(dt)) {
        "padj_pearson"
    } else if (tolower(method) == "spearman" && "padj_spearman" %in% names(dt)) {
        "padj_spearman"
    } else if ("padj" %in% names(dt)) "padj" else stop("Cannot find padj column in CSV")

    # 4) Prepare long format
    df_long <- data.frame(
        pathway = dt$score,
        predictor = paste0("NES_", dt$predictor),
        NES = dt$r,
        padj = dt[[padj_col]],
        stringsAsFactors = FALSE
    )

    # 5) y-axis order (top -> bottom)
    y_order <- c(
        "CIBERSORT_ABSOLUTE_SCORE",
        "XCELL_IMMUNE_SCORE",
        "ESTIMATE_IMMUNE_SCORE",
        "XCELL_STROMAL_SCORE",
        "ESTIMATE_STROMAL_SCORE"
    )

    # 6) Get palette
    pal <- .BRCA_HM_PALETTES[[hm_colorset]]
    if (is.null(pal)) pal <- .BRCA_HM_PALETTES[["GSEA_DEFAULT"]]
    if (is.character(pal) && is.null(names(pal)) && length(pal) >= 3) {
        pal <- c(neg = pal[1], mid = pal[2], pos = pal[3])
    }

    # 7) Plot
    title0 <- paste0("Immune-score vs Predictors (", toupper(method), ")")
    subt0 <- paste(dataset_id, stratum, version, sep = " / ")

    # Write temporary long table CSV for helper to read
    tmp_long_csv <- tempfile(pattern = "immune_long_", fileext = ".csv")
    data.table::fwrite(df_long, tmp_long_csv)

    g_res <- .make_heatmap_plot_with_yorder_colored(
        csv_file = tmp_long_csv,
        y_order  = y_order,
        palette  = pal
    )

    g <- if (is.list(g_res) && !is.null(g_res$plot)) g_res$plot else g_res
    if (is.list(g_res) && !is.null(g_res$width) && !is.null(g_res$height)) {
        width <- g_res$width
        height <- g_res$height
    }

    g <- g + ggplot2::labs(title = title0, subtitle = subt0)

    if (is.null(g) || !inherits(g, "ggplot")) {
        stop("[BRCA-immune][heatmap] Failed to generate heatmap object: g is NULL or not ggplot.")
    }

    # 8) Output
    out_dir <- dirname(csv_file)
    base_bn <- paste0("immune_vs_predictor_heatmap_", toupper(method))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    fmts <- toupper(out_formats)
    fmts <- if (length(fmts)) fmts else c("TIFF")
    any_ok <- FALSE
    for (f in fmts) {
        target <- file.path(out_dir, paste0(base_bn, ".", tolower(f)))
        if (f == "TIFF") {
            tryCatch(
                {
                    ragg::agg_tiff(
                        filename = target, width = width, height = height,
                        units = "in", res = dpi, compression = "lzw"
                    )
                    invisible(print(g))
                    grDevices::dev.off()
                    message(sprintf("[BRCA-immune][heatmap] Output: %s", target))
                },
                error = function(e) {
                    message(sprintf(
                        "[BRCA-immune][heatmap] Failed to save TIFF: %s -> %s",
                        conditionMessage(e), target
                    ))
                    stop(e)
                }
            )
            if (file.exists(target) && file.size(target) > 0) any_ok <- TRUE
        } else if (f %in% c("PNG", "JPG", "JPEG", "PDF")) {
            tryCatch(
                {
                    ggplot2::ggsave(
                        filename = target, plot = g,
                        width = width, height = height, units = "in", dpi = dpi,
                        bg = "white", limitsize = FALSE
                    )
                    message(sprintf("[BRCA-immune][heatmap] Output: %s", target))
                },
                error = function(e) {
                    message(sprintf(
                        "[BRCA-immune][heatmap] Failed to save %s: %s -> %s",
                        f, conditionMessage(e), target
                    ))
                    stop(e)
                }
            )
            if (file.exists(target) && file.size(target) > 0) any_ok <- TRUE
        } else {
            stop("[BRCA-immune][heatmap] Unsupported format: ", f)
        }
    }
    if (!any_ok) {
        stop("[BRCA-immune][heatmap] No output generated, please check if out_dir is writable: ", out_dir)
    }

    invisible(list(plot = g, file_base = file.path(out_dir, base_bn)))
}

#' Create heatmap plot with ordered y-axis and colored tiles
#'
#' @param csv_file Path to CSV file with correlation data
#' @param y_order Character vector specifying y-axis order
#' @param palette Named vector with neg, mid, pos colors
#' @return List with plot, width, and height
.make_heatmap_plot_with_yorder_colored <- function(csv_file, y_order, palette) {
    df_raw <- suppressMessages(.read_csv_safe(csv_file))
    path_candidates <- c("pathway", "Pathway", "term", "Term", "gs_name", "NAME", "set", "Set")
    path_col <- intersect(path_candidates, names(df_raw))[1]
    if (is.na(path_col)) stop("Cannot find pathway column: ", paste(path_candidates, collapse = ", "))

    has_long <- all(c("predictor", "NES", "padj") %in% names(df_raw))
    if (has_long) {
        df_long <- df_raw %>%
            dplyr::rename(pathway = dplyr::all_of(path_col)) %>%
            dplyr::mutate(
                predictor = as.character(.data$predictor),
                NES = as.numeric(.data$NES),
                padj = as.numeric(.data$padj)
            )
    } else {
        nes_cols <- grep("^NES_", names(df_raw), value = TRUE)
        padj_cols <- grep("^padj_", names(df_raw), value = TRUE)
        if (!length(nes_cols)) stop("Wide table format: cannot find NES_* columns.")
        df_wide <- df_raw %>% dplyr::rename(pathway = dplyr::all_of(path_col))
        nes_map <- tibble::tibble(nes_col = nes_cols, key = sub("^NES_", "", nes_cols))
        padj_map <- tibble::tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
        pair_map <- dplyr::left_join(nes_map, padj_map, by = "key")
        df_list <- lapply(seq_len(nrow(pair_map)), function(i) {
            nes_c <- pair_map$nes_col[i]
            padj_c <- pair_map$padj_col[i]
            tibble::tibble(
                pathway   = df_wide$pathway,
                predictor = paste0("NES_", pair_map$key[i]),
                NES       = suppressWarnings(as.numeric(df_wide[[nes_c]])),
                padj      = if (!is.na(padj_c)) suppressWarnings(as.numeric(df_wide[[padj_c]])) else NA_real_
            )
        })
        df_long <- dplyr::bind_rows(df_list)
    }

    present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
    df_long <- df_long %>% dplyr::filter(.data$predictor %in% present_preds)

    # y-axis uses ALL's order
    df_long <- df_long %>% dplyr::filter(.data$pathway %in% y_order)
    y_levels <- y_order[y_order %in% df_long$pathway]
    if (!length(y_levels)) stop("Ordered version: no intersection with ALL's pathway order.")
    df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(y_levels)))

    # gap rules and positions
    gap <- 0.4
    needs_gap1 <- all(c("NES_CSN_SCORE", "NES_GPS1") %in% present_preds)
    needs_gap2 <- "NES_RESIDUAL_GPS1" %in% present_preds
    pos_map <- list()
    pos <- 0
    for (p in present_preds) {
        if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
        if (p == "NES_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
        pos <- pos + 1
        pos_map[[p]] <- pos
    }
    pos_map <- unlist(pos_map)
    df_long <- df_long %>% dplyr::mutate(xpos = unname(pos_map[.data$predictor]))

    # Symmetric color scale range
    L <- max(abs(df_long$NES), na.rm = TRUE)
    if (!is.finite(L) || L == 0) L <- 1

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
        ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
        ggplot2::geom_point(
            data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
            ggplot2::aes(x = xpos, y = pathway),
            shape = 16, size = 1.6, color = "black", inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_gradient2(
            low = palette[["neg"]], mid = palette[["mid"]], high = palette[["pos"]],
            limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "NES"
        ) +
        ggplot2::scale_x_continuous(
            breaks = unname(pos_map[present_preds]),
            labels = present_preds,
            expand = ggplot2::expansion(mult = c(0.01, 0.01)),
            position = "top"
        ) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.x.top = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = ggplot2::margin(b = 8)),
            axis.text.y = ggplot2::element_text(size = 9),
            legend.position = "right",
            plot.margin = ggplot2::margin(6, 12, 6, 6),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.background = ggplot2::element_rect(fill = "white", colour = NA)
        ) +
        ggplot2::coord_cartesian(clip = "off")

    n_path <- nlevels(df_long$pathway)
    W <- max(8, length(present_preds) * 0.45)
    H <- max(6, n_path * 0.22)
    list(plot = p, width = W, height = H)
}
