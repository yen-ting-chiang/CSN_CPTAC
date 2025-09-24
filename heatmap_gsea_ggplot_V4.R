## =========================================================
## 批次將單一 dataset/stratum 的 GSEA Summary CSV 繪製成 heatmap
## - 兼容 RAW 與 BatchAdj（含 limma_t, interaction, spearman）
## - 完全沿用 heatmap_gsea_ggplot.txt 的讀檔、轉長表、配色與打點規則
## - 產生兩套輸出：就地輸出 & 彙整輸出 single_dataset_GSEA_heatmap/<dataset>/
## =========================================================

setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
##setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
})


## ===== Cell 風格：正值橘、負值綠（多組可選）=====
CELL_OG_PALETTES <- list(
  # 深綠/深橘（對比強）
  cell_og1 = c(neg = "#1B9E77", mid = "#FFFFFF", pos = "#D95F02"),
  # 墨綠/琥珀橘（偏 Cell 視覺）
  cell_og2 = c(neg = "#2E7D32", mid = "#FFFFFF", pos = "#FB8C00"),
  # 綠玉/芥末橘（柔和）
  cell_og3 = c(neg = "#009E73", mid = "#FFFFFF", pos = "#E69F00"),
  # 森綠/焦糖橘（高對比）
  cell_og4 = c(neg = "#006D2C", mid = "#FFFFFF", pos = "#E6550D"),
  # 松綠/橘紅（更飽和）
  cell_og5 = c(neg = "#238B45", mid = "#FFFFFF", pos = "#F16913")
)

# 可用 options(csn_interaction_palette = "cell_og3") 選盤；預設用 cell_og2
.get_interaction_palette <- function() {
  key <- getOption("csn_interaction_palette", "cell_og2")
  pal <- CELL_OG_PALETTES[[key]]
  if (is.null(pal)) pal <- CELL_OG_PALETTES[[1]]
  pal
}

# 與 .make_heatmap_plot_with_yorder 規則一致，但允許指定橘/綠配色
.make_heatmap_plot_with_yorder_colored <- function(csv_file, y_order, palette = .get_interaction_palette()) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway","Pathway","term","Term","gs_name","NAME","set","Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("找不到 pathway 欄位：", paste(path_candidates, collapse=", "))
  
  has_long <- all(c("predictor","NES","padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      dplyr::rename(pathway = dplyr::all_of(path_col)) %>%
      dplyr::mutate(predictor = as.character(.data$predictor),
                    NES = as.numeric(.data$NES),
                    padj = as.numeric(.data$padj))
  } else {
    nes_cols  <- grep("^NES_",  names(df_raw), value = TRUE)
    padj_cols <- grep("^padj_", names(df_raw), value = TRUE)
    if (!length(nes_cols)) stop("寬表格式：找不到 NES_* 欄位。")
    df_wide <- df_raw %>% dplyr::rename(pathway = dplyr::all_of(path_col))
    nes_map  <- tibble::tibble(nes_col  = nes_cols,  key = sub("^NES_",  "", nes_cols))
    padj_map <- tibble::tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
    pair_map <- dplyr::left_join(nes_map, padj_map, by = "key")
    df_list <- lapply(seq_len(nrow(pair_map)), function(i){
      nes_c  <- pair_map$nes_col[i]
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
  
  # y 軸用 ALL 的順序；僅取交集
  y_use <- intersect(y_order, unique(df_long$pathway))
  if (!length(y_use)) stop("ordered 版本：與 ALL 的 pathway 順序沒有交集。")
  df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(y_use)))
  
  # gap 規則與位置（與你原本一致）
  gap <- 0.4
  needs_gap1 <- all(c("NES_CSN_SCORE","NES_GPS1") %in% present_preds)
  needs_gap2 <- "NES_RESIDUAL_GPS1" %in% present_preds
  pos_map <- list(); pos <- 0
  for (p in present_preds) {
    if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "NES_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- df_long %>% dplyr::mutate(xpos = unname(pos_map[.data$predictor]))
  
  # 對稱色階範圍
  L <- max(abs(df_long$NES), na.rm = TRUE); if (!is.finite(L) || L == 0) L <- 1
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = NES)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
                        ggplot2::aes(x = xpos, y = pathway),
                        shape = 16, size = 1.6, color = "black", inherit.aes = FALSE) +
    ggplot2::scale_fill_gradient2(low = palette[["neg"]], mid = palette[["mid"]], high = palette[["pos"]],
                                  limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "NES") +
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
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    ) + ggplot2::coord_cartesian(clip = "off")
  
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}


.read_csv_safe <- function(path) {
  p <- path
  
  # Windows: 如果原字串不存在，試著把 / 換成 \
  if (!file.exists(p) && .Platform$OS.type == "windows") {
    p2 <- gsub("/", "\\\\", p, fixed = TRUE)
    if (file.exists(p2)) p <- p2
  }
  if (!file.exists(p)) stop("檔案不存在：", p)
  
  fi <- suppressWarnings(file.info(p))
  if (isTRUE(fi$isdir)) stop("目標不是檔案：", p)
  if (!is.na(fi$size) && fi$size == 0) stop("空檔（size = 0）：", p)
  
  # 1) readr 預設
  out <- try(suppressMessages(readr::read_csv(p, show_col_types = FALSE, progress = FALSE)), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  # 2) readr + UTF-8-BOM
  out <- try(suppressMessages(readr::read_csv(
    p, locale = readr::locale(encoding = "UTF-8-BOM"),
    show_col_types = FALSE, progress = FALSE
  )), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  # 3) 不手動開連線：讀成全文字串再解析
  txt <- try(readr::read_file(p), silent = TRUE)
  if (!inherits(txt, "try-error")) {
    out <- try(suppressMessages(readr::read_csv(
      readr::I(txt), show_col_types = FALSE, progress = FALSE
    )), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  
  # 4) 最後退路：base R
  out <- try(utils::read.csv(p, check.names = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) stop(out)
  tibble::as_tibble(out)
}




.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# 目標 X 軸順序（沿用你原腳本）
.pred_order_all <- c(
  "NES_CSN_SCORE", "NES_GPS1",
  "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
  "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
  "NES_RESIDUAL_GPS1",
  "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
  "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
)

# 單一 CSV → 長表 → ggplot 物件（不輸出檔案）
.make_heatmap_plot <- function(csv_file) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  
  # 自動偵測 pathway 欄位（沿用你原腳本）
  path_candidates <- c("pathway","Pathway","term","Term","gs_name","NAME","set","Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("找不到 pathway 欄位：", paste(path_candidates, collapse=", "))
  
  # 長/寬表判定 + 轉長表（沿用）
  has_long <- all(c("predictor","NES","padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      rename(pathway = all_of(path_col)) %>%
      mutate(predictor = as.character(predictor),
             NES = as.numeric(NES),
             padj = as.numeric(padj))
  } else {
    nes_cols  <- grep("^NES_",  names(df_raw), value = TRUE)
    padj_cols <- grep("^padj_", names(df_raw), value = TRUE)
    if (!length(nes_cols)) stop("寬表格式：找不到 NES_* 欄位。")
    
    df_wide <- df_raw %>% rename(pathway = all_of(path_col))
    nes_map  <- tibble(nes_col  = nes_cols,  key = sub("^NES_",  "", nes_cols))
    padj_map <- tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
    pair_map <- nes_map %>% left_join(padj_map, by = "key")
    
    df_list <- lapply(seq_len(nrow(pair_map)), function(i){
      nes_c  <- pair_map$nes_col[i]
      padj_c <- pair_map$padj_col[i]  # 可能為 NA
      tibble(
        pathway   = df_wide$pathway,
        predictor = paste0("NES_", pair_map$key[i]),
        NES       = suppressWarnings(as.numeric(df_wide[[nes_c]])),
        padj      = if (!is.na(padj_c)) suppressWarnings(as.numeric(df_wide[[padj_c]])) else NA_real_
      )
    })
    df_long <- bind_rows(df_list)
  }
  
  # 僅保留指定順序中的 predictors，並依順序排序
  present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
  if (!length(present_preds)) stop("資料中沒有任何指定的 predictors。")
  df_long <- df_long %>% filter(predictor %in% present_preds)
  
  # Y 軸排序：優先用 NES_CSN_SCORE 由大到小；否則以行均值
  if ("NES_CSN_SCORE" %in% present_preds) {
    nes_csn <- df_long %>% filter(predictor == "NES_CSN_SCORE") %>%
      select(pathway, NES) %>% distinct()
    path_order <- nes_csn %>% arrange(desc(NES)) %>% pull(pathway)
  } else {
    path_order <- df_long %>% group_by(pathway) %>%
      summarise(m = mean(NES, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(m)) %>% pull(pathway)
  }
  df_long <- df_long %>% mutate(pathway = factor(pathway, levels = rev(path_order)))
  
  
  # X 軸不等距位置與兩個 gap
  # 規則調整：即使「沒有 NES_COPS9」，只要有 NES_RESIDUAL_GPS1 仍保留第二個間隔
  gap <- 0.4
  needs_gap1 <- all(c("NES_CSN_SCORE","NES_GPS1") %in% present_preds)
  needs_gap2 <- "NES_RESIDUAL_GPS1" %in% present_preds  # 修改：不再要求同時有 COPS9
  
  pos_map <- list(); pos <- 0
  for (p in present_preds) {
    if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "NES_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- df_long %>% mutate(xpos = unname(pos_map[predictor]))
  
  # 顏色範圍（對稱 0）與 cell 風格色
  L <- max(abs(df_long$NES), na.rm = TRUE); if (!is.finite(L) || L == 0) L <- 1
  col_low <- "#053061"; col_mid <- "#FFFFFF"; col_high <- "#67001F"
  
  # 作圖（黑點：padj < 0.05）
  p <- ggplot(df_long, aes(x = xpos, y = pathway, fill = NES)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(data = df_long %>% filter(is.finite(padj), padj < 0.05),
               aes(x = xpos, y = pathway),
               shape = 16, size = 1.6, color = "black", inherit.aes = FALSE) +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                         limits = c(-L, L), oob = scales::squish, name = "NES") +
    scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = margin(b = 8)),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.margin = margin(6, 12, 6, 6),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) + coord_cartesian(clip = "off")
  
  # 尺寸：隨 pathway 數量調整
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45)   # inches
  H <- max(6, n_path * 0.22)                  # inches
  list(plot = p, width = W, height = H)
}

# 讀取 CSV 並依照 .make_heatmap_plot 的規則計算 pathway 排序（回傳向量）
.compute_pathway_order_from_csv <- function(csv_file) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway","Pathway","term","Term","gs_name","NAME","set","Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("找不到 pathway 欄位：", paste(path_candidates, collapse=", "))
  
  has_long <- all(c("predictor","NES","padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      dplyr::rename(pathway = dplyr::all_of(path_col)) %>%
      dplyr::mutate(predictor = as.character(.data$predictor),
                    NES = as.numeric(.data$NES))
  } else {
    nes_cols  <- grep("^NES_",  names(df_raw), value = TRUE)
    padj_cols <- grep("^padj_", names(df_raw), value = TRUE)  # 只為了與 make 相容，排序不用
    if (!length(nes_cols)) stop("寬表格式：找不到 NES_* 欄位。")
    df_wide <- df_raw %>% dplyr::rename(pathway = dplyr::all_of(path_col))
    nes_map  <- tibble::tibble(nes_col  = nes_cols,  key = sub("^NES_",  "", nes_cols))
    df_list <- lapply(seq_len(nrow(nes_map)), function(i){
      nes_c  <- nes_map$nes_col[i]
      tibble::tibble(
        pathway   = df_wide$pathway,
        predictor = paste0("NES_", nes_map$key[i]),
        NES       = suppressWarnings(as.numeric(df_wide[[nes_c]]))
      )
    })
    df_long <- dplyr::bind_rows(df_list)
  }
  
  present_preds <- intersect(.pred_order_all, unique(df_long$predictor))
  if (!length(present_preds)) stop("資料中沒有任何指定的 predictors。")
  
  if ("NES_CSN_SCORE" %in% present_preds) {
    nes_csn <- df_long %>% dplyr::filter(.data$predictor == "NES_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$NES) %>% dplyr::distinct()
    path_order <- nes_csn %>% dplyr::arrange(dplyr::desc(.data$NES)) %>% dplyr::pull(.data$pathway)
  } else {
    path_order <- df_long %>%
      dplyr::group_by(.data$pathway) %>% dplyr::summarise(m = mean(.data$NES, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$m)) %>% dplyr::pull(.data$pathway)
  }
  path_order
}

# 以既定 y_order（pathway levels）作圖（規則完全仿 .make_heatmap_plot）
.make_heatmap_plot_with_yorder <- function(csv_file, y_order) {
  df_raw <- suppressMessages(.read_csv_safe(csv_file))
  path_candidates <- c("pathway","Pathway","term","Term","gs_name","NAME","set","Set")
  path_col <- intersect(path_candidates, names(df_raw))[1]
  if (is.na(path_col)) stop("找不到 pathway 欄位：", paste(path_candidates, collapse=", "))
  
  has_long <- all(c("predictor","NES","padj") %in% names(df_raw))
  if (has_long) {
    df_long <- df_raw %>%
      dplyr::rename(pathway = dplyr::all_of(path_col)) %>%
      dplyr::mutate(predictor = as.character(.data$predictor),
                    NES = as.numeric(.data$NES),
                    padj = as.numeric(.data$padj))
  } else {
    nes_cols  <- grep("^NES_",  names(df_raw), value = TRUE)
    padj_cols <- grep("^padj_", names(df_raw), value = TRUE)
    if (!length(nes_cols)) stop("寬表格式：找不到 NES_* 欄位。")
    df_wide <- df_raw %>% dplyr::rename(pathway = dplyr::all_of(path_col))
    nes_map  <- tibble::tibble(nes_col  = nes_cols,  key = sub("^NES_",  "", nes_cols))
    padj_map <- tibble::tibble(padj_col = padj_cols, key = sub("^padj_", "", padj_cols))
    pair_map <- dplyr::left_join(nes_map, padj_map, by = "key")
    df_list <- lapply(seq_len(nrow(pair_map)), function(i){
      nes_c  <- pair_map$nes_col[i]
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
  
  # 用 ALL 的順序；僅取交集
  y_use <- intersect(y_order, unique(df_long$pathway))
  if (!length(y_use)) stop("ordered 版本：與 ALL 的 pathway 順序沒有交集。")
  df_long <- df_long %>% dplyr::mutate(pathway = factor(.data$pathway, levels = rev(y_use)))
  
  # gap 規則與配色、打點 — 完全沿用
  gap <- 0.4
  needs_gap1 <- all(c("NES_CSN_SCORE","NES_GPS1") %in% present_preds)
  needs_gap2 <- "NES_RESIDUAL_GPS1" %in% present_preds  # 你先前的修正
  pos_map <- list(); pos <- 0
  for (p in present_preds) {
    if (p == "NES_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "NES_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- df_long %>% dplyr::mutate(xpos = unname(pos_map[.data$predictor]))
  
  L <- max(abs(df_long$NES), na.rm = TRUE); if (!is.finite(L) || L == 0) L <- 1
  col_low <- "#053061"; col_mid <- "#FFFFFF"; col_high <- "#67001F"
  
  p <- ggplot(df_long, aes(x = xpos, y = pathway, fill = NES)) +
    geom_tile(width = 1, height = 0.9, color = NA) +
    geom_point(data = df_long %>% dplyr::filter(is.finite(.data$padj), .data$padj < 0.05),
               aes(x = xpos, y = pathway),
               shape = 16, size = 1.6, color = "black", inherit.aes = FALSE) +
    scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                         limits = c(-L, L), oob = scales::squish, name = "NES") +
    scale_x_continuous(
      breaks = unname(pos_map[present_preds]),
      labels = present_preds,
      expand = expansion(mult = c(0.01, 0.01)),
      position = "top"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 9, margin = margin(b = 8)),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.margin = margin(6, 12, 6, 6),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    ) + coord_cartesian(clip = "off")
  
  n_path <- nlevels(df_long$pathway)
  W <- max(8, length(present_preds) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}

# 只在「就地」輸出四種格式（_ordered 版用）
.save_near_only <- function(p, width, height, csv_file, meta, suffix = "_ordered", dpi = 600) {
  bn <- basename(csv_file)
  dir_near <- dirname(csv_file)
  near_prefix <- paste0("heatmap_", .safe_fs(meta$dataset), "_", .safe_fs(meta$stratum),
                        "_", .safe_fs(meta$variant), "_")
  near_base <- file.path(dir_near, paste0(near_prefix, bn, suffix))
  ggsave(paste0(near_base, ".tiff"), p, width = width, height = height, units = "in",
         dpi = dpi, bg = "white", compression = "lzw")
  #ggsave(paste0(near_base, ".png"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white")
  #ggsave(paste0(near_base, ".jpg"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white", quality = 100)
  #ggsave(paste0(near_base, ".pdf"),  p, width = width, height = height, units = "in",
  #       bg = "white")
}


# 從路徑推斷 dataset / variant / stratum
.parse_meta_from_path <- function(csv_file) {
  parts    <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
  parts_lc <- tolower(parts)
  bn       <- basename(csv_file)
  
  # 變體（RAW / BatchAdj）：抓任一層
  vidx <- which(parts_lc %in% c("raw", "batchadj"))
  variant <- if (length(vidx)) parts[vidx[1]] else NA_character_
  
  # stratum 真正位置：在 variant 的「上一層」
  stratum <- NA_character_
  if (length(vidx) && vidx[1] - 1 >= 1) {
    cand <- parts[vidx[1] - 1]
    if (tolower(cand) %in% c("all", "tp53_mutant", "tp53_wild_type")) {
      stratum <- cand
    }
  }
  # 後備：若上一步沒抓到，就在任何一層找一次
  if (is.na(stratum) || !nzchar(stratum)) {
    hit <- which(parts_lc %in% c("all","tp53_mutant","tp53_wild_type"))
    if (length(hit)) stratum <- parts[hit[1]]
  }
  
  # dataset：若路徑含 csn_gsea_results_TP53，取它的上一層；否則交由呼叫端 fallback
  dataset <- NA_character_
  hit2 <- which(parts_lc == "csn_gsea_results_tp53")
  if (length(hit2) && hit2[1] > 1) dataset <- parts[hit2[1] - 1]
  
  list(dataset = dataset, variant = variant, stratum = stratum)
}


# 實際輸出：四種格式（就地 & 彙整）
.save_both_places <- function(p, width, height, csv_file, meta,
                              collect_root = "single_dataset_GSEA_heatmap", dpi = 600) {
  bn <- basename(csv_file) # e.g. Summary_H_GSEA_limma_t_cont_ALL.csv
  dir_near <- dirname(csv_file)
  
  # 1) 就地檔名
  near_prefix <- paste0("heatmap_", .safe_fs(meta$dataset), "_", .safe_fs(meta$stratum),
                        "_", .safe_fs(meta$variant), "_")
  near_base <- file.path(dir_near, paste0(near_prefix, bn))
  
  # 2) 彙整路徑
  out_dir2 <- file.path(collect_root, .safe_fs(meta$dataset))
  dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
  collect_base <- file.path(out_dir2, paste0("heatmap_", bn))
  
  # 寫四種格式
  ggsave(paste0(near_base, ".tiff"), p, width = width, height = height, units = "in",
         dpi = dpi, bg = "white", compression = "lzw")
  #ggsave(paste0(near_base, ".png"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white")
  #ggsave(paste0(near_base, ".jpg"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white", quality = 100)
  #ggsave(paste0(near_base, ".pdf"),  p, width = width, height = height, units = "in",
  #       bg = "white")
  
  ggsave(paste0(collect_base, ".tiff"), p, width = width, height = height, units = "in",
         dpi = dpi, bg = "white", compression = "lzw")
  #ggsave(paste0(collect_base, ".png"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white")
  #ggsave(paste0(collect_base, ".jpg"),  p, width = width, height = height, units = "in",
  #       dpi = dpi, bg = "white", quality = 100)
  #ggsave(paste0(collect_base, ".pdf"),  p, width = width, height = height, units = "in",
  #       bg = "white")
}

# 掃描單一 dataset 根目錄，抓到指定 Summary CSV 清單
.find_target_csvs_for_dataset <- function(ds_dir) {
  # 搜尋名稱匹配的 Summary 檔
  files <- list.files(
    ds_dir,
    pattern = "^Summary_H_GSEA_(limma_t_cont|limma_interaction|spearman)_(ALL|TP53_mutant|TP53_wild_type)\\.csv$",
    recursive = TRUE, full.names = TRUE
  )
  if (!length(files)) return(files)
  
  # 仍要求路徑中要有 RAW 或 BatchAdj（同時接受 / 或 \）
  sep_pat <- "[/\\\\]"  # 同時匹配 / 與 \
  keep <- grepl(paste0(sep_pat, "raw", sep_pat), files, ignore.case = TRUE) |
    grepl(paste0(sep_pat, "batchadj", sep_pat), files, ignore.case = TRUE)
  files <- files[keep]
  
  
  ## ===== 這裡就是你要新增的區塊（放在 return 前）=====
  # 剔除 size=0（尚未寫完或空檔）
  if (length(files)) {
    finfo <- suppressWarnings(file.info(files))
    files <- files[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  ## ====================================================
  
  files
}


# 主程式：跑所有 dataset
run_gsea_heatmaps_for_all_datasets <- function(dataset_dirs_map = NULL,
                                               collect_root = "single_dataset_GSEA_heatmap") {
  if (is.null(dataset_dirs_map)) {
    if (exists("dataset_dirs_run")) {
      dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
    } else if (exists("dataset_dirs")) {
      dataset_dirs_map <- get("dataset_dirs", inherits = TRUE)
    } else stop("請提供 dataset_dirs_map 或在環境中定義 dataset_dirs_run / dataset_dirs")
  }
  for (ds in names(dataset_dirs_map)) {
    ds_dir <- dataset_dirs_map[[ds]]
    csvs <- .find_target_csvs_for_dataset(ds_dir)
    if (!length(csvs)) {
      message("[heatmap] ", ds, ": 找不到符合規則的 Summary CSV，略過。")
      next
    }
    for (csv in csvs) {
      meta <- .parse_meta_from_path(csv)
      # 若以路徑無法推到 dataset 名稱，退回 dataset_dirs_map 的鍵名
      if (!length(meta$dataset) || is.na(meta$dataset) || meta$dataset == "")
        meta$dataset <- ds
      
      # 如果是 interaction 但 stratum 不是 ALL，就跳過（以防萬一）
      if (grepl("GSEA_limma_interaction", basename(csv), fixed = TRUE) &&
          !identical(meta$stratum, "ALL")) {
        next
      }
      
      h <- try(.make_heatmap_plot(csv), silent = TRUE)
      if (inherits(h, "try-error")) {
        message("[heatmap] 失敗：", csv, " | ", as.character(h))
        next
      }
      .save_both_places(h$plot, h$width, h$height, csv, meta, collect_root = collect_root)
      message("[heatmap] 完成：", csv)
      
      ## === 新增：interaction 的 _ordered 版本（依同 dataset + 同 variant 的 limma_t_cont ALL）===
      # 僅處理 ALL 層的 interaction（你的資料結構也只在 ALL 才有 interaction）
      bn_lc <- tolower(basename(csv))
      if (grepl("gsea_limma_interaction", bn_lc, fixed = TRUE) &&
          identical(meta$stratum, "ALL")) {
        # 目標 limma_t_cont ALL 的路徑：將資料夾從 GSEA_limma_interaction 換成 GSEA_limma_t_cont，
        # 並將檔名中的 limma_interaction 換成 limma_t_cont
        dir_now  <- normalizePath(dirname(csv), winslash = "/", mustWork = FALSE)
        dir_cont <- sub("/GSEA_limma_interaction/?$", "/GSEA_limma_t_cont", dir_now, ignore.case = TRUE)
        bn_cont  <- gsub("GSEA_limma_interaction", "GSEA_limma_t_cont", basename(csv), ignore.case = TRUE)
        cont_csv <- file.path(dir_cont, bn_cont)
        
        if (file.exists(cont_csv)) {
          y_order <- try(.compute_pathway_order_from_csv(cont_csv), silent = TRUE)
          if (!inherits(y_order, "try-error")) {
            h2 <- try(.make_heatmap_plot_with_yorder_colored(csv, y_order), silent = TRUE)
            if (!inherits(h2, "try-error")) {
              .save_near_only(h2$plot, h2$width, h2$height, csv, meta, suffix = "_ordered")
              message("[heatmap][interaction-ordered] 完成（依 limma_t_cont ALL）：", csv)
            } else {
              message("[heatmap][interaction-ordered] 建圖失敗，略過：", csv, " | ", as.character(h2))
            }
          } else {
            message("[heatmap][interaction-ordered] 取 limma_t_cont ALL 順序失敗，略過：", cont_csv, " | ", as.character(y_order))
          }
        } else {
          message("[heatmap][interaction-ordered] 找不到對應 limma_t_cont ALL，略過：", cont_csv)
        }
      }
      
      ## === 新增：TP53_mutant / TP53_wild_type 的 ordered 版本（依 ALL 的 y-axis 順序） ===
      # 僅限 limma_t_cont 與 spearman 兩種方法
      if (tolower(meta$stratum) %in% c("tp53_mutant","tp53_wild_type")) {
        bn_lc <- tolower(basename(csv))
        is_cont     <- grepl("gsea_limma_t_cont", bn_lc, fixed = TRUE)
        is_spearman <- grepl("gsea_spearman",     bn_lc, fixed = TRUE)
        if (is_cont || is_spearman) {
          # 依「同 dataset + 同 variant」構造 ALL 對應路徑：把路徑中 variant 的上一層換成 ALL
          parts    <- strsplit(normalizePath(csv, winslash = "/"), "/")[[1]]
          parts_lc <- tolower(parts)
          vidx <- which(parts_lc %in% c("raw", "batchadj"))
          if (length(vidx) && vidx[1] - 1 >= 1) {
            parts2 <- parts; parts2[vidx[1] - 1] <- "ALL"
            all_csv <- paste(parts2, collapse = "/")
            if (file.exists(all_csv)) {
              # 取 ALL 的 pathway 順序
              y_order <- try(.compute_pathway_order_from_csv(all_csv), silent = TRUE)
              if (!inherits(y_order, "try-error")) {
                # 依 ALL 順序重畫當前 CSV
                h2 <- try(.make_heatmap_plot_with_yorder(csv, y_order), silent = TRUE)
                if (!inherits(h2, "try-error")) {
                  .save_near_only(h2$plot, h2$width, h2$height, csv, meta, suffix = "_ordered")
                  message("[heatmap][ordered] 完成（依 ALL 順序）：", csv)
                } else {
                  message("[heatmap][ordered] 建圖失敗，略過：", csv, " | ", as.character(h2))
                }
              } else {
                message("[heatmap][ordered] 取 ALL 順序失敗，略過：", all_csv, " | ", as.character(y_order))
              }
            } else {
              message("[heatmap][ordered] 找不到對應 ALL 檔案，略過：", all_csv)
            }
          }
        }
      }
      # ---- 每處理完一個檔案就清理一次連線與記憶體（避免 Windows "Too many open files"）----
      try(closeAllConnections(), silent = TRUE)
      invisible(gc(FALSE))
    }
  }
  invisible(TRUE)
}

options(csn_interaction_palette = "cell_og2")

## ---- 執行（例）----
run_gsea_heatmaps_for_all_datasets(
  dataset_dirs_map = if (exists("dataset_dirs_run")) dataset_dirs_run else NULL,
  collect_root = "single_dataset_GSEA_heatmap"
)
