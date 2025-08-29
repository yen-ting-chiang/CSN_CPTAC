setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC")
getwd()

## ===== 路徑設定 =====
project_dir    <- "C:/Users/danny/Documents/R_project/CSN_CPTAC/CSN_CPTAC_cBioPortal"
dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)

# <<< 這裡改成你的 brca_cptac_2020 資料夾實際路徑 >>>
local_brca_dir <- "C:/Users/danny/Documents/R_project/CSN_CPTAC/brca_cptac_2020"

## ===== 套件 =====
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(janitor); library(glue)
  library(SummarizedExperiment); library(MultiAssayExperiment); library(S4Vectors)
  library(limma); library(imputeLCMD); library(msigdbr); library(fgsea)
  library(openxlsx); library(ComplexHeatmap); library(cowplot)
})

## ===== 參數 =====
set.seed(1234)

csn_subunits <- c("GPS1","COPS2","COPS3","COPS4","COPS5","COPS6","COPS7A","COPS7B","COPS8","COPS9")

min_frac_complete   <- 0.75     # 每基因至少 75% 非 NA
hi_lo_quantile      <- 0.25     # 高低各 25%
minSize <- 15; maxSize <- 500; nperm <- 10000

min_per_group      <- 3         # High/Low 各至少多少樣本（不足就跳過 limma）
min_pairs_spearman <- 10        # Spearman 至少配對樣本數

log_msg <- function(text, ..., .envir = parent.frame()) {
  ts  <- format(Sys.time(), "%H:%M:%S")
  msg <- tryCatch(glue::glue(text, ..., .envir = .envir), error = function(.) text)
  cat(sprintf("[%s] %s\n", ts, msg))
}

## ===== 一、讀取本地 BRCA 資料夾 =====
# 1) 讀 sample list（cases_protein_quantification.txt）
read_case_list <- function(path_file){
  if (!file.exists(path_file)) return(character(0))
  x <- readLines(path_file, warn = FALSE, encoding = "UTF-8")
  line <- x[grepl("^case_list_ids:", x)]
  if (length(line) == 0) return(character(0))
  ids <- sub("^case_list_ids:\\s*", "", line[1])
  ids <- unlist(strsplit(ids, "[,\\s]+"))        # 空白或逗號切
  unique(ids[nchar(ids) > 0])
}

# 2) 讀 protein 矩陣（data_protein_quantification.txt 或 zscores）
load_local_protein_matrix <- function(dir){
  f1 <- file.path(dir, "data_protein_quantification.txt")
  f2 <- file.path(dir, "data_protein_quantification_zscores_ref_all_samples.txt")
  if (file.exists(f1))      fp <- f1
  else if (file.exists(f2)) fp <- f2
  else stop("資料夾內找不到 protein 矩陣檔（data_protein_quantification*.txt）")
  
  log_msg("讀取 protein 矩陣：{basename(fp)}")
  dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))
  
  gene_cols <- c("Hugo_Symbol","hugo_symbol","Gene","Gene_Symbol","HugoSymbol",
                 "GENE_SYMBOL","gene","gene_symbol")
  gcol <- intersect(gene_cols, names(dat)); if (length(gcol) == 0) gcol <- names(dat)[1]
  dat  <- dplyr::rename(dat, Gene = !!gcol[1])
  dat$Gene <- sub("\\|.*$", "", dat$Gene)       # 取「|」前的基因名
  
  not_sample <- c("Gene","Entrez_Gene_Id","Entrez_Gene_Id.","ENTREZ_GENE_ID",
                  "Description","Gene_Name","GeneName","Gene_Symbol")
  sample_cols_all <- setdiff(names(dat), not_sample)
  
  case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
  keep_ids  <- read_case_list(case_file)
  if (length(keep_ids)) {
    inter <- intersect(sample_cols_all, keep_ids)
    sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
    if (length(inter) < 10) log_msg("提示：case_list 與欄名交集太小（{length(inter)}），改用全部樣本欄位")
  } else {
    sample_cols <- sample_cols_all
  }
  
  m <- dat %>%
    dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
    janitor::remove_empty("cols")
  
  rn <- m$Gene
  m  <- as.matrix(m[,-1, drop=FALSE]); suppressWarnings(storage.mode(m) <- "double")
  rownames(m) <- rn
  if (ncol(m) == 0) stop("讀到 0 個樣本欄位，請檢查檔案欄名")
  
  if (anyDuplicated(rownames(m))) {
    log_msg("偵測到重複基因，對重複 rows 取平均")
    m <- rowsum(m, group = rownames(m), reorder = FALSE) / as.vector(table(rownames(m)))
  }
  log_msg("矩陣維度：{nrow(m)} genes × {ncol(m)} samples")
  m
}

# 3) 讀 sample metadata（data_clinical_sample.txt）
load_local_sample_meta <- function(dir, sample_ids){
  f <- file.path(dir, "data_clinical_sample.txt")
  sample_ids <- as.character(sample_ids)
  if (!file.exists(f)) {
    return(data.frame(SampleID = sample_ids, row.names = sample_ids, stringsAsFactors = FALSE))
  }
  dat <- suppressMessages(
    readr::read_tsv(f, guess_max = 100000, show_col_types = FALSE, comment = "#")
  )
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  id_cols <- intersect(c("SAMPLE_ID","sample_id","Sample_ID","Sample","sample"), names(dat))
  if (length(id_cols) == 0) {
    return(data.frame(SampleID = sample_ids, row.names = sample_ids, stringsAsFactors = FALSE))
  }
  dat <- dplyr::rename(dat, SampleID = !!id_cols[1])
  dat$SampleID <- as.character(dat$SampleID)
  dat2 <- dat[match(sample_ids, dat$SampleID), , drop = FALSE]
  dat2$SampleID <- sample_ids
  rownames(dat2) <- sample_ids
  dat2
}

# 4) 封成 SE/MAE 以沿用後續流程
load_brca_local_as_mae <- function(dir){
  mat  <- load_local_protein_matrix(dir)
  meta <- load_local_sample_meta(dir, colnames(mat))
  stopifnot(identical(rownames(meta), colnames(mat)))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(protein = mat),
    colData = S4Vectors::DataFrame(meta, row.names = rownames(meta))
  )
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(protein = se),
    colData     = S4Vectors::DataFrame(meta, row.names = rownames(meta))
  )
}

## ===== 二、共用工具 =====
extract_proteome <- function(mae) {
  if (inherits(mae, "MultiAssayExperiment")) {
    exps <- MultiAssayExperiment::experiments(mae)
    exp_names <- names(exps)
    pick_exp <- exp_names[grepl("protein|proteomic", exp_names, ignore.case = TRUE)]
    if (length(pick_exp) == 0) pick_exp <- exp_names[1]
    se <- exps[[pick_exp[1]]]
  } else if (inherits(mae, "SummarizedExperiment")) {
    se <- mae
  } else stop("不支援的物件類型")
  an <- SummarizedExperiment::assayNames(se)
  try_names <- c("protein","proteomics","protein_quantification","protein level","protein_level")
  hit <- intersect(try_names, tolower(an))
  nm  <- if (length(hit) > 0) an[match(hit[1], tolower(an))] else an[1]
  prot <- SummarizedExperiment::assay(se, nm) %>% as.matrix()
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  if (!"SampleID" %in% colnames(meta)) {
    cand <- c("SAMPLE_ID","sampleId","sample_id","Sample","sample")
    hitc <- intersect(cand, colnames(meta))
    meta$SampleID <- if (length(hitc) > 0) meta[[hitc[1]]] else rownames(meta)
  }
  list(mat = prot, meta = meta)
}

impute_and_filter <- function(mat, min_frac=0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop=FALSE]
  if (any(is.na(m))) m <- imputeLCMD::impute.MinProb(m, q = 0.01)
  m
}

# 備用：舊版（不再用於分組）
rank_by_limma_t <- function(mat, meta, subunit, hi_lo=0.25,
                            covars_guess=c("batch","Batch","TMT","TMT_batch","ProteomicsBatch","Sex","Gender","Age","Stage")){
  stopifnot(subunit %in% rownames(mat))
  sub <- mat[subunit, ]
  qs <- quantile(sub, probs=c(hi_lo, 1-hi_lo), na.rm=TRUE)
  grp <- dplyr::case_when(sub <= qs[1] ~ "Low", sub >= qs[2] ~ "High", TRUE ~ "Mid")
  keep_samples <- names(grp)[grp != "Mid"]
  if (length(keep_samples) < 6) stop("可用樣本太少，無法做 High/Low 分組")
  grp2 <- factor(grp[keep_samples], levels=c("Low","High"))
  X <- t(mat[, keep_samples, drop=FALSE])
  meta2 <- meta %>% dplyr::filter(SampleID %in% keep_samples)
  meta2 <- meta2[match(keep_samples, meta2$SampleID), , drop=FALSE]
  covars <- intersect(covars_guess, colnames(meta2))
  if (length(covars) > 0) {
    design <- model.matrix(~ 0 + grp2 + ., data = data.frame(grp2=grp2, meta2[, covars, drop=FALSE]))
    colnames(design)[1:2] <- c("Low","High")
  } else { design <- model.matrix(~ 0 + grp2); colnames(design) <- c("Low","High") }
  fit <- limma::lmFit(t(X), design)
  contr <- limma::makeContrasts(High - Low, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr))
  t_stat <- fit2$t[,1]
  stats <- sort(setNames(t_stat, rownames(mat)), decreasing = TRUE)
  stats[!is.na(stats)]
}

rank_by_spearman <- function(mat, subunit){
  stopifnot(subunit %in% rownames(mat))
  sub <- mat[subunit, ]
  rho <- apply(mat, 1, function(x){
    ok <- stats::complete.cases(x, sub)
    if (sum(ok) >= 5) suppressWarnings(stats::cor(x[ok], sub[ok], method="spearman")) else NA_real_
  })
  sort(rho[!is.na(rho)], decreasing = TRUE)
}

# Spearman：用原始值（pairwise-complete），不補值
rank_by_spearman_raw <- function(mat_raw, subunit, min_pairs = 10){
  stopifnot(subunit %in% rownames(mat_raw))
  sub <- mat_raw[subunit, ]
  keep <- !is.na(sub)
  if (sum(keep) < min_pairs) stop("可用樣本太少，無法做 Spearman")
  sub <- sub[keep]
  m   <- mat_raw[, keep, drop = FALSE]
  rho <- apply(m, 1, function(x){
    ok <- stats::complete.cases(x, sub)
    if (sum(ok) >= min_pairs) suppressWarnings(stats::cor(x[ok], sub[ok], method = "spearman")) else NA_real_
  })
  sort(rho[!is.na(rho)], decreasing = TRUE)
}

## ---- 檔名清理（避免 : / 空白 等）----
safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

## ---- 通用 GSEA 跑法（依任何基因集清單）----
run_fgsea_save_generic <- function(stats, genesets, out_prefix, top_plot_n = 6, plot_title = "") {
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  res <- fgsea(pathways = genesets, stats = stats, minSize = minSize, maxSize = maxSize, nperm = nperm) |>
    dplyr::arrange(padj, dplyr::desc(NES))
  data.table::fwrite(as.data.frame(res), paste0(out_prefix, ".csv"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "GSEA")
  openxlsx::writeData(wb, 1, as.data.frame(res))
  openxlsx::saveWorkbook(wb, paste0(out_prefix, ".xlsx"), overwrite = TRUE)
  if (nrow(res) > 0) {
    to_plot <- head(res$pathway, min(top_plot_n, nrow(res)))
    for (pw in to_plot) {
      fn <- paste0(out_prefix, "_", safe_name(pw), ".pdf")
      try({
        pdf(fn, width = 6, height = 4)
        print(plotEnrichment(genesets[[pw]], stats) + ggtitle(glue::glue("{plot_title}\n{pw}")))
        dev.off()
      }, silent = TRUE)
    }
  }
  res
}

## ===== 三、基因集（H + C1–C8，依 category / subcategory 分組；不依賴欄位名稱） =====
log_msg("下載 MSigDB 全人類 collections（H, C1–C8）並依組別分開")

cats <- c("H", paste0("C", 1:8))

# 產生一個清單：names 為 group（如 "H"、"C2__REACTOME"），值為 {pathway -> gene vector}
genesets_by_group <- list()

for (cat in cats) {
  df <- msigdbr(species = "Homo sapiens", category = cat)
  
  # 相容不同版本欄位名稱，找「子集合」欄位
  subcat_candidates <- c("gs_subcat", "subcategory", "sub_category", "gs_subcategory")
  subcat_col <- subcat_candidates[subcat_candidates %in% names(df)]
  if (length(subcat_col) == 0) {
    df$.__subcat__ <- NA_character_
  } else {
    df$.__subcat__ <- df[[subcat_col[1]]]
  }
  
  # 建立 group 名稱（有 subcat 就 "C2__REACTOME"，否則就是 "H" 或 "C3"）
  df$.__group__ <- ifelse(is.na(df$.__subcat__) | df$.__subcat__ == "",
                          cat, paste0(cat, "__", df$.__subcat__))
  
  # 依 group 分拆，並轉成 pathway -> genes 的清單
  split_by_group <- split(df, df$.__group__)
  for (grp in names(split_by_group)) {
    dfg <- split_by_group[[grp]]
    # 注意 pathway 欄位固定是 gs_name；gene 欄位固定是 gene_symbol
    gs_list <- lapply(split(dfg$gene_symbol, dfg$gs_name), unique)
    genesets_by_group[[grp]] <- gs_list
  }
}

log_msg("完成載入 {length(genesets_by_group)} 個 gene set groups：{paste(head(names(genesets_by_group), 10), collapse=', ')}{if (length(genesets_by_group)>10) ' ...' else ''}")


## ===== 四、執行：只跑本地 BRCA =====
study_id   <- "brca_cptac_2020"
study_name <- "Proteogenomic landscape of breast cancer (CPTAC, Cell 2020)"
log_msg("== 研究 {study_id} : {study_name} ==")

out_root <- file.path(project_dir, "results", study_id)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# 1) 載入本地 MAE
mae <- load_brca_local_as_mae(local_brca_dir)

# 2) 取出 protein 與 metadata
ex   <- extract_proteome(mae)
mat0 <- ex$mat; meta <- ex$meta

# 3) 對齊樣本（保守檢查）
common <- intersect(colnames(mat0), meta$SampleID)
if (length(common) < 10) stop("可對齊樣本過少")
mat0 <- mat0[, common, drop=FALSE]
meta <- meta %>% dplyr::filter(SampleID %in% common)

# 4) 若值很大，視為未 log2，做 log2(x+1)
mx <- suppressWarnings(max(mat0, na.rm=TRUE))
if (is.finite(mx) && mx > 100) {
  log_msg("偵測到非 log2 標度，進行 log2(x+1) 轉換")
  mat0 <- log2(mat0 + 1)
}

# 5) 缺失過濾 + MNAR 填補（僅供 limma 使用；分組與 Spearman 用 mat0 原值）
log_msg("缺失過濾與填補")
mat <- impute_and_filter(mat0, min_frac = min_frac_complete)

# 6) 逐 CSN subunit 跑 GSEA（分組用 mat0 原值；limma 用補值矩陣但排除該 subunit；Spearman 用 raw）
present_sub <- intersect(csn_subunits, rownames(mat0))
missing_sub <- setdiff(csn_subunits, rownames(mat0))
if (length(missing_sub)) log_msg("以下次單元在原始矩陣不存在、此研究將跳過：{paste(missing_sub, collapse=', ')}")
if (length(present_sub) == 0) stop("找不到任何 CSN 次單元於原始矩陣 rownames")

for (su in csn_subunits) {
  if (!(su %in% present_sub)) next  # 原始檔就沒有 → 跳過（跨研究彙整會是 NA）
  
  log_msg("  >> 次單元：{su}")
  out_dir <- file.path(out_root, su)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## --- 1) 用原始值分組（不補值；保留樣本名稱） ---
  sub_raw <- mat0[su, ]
  v <- sub_raw[!is.na(sub_raw)]                    # 非 NA 值（帶 names=樣本ID）
  if (length(v) < (2 * min_per_group)) {
    log_msg("    可用樣本過少（非NA={length(v)}），略過 limma")
  } else {
    qs <- quantile(v, probs = c(hi_lo_quantile, 1 - hi_lo_quantile), names = FALSE)
    grp <- setNames(rep("Mid", length(sub_raw)), names(sub_raw))   # 關鍵：保留 names
    grp[names(v)[v <= qs[1]]] <- "Low"
    grp[names(v)[v >= qs[2]]] <- "High"
    
    keep_samples <- names(sub_raw)[grp != "Mid" & !is.na(sub_raw)]
    n_high <- sum(grp[keep_samples] == "High")
    n_low  <- sum(grp[keep_samples] == "Low")
    
    ## 若分位因 ties 造成 High/Low=0，改用排名 fallback 取前/後 k 名
    if (n_high == 0 || n_low == 0) {
      k <- max(min_per_group, ceiling(hi_lo_quantile * length(v)))
      ord_asc  <- names(sort(v, decreasing = FALSE))
      ord_desc <- names(sort(v, decreasing = TRUE))
      idx_low  <- head(ord_asc,  k)
      idx_high <- head(ord_desc, k)
      grp <- setNames(rep("Mid", length(sub_raw)), names(sub_raw))
      grp[idx_low]  <- "Low"
      grp[idx_high] <- "High"
      keep_samples <- union(idx_low, idx_high)
      n_high <- sum(grp[keep_samples] == "High")
      n_low  <- sum(grp[keep_samples] == "Low")
      log_msg("    量化分位遇到 ties，已改用排名 fallback：High={n_high}, Low={n_low}")
    }
    
    ## --- 2) limma：樣本足夠 → 用補值矩陣且移除該 subunit 列；並對每個 gene set group 各自輸出 ---
    if (length(keep_samples) < (2 * min_per_group) || n_high < min_per_group || n_low < min_per_group) {
      log_msg("    樣本不足（High={n_high}, Low={n_low}，門檻各≥{min_per_group}），略過 limma")
    } else {
      mat_limma <- mat[setdiff(rownames(mat), su), keep_samples, drop = FALSE]
      meta2     <- meta %>% dplyr::filter(SampleID %in% keep_samples)
      meta2     <- meta2[match(keep_samples, meta2$SampleID), , drop = FALSE]
      
      stats_t <- tryCatch({
        grp2 <- factor(ifelse(grp[keep_samples] == "High", "High", "Low"), levels = c("Low","High"))
        covars_guess <- c("batch","Batch","TMT","TMT_batch","ProteomicsBatch","Sex","Gender","Age","Stage")
        covars <- intersect(covars_guess, colnames(meta2))
        design <- if (length(covars) > 0) {
          model.matrix(~ 0 + grp2 + ., data = data.frame(grp2 = grp2, meta2[, covars, drop = FALSE]))
        } else model.matrix(~ 0 + grp2)
        colnames(design)[1:2] <- c("Low","High")
        fit  <- limma::lmFit(mat_limma, design)
        fit2 <- limma::eBayes(limma::contrasts.fit(fit, limma::makeContrasts(High - Low, levels = design)))
        t_stat <- fit2$t[, 1]
        sort(t_stat[!is.na(t_stat)], decreasing = TRUE)
      }, error = function(e){ log_msg("    limma 排名失敗：{e$message}"); NULL })
      
      if (!is.null(stats_t)) {
        for (grp_label in names(genesets_by_group)) {
          gs <- genesets_by_group[[grp_label]]
          subgrp_dir  <- file.path(out_dir, safe_name(grp_label))
          dir.create(subgrp_dir, showWarnings = FALSE, recursive = TRUE)
          out_prefix  <- file.path(subgrp_dir, "GSEA_limma_t")
          run_fgsea_save_generic(
            stats = stats_t,
            genesets = gs,
            out_prefix = out_prefix,
            top_plot_n = 6,
            plot_title = glue::glue("{study_id} | {su} | High vs Low (Q{round(hi_lo_quantile*100)} vs Q{round((1-hi_lo_quantile)*100)}) | {grp_label}")
          )
        }
        saveRDS(stats_t, file.path(out_dir, "rank_stats_limma_t.rds"))
      }
    }
  }
  
  ## --- 3) Spearman：用原始值（pairwise-complete），不補值；對每個 group 各自輸出 ---
  stats_rho <- tryCatch(
    rank_by_spearman_raw(mat0, su, min_pairs = min_pairs_spearman),
    error = function(e){ log_msg("    Spearman 無法進行：{e$message}"); NULL }
  )
  if (!is.null(stats_rho)) {
    for (grp_label in names(genesets_by_group)) {
      gs <- genesets_by_group[[grp_label]]
      subgrp_dir  <- file.path(out_dir, safe_name(grp_label))
      dir.create(subgrp_dir, showWarnings = FALSE, recursive = TRUE)
      out_prefix  <- file.path(subgrp_dir, "GSEA_spearman_raw")
      run_fgsea_save_generic(
        stats = stats_rho,
        genesets = gs,
        out_prefix = out_prefix,
        top_plot_n = 4,
        plot_title = glue::glue("{study_id} | {su} | Spearman (raw, pairwise-complete) | {grp_label}")
      )
    }
    saveRDS(stats_rho, file.path(out_dir, "rank_stats_spearman_raw.rds"))
  }
}

log_msg("BRCA 本地流程完成，輸出位於：{file.path(out_root)}")
