## =========================================================
## CSN subunits → proteomic GSEA (CPTAC/TCGA, 加入 TP53 分層)
##  - 依 strata = c("ALL","TP53_mutant","TP53_wild_type") 輸出到
##    <dataset>/csn_gsea_results_TP53/<STRATUM>/...
##  - 另做 pan-cancer 統整與兩種整理表格到 csn_gsea_pan_summary_TP53/
## =========================================================

setwd("C:/Users/danny/Documents/R_project/CSN_CPTAC") ## YTC laptop
##setwd("C:/Users/cmuh/Documents/YenTing_document/CSN_CPTAC") ## lab computer
getwd()

## ===== 在「套件」區塊多加 yaml（1）=====
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(janitor); library(glue)
  library(SummarizedExperiment); library(MultiAssayExperiment); library(S4Vectors)
  library(limma); library(imputeLCMD); library(msigdbr); library(fgsea)
  library(openxlsx); library(ComplexHeatmap); library(cowplot); library(matrixStats)
  library(yaml)   # <— 新增：用來寫 run_manifest.yml
})

## ===== global helper: opt (must exist before any writer uses it) =====
if (!exists("opt", mode = "function")) {
  opt <- function(nm, default) {
    if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  }
}


## ===== helper: coerce_covariates_safely (must be global, defined before use) =====
if (!exists("coerce_covariates_safely", mode = "function")) {
  coerce_covariates_safely <- function(df){
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)
    
    for (cn in colnames(df)) {
      v <- df[[cn]]
      if (is.factor(v) || is.character(v) || is.logical(v)) {
        v <- factor(v)
        lv <- levels(droplevels(v[!is.na(v)]))
        if (length(lv) <= 1) {            # 單一層級因子 → 丟掉，避免 contrasts 錯誤
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
      } else {
        df[[cn]] <- suppressWarnings(as.numeric(v))  # 數值保持 numeric
        # 如需一併濾掉「完全常數」數值欄，解除下列註解：
        # vv <- df[[cn]]
        # if (sum(!is.na(vv)) >= 2 && isTRUE(all(vv[!is.na(vv)] == vv[which(!is.na(vv))[1]]))) {
        #   keep[cn] <- FALSE
        #   if (exists("logf")) try(logf("  [covars] drop constant numeric: %s", cn), silent = TRUE)
        # }
      }
    }
    df <- df[, keep, drop = FALSE]
    df
  }
}


## ==== 強制單執行緒、關閉所有平行後端（跨常見框架）====
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
  if (requireNamespace("doParallel", quietly = TRUE)) {
    # 嘗試關閉殘留 cluster（若你有把 cluster 存在 .GlobalEnv 可自行 stopCluster）
    # 這裡先不硬停未知物件，改用 DoSEQ 即可。
  }
  
  # future
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  Sys.setenv("R_FUTURE_FORK_ENABLE" = "FALSE")
  
  # fgsea：有些版本支援 nproc，這裡只設定一個選項供 wrapper 讀取
  options(.fgsea_nproc = 1L)
}
.force_serial_execution()

## 建立 run_info 目錄（寫任何 run_info/* 檔前）
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)


## 附上多重假設層級說明（更新：加入 meta-FDR）
cat(paste(
  "Multiple-testing policy:",
  " - Within each gene-set group and per statistic, we control FDR using fgsea padj (per-dataset/stratum).",
  " - Pan-cancer summaries aggregate directions and counts across datasets/strata (descriptive).",
  " - Additionally, we provide a simple meta-analysis across datasets using Stouffer's z to obtain meta p-values,",
  "   followed by Benjamini-Hochberg to report meta-level FDR (see csn_gsea_pan_summary_TP53/meta_fdr).",
  sep = "\n"),
  file = file.path("run_info","analysis_notes.txt"))


## ===== 參數 =====
set.seed(1234)

csn_subunits <- c("GPS1","COPS2","COPS3","COPS4","COPS5","COPS6","COPS7A","COPS7B","COPS8","COPS9")

min_frac_complete   <- 0.75     # 每基因至少 75% 非 NA
hi_lo_quantile      <- 0.25     # 高低各 25%
minSize <- 5; maxSize <- 1000
fgsea_eps <- 1e-10              # fgseaMultilevel() 精度

min_per_group      <- 8         # High/Low 各至少多少樣本（不足就跳過 limma）
min_pairs_spearman <- 10        # Spearman 至少配對樣本數
MAKE_PLOTS <- FALSE

## ---- 全域開關：是否跑 High/Low 敏感度（預設 FALSE；主分析用連續 predictor） ----
RUN_HILO_SENSITIVITY <- FALSE

## ---- [USER CONFIG | 內嵌設定，不用外部 options()] ----
## 依需求修改下兩行：
## 1) 要跑哪些 pass：
##    - 只跑 BatchAdj：RUN_PASSES <- c("BatchAdj")    (預設)
##    - 只跑 RAW：     RUN_PASSES <- "RAW"
##    - 兩者都跑：     RUN_PASSES <- c("RAW","BatchAdj")
RUN_PASSES <- c("BatchAdj")

## 2) 是否跑 camera / mroast 輸出（GSEA_camera*.csv、GSEA_mroast*.csv）：
##    - 關閉（預設）：FALSE
##    - 開啟：        TRUE
RUN_CAMERA_MROAST <- FALSE

## 將上面的內嵌設定寫入 options，供程式其他地方的 getOption(...) 讀取
options(csn.run_passes = RUN_PASSES)
options(csn.run_camera_mroast = RUN_CAMERA_MROAST)
## ---- [END USER CONFIG] ----

## ---- [USER CONFIG | heatmap top/bottom limits] ----
## 非 H 的 collection 熱圖篩選門檻（單一 dataset 用 NES，PAN 用 Z；皆以 CSN_SCORE 為篩選依據）
## - 單一 dataset：只保留 padj<0.05 的 NES 前 N 大與前 M 小（依 CSN_SCORE）
## - PAN（meta-FDR）：只保留 padj_meta<0.05 的 Z 前 N 大與前 M 小（依 CSN_SCORE）
DATASET_HEATMAP_TOP_N    <- 25   # 預設前 25 大（NES）
DATASET_HEATMAP_BOTTOM_N <- 25   # 預設前 25 小（NES）
PAN_HEATMAP_TOP_N        <- 25   # 預設前 25 大（Z）
PAN_HEATMAP_BOTTOM_N     <- 25   # 預設前 25 小（Z）
## ---- [END USER CONFIG | heatmap top/bottom limits] ----


## ---- [USER CONFIG | heatmap collections] ----
## 想畫哪些 collections 的熱圖（字串向量；例：c("H", "C6", "C2:CP:BIOCARTA")）
## - 單一 dataset 熱圖：留空或 NULL 表示「全部能抓到的 collection 都畫」
## - PAN 熱圖：留空或 NULL 表示「全部能抓到的 collection 都畫」
PLOT_DATASET_COLLECTIONS <- NULL
PLOT_PAN_COLLECTIONS     <- NULL
## ---- [END USER CONFIG | heatmap collections] ----


## ---- Pipeline selection（limma t / spearman / both；預設 limma t）----
## 用法：
PIPELINES_TO_RUN <- c("limma_t")               # 只跑 limma t（預設）
##   PIPELINES_TO_RUN <- c("spearman")              # 只跑 spearman（含 diffcorr）
##   PIPELINES_TO_RUN <- c("limma_t","spearman")    # 兩者都跑
PIPELINES_TO_RUN <- get0("PIPELINES_TO_RUN", ifnotfound = c("limma_t"))
.RUN_LIMMA     <- any(tolower(PIPELINES_TO_RUN) %in% c("limma_t","limma","both","all"))
.RUN_SPEARMAN  <- any(tolower(PIPELINES_TO_RUN) %in% c("spearman","rho","both","all"))

log_msg <- function(text, ..., .envir = parent.frame()) {
  ts <- format(Sys.time(), "%H:%M:%S")
  msg <- tryCatch({
    if (grepl("\\{[^}]+\\}", text)) {
      # 有 { } → 用 glue
      glue::glue(text, ..., .envir = .envir)
    } else if (grepl("%", text)) {
      # 有 % → 用 sprintf
      do.call(sprintf, c(list(fmt = text), list(...)))
    } else {
      # 純文字；若有額外引數也當 sprintf 格式用
      if (nargs() > 1L) do.call(sprintf, c(list(fmt = text), list(...))) else text
    }
  }, error = function(e) text)
  cat(sprintf("[%s] %s\n", ts, msg))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_fs_name <- function(s) {
  s <- gsub('[<>:"/\\\\|?*]', "_", s)
  s <- gsub("\\s+", "_", s)
  s
}

## ===== Gene set 選擇（支援所有 MSigDB collections/subcollections）=====
## 可選 token（不分大小寫）：
##   "H",
##   "C1",
##   "C2", "C2:REACTOME", "C2:BIOCARTA", "C2:CGP", "C2:KEGG"  # 註：新版多數無 KEGG
##   "C3", "C3:TFT", "C3:MIR",
##   "C4", "C4:CGN", "C4:CM",
##   "C5", "C5:BP", "C5:CC", "C5:MF",
##   "C6", "C7", "C8"



#GENESET_GROUPS_TO_RUN <- c("H")  # ← 需要時自行改，例如 c("H","C5:BP","C2:REACTOME","C7")
#GENESET_GROUPS_TO_RUN <- c("H", "C5:BP", "C2:REACTOME", "C3:TFT", "C7")
#GENESET_GROUPS_TO_RUN <- c("H","C1","C2","C3","C4","C5","C6","C7","C8")
#GENESET_GROUPS_TO_RUN <- c("C5:BP","C5:CC","C5:MF")
#GENESET_GROUPS_TO_RUN <- c("H","C3:MIR","C2:CP:REACTOME","C5:GO:BP")

GENESET_GROUPS_TO_RUN <- c("H")



## ===== Gene sets（只建置被選到的群組）=====
log_msg("準備 MSigDB gene sets；將依 GENESET_GROUPS_TO_RUN 建置：%s",
        paste(GENESET_GROUPS_TO_RUN, collapse = ", "))
# (Optional) 檢視並輸出目前 msigdbr 可用的 collections / subcollections（含 geneset 計數）
# 放置點：上方 log_msg(...) 之後、genesets_by_group <- list() 之前
df0 <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)

if (!is.null(df0) && NROW(df0) > 0) {
  df0 <- as.data.frame(df0)  # 避免 tibble 子集嚴格行為
  
  ## 1) 偵測實際欄位名稱（2025.1 為 gs_collection / gs_subcollection）
  cat_candidates <- c("gs_collection", "gs_cat", "gs_category", "category", "collection")
  sub_candidates <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
  name_candidates <- c("gs_name", "geneset_name", "set_name")
  
  cat_hits  <- intersect(cat_candidates,  names(df0))
  sub_hits  <- intersect(sub_candidates,  names(df0))
  name_hits <- intersect(name_candidates, names(df0))
  
  if (length(cat_hits) == 0 || length(name_hits) == 0) {
    log_msg("（略過）找不到必要欄位（collection 或 gs_name）；未輸出 collections 清單。")
  } else {
    catcol  <- cat_hits[1]
    subcol  <- if (length(sub_hits)  >= 1) sub_hits[1]  else NULL
    namecol <- name_hits[1]
    
    ## 2) 組資料框並計數「唯一 geneset（gs_name）數量」
    tmp <- data.frame(
      gs_cat    = as.character(df0[[catcol]]),
      gs_subcat = if (!is.null(subcol)) as.character(df0[[subcol]]) else NA_character_,
      gs_name   = as.character(df0[[namecol]]),
      stringsAsFactors = FALSE
    )
    
    counts <- stats::aggregate(gs_name ~ gs_cat + gs_subcat, data = tmp,
                               FUN = function(x) length(unique(x)))
    names(counts)[3] <- "n_genesets"
    
    ## 3) 排序（把 NA 子類別視為空字串避免 order 遇 NA）
    ord <- order(counts$gs_cat, ifelse(is.na(counts$gs_subcat), "", counts$gs_subcat))
    counts <- counts[ord, , drop = FALSE]
    
    ## 4) 輸出清單到 run_info
    dir.create("run_info", showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(counts, file.path("run_info", "msigdb_available_collections.csv"))
    log_msg("已輸出 msigdb 可用 collections（含 geneset 計數）至 run_info/msigdb_available_collections.csv")
  }
} else {
  log_msg("（略過）msigdbr() 回傳為空或呼叫失敗；未輸出 collections 清單。")
}

genesets_by_group <- list()

.get_subcat_col <- function(df){
  cand <- c("gs_subcollection", "gs_subcat", "gs_subcategory", "subcategory", "sub_category")
  hit  <- cand[cand %in% names(df)]
  if (length(hit)) hit[1] else NULL
}


.fetch_msig <- function(cat) {
  try(msigdbr(species = "Homo sapiens", category = toupper(cat)), silent = TRUE)
}

.parse_group_token <- function(tok){
  t <- toupper(trimws(tok))
  sp <- unlist(strsplit(t, ":", fixed = TRUE))
  list(cat = sp[1], sub = if (length(sp) >= 2) paste(sp[-1], collapse=":") else NULL, raw = t)
}

.make_group_label <- function(cat, df, want_sub){
  subcol <- .get_subcat_col(df)
  if (!is.null(want_sub)) {
    if (!is.null(subcol)) {
      uniq <- unique(df[[subcol]])
      hit  <- uniq[grepl(want_sub, uniq, ignore.case = TRUE)]
      if (length(hit) >= 1) return(sprintf("%s__%s", cat, hit[1]))
    }
    return(sprintf("%s__%s", cat, safe_fs_name(want_sub)))
  }
  return(cat)
}

for (tok in GENESET_GROUPS_TO_RUN) {
  pp  <- .parse_group_token(tok)
  cat <- pp$cat
  df  <- suppressMessages(.fetch_msig(cat))
  if (inherits(df, "try-error") || is.null(df) || !nrow(df)) {
    log_msg("  注意：msigdbr 無法提供分類 %s；略過 %s", cat, pp$raw)
    next
  }
  
  if (!is.null(pp$sub)) {
    subcol <- .get_subcat_col(df)
    subtxt <- if (!is.null(subcol)) as.character(df[[subcol]]) else ""
    # 只在子類別字串的「開頭」匹配（例：MIR、CP:REACTOME、GO:BP）
    # 並且允許後面接 ":" 或結束
    esc  <- function(s) gsub("([\\^\\$\\.|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])","\\\\\\1", s)
    pattern <- paste0("^", esc(pp$sub), "($|:)")
    keep <- grepl(pattern, subtxt, ignore.case = TRUE)
    df   <- df[keep, , drop = FALSE]
    if (!nrow(df)) { 
      log_msg("  注意：%s 找不到任何 %s 的子集合；略過", cat, pp$sub)
      next 
    }
    grp <- .make_group_label(cat, df, pp$sub)
  } else {
    grp <- cat
  }
  
  genesets_by_group[[grp]] <- lapply(split(df$gene_symbol, df$gs_name), unique)
  log_msg("  取得 %s：%d 個集合", grp, length(genesets_by_group[[grp]]))
}

if (!length(genesets_by_group)) {
  stop("沒有可用的 gene-set。請檢查 GENESET_GROUPS_TO_RUN 設定。")
}
log_msg("最終 gene-set groups：%s", paste(names(genesets_by_group), collapse = ", "))


## ===== 讀檔與共用工具 =====
read_case_list <- function(path_file){
  if (!file.exists(path_file)) return(character(0))
  x <- readLines(path_file, warn=FALSE, encoding="UTF-8")
  line <- x[grepl("^case_list_ids:", x)]
  if (!length(line)) return(character(0))
  ids <- sub("^case_list_ids:\\s*", "", line[1])
  ids <- unlist(strsplit(ids, "[,\\s]+"))
  unique(ids[nchar(ids)>0])
}

load_matrix_from_dataset_dir <- function(dir){
  fp <- file.path(dir, "data_protein_quantification.txt")
  if (!file.exists(fp)) stop(glue::glue("找不到檔案：{fp}"))
  log_msg("讀取 protein 矩陣：{basename(fp)}")
  dat <- suppressMessages(readr::read_tsv(fp, guess_max=200000, show_col_types=FALSE))
  gene_cols <- c("Hugo_Symbol","hugo_symbol","Gene","Gene_Symbol","HugoSymbol","GENE_SYMBOL","gene","gene_symbol")
  gcol <- intersect(gene_cols, names(dat)); if (!length(gcol)) gcol <- names(dat)[1]
  dat  <- dplyr::rename(dat, Gene = !!gcol[1]); dat$Gene <- sub("\\|.*$","",dat$Gene)
  not_sample <- c("Gene","Entrez_Gene_Id","Entrez_Gene_Id.","ENTREZ_GENE_ID","Description","Gene_Name","GeneName","Gene_Symbol")
  sample_cols_all <- setdiff(names(dat), not_sample)
  # 若資料夾有 case_list 就盡量比對
  case_file <- file.path(dir, "case_lists", "cases_protein_quantification.txt")
  keep_ids  <- read_case_list(case_file)
  if (length(keep_ids)) {
    inter <- intersect(sample_cols_all, keep_ids)
    sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
    if (length(inter) < 10) log_msg("提示：case_list 交集太小（{length(inter)}），改用全部樣本欄位")
  } else sample_cols <- sample_cols_all
  m <- dat %>% dplyr::select(Gene, dplyr::all_of(sample_cols)) %>% janitor::remove_empty("cols")
  rn <- m$Gene; m <- as.matrix(m[,-1, drop=FALSE]); storage.mode(m) <- "double"; rownames(m) <- rn
  if (ncol(m)==0) stop("讀到 0 個樣本欄位")
  if (anyDuplicated(rownames(m))) {
    log_msg("偵測到重複基因，對重複 rows 取平均")
    m <- rowsum(m, group=rownames(m), reorder=FALSE) / as.vector(table(rownames(m)))
  }
  log_msg("矩陣維度：{nrow(m)} genes × {ncol(m)} samples")
  m
}

write_geneset_manifest <- function(genesets_by_group, out_csv = file.path("run_info","geneset_manifest.csv")){
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  ver <- tryCatch(as.character(utils::packageVersion("msigdbr")), error=function(e) NA_character_)
  rows <- lapply(names(genesets_by_group), function(g){
    gs <- genesets_by_group[[g]]
    if (is.null(gs) || !length(gs)) return(NULL)
    data.frame(group = g,
               pathway = names(gs),
               genes_n = vapply(gs, function(v) length(unique(v)), integer(1)),
               stringsAsFactors = FALSE)
  })
  df <- dplyr::bind_rows(rows)
  if (!is.null(df)) data.table::fwrite(df, out_csv)
  return(ver)
}

spearman_prerank <- function(mat, predictor, min_non_na = 16L) {
  stopifnot(!is.null(colnames(mat)))
  sam <- colnames(mat)
  if (is.null(names(predictor))) names(predictor) <- sam
  predictor <- predictor[sam]
  
  ok <- is.finite(predictor)
  if (sum(ok) < min_non_na) {
    stop(sprintf("predictor 非NA 樣本不足（%d < %d）", sum(ok), min_non_na))
  }
  mat <- mat[, ok, drop = FALSE]
  predictor <- predictor[ok]
  
  # gene-wise Spearman rho
  rho <- apply(mat, 1L, function(v) suppressWarnings(
    suppressWarnings(stats::cor(v, predictor, method = "spearman",
                                use = "pairwise.complete.obs"))
  ))
  rho <- rho[is.finite(rho)]
  sort(rho, decreasing = TRUE)   # 給 fgsea 的 named vector
}

## --- Partial Spearman（殘差化→Spearman）---
## 讓 Spearman-BatchAdj 與 limma「設計式內納入 batch/covars」的理念一致
spearman_prerank_partial <- function(M_ba, DF_ba,
                                     min_non_na = 16L,
                                     residualize_predictor = TRUE) {
  stopifnot(is.matrix(M_ba), !is.null(colnames(M_ba)),
            is.data.frame(DF_ba), "predictor" %in% colnames(DF_ba))
  # 僅保留 DF_ba 與 M_ba 的共同樣本（已在上游對齊，這裡再保險一次）
  sam <- intersect(colnames(M_ba), rownames(DF_ba))
  if (length(sam) < min_non_na) {
    stop(sprintf("[partialSpearman] 可用樣本不足（%d < %d）", length(sam), min_non_na))
  }
  M_ba <- M_ba[, sam, drop = FALSE]
  DF_ba <- DF_ba[sam, , drop = FALSE]
  
  # nuisance 設計：把 predictor 拿掉，只留 batch + 其它 covariates
  nui_df <- DF_ba[, setdiff(colnames(DF_ba), "predictor"), drop = FALSE]
  # 若沒有任何 nuisance，退回原 spearman_prerank（理論上都有 batch/covars）
  if (!ncol(nui_df)) {
    pred <- DF_ba$predictor; names(pred) <- rownames(DF_ba)
    return(spearman_prerank(M_ba, predictor = pred, min_non_na = min_non_na))
  }
  # 0+ 去掉截距，讓下游殘差化工具自行加截距（避免重複）
  des_nui <- stats::model.matrix(
    as.formula(paste("~ 0 +", paste(colnames(nui_df), collapse = " + "))),
    data = nui_df
  )
  
  
  # --- 殘差化表達矩陣：Y = samples × genes，用 orthogonalize_to 一次性殘差化 ---
  #   你的程式已經內建 orthogonalize_to(Y, nuisance) 可用（會自動加截距）:contentReference[oaicite:1]{index=1}
  Y <- t(M_ba)                                             # samples × genes
  Y_res <- orthogonalize_to(Y, des_nui)                    # samples × genes 的殘差
  M_res <- t(Y_res)                                        # 回到 genes × samples
  
  # --- 殘差化 predictor（可選；預設 TRUE）---
  p <- DF_ba$predictor; names(p) <- rownames(DF_ba)
  if (isTRUE(residualize_predictor)) {
    Xp <- cbind("(Intercept)" = 1, des_nui)
    fitp <- lm.fit(x = Xp, y = as.numeric(p))
    p_res <- rep(NA_real_, nrow(Xp)); p_res[] <- fitp$residuals
    names(p_res) <- rownames(Xp)
    p_use <- p_res
  } else {
    p_use <- as.numeric(p)
    names(p_use) <- rownames(DF_ba)
  }
  
  # --- gene-wise Spearman（pairwise.complete.obs）---
  rho <- apply(M_res, 1L, function(v)
    suppressWarnings(stats::cor(as.numeric(v), p_use,
                                method = "spearman",
                                use    = "pairwise.complete.obs")))
  rho <- rho[is.finite(rho)]
  sort(rho, decreasing = TRUE)  # preranked vector for fgsea
}
run_and_write_gsea_spearman <- function(stats_named_vec, pathways,
                                        group_tag, out_dir_root, subunit, pass_label) {
  stopifnot(is.numeric(stats_named_vec), !is.null(names(stats_named_vec)))
  out_dir <- file.path(out_dir_root, subunit)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  minSize <- opt("minSize", 15L)
  maxSize <- opt("maxSize", 500L)
  
  orig_n <- length(pathways)
  pw_use <- ._intersect_and_filter_pathways(pathways, names(stats_named_vec),
                                            minSize = minSize, maxSize = maxSize)
  message(sprintf("[gsea:%s] %s |%s| 原=%d → 用=%d",
                  pass_label, subunit, group_tag, orig_n, length(pw_use)))
  
  if (!length(pw_use)) {
    out_csv <- file.path(out_dir, "GSEA_spearman.csv")
    data.table::fwrite(data.frame(), out_csv)
    message(sprintf("[spearman:%s] %s -> %s (n_pathways=%d)", pass_label, subunit, out_csv, 0L))
    return(invisible(NULL))
  }
  
  set.seed(1L)
  res <- tryCatch({
    suppressWarnings(fgsea::fgseaMultilevel(
      pathways = pw_use, stats = stats_named_vec,
      minSize = minSize, maxSize = maxSize,
      eps = 1e-10, nproc = getOption(".fgsea_nproc", 1L)
    ))
  }, error = function(e) {
    message(sprintf("[spearman:%s] fgseaMultilevel failed: %s; fallback to fgseaSimple",
                    pass_label, conditionMessage(e)))
    suppressWarnings(fgsea::fgseaSimple(
      pathways = pw_use, stats = stats_named_vec,
      nperm = 10000, minSize = minSize, maxSize = maxSize
    ))
  })
  res <- as.data.frame(res)
  out_csv <- file.path(out_dir, "GSEA_spearman.csv")
  data.table::fwrite(res, out_csv)
  message(sprintf("[spearman:%s] %s -> %s (n_pathways=%d)", pass_label, subunit, out_csv, nrow(res)))
  invisible(res)
}


# 供 PCA 前清理：移除 Inf→NA、保留至少 min_samples 個有限值的列；行內中位數插補 NA
.clean_for_pca <- function(X, min_samples = 10L, min_genes = 5L){
  X <- as.matrix(X); X[!is.finite(X)] <- NA
  keep_rows <- rowSums(is.finite(X)) >= min_samples
  if (!any(keep_rows)) return(NULL)
  X <- X[keep_rows, , drop = FALSE]
  if (nrow(X) < min_genes) return(NULL)
  if (anyNA(X)) {
    med <- apply(X, 1, function(r) median(r[is.finite(r)], na.rm = TRUE))
    for (i in seq_len(nrow(X))) {
      xi <- X[i, ]; xi[!is.finite(xi)] <- med[i]; X[i, ] <- xi
    }
  }
  X
}

# 安全版 CSN SCORE（先試原法；失敗再 fallback）
build_csn_score_safe <- function(mat0, subunits, combine_7AB = TRUE,
                                 min_members = 5L, pca_min_samples = 10L){
  sub <- intersect(subunits, rownames(mat0))
  out_na <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (length(sub) < min_members) return(out_na)
  
  cs_try <- try({
    build_csn_score(mat0, subunits = sub, combine_7AB = combine_7AB, min_members = min_members)
  }, silent = TRUE)
  
  if (!inherits(cs_try, "try-error") && sum(is.finite(cs_try)) >= pca_min_samples) return(cs_try)
  
  X <- .clean_for_pca(mat0[sub, , drop = FALSE], min_samples = pca_min_samples, min_genes = min_members)
  if (is.null(X)) {
    if (exists("log_msg", mode="function")) log_msg("[CSN_SCORE-safe] 可用子單元或樣本不足 → 全 NA")
    return(out_na)
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode="function")) log_msg("[CSN_SCORE-safe] prcomp 失敗 → 全 NA")
    return(out_na)
  }
  sc <- pc$x[,1]; names(sc) <- rownames(pc$x)
  ref <- colMeans(X, na.rm = TRUE)
  rr  <- suppressWarnings(stats::cor(sc, ref, use = "pairwise.complete.obs"))
  if (is.finite(rr) && rr < 0) sc <- -sc
  out <- out_na; out[names(sc)] <- as.numeric(sc)
  if (exists("log_msg", mode="function")) {
    varpc1 <- if (!is.null(pc$sdev)) (pc$sdev[1]^2) / sum(pc$sdev^2) else NA_real_
    log_msg("[CSN_SCORE-safe] fallback：genes=%d；PC1%%=%.1f；nonNA=%d/%d",
            nrow(X), 100*varpc1, sum(is.finite(out)), length(out))
  }
  out
}

# 安全版審核：先試原本 audit；失敗再用乾淨矩陣跑 prcomp 只做紀錄，不中止流程
audit_csn_score_feasibility_safe <- function(ds_id, stratum, mat0, present_sub,
                                             min_members = 5L, pca_min_samples = 10L,
                                             min_per_group = get0("min_per_group", ifnotfound = 8L),
                                             out_dir = file.path("run_info","csn_score_audit")) {
  ok <- try({
    audit_csn_score_feasibility(
      ds_id    = ds_id,
      stratum  = stratum,
      mat0     = mat0,
      present_sub = present_sub,
      min_members     = min_members,
      pca_min_samples = pca_min_samples,
      min_per_group   = min_per_group,
      out_dir = out_dir
    )
    TRUE
  }, silent = TRUE)
  if (!inherits(ok, "try-error")) return(invisible(TRUE))
  
  X <- .clean_for_pca(mat0[intersect(present_sub, rownames(mat0)), , drop = FALSE],
                      min_samples = pca_min_samples, min_genes = min_members)
  if (is.null(X)) {
    if (exists("log_msg", mode="function"))
      log_msg("[CSN-audit-safe] %s | %s：可用子單元或樣本不足，略過 audit（不中止）", ds_id, stratum)
    return(invisible(FALSE))
  }
  pc <- try(stats::prcomp(t(X), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pc, "try-error")) {
    if (exists("log_msg", mode="function"))
      log_msg("[CSN-audit-safe] %s | %s：fallback prcomp 仍失敗，略過（不中止）", ds_id, stratum)
    return(invisible(FALSE))
  }
  varpc1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
  if (exists("log_msg", mode="function"))
    log_msg("[CSN-audit-safe] %s | %s：fallback OK；genes=%d；PC1%%=%.1f",
            ds_id, stratum, nrow(X), 100*varpc1)
  invisible(TRUE)
}

# --- helper: 確保 limma t-stat 有正確 gene names ---
._ensure_stats_names <- function(stats, gene_names, label = NULL) {
  # stats: numeric 向量（limma t 等）
  # gene_names: 對應矩陣的 rownames（基因 ID）
  if (is.null(stats)) return(NULL)
  v <- suppressWarnings(as.numeric(stats))
  if (length(v) != length(gene_names)) {
    stop(sprintf("[ensure-names%s] stats(%d) 與 gene_names(%d) 長度不符",
                 if (!is.null(label)) paste0("-", label) else "",
                 length(v), length(gene_names)))
  }
  names(v) <- as.character(gene_names)
  # 先不要過濾，交給下一步清理（會記錄 log）
  v
}


# 從已排序的統計量做 GSEA（優先 fgseaMultilevel，失敗時退回 fgseaSimple）
.gsea_from_ranks <- function(pathways, stats, minSize, maxSize, gsea_eps = 0, label = NULL) {
  if (is.null(pathways) || !length(pathways) || is.null(stats) || !length(stats)) return(NULL)
  orig_n <- length(pathways)
  pw_use <- ._intersect_and_filter_pathways(pathways, names(stats), minSize = minSize, maxSize = maxSize)
  message(sprintf("[gsea:%s] 原=%d → 用=%d",
                  ifelse(is.null(label), "NA", label),
                  orig_n, length(pw_use)))
  if (!length(pw_use)) return(NULL)
  
  nproc <- getOption(".fgsea_nproc", 1L)
  set.seed(1L)
  res <- tryCatch({
    suppressWarnings(fgsea::fgseaMultilevel(
      pathways = pw_use, stats = stats,
      minSize = minSize, maxSize = maxSize,
      eps = gsea_eps, nproc = nproc
    ))
  },
  error = function(e) {
    message(sprintf("[gsea_from_ranks] fgseaMultilevel failed: %s; fallback to fgseaSimple", conditionMessage(e)))
    suppressWarnings(fgsea::fgseaSimple(
      pathways = pw_use, stats = stats,
      nperm = 10000, minSize = minSize, maxSize = maxSize
    ))
  }
  )
  if (!is.null(label)) attr(res, "label") <- label
  res
}


# --- 工具：清理並排序 GSEA 的 stats（保留有限值、補名稱、去重、排序） ---
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
      log_msg("[gsea-%s] drop non-finite stats: %d", label, sum(!ok))
    }
  }
  stats <- stats[ok]
  nm <- nm[ok]
  names(stats) <- nm
  # 名稱必須存在
  if (is.null(nm)) {
    log_msg("[gsea-%s] stats has no names after filtering → skip", label); 
    return(NULL)
  }
  # 去除重複名稱（保留首次）
  dup <- duplicated(nm)
  if (any(dup)) {
    if (!is.null(label) && nzchar(label)) {
      log_msg("[gsea-%s] drop duplicated gene names in stats: %d", label, sum(dup))
    }
    stats <- stats[!dup]; nm <- nm[!dup]
  }
  names(stats) <- nm
  if (length(stats) < min_n) {
    log_msg("[gsea-%s] too few finite stats after filtering: %d < %d → skip", label, length(stats), min_n)
    return(NULL)
  }
  # 依照分數由大到小排序（preranked）
  stats[order(stats, decreasing = TRUE)]
}

# --- 工具：把 pathways 映射到 stats 的宇宙並過濾大小 ---
._intersect_and_filter_pathways <- function(pathways, universe, minSize, maxSize, label = "") {
  pw <- lapply(pathways, function(gs) unique(intersect(gs, universe)))
  lens <- vapply(pw, length, integer(1))
  keep <- which(lens >= minSize & lens <= maxSize)
  if (length(keep) == 0L) {
    log_msg("[gsea-%s] no pathways within size bounds after intersect (min=%d, max=%d)", label, minSize, maxSize)
    return(NULL)
  }
  pw[keep]
}

.align_limma_inputs <- function(M, design_df, rhs_terms, label = "limma") {
  stopifnot(!is.null(M), is.matrix(M), ncol(M) > 0, !is.null(colnames(M)))
  # 先把設計表轉 data.frame，確保 rownames 是樣本 ID
  design_df <- as.data.frame(design_df, check.names = FALSE)
  stopifnot(!is.null(rownames(design_df)))
  samp_M   <- as.character(colnames(M))
  samp_des <- as.character(rownames(design_df))
  
  # 先把設計表中的 NA 样本踢掉（predictor / 數值共變項 / 因子 NA 都視為不可用）
  na_mask <- rep(FALSE, nrow(design_df))
  rn <- rownames(design_df)
  
  # predictor 必須存在
  if (!("predictor" %in% colnames(design_df)))
    stop(sprintf("[%s] design_df 缺少 predictor 欄位", label))
  
  # predictor 非有限
  v <- design_df$predictor
  na_mask <- na_mask | is.na(v) | !is.finite(v)
  
  # 其它數值欄
  num_cols <- setdiff(colnames(design_df), c("predictor"))
  for (cn in num_cols) {
    v <- design_df[[cn]]
    if (is.numeric(v)) {
      na_mask <- na_mask | is.na(v) | !is.finite(v)
    } else {
      # 因子/字串/邏輯：只要 NA 就剔除
      na_mask <- na_mask | is.na(v)
    }
  }
  design_df <- design_df[!na_mask, , drop = FALSE]
  
  # 取交集並同序
  common <- intersect(as.character(colnames(M)), as.character(rownames(design_df)))
  if (!length(common)) return(NULL)
  common <- sort(common)
  
  M2         <- M[, common, drop = FALSE]
  design_df2 <- design_df[common, , drop = FALSE]
  
  # 建立設計矩陣
  rhs <- unique(rhs_terms)
  des2 <- stats::model.matrix(stats::reformulate(rhs), data = design_df2)
  
  # 最後再檢查一次
  if (nrow(des2) != ncol(M2)) {
    msg <- sprintf("[%s] nrow(des)=%d != ncol(M)=%d（對齊後仍不一致）",
                   label, nrow(des2), ncol(M2))
    stop(msg)
  }
  list(M = M2, des = des2, design_df = design_df2)
}

audit_spearman_pairs <- function(predictor_name, predictor_vec, mat0, out_dir, min_pairs_spearman = 10L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pv <- as.numeric(predictor_vec); names(pv) <- names(predictor_vec)
  fin_p <- is.finite(pv)
  n_pairs <- vapply(rownames(mat0), function(g){
    x <- as.numeric(mat0[g, ])
    sum(fin_p & is.finite(x))
  }, FUN.VALUE = integer(1))
  df <- data.frame(gene = rownames(mat0), n_pairs = as.integer(n_pairs),
                   below_threshold = as.integer(n_pairs < min_pairs_spearman))
  data.table::fwrite(df, file.path(out_dir, "spearman_pairs_audit.csv"))
  qs <- stats::quantile(df$n_pairs, probs = c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE, names = FALSE)
  cat(sprintf("[Spearman-audit] %s | n_genes=%d | <thr(%d)=%d | min/5%%/25%%/50%%/75%%/95%%/max=%s\n",
              predictor_name, nrow(df), min_pairs_spearman, sum(df$below_threshold, na.rm=TRUE),
              paste(qs, collapse = "/")),
      file = file.path(out_dir, "spearman_pairs_summary.txt"))
}

audit_design_alignment <- function(tag, samples, mod_interest, mod_nuisance, out_dir = NULL) {
  safe_rn <- function(X) if (is.null(X)) rep(NA_character_, length(samples)) else rownames(as.matrix(X))
  df <- data.frame(
    pos          = seq_along(samples),
    M_sample     = samples,
    interest_rn  = safe_rn(mod_interest),
    nuisance_rn  = safe_rn(mod_nuisance),
    interest_ok  = if (!is.null(mod_interest)) safe_rn(mod_interest) == samples else NA,
    nuisance_ok  = if (!is.null(mod_nuisance)) safe_rn(mod_nuisance) == samples else NA,
    stringsAsFactors = FALSE
  )
  n_i <- sum(df$interest_ok,  na.rm = TRUE)
  n_n <- sum(df$nuisance_ok,  na.rm = TRUE)
  log_msg("[align:%s] interest match = %d/%d, nuisance match = %d/%d, both = %d",
          tag, n_i, length(samples), n_n, length(samples), sum(df$interest_ok & df$nuisance_ok, na.rm = TRUE))
  if (!is.null(out_dir)) {
    fn <- file.path(out_dir, sprintf("ALIGN_%s.csv", gsub("[^A-Za-z0-9]+","_", tag)))
    utils::write.csv(df, fn, row.names = FALSE)
    log_msg("[align:%s] wrote %s", tag, fn)
  }
  # 額外：列出前 6 筆對照，方便你在 log 快速看
  head_show <- utils::head(df[, c("pos","M_sample","interest_rn","nuisance_rn","interest_ok","nuisance_ok")], 6)
  log_msg("[align:%s] head:\n%s", tag, utils::capture.output(print(head_show)) |> paste(collapse="\n"))
  invisible(df)
}


## ===== impute_and_filter：插補前固定 seed（5）=====
impute_and_filter <- function(mat, min_frac=0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop=FALSE]
  if (any(is.na(m))) { set.seed(1234); m <- imputeLCMD::impute.MinProb(m, q=0.01) }
  m
}

## ===== SVA sanity check（不依賴 keep_s）=====
audit_sva_sanity <- function(
    M,                        # 基因 x 樣本 的表達矩陣（或子矩陣）
    sv = NULL,                # SV 矩陣（樣本 x K）
    SV = NULL,                # ← 與呼叫端相容的別名（會自動轉給 sv）
    predictor = NULL,         # 連續生物學訊號（向量；與樣本對齊）
    tp53 = NULL,              # 可選：TP53 狀態（因子/字串；與樣本對齊）
    out_dir = file.path("run_info","covars_audit"),  # 報表輸出目錄
    tag = NULL,               # 顯示在檔名/日誌的標籤（如 "GPS1|BatchAdj"）
    sample_order = NULL,      # 可選：欲使用的樣本順序
    purity = NULL, sex = NULL, age = NULL            # 保留參數（目前不強制用）
){
  logf <- function(...) if (exists("log_msg", mode="function"))
    try(log_msg(...), silent = TRUE)
  
  ## 參數別名對齊：讓呼叫端傳 SV / out_dir / tag 都能吃
  if (is.null(sv) && !is.null(SV)) sv <- SV
  
  ## 基本防禦
  if (is.null(sv) || !NCOL(as.matrix(sv))) {
    logf("[SVA-sanity%s] no SV provided → skip",
         ifelse(is.null(tag), "", paste0("|", tag)))
    return(invisible(NULL))
  }
  
  ## 樣本對齊（更穩健）：先決定 nm，再必要時補上 sv 的 rownames
  M  <- as.matrix(M)
  sv <- as.matrix(sv)
  
  ## NEW: 把 "1..n" 這種無意義 rownames 視為未命名
  rn <- rownames(sv)
  if (!is.null(rn) && identical(rn, as.character(seq_len(nrow(sv))))) {
    rownames(sv) <- NULL
  }
  
  nm <- if (!is.null(sample_order)) as.character(sample_order) else colnames(M)
  
  ## 若 sv 沒 rownames，但行數剛好對應 nm，直接補上；若偵測到 t(sv) 也自動轉回
  if (is.null(rownames(sv)) || anyNA(rownames(sv)) || any(rownames(sv) == "")) {
    if (nrow(sv) == length(nm)) {
      rownames(sv) <- nm
    } else if (ncol(sv) == length(nm) && nrow(sv) < ncol(sv)) {
      sv <- t(sv); rownames(sv) <- nm
    }
  }
  
  ## 交集樣本（維持 nm 的順序）
  samp <- intersect(nm, rownames(sv))
  if (!length(samp)) {
    logf("[SVA-sanity%s] no overlapping samples between M and SV → skip",
         ifelse(is.null(tag) || !nzchar(tag), "", paste0("|", tag)))
    return(invisible(NULL))
  }
  
  M  <- M[, samp, drop = FALSE]
  sv <- sv[samp, , drop = FALSE]
  
  if (!is.null(predictor)) predictor <- predictor[samp]
  if (!is.null(tp53))      tp53      <- tp53[samp]
  
  ## 計算每個 SV 與 predictor / TP53 的關聯（Spearman；與差均值）
  out_list <- lapply(seq_len(ncol(sv)), function(j){
    v <- sv[, j]
    df1 <- data.frame(SV = colnames(sv)[j], stringsAsFactors = FALSE)
    
    # 與連續 predictor 的相關
    if (!is.null(predictor)) {
      r  <- suppressWarnings(stats::cor(v, predictor, method = "spearman",
                                        use = "pairwise.complete.obs"))
      df1$rho_predictor <- as.numeric(r)
    } else {
      df1$rho_predictor <- NA_real_
    }
    
    # 與 TP53（兩群）之均值差
    if (!is.null(tp53)) {
      tp <- factor(tp53)
      if (nlevels(tp) == 2L) {
        m1 <- mean(v[tp == levels(tp)[1]], na.rm = TRUE)
        m2 <- mean(v[tp == levels(tp)[2]], na.rm = TRUE)
        df1$delta_mean_TP53 <- m2 - m1
      } else {
        df1$delta_mean_TP53 <- NA_real_
      }
    } else {
      df1$delta_mean_TP53 <- NA_real_
    }
    
    df1
  })
  out <- do.call(rbind, out_list)
  
  ## 輸出 CSV 報表
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  safe_tag <- gsub("[^A-Za-z0-9._-]+", "_", ifelse(is.null(tag), "SVA", tag))
  fp <- file.path(out_dir, paste0("SVA_sanity_", safe_tag, ".csv"))
  data.table::fwrite(out, fp)
  
  logf("[SVA-sanity|%s] wrote %s", safe_tag, fp)
  invisible(out)
}


## ===== 覆蓋率與共變項稽核：audit_covars_coverage() =====
## 建議放在 run_one_stratum() 之前
## 會把摘要附加到 run_info/covars_audit/audit_rows.csv，並印一行 log
if (!exists(".ensure_mat_or_null", mode = "function")) {
  .ensure_mat_or_null <- function(x) {
    if (is.null(x)) return(NULL)
    m <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (is.null(dim(m))) {               # 向量 → 單欄矩陣
      m <- matrix(m, ncol = 1)
    }
    if (ncol(m) == 0) return(NULL)
    m
  }
}



rank_by_spearman_raw <- function(mat_raw, subunit, min_pairs = 10){
  stopifnot(subunit %in% rownames(mat_raw))
  sub <- mat_raw[subunit, ]
  keep <- !is.na(sub)
  if (sum(keep) < min_pairs) stop("可用樣本太少，無法做 Spearman")
  sub <- sub[keep]; m <- mat_raw[, keep, drop = FALSE]
  rho <- apply(m, 1, function(x){
    ok <- stats::complete.cases(x, sub)
    if (sum(ok) >= min_pairs) suppressWarnings(stats::cor(x[ok], sub[ok], method = "spearman")) else NA_real_
  })
  sort(rho[!is.na(rho)], decreasing = TRUE)
}
gsea_from_diffcorr <- function(pathways, z_stats, label = "diff-corr", out_prefix = NULL) {
  stopifnot(is.numeric(z_stats), !is.null(names(z_stats)))
  minSize <- opt("minSize", 15L); maxSize <- opt("maxSize", 500L)
  orig_n <- length(pathways)
  pw_use <- ._intersect_and_filter_pathways(pathways, names(z_stats), minSize = minSize, maxSize = maxSize)
  if (!is.null(out_prefix)) {
    message(sprintf("[gsea:%s] %s 原=%d → 用=%d", label, out_prefix, orig_n, length(pw_use)))
  } else {
    message(sprintf("[gsea:%s] 原=%d → 用=%d", label, orig_n, length(pw_use)))
  }
  if (!length(pw_use)) return(invisible(NULL))
  
  set.seed(1L)
  res <- tryCatch(
    {
      suppressWarnings(
        fgsea::fgseaMultilevel(
          pathways = pw_use, stats = z_stats,
          minSize = minSize, maxSize = maxSize,
          eps = 1e-10, nproc = getOption(".fgsea_nproc", 1L)
        )
      )
    },
    error = function(e) {
      message(sprintf("[gsea_from_diffcorr] Multilevel 失敗：%s → 改用 fgseaSimple", conditionMessage(e)))
      suppressWarnings(
        fgsea::fgseaSimple(
          pathways = pw_use, stats = z_stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        )
      )
    }
  )
  invisible(as.data.frame(res))
}

run_fgsea_save <- function(stats, genesets, out_prefix, top_plot_n = 0, plot_title = "") {
  set.seed(1234)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  
  minSize <- opt("minSize", 15L); maxSize <- opt("maxSize", 500L)
  
  orig_n <- length(genesets)
  pw_use <- ._intersect_and_filter_pathways(genesets, names(stats), minSize = minSize, maxSize = maxSize)
  message(sprintf("[gsea] %s 原=%d → 用=%d", out_prefix, orig_n, length(pw_use)))
  if (!length(pw_use)) {
    data.table::fwrite(data.frame(), paste0(out_prefix, ".csv")); return(invisible(NULL))
  }
  
  bp <- tryCatch({
    if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::SerialParam() else NULL
  }, error = function(e) NULL)
  
  res <- tryCatch({
    suppressWarnings(
      fgsea::fgseaMultilevel(
        pathways = pw_use,
        stats    = stats,
        minSize  = minSize,
        maxSize  = maxSize,
        eps      = 1e-10,
        BPPARAM  = bp
      )
    )
  },
  error = function(e) {
    log_msg("    [fgseaMultilevel] 失敗：%s → 改用 fgsea()", e$message)
    suppressWarnings(
      fgsea::fgsea(
        pathways = pw_use,
        stats    = stats,
        minSize  = minSize,
        maxSize  = maxSize,
        nperm    = 1000,
        BPPARAM  = bp
      )
    )
  }
  )
  res <- as.data.frame(res)
  data.table::fwrite(res, paste0(out_prefix, ".csv"))
  if (top_plot_n > 0) {
    # 若需保留原繪圖，可在此加入
  }
  invisible(res)
}


read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
  grp <- safe_fs_name(group_name)
  su  <- subunit
  f   <- paste0(stat_tag, ".csv")
  
  # 依你實際寫檔結構（RAW/BatchAdj）與舊假設，列候選路徑
  candidates <- c(
    file.path(out_root, su, f),            # [NEW] collection-first: ./<subunit>/STAT.csv
    file.path(out_root, grp, su, "H", f),  # 舊：./<group>/<subunit>/H/STAT.csv
    file.path(out_root, grp, su, f),       # 舊：無 H 子夾
    file.path(out_root, su, grp, "H", f),  # 舊：./<subunit>/<group>/H/STAT.csv
    file.path(out_root, su, grp, f)        # 舊：無 H 子夾
  )
  
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    if (exists("log_msg", mode="function")) try(log_msg("    檔案不存在：%s", paste(candidates, collapse=" | ")), silent=TRUE)
    return(NULL)
  }
  fp <- hit[1L]
  
  dt <- tryCatch(data.table::fread(fp, na.strings = c("NA","NaN","")), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  
  need_cols <- c("pathway","NES","padj")
  if (!all(need_cols %in% names(dt))) {
    # 寬鬆對應欄名
    nms <- tolower(names(dt))
    if (!"pathway" %in% names(dt) && "pathway" %in% nms) names(dt)[match("pathway", nms)] <- "pathway"
    if (!"NES" %in% names(dt) && "nes" %in% nms)       names(dt)[match("nes", nms)]       <- "NES"
    if (!"padj" %in% names(dt)) dt$padj <- NA_real_
  }
  dt <- as.data.frame(dt[, c("pathway","NES","padj")])
  names(dt)[names(dt) == "NES"]  <- paste0("NES_", subunit)
  names(dt)[names(dt) == "padj"] <- paste0("padj_", subunit)
  dt
}

merge_subunit_tables <- function(tbl_list) {
  keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
  if (!length(keep)) return(NULL)
  out <- Reduce(function(x,y) dplyr::full_join(x,y,by="pathway"), keep)
  num_cols <- setdiff(names(out), "pathway")
  out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
  out
}

add_sig_counts <- function(df, alphas = c(0.05, 0.25)) {
  if (is.null(df) || !nrow(df)) return(df)
  padj_cols <- grep("^padj_", names(df), value = TRUE)
  subunits  <- sub("^padj_", "", padj_cols)
  nes_cols  <- paste0("NES_", subunits)
  keep_idx  <- nes_cols %in% names(df)
  padj_cols <- padj_cols[keep_idx]; nes_cols <- nes_cols[keep_idx]
  if (!length(padj_cols)) return(df)
  for (a in alphas) {
    a_tag <- gsub("\\.", "_", as.character(a))
    sig_mat <- mapply(function(p, n){ pv <- df[[p]]; as.integer(!is.na(pv) & pv < a) }, padj_cols, nes_cols, SIMPLIFY = TRUE)
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    df[[sprintf("sig_n_padj_%s", a_tag)]] <- rowSums(sig_mat, na.rm=TRUE)
    pos_mat <- mapply(function(p, n){ pv <- df[[p]]; nv <- df[[n]]; as.integer(!is.na(pv) & pv < a & !is.na(nv) & nv > 0) }, padj_cols, nes_cols, SIMPLIFY = TRUE)
    if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
    df[[sprintf("pos_n_padj_%s", a_tag)]] <- rowSums(pos_mat, na.rm=TRUE)
    neg_mat <- mapply(function(p, n){ pv <- df[[p]]; nv <- df[[n]]; as.integer(!is.na(pv) & pv < a & !is.na(nv) & nv < 0) }, padj_cols, nes_cols, SIMPLIFY = TRUE)
    if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
    df[[sprintf("neg_n_padj_%s", a_tag)]] <- rowSums(neg_mat, na.rm=TRUE)
  }
  df
}

write_summary_outputs <- function(df, out_dir, group_name, stat_tag) {
  if (is.null(df) || !nrow(df)) { log_msg("  [略過輸出] {group_name} | {stat_tag} 無可用結果"); return(invisible(NULL)) }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))
  ord_keys <- c("sig_n_padj_0_05","pos_n_padj_0_05","neg_n_padj_0_05"); ord_keys <- intersect(ord_keys, names(df))
  if (length(ord_keys)) {
    df <- df |>
      dplyr::arrange(dplyr::desc(.data[[ord_keys[1]]]),
                     dplyr::desc(ifelse(length(ord_keys) > 1, .data[[ord_keys[2]]], 0)),
                     dplyr::desc(ifelse(length(ord_keys) > 2, .data[[ord_keys[3]]], 0)),
                     .data[["pathway"]])
  } else df <- df |> dplyr::arrange(.data[["pathway"]])
  data.table::fwrite(df, paste0(base, "_ALL.csv"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "ALL"); openxlsx::writeData(wb, "ALL", df)
  if ("sig_n_padj_0_05" %in% names(df)) {
    df005 <- df |> dplyr::filter(.data[["sig_n_padj_0_05"]] > 0)
    openxlsx::addWorksheet(wb, "padj_lt_0.05"); openxlsx::writeData(wb, "padj_lt_0.05", df005)
    data.table::fwrite(df005, paste0(base, "_padjLT0.05.csv"))
  }
  if ("sig_n_padj_0_25" %in% names(df)) {
    df025 <- df |> dplyr::filter(.data[["sig_n_padj_0_25"]] > 0)
    openxlsx::addWorksheet(wb, "padj_lt_0.25"); openxlsx::writeData(wb, "padj_lt_0.25", df025)
    data.table::fwrite(df025, paste0(base, "_padjLT0.25.csv"))
  }
  openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite=TRUE)
  log_msg("  [完成輸出] {group_name} | {stat_tag} -> {dirname(base)}")
}

summarize_all_groups <- function(out_root, csn_subunits, genesets_by_group,
                                 stat_tags = c("GSEA_limma_t","GSEA_spearman")) {
  sum_root <- file.path(out_root, "summary")
  dir.create(sum_root, recursive = TRUE, showWarnings = FALSE)
  for (grp_name in names(genesets_by_group)) {
    grp_safe <- safe_fs_name(grp_name)
    for (stat_tag in stat_tags) {
      log_msg("== 彙整：group={grp_name} | stat={stat_tag} ==")
      tbl_list <- setNames(vector("list", length(csn_subunits)), csn_subunits)
      for (su in csn_subunits) tbl_list[[su]] <- read_gsea_table(out_root, su, grp_name, stat_tag)
      wide <- merge_subunit_tables(tbl_list)
      wide <- add_sig_counts(wide, alphas = c(0.05, 0.25))
      out_dir <- file.path(sum_root, stat_tag)
      write_summary_outputs(wide, out_dir, grp_name, stat_tag)
    }
  }
}

## ===== Meta-FDR via Stouffer (per stratum × predictor × group × stat) =====
safe_read_gsea <- function(fp){
  if (!file.exists(fp)) return(NULL)
  dt <- tryCatch(data.table::fread(fp, na.strings=c("NA","NaN","")), error=function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  need <- c("pathway","NES","pval")
  if (!all(need %in% names(dt))) return(NULL)
  as.data.frame(dt[, need])
}

.z_from_p_signed <- function(p, sign){
  p <- pmax(pmin(p, 1-1e-16), 1e-16)
  z <- stats::qnorm(1 - p/2)
  z * sign
}

meta_fdr_stouffer <- function(dataset_dirs,
                              strata = c("ALL","TP53_mutant","TP53_wild_type"),
                              stat_tags = c("GSEA_limma_t_cont","GSEA_spearman"),
                              groups = names(genesets_by_group),
                              out_root     = if (is.null(COMBO_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")) {
  stopifnot(length(dataset_dirs) > 0)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("請先安裝 data.table")
  if (!requireNamespace("stats", quietly = TRUE)) stop("缺 stats 套件")
  data.table::setDTthreads(1L)
  
  # ---- helpers --------------------------------------------------------------
  normalize_cols <- function(DT) {
    # 將常見欄名映射為通用名稱：pathway, NES, pval, padj, size
    nm <- names(DT)
    # 先全轉小寫做比對
    low <- tolower(nm)
    map <- c(
      "pathway"      = "pathway",
      "term"         = "pathway",
      "pathway_name" = "pathway",
      "nes"          = "NES",
      "enrichment_score" = "NES",
      "pval"         = "pval",
      "p.value"      = "pval",
      "pvalue"       = "pval",
      "p"            = "pval",
      "padj"         = "padj",
      "fdr"          = "padj",
      "qval"         = "padj",
      "size"         = "size",
      "setsize"      = "size",
      "n"            = "size"
    )
    # 套用映射：若 low 名稱存在 map，就改名成 map 目標
    for (i in seq_along(nm)) {
      key <- low[i]
      if (key %in% names(map)) data.table::setnames(DT, nm[i], map[[key]])
    }
    DT
  }
  
  pick_cols <- function(DT, need) {
    need <- unique(need)
    if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
    # 補齊缺欄，避免 rbind / 選欄失敗
    miss <- setdiff(need, names(DT))
    for (m in miss) DT[, (m) := NA]
    # data.table NSE 正確選欄法：..need 或 .SDcols
    DT[, ..need]
  }
  
  stouffer_df <- function(D) {
    # D 必須至少有 pathway, NES, pval
    D <- D[is.finite(pval) & is.finite(NES)]
    if (nrow(D) == 0) return(NULL)
    # 以「兩尾 p 值 + NES 方向」構成有號 Z：z = sign(NES) * qnorm(1 - p/2)
    D[, z := sign(NES) * stats::qnorm(p = pval/2, lower.tail = FALSE)]
    # 無權重 Stouffer（可自行改成加權）
    Out <- D[, .(k = .N, Z = sum(z) / sqrt(.N)), by = .(pathway)]
    Out[, p_meta := 2 * stats::pnorm(abs(Z), lower.tail = FALSE)]
    Out[, padj_meta := p.adjust(p_meta, method = "BH")]
    data.table::setorder(Out, p_meta)
    Out[]
  }
  
  # ---- main loop ------------------------------------------------------------
  passes <- c("RAW","BatchAdj")
  base_out <- out_root
  dir.create(base_out, recursive = TRUE, showWarnings = FALSE)
  
  for (stratum in strata) {
    for (pass_label in passes) {
      for (stat_tag in stat_tags) {
        for (grp in groups) {
          
          # 1) 找出所有 subunit（以各 dataset 的資料夾 union）
          subunits <- unique(unlist(lapply(dataset_dirs, function(dsdir) {
            base <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection",
                              safe_fs_name(grp), pass_label, basename(dsdir), stratum)
            if (!dir.exists(base)) return(character(0))
            # 只列最上層子目錄名（每一個 subunit）
            subs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
            subs[subs != ""]
          })))
          
          if (length(subunits) == 0) {
            message(sprintf("[meta] %s/%s/%s/%s 沒有可用 subunit，略過",
                            stratum, pass_label, stat_tag, grp))
            next
          }
          
          # 2) 對每個 subunit，收集各 dataset 的 CSV、做 Stouffer 合併並落檔
          for (su in subunits) {
            parts <- list()
            for (dsdir in dataset_dirs) {
              f <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection",
                             safe_fs_name(grp), pass_label, basename(dsdir), stratum, su,
                             sprintf("%s.csv", stat_tag))
              if (!file.exists(f)) next
              dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
              if (is.null(dt) || nrow(dt) == 0) next
              dt <- normalize_cols(dt)
              dt <- pick_cols(dt, need = c("pathway","NES","pval","padj","size"))
              dt[, dataset := basename(dsdir)]
              parts[[length(parts) + 1L]] <- dt
            }
            
            if (length(parts) == 0) {
              # 沒有任何隊列有這個 subunit 的輸出 → 不寫檔也不報錯
              next
            }
            
            D <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
            # 去掉 pathway 缺失
            D <- D[!is.na(pathway)]
            res <- stouffer_df(D)
            if (is.null(res)) next
            
            # 3) 寫檔
            out_dir <- file.path(base_out, stratum, pass_label, su, grp)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_file <- file.path(out_dir, sprintf("%s_meta_fdr.csv", stat_tag))
            data.table::fwrite(res, out_file)
            message(sprintf("[meta] %s/%s/%s/%s -> %s (n=%d pathways, k>=1 datasets)",
                            stratum, pass_label, su, stat_tag, out_file, nrow(res)))
          } # su
        } # grp
      } # stat_tag
    } # pass
  } # stratum
  
  invisible(TRUE)
}

## ==== PATCH: helpers ====

## === NEW: CSN complex score (PC1) & residualization helpers ==================

# 以 subunits 的「跨樣本 z 值」做 PCA，取 PC1 當 CSN score；方向校正成與平均 z 同號
build_csn_score <- function(mat0,
                            subunits = csn_subunits,
                            combine_7AB = TRUE,
                            min_members = 5L) {
  present <- intersect(subunits, rownames(mat0))
  # 預先建立回傳骨架（一定有名字）
  s <- setNames(rep(NA_real_, ncol(mat0)), colnames(mat0))
  if (!length(present)) return(s)
  
  # z-score（保留樣本名稱）
  get_z <- function(v) {
    nm <- names(v)                  # 先存下樣本名
    v  <- as.numeric(v)
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    out <- (v - mu) / sdv
    names(out) <- nm                # 放回樣本名
    out
  }
  
  X <- do.call(rbind, lapply(present, function(g) get_z(mat0[g, ])))
  rownames(X)  <- present
  colnames(X)  <- colnames(mat0)    # 這行很重要：確保有樣本名
  
  # 可選：合併 COPS7A/7B
  if (combine_7AB && all(c("COPS7A","COPS7B") %in% rownames(X))) {
    Z7 <- colMeans(X[c("COPS7A","COPS7B"), , drop = FALSE], na.rm = TRUE)
    X  <- rbind(X[setdiff(rownames(X), c("COPS7A","COPS7B")), , drop = FALSE],
                "COPS7*" = Z7)
  }
  
  # 用「原始 mat0」的非 NA 數來判斷覆蓋
  enough <- colSums(is.finite(mat0[present, , drop = FALSE])) >= min_members
  keep_sam <- names(s)[enough]
  
  if (length(keep_sam) >= 10) {
    pc <- stats::prcomp(t(X[, keep_sam, drop = FALSE]), center = TRUE, scale. = FALSE)
    sc <- pc$x[, 1]
    # 方向校正：與 subunit 平均 z 同號
    mu <- colMeans(X[, keep_sam, drop = FALSE], na.rm = TRUE)
    if (suppressWarnings(cor(sc, mu, use = "pairwise.complete.obs")) < 0) sc <- -sc
    s[keep_sam] <- sc               # 這裡一定能以樣本名對齊成功
  }
  
  s
}

# 將向量 y（某 subunit）對 CSN score + batch + 其他共變項做線回歸後取殘差
# 取代原 residualize_vector()
residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8L) {
  # 1) 對齊樣本名稱
  if (is.null(names(y))) stop("[residualize_vector] y 必須是具樣本名稱的數值向量")
  sam <- names(y)
  common <- sam
  if (!is.null(csn_score)) common <- intersect(common, names(csn_score))
  if (!is.null(batch))     common <- intersect(common, names(batch))
  if (!is.null(covars)) {
    rn <- rownames(as.data.frame(covars, check.names = FALSE))
    if (is.null(rn)) stop("[residualize_vector] covars 必須 rownames=樣本")
    common <- intersect(common, rn)
  }
  if (length(common) < min_n) return(setNames(rep(NA_real_, length(y)), names(y)))
  
  # 2) 建立資料框 → 單一 model.matrix 展開
  DF <- data.frame(
    csn = as.numeric(csn_score[common]),
    row.names = common, check.names = FALSE
  )
  if (!is.null(batch)) {
    bf <- droplevels(batch[common])
    DF[["batch"]] <- bf
  }
  if (!is.null(covars)) {
    C <- as.data.frame(covars[common, , drop = FALSE], check.names = FALSE)
    C <- coerce_covariates_safely(C)                # ← 新增
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]]     # ← 直接放回 DF
  }
  DF_y <- suppressWarnings(as.numeric(y[common]))
  
  ## === NEW: 先把會讓 complete.cases 全滅的欄位清掉 ===
  # 2a) 丟全 NA 欄
  all_na_col <- vapply(DF, function(v) all(is.na(v)), logical(1))
  if (any(all_na_col)) {
    if (exists("log_msg", mode="function")) try(log_msg("  [resid] drop all-NA columns: %s",
                                                        paste(names(all_na_col)[all_na_col], collapse=",")), silent=TRUE)
    DF <- DF[, !all_na_col, drop = FALSE]
  }
  
  # 2b) 丟單一水準的因子或常數數值欄（避免設計矩陣奇異）
  is_const <- vapply(DF, function(v){
    vv <- v[!is.na(v)]
    if (!length(vv)) TRUE else {
      if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
    }
  }, logical(1))
  if (any(is_const)) {
    if (exists("log_msg", mode="function")) try(log_msg("  [resid] drop constant/single-level columns: %s",
                                                        paste(names(is_const)[is_const], collapse=",")), silent=TRUE)
    DF <- DF[, !is_const, drop = FALSE]
  }
  ## === /NEW ===
  
  ok <- is.finite(DF_y) & stats::complete.cases(DF)
  if (sum(ok) < min_n) {
    out <- setNames(rep(NA_real_, length(y)), names(y))
    return(out)
  }
  DF   <- DF[ok, , drop = FALSE]
  y_ok <- DF_y[ok]
  
  # 3) 展開設計矩陣
  # - batch 若存在且為因子，自動 dummy encode
  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == length(y_ok))
  
  # 4) 線性回歸 y ~ csn + batch + covars
  fit <- lm.fit(x = des, y = y_ok)
  res <- rep(NA_real_, length(DF_y))
  res[ok] <- fit$residuals
  
  # 5) 對應回完整樣本並（可選）z-score
  out <- setNames(rep(NA_real_, length(y)), names(y))
  out[common] <- res
  fin <- is.finite(out)
  if (sum(fin) >= 3L) {
    mu <- mean(out[fin]); sdv <- stats::sd(out[fin]); if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out[fin] <- (out[fin] - mu)/sdv
  }
  out
}

audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL,
                                  sv = NULL,
                                  tech = NULL) {
  dir.create(file.path("run_info","covars_audit"), recursive = TRUE, showWarnings = FALSE)
  
  cov_df  <- if (is.null(covars)) NULL else as.data.frame(covars, check.names = FALSE)
  cov_nms <- if (!is.null(cov_df)) colnames(cov_df) else character(0)
  
  get_cov <- function(k){
    if (is.null(cov_df) || !(k %in% names(cov_df))) return(NA_real_)
    v <- cov_df[[k]]
    if (is.numeric(v)) {
      mean(is.finite(v)) * 100
    } else {
      mean(!is.na(v)) * 100   # 因子/字串 → 以非 NA 比例計
    }
  }
  cov_purity <- get_cov("purity"); cov_sex <- get_cov("sex"); cov_age <- get_cov("age")
  
  sv_nms   <- if (!is.null(sv))   colnames(.ensure_mat_or_null(sv))   else NULL
  tech_nms <- if (!is.null(tech)) colnames(.ensure_mat_or_null(tech)) else NULL
  
  line <- sprintf("  [audit:%s] covars={%s} | coverage: purity=%.1f%%, sex=%.1f%%, age=%.1f%% | batch_levels=%s | SV={%s} | tech={%s}",
                  tag,
                  if (length(cov_nms)) paste(cov_nms, collapse=",") else "NULL",
                  cov_purity %||% NaN, cov_sex %||% NaN, cov_age %||% NaN,
                  if (!is.null(batch)) nlevels(batch) else 0L,
                  if (length(sv_nms)) paste(sv_nms, collapse=",") else "none",
                  if (length(tech_nms)) paste(tech_nms, collapse=",") else "none")
  log_msg(line)
  
  # 追加一行到總表
  df <- data.frame(
    dataset = ds_id, stratum = stratum, subunit = su, tag = tag,
    batch_levels = if (!is.null(batch)) nlevels(batch) else 0L,
    batch_sizes = if (!is.null(batch)) paste(sprintf("%s=%d", names(sort(table(batch), decreasing=TRUE)),
                                                     as.integer(sort(table(batch), decreasing=TRUE))), collapse="; ") else NA_character_,
    covars_cols = if (length(cov_nms)) paste(cov_nms, collapse=";") else "NULL",
    sv_cols     = if (length(sv_nms)) paste(sv_nms, collapse=";") else "NULL",
    tech_cols   = if (length(tech_nms)) paste(tech_nms, collapse=";") else "NULL",
    purity_cov  = cov_purity, sex_cov = cov_sex, age_cov = cov_age,
    stringsAsFactors = FALSE
  )
  fp <- file.path("run_info","covars_audit","audit_rows.csv")
  data.table::fwrite(df, fp, append = file.exists(fp))
}

orthogonalize_to <- function(mat, nuisance){
  if (is.null(mat) || is.null(nuisance)) return(mat)
  Y <- as.matrix(mat)
  X <- as.matrix(nuisance)
  
  # 1) 加截距
  X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
  
  # 2) 非有限值處理 + 移除零變異欄
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (any(!is.finite(v))) {
      mu <- mean(v[is.finite(v)], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[!is.finite(v)] <- mu
      X[, j] <- v
    }
  }
  var_ok <- apply(X, 2, function(z) { z <- as.numeric(z); is.finite(var(z)) && var(z) > 0 })
  if (!all(var_ok)) X <- X[, var_ok, drop = FALSE]
  
  beta <- tryCatch(qr.coef(qr(X), Y), error = function(e) NULL)
  if (is.null(beta)) return(as.matrix(mat))
  
  resid <- Y - X %*% beta
  rownames(resid) <- rownames(Y)
  colnames(resid) <- colnames(Y)
  return(resid)
}

gate_tech_pcs <- function(cv, v, thr = 0.3, alpha = 0.01){
  # 只 gate 技術 PC（PC1/PC2…），median/missrate 不動（但我們下游只會取 PC1/PC2 當 tech covars）
  if (is.null(cv) || ncol(cv) == 0) return(cv)
  pc_cols <- grep("^PC[0-9]+$", colnames(cv), value = TRUE)
  if (!length(pc_cols)) return(cv)
  keep <- sapply(pc_cols, function(cn){
    r <- suppressWarnings(stats::cor(cv[, cn], v, use = "pairwise.complete.obs"))
    if (!is.finite(r)) return(TRUE)
    n <- sum(stats::complete.cases(cv[, cn], v))
    if (n < 15) return(TRUE)
    t <- r * sqrt((n - 2)/(1 - r^2)); p <- 2 * stats::pt(-abs(t), df = n - 2)
    !(abs(r) >= thr && p < alpha)  # 太像生物訊號 → 不放
  })
  cbind(cv[, setdiff(colnames(cv), pc_cols), drop = FALSE],
        cv[, pc_cols[keep], drop = FALSE])
}

diagnose_MNAR <- function(mat_raw){
  gene_miss <- rowMeans(!is.finite(mat_raw))
  gene_med  <- apply(mat_raw, 1, function(x) stats::median(x[is.finite(x)], na.rm=TRUE))
  rho <- suppressWarnings(stats::cor(gene_miss, gene_med, method="spearman", use="pairwise.complete.obs"))
  
  samp_med <- apply(mat_raw, 2, function(x) stats::median(x[is.finite(x)], na.rm=TRUE))
  idx <- which(is.finite(gene_med))
  mat_sub <- mat_raw[idx,, drop=FALSE]
  miss_vec <- as.integer(is.na(as.vector(mat_sub)))
  df <- data.frame(
    miss = miss_vec,
    gmed = rep(gene_med[idx], times=ncol(mat_sub)),
    smed = rep(samp_med,     each = nrow(mat_sub))
  )
  fit  <- suppressWarnings(stats::glm(miss ~ gmed + smed, data=df, family=stats::binomial()))
  summ <- suppressWarnings(summary(fit)$coefficients)
  list(
    spearman_gene_missing_vs_abundance = rho,
    glm_coef = summ[c("gmed","smed"), , drop=FALSE]
  )
}

## ---- 自動偵測 batch 欄位 & 讀取（擴充候選欄位）----
## 小工具（只在此檔內用）
.take_first <- function(x) sub("\\|.*$", "", as.character(x))
.sanitize_levels <- function(f) {
  if (is.null(f)) return(NULL)
  f <- droplevels(factor(f))
  lv <- levels(f)
  lv2 <- make.names(lv)
  lv2 <- paste0("b_", lv2)     # 避免數字開頭
  levels(f) <- lv2
  f
}
.align_by_colnames <- function(vec_or_df, target_names) {
  # 若提供的是 named vector/factor/data.frame，依 colnames(mat) 對齊；否則當作已對齊
  if (is.null(vec_or_df)) return(NULL)
  if (is.vector(vec_or_df) || is.factor(vec_or_df)) {
    nm <- names(vec_or_df)
    if (!is.null(nm) && length(nm) == length(vec_or_df)) {
      return(vec_or_df[match(target_names, nm)])
    } else return(vec_or_df)
  } else {
    rn <- rownames(vec_or_df)
    if (!is.null(rn)) return(vec_or_df[match(target_names, rn), , drop=FALSE])
    return(vec_or_df)
  }
}
.mk_batch_factor <- function(batch, sample_order, min_count_collapse = 1L) {
  if (is.null(batch)) return(NULL)
  b <- .align_by_colnames(batch, sample_order)
  b <- .take_first(b)
  b[!nzchar(b) | is.na(b)] <- "unknown"
  b <- .sanitize_levels(b)
  # 可選：把極小批次合併（預設單例就合併，設 0 代表不合併）
  if (min_count_collapse > 0L) {
    tb <- table(b)
    small <- names(tb)[tb <= min_count_collapse]
    if (length(small)) {
      b[b %in% small] <- factor("b_small")
      b <- droplevels(factor(b))
    }
  }
  b
}

## ===== 批次清理設定（你可以自行調整）=====
BATCH_PIPE_POLICY    <- "NA"   # "NA" or "b_small"：含 '|' 的值要設為 NA 還是併到 b_small
BATCH_MIN_PER_LEVEL  <- 2      # 小於這個門檻的 level 會併到 b_small（避免不可估係數）

## ---- 批次值清理：處理含 '|', 合法化名稱、合併稀疏 level ----
sanitize_batch_levels <- function(x,
                                  pipe_policy   = BATCH_PIPE_POLICY,
                                  min_per_level = BATCH_MIN_PER_LEVEL){
  x0 <- as.character(x)
  has_pipe <- grepl("\\|", x0 %||% "")
  if (any(has_pipe)) {
    log_msg("  [batch] 偵測到含 '|' 的值 {sum(has_pipe)} 個 → 依政策 {pipe_policy} 處理")
    x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
  }
  # 合法化名稱（避免 eBayes/設計矩陣列名問題）
  fac <- factor(make.names(x0))
  fac <- droplevels(fac)
  
  # 合併稀疏 level（例如只有 1 個樣本）
  if (!is.null(min_per_level) && min_per_level > 1) {
    tab <- table(fac, useNA = "no")
    small <- names(tab)[tab < min_per_level]
    if (length(small)) {
      log_msg("  [batch] 合併稀疏 level 至 'b_small'：%s",
              paste(sprintf("%s(n=%d)", small, as.integer(tab[small])), collapse=", "))
      fac_chr <- as.character(fac)
      fac_chr[fac_chr %in% small] <- "b_small"
      fac <- droplevels(factor(fac_chr))
    }
  }
  fac
}

## ---- 自動偵測 batch 欄位（加入清理管線）----
detect_batch_column <- function(meta,
                                pipe_policy   = BATCH_PIPE_POLICY,
                                min_per_level = BATCH_MIN_PER_LEVEL){
  cand <- c("TMT_PLEX", "EXPERIMENT", "PROTEOMICS_TMT_BATCH")
  hit  <- intersect(cand, colnames(meta))
  if (!length(hit)) return(NULL)
  
  for (cn in hit) {
    fac <- sanitize_batch_levels(meta[[cn]],
                                 pipe_policy   = pipe_policy,
                                 min_per_level = min_per_level)
    # 至少要有 2 個有效 level 且有效樣本數要 >= 3
    if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
      return(list(name = cn, fac = fac))
    }
  }
  NULL
}

## ---- 依樣本順序取出 batch factor（已清理）----
get_batch_factor <- function(ds_dir, sample_ids,
                             pipe_policy   = BATCH_PIPE_POLICY,
                             min_per_level = BATCH_MIN_PER_LEVEL){
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  if (!file.exists(meta_fp)) return(NULL)
  
  meta <- suppressMessages(
    readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
  ) |> as.data.frame()
  
  id_cols <- intersect(c("SAMPLE_ID","sample_id","Sample_ID","Sample","sample"), names(meta))
  if (!length(id_cols)) return(NULL)
  
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)
  
  # 先按照 sample_ids 對齊，確保偵測與回傳的向量長度一致
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  rownames(meta) <- sample_ids
  
  det <- detect_batch_column(meta,
                             pipe_policy   = pipe_policy,
                             min_per_level = min_per_level)
  if (is.null(det)) return(NULL)
  
  fac <- det$fac   # 已與 sample_ids 對齊
  names(fac) <- sample_ids
  list(name = det$name, fac = fac)
}

## ---- 批次需求快篩（沿用你的邏輯）----
screen_batch_need <- function(ds_dir, min_frac_complete = 0.75){
  log_msg("== Batch 檢查：%s ==", basename(ds_dir))
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  mx <- suppressWarnings(max(mat0, na.rm=TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  m <- impute_and_filter(mat0, min_frac = min_frac_complete)
  
  bi <- get_batch_factor(ds_dir, colnames(m))
  if (is.null(bi)) { log_msg("  [batch] 找不到明確批次欄位 → 先不校正"); 
    return(invisible(list(dataset = basename(ds_dir),
                          batch_col = NA, n_levels = 0, sizes = NA,
                          pc_R2 = NA, pc_p = NA, frac_FDR05 = NA, frac_FDR25 = NA,
                          recommend = FALSE))) }
  
  batch <- droplevels(bi$fac[colnames(m)])
  tab   <- sort(table(batch), decreasing=TRUE)
  log_msg("  [batch] 欄位：%s；levels=%d；每層 n：%s",
          bi$name, nlevels(batch),
          paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse=", "))
  
  X  <- t(scale(t(m), center = TRUE, scale = TRUE))
  pc <- prcomp(t(X), scale. = FALSE)
  K  <- min(5, ncol(pc$x))
  r2 <- vapply(seq_len(K), function(i) summary(lm(pc$x[,i] ~ batch))$r.squared, numeric(1))
  pv <- vapply(seq_len(K), function(i) { a <- anova(lm(pc$x[,i] ~ batch)); as.numeric(a$`Pr(>F)`[1]) }, numeric(1))
  log_msg("  PCA (by batch) R²：%s ; p：%s", paste(round(r2,3), collapse=", "),
          paste(signif(pv,3), collapse=", "))
  
  design <- model.matrix(~ batch)
  fit <- limma::lmFit(m, design); fit <- limma::eBayes(fit)
  padj <- p.adjust(fit$F.p.value, "BH")
  prop05 <- mean(padj < 0.05, na.rm=TRUE)
  prop25 <- mean(padj < 0.25, na.rm=TRUE)
  log_msg("  基因層級 F 檢定：FDR<0.05 比例 = %.1f%%；FDR<0.25 = %.1f%%",
          100*prop05, 100*prop25)
  
  recommend <- (any(r2 >= 0.10 & pv[seq_along(r2)] < 0.01) || prop05 >= 0.05)
  if (recommend) log_msg("  **建議校正**：R²或基因比例達門檻（≥10%% R² 且 p<0.01，或 FDR<0.05 基因 ≥5%%）")
  else log_msg("  **暫可不校正**：未見明顯 batch 影響（保留紀錄）")
  
  invisible(list(dataset = basename(ds_dir),
                 batch_col = bi$name, n_levels = nlevels(batch),
                 sizes = tab, pc_R2 = r2, pc_p = pv,
                 frac_FDR05 = prop05, frac_FDR25 = prop25,
                 recommend = recommend))
}

datasets_root <- getwd()

dataset_ids <- c(
  "brca_cptac_2020",
  "brca_tcga_pan_can_atlas_2018",
  "coad_cptac_2019",
  "coadread_tcga_pan_can_atlas_2018",
  "gbm_cptac_2021",
  "luad_cptac_2020",
  "lusc_cptac_2021",
  "ov_tcga_pan_can_atlas_2018",
  "paad_cptac_2021",
  "ucec_cptac_2020"
)


COMBO_MODE <- "combo_1"
##COMBO_MODE <- "combo_2"
##COMBO_MODE <- "combo_3"


## ---- [NEW | COMBO PRESETS: choose one] ----
## 選擇其中一個： "combo_1" / "combo_2" / "combo_3"；若不設定 (NULL) 則維持原行為
COMBO_MODE   <- get0("COMBO_MODE", ifnotfound = NULL)
COMBO_PREFIX <- if (!is.null(COMBO_MODE)) COMBO_MODE else NULL

.combo_ds_1 <- c("brca_cptac_2020","luad_cptac_2020","lusc_cptac_2021","ucec_cptac_2020")
.combo_ds_7 <- c("brca_cptac_2020","luad_cptac_2020","lusc_cptac_2021","ucec_cptac_2020",
                 "coad_cptac_2019","gbm_cptac_2021","paad_cptac_2021")
.combo2_raw_ds <- c("coad_cptac_2019","gbm_cptac_2021","paad_cptac_2021")

if (!is.null(COMBO_MODE)) {
  if (COMBO_MODE == "combo_1") {
    dataset_ids <- .combo_ds_1
  } else if (COMBO_MODE == "combo_2") {
    dataset_ids <- .combo_ds_7
  } else if (COMBO_MODE == "combo_3") {
    dataset_ids <- .combo_ds_7
  } else {
    stop(sprintf("未知的 COMBO_MODE: %s", COMBO_MODE))
  }
  ## 強制本輪只跑 limma_t，關閉 camera/mroast（與你要求一致）
  PIPELINES_TO_RUN <<- c("limma_t"); .RUN_LIMMA <<- TRUE; .RUN_SPEARMAN <<- FALSE
  RUN_CAMERA_MROAST <<- FALSE; options(csn.run_camera_mroast = FALSE)
}
## ---- [END NEW | COMBO PRESETS] ----



dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)
strata <- c("ALL","TP53_mutant","TP53_wild_type")
message("datasets_root = ", datasets_root)

## 舊版是 stopifnot(all(dir.exists(dataset_dirs)))，會在某些目錄缺失時整批中斷。
## 改為：列出缺失者、以 dataset_dirs_run 篩出「本輪實際可跑」的集合。
missing_dirs <- names(dataset_dirs)[!dir.exists(dataset_dirs)]
if (length(missing_dirs)) {
  log_msg("偵測到 %d 個缺失資料夾，將略過：%s", length(missing_dirs), paste(missing_dirs, collapse=", "))
}


## 本輪實際要跑與可用的資料集清單（存在資料夾且含 data_protein_quantification.txt）
dataset_dirs_run <- dataset_dirs[
  dir.exists(dataset_dirs) &
    file.exists(file.path(dataset_dirs, "data_protein_quantification.txt"))
]

if (!length(dataset_dirs_run)) stop("dataset_dirs_run 為空，請確認資料夾與 data_protein_quantification.txt 是否存在")

log_msg("本輪可用資料集（%d）：%s", length(dataset_dirs_run), paste(names(dataset_dirs_run), collapse=", "))

## ==== PATCH: SVA (full-rank + nested mods + robust alignment v3; non-count safe) ====
estimate_svs <- function(
    M,                      # genes x samples (numeric matrix)
    mod_interest,           # model matrix (含 predictor 與 covariates)
    mod_nuisance = NULL,    # nuisance-only design（不含截距；本函式會自己加）
    max_k = 5L,
    label = NULL
){
  tag <- ifelse(is.null(label) || !nzchar(label), "", paste0(label))
  
  # ---- 型別/對齊 ----
  M <- as.matrix(M); storage.mode(M) <- "double"
  mod_interest <- as.matrix(mod_interest)
  if (!is.null(mod_nuisance)) mod_nuisance <- as.matrix(mod_nuisance)
  
  # 共同樣本順序：以 M 的樣本名為準
  sam <- intersect(colnames(M), rownames(mod_interest))
  if (!length(sam)) { log_msg("[SVA%s] no overlapping samples -> SV=NULL", tag); return(NULL) }
  sam <- sam[order(match(sam, colnames(M)))]
  M <- M[, sam, drop=FALSE]
  mod_interest <- mod_interest[sam, , drop=FALSE]
  if (!is.null(mod_nuisance)) mod_nuisance <- mod_nuisance[sam, , drop=FALSE]
  
  # ---- 建 mod0 / 合併 full mod，補齊秩 ----
  fix_design <- function(X){
    X <- as.matrix(X); rownames(X) <- sam
    X[!is.finite(X)] <- 0
    # 去掉常數/全 NA 欄
    keep <- vapply(seq_len(ncol(X)), function(j){
      v <- X[, j]; v <- v[is.finite(v)]
      if (!length(v)) return(FALSE)
      stats::var(as.numeric(v)) > 0
    }, logical(1))
    if (any(!keep)) X <- X[, keep, drop=FALSE]
    X
  }
  mod0 <- if (!is.null(mod_nuisance) && ncol(mod_nuisance)>0) fix_design(mod_nuisance) else NULL
  mod0 <- cbind("(Intercept)" = 1, mod0)  # 本函式內自加截距
  mod  <- {
    add_full_rank <- function(B, Xadd){
      B <- as.matrix(B)
      if (is.null(Xadd) || ncol(as.matrix(Xadd))==0) return(B)
      Xadd <- as.matrix(Xadd)
      # 去掉與 B 完全共線/常數的欄位
      keep_nc <- apply(Xadd, 2, function(z) sd(z) > 0)
      if (any(!keep_nc)) Xadd <- Xadd[, keep_nc, drop = FALSE]
      if (ncol(Xadd)==0) return(B)
      rB <- qr(B)$rank
      for (j in seq_len(ncol(Xadd))) {
        Tst <- cbind(B, Xadd[, j])
        if (qr(Tst)$rank > rB) { B <- Tst; rB <- rB + 1 }
      }
      B
    }
    add_full_rank(mod0, fix_design(mod_interest))
  }
  
  if (ncol(M) != nrow(mod)) {
    log_msg("[SVA%s] DIM MISMATCH: ncol(M)=%d vs nrow(mod)=%d -> SV=NULL", tag, ncol(M), nrow(mod))
    return(NULL)
  }
  
  r_mod  <- qr(mod)$rank
  r_mod0 <- qr(mod0)$rank
  res_df <- ncol(M) - r_mod
  log_msg("[SVA%s] dims: genes=%d, samples=%d | rank(mod)=%s, rank(mod0)=%s | residual_df=%s",
          tag, nrow(M), ncol(M), r_mod, r_mod0, res_df)
  if (!is.finite(res_df) || res_df < 5) { log_msg("[SVA%s] residual df too small -> SV=NULL", tag); return(NULL) }
  
  # ---- 清理非有限值 → 0；濾掉「零變異列」 ----
  M_finite <- M
  M_finite[!is.finite(M_finite)] <- 0
  row_sd <- tryCatch(
    if (requireNamespace("matrixStats", quietly=TRUE)) matrixStats::rowSds(M_finite) else apply(M_finite,1,sd),
    error=function(e) rep(NA_real_, nrow(M_finite))
  )
  keep_var <- is.finite(row_sd) & row_sd > 0
  if (any(!keep_var)) {
    log_msg("[SVA%s] filter: drop %d/%d zero-variance genes (raw)", tag, sum(!keep_var), length(keep_var))
    M_finite <- M_finite[keep_var, , drop=FALSE]
  }
  if (nrow(M_finite) < 50) { log_msg("[SVA%s] too few genes after filter -> SV=NULL", tag); return(NULL) }
  
  # ---- 以與 SVA 相同設計先做「殘差」，再濾掉「零殘差變異列」 ----
  R_pre <- tryCatch({ fit <- stats::lm.fit(x = mod, y = t(M_finite)); t(fit$residuals) }, error=function(e) NULL)
  if (!is.null(R_pre)) {
    rsd <- tryCatch(
      if (requireNamespace("matrixStats", quietly=TRUE)) matrixStats::rowSds(R_pre) else apply(R_pre,1,sd),
      error=function(e) rep(NA_real_, nrow(R_pre))
    )
    keep_r <- is.finite(rsd) & rsd > 0
    if (any(!keep_r)) {
      log_msg("[SVA%s] filter: drop %d/%d zero-variance genes (residual)", tag, sum(!keep_r), length(keep_r))
      M_finite <- M_finite[keep_r, , drop=FALSE]
    }
  }
  
  # ---- 估 k（be/leek）並裁到 res_df-1 與 max_k ----
  k_be <- tryCatch(sva::num.sv(dat = M_finite, mod = mod, method = "be"),
                   error=function(e){ log_msg("[SVA%s] num.sv(be) ERROR: %s", tag, e$message); NA_integer_ })
  k_le <- tryCatch(sva::num.sv(dat = M_finite, mod = mod, method = "leek"),
                   error=function(e){ log_msg("[SVA%s] num.sv(leek) ERROR: %s", tag, e$message); NA_integer_ })
  log_msg("[SVA%s] num.sv -> be=%s, leek=%s", tag, k_be, k_le)
  k <- if (is.finite(k_be) && k_be > 0) k_be else if (is.finite(k_le) && k_le > 0) k_le else 0L
  k <- as.integer(max(0, min(max_k, k, res_df - 1L)))
  
  # 若 k=0，看看殘差 PCA 是否強烈，若 PC1 解釋 >= 20% 就硬給 k=1
  if (k == 0) {
    ve1 <- 0
    R0 <- tryCatch({ fit <- stats::lm.fit(x = mod, y = t(M_finite)); fit$residuals }, error=function(e) NULL)
    if (!is.null(R0)) {
      pc0 <- tryCatch(stats::prcomp(R0, center=TRUE, scale.=FALSE), error=function(e) NULL)
      if (!is.null(pc0)) {
        ve <- 100 * pc0$sdev^2 / sum(pc0$sdev^2); ve1 <- if (length(ve)) ve[1] else 0
      }
    }
    if (is.finite(ve1) && ve1 >= 20) {
      k <- 1L
      log_msg("[SVA%s] k=0 but residual PC1=%.1f%% -> force k=1", tag, ve1)
    } else {
      log_msg("[SVA%s] k=0 and residual PC1=%.1f%% -> SV=NULL", tag, ve1)
      return(NULL)
    }
  }
  
  # ---- 引擎選擇：非 counts → sva；count-like 優先 svaseq，失敗回退 ----
  is_counts <- all(is.finite(M_finite)) && min(M_finite) >= 0 && mean(abs(M_finite - round(M_finite))) < 0.05
  engine <- if (is_counts) "svaseq" else "sva"
  sv <- NULL
  if (engine == "svaseq") {
    sv <- tryCatch(sva::svaseq(dat=M_finite, mod=as.matrix(mod), mod0=as.matrix(mod0), n.sv=k)$sv,
                   error=function(e){ log_msg("[SVA%s] svaseq error: %s -> fallback to sva", tag, e$message); NULL })
    if (is.null(sv)) {
      sv <- tryCatch(sva::sva(dat=M_finite, mod=as.matrix(mod), mod0=as.matrix(mod0), n.sv=k)$sv,
                     error=function(e){ log_msg("[SVA%s] sva fallback error: %s", tag, e$message); NULL })
      engine <- "sva(fallback)"
    }
  } else {
    sv <- tryCatch(sva::sva(dat=M_finite, mod=as.matrix(mod), mod0=as.matrix(mod0), n.sv=k)$sv,
                   error=function(e){ log_msg("[SVA%s] sva error: %s", tag, e$message); NULL })
  }
  
  # ---- 最後防線：若 sva 仍失敗，改用「殘差 PCA」產生 k 個 SV（不中斷）----
  if (is.null(sv)) {
    log_msg("[SVA%s] returned NULL with k=%d -> fallback to residual PCA", tag, k)
    R <- tryCatch({ fit <- stats::lm.fit(x = mod, y = t(M_finite)); fit$residuals }, error=function(e) NULL)
    if (!is.null(R)) {
      R[!is.finite(R)] <- 0
      pc <- tryCatch(stats::prcomp(R, center=TRUE, scale.=FALSE), error=function(e) NULL)
      if (!is.null(pc) && ncol(pc$x) >= k) {
        sv <- pc$x[, seq_len(k), drop=FALSE]
      }
    }
    if (is.null(sv)) { log_msg("[SVA%s] residual PCA fallback also failed -> SV=NULL", tag); return(NULL) }
    engine <- "residual-PCA"
  }
  
  sv <- as.matrix(sv)
  colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
  rownames(sv) <- colnames(M)   # **保證樣本名對齊**
  log_msg("[SVA%s] %s produced %d SV(s)", tag, engine, ncol(sv))
  sv
}


## thin wrapper so old call-sites keep working
do_sva <- function(M, mod, mod0 = NULL, tag = "") {
  # estimate_svs() 會自己加 (Intercept) 到 mod0；為避免重複，這裡拿掉既有截距欄
  if (!is.null(mod0)) {
    mod0 <- as.matrix(mod0)
    if ("(Intercept)" %in% colnames(mod0)) {
      mod0 <- mod0[, setdiff(colnames(mod0), "(Intercept)"), drop = FALSE]
    }
  }
  estimate_svs(
    M             = M,
    mod_interest  = as.matrix(mod),
    mod_nuisance  = mod0,
    max_k         = 5,
    label         = tag
  )
}

## ---- 技術代理變數（穩健版：自動移除零變異欄，PCA 不 scale；不足時回傳 0）----
build_tech_covars <- function(M){
  # 1) 每樣本粗略技術指標（僅做 QC/探索用；TMT 模型不直接納入 median/missrate）
  med  <- if (requireNamespace("matrixStats", quietly = TRUE))
    matrixStats::colMedians(M, na.rm = TRUE) else apply(M, 2, stats::median, na.rm = TRUE)
  miss <- colMeans(!is.finite(M))
  
  # 2) gene-wise z-score 後取樣本 PCA（最多取到 PC2；不足就 0）
  Z <- t(scale(t(M), center = TRUE, scale = TRUE))
  Y <- t(Z)
  sds <- if (requireNamespace("matrixStats", quietly = TRUE))
    matrixStats::colSds(Y, na.rm = TRUE) else apply(Y, 2, stats::sd, na.rm = TRUE)
  keep <- which(is.finite(sds) & sds > 0)
  
  pc1 <- rep(0, nrow(Y)); pc2 <- rep(0, nrow(Y))
  if (length(keep) >= 2) {
    Yk <- Y[, keep, drop = FALSE]
    pc <- tryCatch({
      if (requireNamespace("irlba", quietly = TRUE)) irlba::prcomp_irlba(Yk, n = 2, center = TRUE, scale. = FALSE)
      else stats::prcomp(Yk, rank. = 2, center = TRUE, scale. = FALSE)
    }, error = function(e) NULL)
    if (!is.null(pc)) {
      pc1 <- pc$x[, 1]
      if (ncol(pc$x) >= 2) pc2 <- pc$x[, 2]
    }
  }
  out <- cbind(median = med, missrate = miss, PC1 = pc1, PC2 = pc2)
  rownames(out) <- colnames(M)
  as.matrix(out)
}

## ===== Imputation sensitivity demo: MinProb vs QRILC (one dataset/stratum/predictor) =====
run_imputation_sensitivity_demo <- function(ds_id, ds_dir,
                                            stratum = "ALL",
                                            predictor_name = "CSN_SCORE",
                                            geneset_group = "H") {
  out_dir <- file.path(ds_dir, "csn_gsea_results_TP53", stratum, "impute_sensitivity", predictor_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  mat0_full <- load_matrix_from_dataset_dir(ds_dir)
  keep <- colnames(mat0_full)
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  
  present_sub <- intersect(csn_subunits, rownames(mat0))
  if (predictor_name == "CSN_SCORE") {
    x_vec <- build_csn_score(mat0, subunits = present_sub, combine_7AB = TRUE, min_members = 5L)
    exclude_genes <- present_sub
  } else if (grepl("^RESIDUAL_", predictor_name)) {
    su <- sub("^RESIDUAL_", "", predictor_name)
    base_covars_all <- data.frame(
      purity = as.numeric(get_purity_covariate(ds_id, ds_dir, colnames(mat0))),
      sa_get = get_sex_age_covariates(ds_dir, colnames(mat0))
    )
    base_covars_all$sex <- as.numeric(base_covars_all$sa_get[, "sex"])
    base_covars_all$age <- as.numeric(base_covars_all$sa_get[, "age"])
    base_covars_all$sa_get <- NULL
    rownames(base_covars_all) <- colnames(mat0)
    
    x_vec <- residualize_vector(mat0[su, ], build_csn_score(mat0, present_sub),
                                batch = {bi <- get_batch_factor(ds_dir, colnames(mat0)); if (is.null(bi)) NULL else bi$fac},
                                covars = base_covars_all, min_n = 8L)
    exclude_genes <- su
  } else {
    stopifnot(predictor_name %in% rownames(mat0))
    x_vec <- mat0[predictor_name, ]
    exclude_genes <- predictor_name
  }
  
  keep_rows <- rowMeans(is.finite(mat0)) >= min_frac_complete
  M_minprob <- mat0[keep_rows, , drop = FALSE]
  M_qrilc   <- mat0[keep_rows, , drop = FALSE]
  set.seed(1234); M_minprob <- imputeLCMD::impute.MinProb(M_minprob, q = 0.01)
  set.seed(1234); M_qrilc   <- imputeLCMD::impute.QRILC(M_qrilc)
  
  batch_all  <- { bi <- get_batch_factor(ds_dir, colnames(M_minprob)); if (is.null(bi)) NULL else droplevels(bi$fac) }
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(M_minprob))
  sa_all     <- get_sex_age_covariates(ds_dir, colnames(M_minprob))
  covars_all <- data.frame(
    purity = as.numeric(purity_all),
    sex    = as.numeric(sa_all[, "sex"]),
    age    = as.numeric(sa_all[, "age"]),
    row.names = colnames(M_minprob), check.names = FALSE
  )
  covars_all <- select_covars_safely(covars_all, colnames(M_minprob), label = "impute-demo",
                                     y = as.numeric(x_vec[colnames(M_minprob)]),
                                     min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
                                     max_abs_cor = 0.30)
  
  M1 <- M_minprob[setdiff(rownames(M_minprob), exclude_genes), , drop = FALSE]
  M2 <- M_qrilc  [setdiff(rownames(M_qrilc),   exclude_genes), , drop = FALSE]
  x  <- as.numeric(x_vec[colnames(M1)])
  
  t1 <- tryCatch(limma_cont_with_covars(M1, x, batch = batch_all, covars = covars_all), error = function(e) NULL)
  t2 <- tryCatch(limma_cont_with_covars(M2, x, batch = batch_all, covars = covars_all), error = function(e) NULL)
  
  if (is.null(genesets_by_group[[geneset_group]])) stop("geneset_group 不存在於目前設定")
  gs <- genesets_by_group[[geneset_group]]
  res1 <- run_fgsea_save(t1, gs, file.path(out_dir, "GSEA_limma_t_cont__MinProb"), 0,
                         sprintf("%s | %s | MinProb", ds_id, predictor_name))
  res2 <- run_fgsea_save(t2, gs, file.path(out_dir, "GSEA_limma_t_cont__QRILC"),   0,
                         sprintf("%s | %s | QRILC", ds_id, predictor_name))
  
  pick <- function(dt, tag){
    if (is.null(dt)) return(NULL)
    data.frame(pathway = dt$pathway,
               NES  = dt$NES,
               pval = dt$pval,
               padj = dt$padj,
               stringsAsFactors = FALSE, check.names = FALSE) |>
      dplyr::rename_with(~ paste0(., "_", tag), -pathway)
  }
  tbl <- list(pick(res1, "MinProb"), pick(res2, "QRILC")) |>
    purrr::compact() |> Reduce(function(a,b) dplyr::full_join(a,b,by="pathway"), x = _)
  
  if (!is.null(tbl) && nrow(tbl)) {
    cor_NES <- suppressWarnings(stats::cor(tbl$NES_MinProb, tbl$NES_QRILC, use = "pairwise.complete.obs"))
    n_conc  <- sum(sign(tbl$NES_MinProb) == sign(tbl$NES_QRILC), na.rm = TRUE)
    note <- sprintf("NES correlation (MinProb vs QRILC) = %.3f; sign-concordance = %d/%d",
                    cor_NES, n_conc, nrow(tbl))
    data.table::fwrite(tbl, file.path(out_dir, sprintf("impute_sensitivity_%s_%s.csv", geneset_group, predictor_name)))
    wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "compare")
    openxlsx::writeData(wb, 1, tbl)
    openxlsx::addWorksheet(wb, "notes"); openxlsx::writeData(wb, 2, data.frame(note = note))
    openxlsx::saveWorkbook(wb, file.path(out_dir, sprintf("impute_sensitivity_%s_%s.xlsx", geneset_group, predictor_name)), overwrite = TRUE)
    log_msg("[impute-demo] %s", note)
  } else {
    log_msg("[impute-demo] 無法產生比較表")
  }
}

## ---- 小工具：把任何向量/資料框對齊到指定樣本順序 ----
.align_to_samples <- function(x, sam, what = "covariate") {
  if (is.null(x)) return(NULL)
  if (is.vector(x) || is.factor(x)) {
    # 允許長度剛好等於 sam 而且沒有名字；否則優先用名字對齊
    if (is.null(names(x))) {
      if (length(x) != length(sam)) stop(sprintf("[%s] 長度不符：%d vs %d", what, length(x), length(sam)))
      names(x) <- sam
    }
    if (!all(sam %in% names(x))) stop(sprintf("[%s] 缺少樣本：%s", what, paste(setdiff(sam, names(x)), collapse=", ")))
    out <- x[sam]
    return(out)
  } else {
    x <- as.data.frame(x)
    # 期望 rownames(x) 是樣本
    if (is.null(rownames(x))) {
      if (nrow(x) != length(sam)) stop(sprintf("[%s] 列數不符：%d vs %d 且無 rownames 可對齊", what, nrow(x), length(sam)))
      rownames(x) <- sam
    }
    if (!all(sam %in% rownames(x))) stop(sprintf("[%s] 缺少樣本：%s", what, paste(setdiff(sam, rownames(x)), collapse=", ")))
    out <- x[sam, , drop = FALSE]
    return(out)
  }
}

# 任何物件 → 矩陣；強制設 rownames；若 0 欄或列數不符指定 rows 就回 NULL
.mk_mat_or_null <- function(x, rows){
  if (is.null(x)) return(NULL)
  m <- as.matrix(x)
  if (is.null(rownames(m))) rownames(m) <- rows
  if (nrow(m) != length(rows) || ncol(m) == 0) return(NULL)
  m
}
# 更通用的安全轉矩陣（允許向量、1 欄、或零長度；不指定 rows）
.ensure_mat_or_null <- function(x){
  if (is.null(x)) return(NULL)
  m <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(m) || length(m) == 0) return(NULL)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  if (nrow(m) == 0 || ncol(m) == 0) return(NULL)
  m[!is.finite(m)] <- NA_real_
  m
}


## 把向量安全轉成 factor，並把 NA 显式化成 "NA"（之後再加上前綴）
.factorize_with_explicit_NA <- function(x) {
  x_chr <- as.character(x)
  x_chr[!nzchar(x_chr) | is.na(x_chr)] <- "NA"   # 空字串/NA 都視為 "NA"
  factor(x_chr)
}

## ---- limma t-rank（可併入 batch / 其他協變數；強制對齊）----
## ==== PATCH: limma t with covariates (guard p≈n) ====
limma_t_with_covars <- function(mat, grp2, batch = NULL, covars = NULL) {
  stopifnot(!is.null(colnames(mat)))
  sam <- colnames(mat)
  
  if (is.null(names(grp2))) names(grp2) <- sam
  grp2 <- factor(as.character(grp2[sam]), levels = c("Low","High"))
  if (anyNA(grp2)) stop("grp2 仍有 NA，請先處理分組樣本")
  
  des_grp <- model.matrix(~ 0 + grp2); colnames(des_grp) <- c("Low","High")
  
  Xb <- NULL
  if (!is.null(batch)) {
    b <- .align_to_samples(batch, sam, what = "batch")
    b <- .factorize_with_explicit_NA(b); levels(b) <- paste0("b", make.names(levels(b)))
    if (nlevels(b) >= 1) {
      Xb <- model.matrix(~ 0 + b); colnames(Xb) <- levels(b)
    }
  }
  
  Xc <- NULL
  if (!is.null(covars)) {
    cv <- .align_to_samples(covars, sam, what = "covars")
    cv <- as.matrix(cv)
    if (!is.null(cv) && nrow(cv) == length(sam) && ncol(cv) > 0) Xc <- cv
  }
  if (!is.null(Xc)) {
    Xc[!is.finite(Xc)] <- 0
  }
  
  parts <- list(des_grp)
  if (!is.null(Xb) && nrow(Xb) == length(sam) && ncol(Xb) > 0) parts <- c(parts, list(Xb))
  if (!is.null(Xc) && nrow(Xc) == length(sam) && ncol(Xc) > 0) parts <- c(parts, list(Xc))
  design <- do.call(cbind, parts)
  storage.mode(design) <- "double"
  
  log_msg("  [design] 初始欄位：%s", paste(colnames(design), collapse=", "))
  log_msg("  [design] dim=%dx%d", nrow(design), ncol(design))
  
  # 允許 NA；至少有 3 個有限值且變異>0
  var_ok <- vapply(seq_len(ncol(design)), function(j) {
    z <- design[, j]; fin <- is.finite(z)
    if (sum(fin) < 3) return(FALSE)
    v <- var(z[fin]); is.finite(v) && v > 0
  }, logical(1))
  if (!all(var_ok)) {
    drop <- setdiff(colnames(design)[!var_ok], c("Low","High"))
    if (length(drop)) log_msg("  [design] 零/近零變異 → 刪除：%s", paste(drop, collapse=", "))
    design <- design[, var_ok, drop = FALSE]
  }
  
  if (!all(c("Low","High") %in% colnames(design))) {
    log_msg("  [design] Low/High 不完整 → 略過本次 limma")
    return(NULL)
  }
  
  # 去共線
  q <- qr(design)
  keep_idx <- q$pivot[seq_len(q$rank)]
  keep_nms <- colnames(design)[keep_idx]
  removed  <- setdiff(colnames(design), keep_nms)
  removed  <- setdiff(removed, c("Low","High"))
  if (length(removed)) log_msg("  [design] QR 去共線：移除 %s", paste(removed, collapse=", "))
  design <- design[, keep_nms, drop = FALSE]
  
  # p≈n 防呆
  if (ncol(design) > (ncol(mat) - 2)) {
    keep <- intersect(colnames(design),
                      c("Low","High",
                        grep("^b",  colnames(design), value=TRUE),
                        grep("^SV", colnames(design), value=TRUE),   # ← 保留 SV
                        "purity","sex","age","PC1","PC2"))
    if (length(keep) >= 2) design <- design[, keep, drop=FALSE]
  }
  
  fit  <- limma::lmFit(mat, design)
  fit2 <- limma::contrasts.fit(fit, contrasts = limma::makeContrasts(High - Low, levels = design))
  fit2 <- limma::eBayes(fit2)
  t <- fit2$t[, 1]
  t <- sort(t[is.finite(t)], decreasing = TRUE)
  attr(t, "design_cols") <- colnames(design)
  t
}

## ---- NEW: limma continuous predictor with covariates ----
limma_cont_with_covars <- function(mat, x, batch = NULL, covars = NULL, min_per_group = 8L) {
  stopifnot(!is.null(colnames(mat)))
  sam <- colnames(mat)
  
  if (is.null(names(x))) names(x) <- sam
  x <- as.numeric(x[sam]); names(x) <- sam
  if (sum(is.finite(x)) < 4) stop("連續 predictor 有效樣本太少")
  
  # 1) 組成單一 DF → model.matrix
  DF <- data.frame(x = x, row.names = sam, check.names = FALSE)
  if (!is.null(batch)) {
    b <- .align_to_samples(batch, sam, what = "batch")
    DF[["batch"]] <- droplevels(b)
  }
  if (!is.null(covars)) {
    cv <- .align_to_samples(covars, sam, what = "covars")
    cv <- as.data.frame(cv, check.names = FALSE)
    cv <- coerce_covariates_safely(cv)              # ← 新增
    for (cn in colnames(cv)) DF[[cn]] <- cv[[cn]]   # ← 保留原本的填入 DF
  }
  ok <- stats::complete.cases(DF)
  if (sum(ok) < (2L * min_per_group)) stop(sprintf("可用樣本不足（%d < %d）", sum(ok), 2L * min_per_group))
  
  DF  <- DF[ok, , drop = FALSE]
  M   <- as.matrix(mat[, rownames(DF), drop = FALSE])
  des <- stats::model.matrix(~ 1 + ., data = DF, na.action = stats::na.fail)
  stopifnot(nrow(des) == ncol(M))
  
  fit <- limma::lmFit(M, des)
  # 係數名稱「x」一定存在
  coef_idx <- which(colnames(des) == "x"); if (length(coef_idx) != 1) stop("找不到 x 係數")
  # contrasts.fit 要求對比矩陣；建一個單位向量取出 x 的 t 值
  C <- matrix(0, nrow = ncol(des), ncol = 1, dimnames = list(colnames(des), "x"))
  C["x", 1] <- 1
  fit2 <- limma::contrasts.fit(fit, C); fit2 <- limma::eBayes(fit2)
  
  t <- fit2$t[, 1]
  t <- sort(t[is.finite(t)], decreasing = TRUE)
  attr(t, "design_cols") <- colnames(des)
  t
}

## ---- Spearman 排名（removeBatchEffect；強制對齊）----
## ==== PATCH: Spearman with covariates (impute-for-RBE-then-restore) ====
spearman_rank_with_covars <- function(mat_raw, sub_raw, batch=NULL, covars=NULL,
                                      min_pairs=10, impute_for_rbe=TRUE){
  v <- sub_raw; keep <- !is.na(v)
  if (sum(keep) < min_pairs) stop("可用樣本太少，無法做 Spearman")
  v <- v[keep]; M_raw <- mat_raw[, keep, drop=FALSE]
  des <- stats::model.matrix(~ v)
  
  b <- if (!is.null(batch)) .align_to_samples(batch, colnames(mat_raw), "batch")[keep] else NULL
  if (!is.null(b)) { b <- .factorize_with_explicit_NA(b); levels(b) <- paste0("b", make.names(levels(b))) }
  
  cv <- if (!is.null(covars)) .align_to_samples(covars, colnames(mat_raw), "covars")[keep,,drop=FALSE] else NULL
  if (!is.null(cv)) {
    cv <- as.matrix(as.data.frame(cv, check.names=FALSE)); cv[!is.finite(cv)] <- 0
    if (ncol(cv) > 0) {
      var_ok <- apply(cv, 2, function(z) is.finite(stats::var(z)) && stats::var(z) > 0)
      if (any(!var_ok)) cv <- cv[, var_ok, drop=FALSE]
      if (ncol(cv) > 1) { qrk <- qr(cv); cv <- cv[, qrk$pivot[seq_len(qrk$rank)], drop=FALSE] }
    }
  }
  
  M_adj <- M_raw
  need_adj <- (!is.null(b) && nlevels(b) >= 2) || (!is.null(cv) && ncol(cv) > 0)
  if (need_adj) {
    if (anyNA(M_raw)) {
      if (isTRUE(impute_for_rbe)) {
        na_mask <- is.na(M_raw)
        M_imp <- imputeLCMD::impute.MinProb(M_raw, q=0.01)
        M_tmp <- limma::removeBatchEffect(M_imp, batch=b, covariates=cv, design=des)
        M_adj <- M_tmp; M_adj[na_mask] <- NA  # 關鍵：不增加原本配對數
      } else {
        log_msg("  [spearman] NA 存在且未允許插補 → 回 RAW")
      }
    } else {
      M_adj <- limma::removeBatchEffect(M_raw, batch=b, covariates=cv, design=des)
    }
  }
  
  rho <- apply(M_adj, 1, function(x){
    ok <- stats::complete.cases(x, v)
    if (sum(ok) >= min_pairs) suppressWarnings(stats::cor(x[ok], v[ok], method="spearman")) else NA_real_
  })
  rho <- sort(rho[is.finite(rho)], decreasing=TRUE)
  attr(rho, "batch_levels") <- if (!is.null(b)) levels(b) else NULL
  attr(rho, "covar_cols")   <- if (!is.null(cv) && ncol(cv)>0) colnames(cv) else NULL
  attr(rho, "imputed_for_rbe") <- isTRUE(impute_for_rbe)
  rho
}

## ===== 取得 TP53 mutation status =====
## ===== TP53 狀態（protein-altering baseline）=====
## 只把會改變蛋白的變異算 TP53-mutant；其餘（Silent、UTR、Intron、IGR、RNA、lincRNA、Flank…）皆視為 wild type

TP53_KEEP_CLASSES <- c(
  "MISSENSE_MUTATION", "NONSENSE_MUTATION",
  "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS",
  "IN_FRAME_DEL", "IN_FRAME_INS",
  "SPLICE_SITE", "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
)

normalize_vc <- function(x){
  # 標準化 Variant_Classification：轉大寫、把各種分隔符號統一成底線
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]+", "_", x)
}

get_tp53_status <- function(ds_dir, sample_ids){
  # 預設全部 wild-type
  status <- setNames(rep("TP53_wild_type", length(sample_ids)), sample_ids)
  
  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  if (!file.exists(mut_fp)) {
    log_msg("  [TP53] 找不到 data_mutations.txt，整批視為 wild-type/ALL 可用")
    return(status)
  }
  
  mutation_df <- tryCatch(
    readr::read_tsv(mut_fp, comment = "#", show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(mutation_df)) return(status)
  
  req <- c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")
  miss <- setdiff(req, colnames(mutation_df))
  if (length(miss)) {
    log_msg("  [TP53] data_mutations.txt 缺欄位：%s → 視為 wild-type", paste(miss, collapse=", "))
    return(status)
  }
  
  tp53_df <- subset(mutation_df, Hugo_Symbol == "TP53",
                    select = c("Variant_Classification","Tumor_Sample_Barcode"))
  if (!nrow(tp53_df)) return(status)
  
  vc_norm <- normalize_vc(tp53_df$Variant_Classification)
  
  # 嚴格採 protein-altering 類別；同時容忍 FRAME_SHIFT_* / IN_FRAME_* 這類前綴
  keep <- vc_norm %in% TP53_KEEP_CLASSES |
    grepl("^FRAME_SHIFT_", vc_norm) |
    grepl("^IN_FRAME_",   vc_norm)
  
  if (!any(keep)) return(status)
  
  tp53_samples <- toupper(unique(tp53_df$Tumor_Sample_Barcode[keep]))
  sid_up <- toupper(sample_ids)
  status[sid_up %in% tp53_samples] <- "TP53_mutant"
  status
}

## ===== TP53 status auditing: per-dataset sample counts =====
## 產出：
##   - run_info/tp53_status/<dataset>_tp53_class_sample_counts.csv
##      （原始 Variant_Classification → 「樣本數」；另含 Any_TP53_mutation 與 protein_altering）
##   - run_info/tp53_status/tp53_class_sample_counts_long.csv（彙整 long 格式）
##   - run_info/tp53_status/tp53_binary_counts_by_dataset.csv（wild type / mutant 總表）

dir.create(file.path("run_info","tp53_status"), recursive = TRUE, showWarnings = FALSE)

# 保險：若前面尚未定義，這裡補一份
if (!exists("normalize_vc", mode = "function")) {
  normalize_vc <- function(x){
    x <- toupper(trimws(as.character(x)))
    gsub("[^A-Z0-9]+", "_", x)
  }
}
if (!exists("TP53_KEEP_CLASSES")) {
  TP53_KEEP_CLASSES <- c(
    "MISSENSE_MUTATION","NONSENSE_MUTATION",
    "FRAME_SHIFT_DEL","FRAME_SHIFT_INS",
    "IN_FRAME_DEL","IN_FRAME_INS",
    "SPLICE_SITE","TRANSLATION_START_SITE","NONSTOP_MUTATION"
  )
}

# 單一 dataset 的 TP53 樣本數摘要
summarize_tp53_counts_for_dataset <- function(ds_dir){
  ds_id <- basename(ds_dir)
  
  # 以蛋白矩陣樣本為「母體」
  M <- try(load_matrix_from_dataset_dir(ds_dir), silent = TRUE)
  if (inherits(M, "try-error")) {
    log_msg("[TP53-audit] %s: 無法讀取矩陣，略過", ds_id)
    return(NULL)
  }
  sample_ids <- colnames(M)
  sid_up <- toupper(sample_ids)
  n_all <- length(sample_ids)
  
  # 取得 mutant/wild type（protein-altering 定義）
  status <- get_tp53_status(ds_dir, sample_ids)
  tb_bin <- table(factor(status, levels = c("TP53_wild_type","TP53_mutant")))
  bin_row <- data.frame(
    dataset = ds_id,
    in_matrix_n = n_all,
    WT_n = as.integer(tb_bin["TP53_wild_type"]),
    MUT_n = as.integer(tb_bin["TP53_mutant"]),
    stringsAsFactors = FALSE
  )
  
  # 讀 MAF（原始 Variant_Classification 分布 → 「樣本數」）
  mut_fp <- file.path(ds_dir, "data_mutations.txt")
  class_df <- NULL
  if (file.exists(mut_fp)) {
    mutation_df <- try(mutation_df <- readr::read_tsv(mut_fp, comment = "#", show_col_types = FALSE), silent = TRUE)
    if (!inherits(mutation_df, "try-error") &&
        all(c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode") %in% colnames(mutation_df))) {
      
      tp53 <- subset(mutation_df, Hugo_Symbol == "TP53",
                     select = c("Variant_Classification","Tumor_Sample_Barcode"))
      
      if (nrow(tp53) > 0) {
        tp53$Variant_Classification <- normalize_vc(tp53$Variant_Classification)
        tp53$Tumor_Sample_Barcode   <- toupper(tp53$Tumor_Sample_Barcode)
        
        # 只統計「在蛋白矩陣中」的樣本
        tp53 <- tp53[tp53$Tumor_Sample_Barcode %in% sid_up, , drop = FALSE]
        
        # 每個 classification → 有幾個「獨立樣本」命中
        class_counts <- tp53 |>
          dplyr::group_by(Variant_Classification) |>
          dplyr::summarise(sample_n = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") |>
          dplyr::arrange(dplyr::desc(sample_n), Variant_Classification)
        
        # 加上「Any_TP53_mutation」與「protein_altering」匯總列（樣本層級）
        any_samples <- unique(tp53$Tumor_Sample_Barcode)
        
        keep_flag <- tp53$Variant_Classification %in% TP53_KEEP_CLASSES |
          grepl("^FRAME_SHIFT_", tp53$Variant_Classification) |
          grepl("^IN_FRAME_",   tp53$Variant_Classification)
        prot_alt_samples <- unique(tp53$Tumor_Sample_Barcode[keep_flag])
        
        add_rows <- data.frame(
          Variant_Classification = c("Any_TP53_mutation", "protein_altering"),
          sample_n = c(length(any_samples), length(prot_alt_samples)),
          stringsAsFactors = FALSE
        )
        
        class_df <- dplyr::bind_rows(class_counts, add_rows)
        class_df$dataset <- ds_id
        class_df <- class_df[, c("dataset","Variant_Classification","sample_n")]
      }
    }
  }
  
  # 各自寫檔
  if (!is.null(class_df)) {
    data.table::fwrite(
      class_df,
      file.path("run_info","tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  } else {
    # 若無 MAF 或無 TP53 記錄，輸出空殼
    data.table::fwrite(
      data.frame(dataset = ds_id, Variant_Classification = NA_character_, sample_n = NA_integer_),
      file.path("run_info","tp53_status", paste0(ds_id, "_tp53_class_sample_counts.csv"))
    )
  }
  
  list(binary = bin_row, class_long = class_df)
}

# 跑所有 datasets，彙整總表
summarize_tp53_counts_all_datasets <- function(dataset_dirs){
  all_bin <- list(); all_class <- list(); k <- 1L; j <- 1L
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) next
    log_msg("[TP53-audit] 開始：%s", ds)
    res <- summarize_tp53_counts_for_dataset(ds_dir)
    if (is.null(res)) next
    all_bin[[k]]   <- res$binary; k <- k + 1L
    if (!is.null(res$class_long)) { all_class[[j]] <- res$class_long; j <- j + 1L }
  }
  
  if (length(all_bin)) {
    bin_df <- dplyr::bind_rows(all_bin) |>
      dplyr::mutate(Any_TP53_mutation_n = in_matrix_n - WT_n)
    data.table::fwrite(bin_df, file.path("run_info","tp53_status","tp53_binary_counts_by_dataset.csv"))
    log_msg("[TP53-audit] 寫出：tp53_binary_counts_by_dataset.csv")
  }
  if (length(all_class)) {
    class_df <- dplyr::bind_rows(all_class)
    data.table::fwrite(class_df, file.path("run_info","tp53_status","tp53_class_sample_counts_long.csv"))
    log_msg("[TP53-audit] 寫出：tp53_class_sample_counts_long.csv")
  }
}

## === 呼叫（放在主流程前或後皆可；不影響下游分析）===
##summarize_tp53_counts_all_datasets(dataset_dirs)

## === robust run_predictor_analyses (center predictors for all limma models; keep Spearman unchanged) ===
run_predictor_analyses <- function(
    predictor_name,
    predictor_vec,             # named numeric by sample (可含 NA)
    exclude_genes = NULL,      # 在排名前要排除的基因（例如自體/整個 complex）
    ds_id, ds_dir,
    mat0,                      # RAW（本函式處理 limma 連續；Spearman 另處理）
    mat,                       # impute+filter 後（limma 用）
    out_root,
    genesets_by_group,         # list(group -> pathways)
    batch_all = NULL,          # factor by sample 或 NULL
    purity_all = NULL,         # named numeric by sample 或 NULL
    sa_all = NULL,             # data.frame(sex, age[, age_missing, age_z_imputed]); rownames=sample
    tp53_num_all = NULL,       # named numeric (0/1) 或 NULL
    is_ALL = FALSE             # TRUE 時才可能加 TP53_mutant
){
  ## -------- 安全參數 --------
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
  min_per_group <- opt("min_per_group", 8L)
  minSize       <- opt("minSize", 15L)
  maxSize       <- opt("maxSize", 500L)
  fgsea_eps     <- opt("fgsea_eps", 0)
  USE_AGE_MISSING_INDICATOR <- isTRUE(opt("USE_AGE_MISSING_INDICATOR", FALSE))
  
  logf <- function(...) if (exists("log_msg", mode="function")) try(log_msg(...), silent=TRUE)
  sfn  <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  
  stopifnot(!is.null(genesets_by_group), length(genesets_by_group) > 0)
  
  ## -------- 對齊樣本（以 predictor 為主）--------
  if (is.null(names(predictor_vec))) names(predictor_vec) <- colnames(mat)
  keep <- intersect(colnames(mat), names(predictor_vec))
  pred_all <- suppressWarnings(as.numeric(predictor_vec[keep])); names(pred_all) <- keep
  fin <- is.finite(pred_all)
  if (sum(fin) < (2L * min_per_group)) {
    logf("  [%s] predictor 非NA 樣本不足（%d < %d）→ skip", predictor_name, sum(fin), 2L * min_per_group)
    return(invisible(NULL))
  }
  sample_order <- keep[fin]
  pred <- pred_all[fin]                  # ← 原始（未中心化）predictor；Spearman 會用這個
  stopifnot(identical(names(pred), sample_order))
  
  ## -------- 建共變項資料框（未強補 NA；僅對齊樣本順序）--------
  build_covars_df <- function(so){
    df <- data.frame(row.names = so, check.names = FALSE)
    if (!is.null(purity_all)) df$purity <- suppressWarnings(as.numeric(purity_all[so]))
    if (!is.null(sa_all)) {
      if ("sex" %in% colnames(sa_all)) df$sex <- factor(sa_all[so, "sex"])
      if ("age" %in% colnames(sa_all)) df$age <- suppressWarnings(as.numeric(sa_all[so, "age"]))
      if (USE_AGE_MISSING_INDICATOR) {
        if ("age_missing"   %in% colnames(sa_all)) df$age_missing   <- suppressWarnings(as.numeric(sa_all[so, "age_missing"]))
        if ("age_z_imputed" %in% colnames(sa_all)) df$age_z_imputed <- suppressWarnings(as.numeric(sa_all[so, "age_z_imputed"]))
      }
    }
    # [POLICY] TP53 covariate disabled in ALL strata
    df
  }
  df_covars0 <- build_covars_df(sample_order)
  
  coerce_covariates_safely <- function(df){
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)
    
    for (cn in colnames(df)) {
      v <- df[[cn]]
      
      if (is.factor(v) || is.character(v) || is.logical(v)) {
        v <- factor(v)  # 保留為因子
        # 只用非 NA 樣本檢查實際層級數
        lv <- levels(droplevels(v[!is.na(v)]))
        if (length(lv) <= 1) {
          # 單層級因子 → 丟掉，避免 contrasts 錯誤
          keep[cn] <- FALSE
          if (exists("logf")) try(logf("  [covars] drop single-level factor: %s", cn), silent = TRUE)
        } else {
          df[[cn]] <- v
        }
        
      } else {
        # 數值就保持數值；不要硬轉因子
        df[[cn]] <- suppressWarnings(as.numeric(v))
        # （可選）若你想一併去掉「完全常數」數值欄，解除下面註解：
        # vv <- df[[cn]]; if (sum(!is.na(vv)) >= 2 && isTRUE(all(vv[!is.na(vv)] == vv[which(!is.na(vv))[1]]))) {
        #   keep[cn] <- FALSE
        #   if (exists("logf")) try(logf("  [covars] drop constant numeric: %s", cn), silent = TRUE)
        # }
      }
    }
    
    df <- df[, keep, drop = FALSE]
    df
  }
  
  
  ## -------- coverage 門檻 --------
  base_thr <- c(purity = 0.60, sex = 0.80, age = 0.80)
  present  <- colnames(df_covars0)
  extra    <- setdiff(present, names(base_thr))
  if (length(extra)) base_thr <- c(base_thr, stats::setNames(rep(min(base_thr), length(extra)), extra))
  
  ## -------- 擇共變項（回傳與 sample_order 對齊的 numeric data.frame）--------
  pick_covars_df <- function(label){
    ## NEW: 決定共變數來源（covars_all 不在此 scope 就回退到 df_covars0）
    cov_src <- tryCatch(get("covars_all", inherits = TRUE), error = function(e) NULL)
    if (is.null(cov_src)) {
      cov_src <- tryCatch(get("df_covars0", inherits = TRUE), error = function(e) NULL)
    }
    if (is.null(cov_src)) {
      logf("  [covars-%s] NO covariate table in scope → skip", label)
      return(NULL)
    }
    sel <- tryCatch(
      select_covars_safely(
        df           = cov_src,
        sample_order = sample_order,
        label        = label,
        y            = pred
      ),
      error = function(e) {
        logf("  [covars-%s] SELECT-FAIL: %s", label, conditionMessage(e))
        ## 診斷資訊：同樣用 cov_src，避免再因不存在物件出錯
        rns   <- tryCatch(rownames(cov_src), error = function(e) NULL)
        nrowC <- tryCatch(NROW(cov_src),     error = function(e) -1L)
        ncolC <- tryCatch(NCOL(cov_src),     error = function(e) -1L)
        hasrn <- is.character(rns) || is.factor(rns)
        logf("  [covars-%s] diag: C_dim = %d x %d; has_rownames=%s", label, nrowC, ncolC, hasrn)
        miss <- tryCatch(setdiff(sample_order, rns), error = function(e) character(0))
        logf("  [covars-%s] samples_not_in_covars = %d / %d (eg: %s)",
             label, length(miss), length(sample_order), paste(head(miss, 5), collapse = ","))
        NULL
      }
    )
    if (is.null(sel)) return(NULL)
    # === NEW: 先把 select_covars_safely() 回傳的 drop 原因抓出來（先存變數） ===
    dl_sel <- tryCatch(attr(sel, "drop_log"), error = function(e) NULL)
    ## === /NEW ===
    X <- NULL
    if (is.character(sel) && is.null(dim(sel))) {
      nm <- intersect(sel, colnames(cov_src))
      if (length(nm)) X <- cov_src[, nm, drop = FALSE]
    } else {
      X <- as.data.frame(sel, stringsAsFactors = FALSE)
      if (is.null(rownames(X))) {
        if (nrow(X) == length(sample_order)) rownames(X) <- sample_order
        else if (!is.null(colnames(X)) && ncol(X) == length(sample_order) && all(colnames(X) == sample_order)) {
          X <- as.data.frame(t(as.matrix(X)))
        } else {
          logf("  [covars-%s] 回傳形狀/對齊不明，捨棄", label); return(NULL)
        }
      }
      X <- X[sample_order, , drop = FALSE]
    }
    if (is.null(X) || !ncol(X)) return(NULL)
    X <- coerce_covariates_safely(X)  # 讓 character/logical → factor；數值保持 numeric
    # 因子友善的最小品質門檻：至少有 3 個非 NA
    good <- vapply(X, function(v) sum(!is.na(v)) >= 3, logical(1))
    if (any(!good)) {
      logf("  [covars-%s] drop low-coverage columns: %s", label, paste(colnames(X)[!good], collapse=","))
      X <- X[, good, drop = FALSE]
    }
    if (!ncol(X)) return(NULL)
    logf("  [covars-picked:%s] kept = {%s}", label, if (!is.null(X)) paste(colnames(X), collapse=",") else "NULL")
    ## === NEW: 把被丟掉的 covariate 與原因落地到 CSV（追加模式） ===
    if (!is.null(dl_sel) && nrow(dl_sel)) {
      dl_sel$dataset <- ds_id
      dl_sel$stratum <- basename(out_root)
      dl_sel$subunit <- predictor_name
      dl_sel$pass    <- label                          # "limma-cont:RAW" 或 "limma-cont:base"
      dl_sel$samples <- length(sample_order)
      
      out_dir <- file.path("run_info", "covars_audit")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fp <- file.path(out_dir, "covariate_drop_reasons.csv")
      
      data.table::fwrite(dl_sel, fp, append = file.exists(fp))
    }
    ## === /NEW ===
    X
  }
  X_raw_cov <- pick_covars_df("limma-cont:RAW")
  if (!is.null(X_raw_cov)) X_raw_cov <- coerce_covariates_safely(X_raw_cov)
  X_ba_cov  <- pick_covars_df("limma-cont:base")
  if (!is.null(X_ba_cov))  X_ba_cov  <- coerce_covariates_safely(X_ba_cov)
  ## -------- RAW：建 DF→剔 NA→**中心化 predictor**→model.matrix --------
  DF_raw <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_raw_cov)) for (nm in colnames(X_raw_cov)) DF_raw[[nm]] <- X_raw_cov[[nm]]
  ok_raw <- stats::complete.cases(DF_raw)
  if (sum(ok_raw) < (2L * min_per_group)) {
    logf("  [%s|RAW] 可用樣本不足（%d < %d）→ skip", predictor_name, sum(ok_raw), 2L * min_per_group)
    return(invisible(NULL))
  }
  DF_raw <- DF_raw[ok_raw, , drop = FALSE]
  ## ——中心化（只移動平均，不縮放）：不影響係數/檢定，僅截距解釋不同
  DF_raw$predictor <- DF_raw$predictor - mean(DF_raw$predictor, na.rm = TRUE)
  M_raw   <- as.matrix(mat[, rownames(DF_raw), drop = FALSE])
  des_raw <- stats::model.matrix(~ 1 + ., data = DF_raw, na.action = stats::na.fail)
  stopifnot(nrow(des_raw) == ncol(M_raw))
  # 稽核：覆蓋率/共變項摘要 + 對齊
  audit_covars_coverage(
    tag       = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    ds_id     = ds_id,
    stratum   = basename(out_root),
    su        = predictor_name,
    sample_ids= rownames(DF_raw),
    batch     = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_raw)])) else NULL,
    covars    = DF_raw[, setdiff(colnames(DF_raw), "predictor"), drop = FALSE],
    sv        = NULL,
    tech      = NULL
  )  # → 會 append 到 run_info/covars_audit/audit_rows.csv  :contentReference[oaicite:17]{index=17}
  
  audit_design_alignment(
    tag         = sprintf("%s_%s_%s|RAW", ds_id, basename(out_root), predictor_name),
    samples     = colnames(M_raw),
    mod_interest= des_raw,
    mod_nuisance= NULL,
    out_dir     = file.path("run_info","covars_audit")
  )  # → 會寫 ALIGN_RAW.csv 類檔名  :contentReference[oaicite:18]{index=18}
  
  ## -------- BatchAdj：同上；batch 進設計表 → **中心化 predictor** --------
  DF_ba <- data.frame(predictor = as.numeric(pred), row.names = sample_order, check.names = FALSE)
  if (!is.null(X_ba_cov)) for (nm in colnames(X_ba_cov)) DF_ba[[nm]] <- X_ba_cov[[nm]]
  if (!is.null(batch_all)) DF_ba$batch <- droplevels(as.factor(batch_all[sample_order]))
  ## NEW: 若沒有 batch，嘗試 SVA；失敗再退 tech PCs
  sv_ba   <- NULL
  tech_ba <- NULL
  if (is.null(batch_all)) {
    logf("[SVA|pre] DF_ba cols = {%s}", paste(colnames(DF_ba), collapse = ","))
    ## 1) SVA 用的設計（允許 NA，避免丟樣本 → 與 M 對齊）
    ## --- 組 des_for_sva / des0_for_sva：用欄名動態建式（0 欄退成 ~ 1） ---
    ## --- 組 des_for_sva / des0_for_sva（動態公式已就位）---
    DF0 <- DF_ba[, setdiff(colnames(DF_ba), "predictor"), drop = FALSE]
    
    form_all <- if (ncol(DF_ba) > 0) {
      as.formula(paste("~ 1 +", paste(colnames(DF_ba), collapse = " + ")))
    } else { ~ 1 }
    
    form_nui <- if (ncol(DF0) > 0) {
      as.formula(paste("~ 1 +", paste(colnames(DF0), collapse = " + ")))
    } else { ~ 1 }
    
    ## ====== NEW: 先決定 SVA 要用的共同樣本集合，避免 95 vs 97 不一致 ======
    samp_sva <- intersect(colnames(mat), rownames(DF_ba))
    if (!length(samp_sva)) {
      logf("[SVA%s] no overlapping samples between M and DF_ba; skip SVA", sprintf("|%s|BatchAdj", predictor_name))
    } else {
      ## 以共同樣本建立 SVA 專用的資料與設計
      DF_ba_sva <- DF_ba[samp_sva, , drop = FALSE]
      
      des_for_sva <- stats::model.matrix(form_all, data = DF_ba_sva, na.action = stats::na.pass)
      des0_for_sva <- stats::model.matrix(
        form_nui,
        data = if (ncol(DF0) > 0) DF_ba_sva[, setdiff(colnames(DF_ba_sva), "predictor"), drop = FALSE]
        else DF_ba_sva[, 0, drop = FALSE],
        na.action = stats::na.pass
      )
      
      ## 若 model.matrix 仍產生了列名不一致（極少見），再做一次交集對齊
      common <- Reduce(intersect, list(samp_sva, rownames(des_for_sva), rownames(des0_for_sva)))
      if (length(common) < length(samp_sva)) {
        des_for_sva  <- des_for_sva [common, , drop = FALSE]
        des0_for_sva <- des0_for_sva[common, , drop = FALSE]
        samp_sva     <- common
      }
      
      ## 用完全一致的樣本子集出 M
      M_sva <- as.matrix(mat[, samp_sva, drop = FALSE])
      
      # --- preflight: 可視化列出非有限值數量（只記錄 log，不改控制流）
      nf_cnt <- sum(!is.finite(M_sva), na.rm = TRUE)
      if (nf_cnt > 0) {
        logf(sprintf("[SVA|pre] M_sva 非有限值（含NA/NaN/Inf）個數 = %d；將在 estimate_svs() 內轉為 0", nf_cnt))
      }
      
      ## 呼叫 SVA（其餘程式不動）
      sv_ba <- tryCatch(
        {
          do_sva <- get("do_sva", mode = "function")
          do_sva(
            M   = M_sva,
            mod = des_for_sva,
            mod0 = des0_for_sva,
            tag  = sprintf("|%s|BatchAdj", predictor_name)
          )
        },
        error = function(e) NULL
      )
    }
    
    ## --- NEW: 給 SVA 產出的矩陣補齊樣本列名（保證與 M_sva 一致） ---
    if (!is.null(sv_ba)) {
      rn <- rownames(sv_ba)
      need_fix <- is.null(rn) || !length(rn) || all(rn == seq_len(nrow(sv_ba)))
      if (need_fix) rownames(sv_ba) <- colnames(M_sva)
    }
    
    ## --- NEW: 用共同樣本交集呼叫 audit（避免「no overlap」） ---
    try({
      common <- intersect(colnames(M_sva), rownames(sv_ba))
      if (length(common)) {
        audit_sva_sanity(
          ## 與 SVA 同一批樣本
          M   = as.matrix(M_sva[, common, drop = FALSE]),
          sv  = as.matrix(sv_ba[common, , drop = FALSE]),
          
          ## 關鍵：給「帶樣本名」的向量，名稱須對齊 sample_order = common
          predictor = setNames(as.numeric(DF_ba[common, "predictor", drop = TRUE]), common),
          
          ## 可選：TP53 也給「帶樣本名」的向量
          tp53 = {
            if (is_ALL && exists("tp53_num_all") && !is.null(tp53_num_all)) {
              setNames(suppressWarnings(as.numeric(tp53_num_all[sample_order])), sample_order)
            } else {
              NULL
            }
          },
          
          ## 正確參數名是 tag（不是 label）
          tag          = sprintf("%s_%s_%s|BatchAdj", ds_id, basename(out_root), predictor_name),
          out_dir      = file.path("run_info","covars_audit"),
          sample_order = common
        )
      } else {
        if (exists("log_msg", mode="function"))
          log_msg("[SVA-sanity|%s|BatchAdj] still no overlap even after renaming; head(M)=%s; head(SV)=%s",
                  predictor_name,
                  paste(head(colnames(M_sva), 3), collapse=","),
                  paste(head(rownames(sv_ba), 3), collapse=","))
      }
    }, silent = FALSE)
    
    ## ====== /NEW ======
    ## --- [NEW] 稽核 SVA：只要抓到 SV，就把與 predictor/TP53 的關係寫到 run_info/covars_audit ---
    if (!is.null(sv_ba) && getOption("csn.audit_sva_sanity", FALSE)) {
      if (exists("audit_sva_sanity", mode = "function")) {
        tp53_col <- intersect(c("TP53_mutant","TP53_status","TP53"), colnames(DF_ba))
        tp53_vec <- if (length(tp53_col)) DF_ba[[tp53_col[1]]] else NULL
        try({
          sva_samps <- rownames(sv_ba)
          audit_sva_sanity(
            M          = as.matrix(mat[, sva_samps, drop = FALSE]),
            sv         = as.matrix(sv_ba),
            predictor  = setNames(as.numeric(DF_ba[sva_samps, "predictor"]), sva_samps),
            tp53 = {
              if (is_ALL && exists("tp53_num_all") && !is.null(tp53_num_all)) {
                setNames(suppressWarnings(as.numeric(tp53_num_all[sample_order])), sample_order)
              } else {
                NULL
              }
            },
            tag          = sprintf("%s_%s_%s|BatchAdj", ds_id, basename(out_root), predictor_name),
            out_dir     = file.path("run_info","covars_audit"),
            sample_order= sva_samps
          )
        }, silent = FALSE)
      }
    }
    ## 2) 若 SVA 仍失敗（例如殘差 df 過小），改用 tech PCs（這裡需要沒有 NA 的設計）
    if (is.null(sv_ba) || ncol(as.matrix(sv_ba)) == 0) {
      form_lm <- if (ncol(DF_ba) > 0) {
        as.formula(paste("~ 1 +", paste(colnames(DF_ba), collapse = " + ")))
      } else {
        ~ 1
      }
      des_for_lm <- tryCatch(
        stats::model.matrix(form_lm, data = DF_ba, na.action = stats::na.fail),
        error = function(e) NULL
      )
      if (!is.null(des_for_lm)) {
        samp_ok <- rownames(des_for_lm)
        R <- tryCatch({
          fit <- stats::lm.fit(
            x = des_for_lm,
            y = t(as.matrix(mat[, samp_ok, drop = FALSE]))
          )
          fit$residuals
        }, error = function(e) NULL)
        
        if (!is.null(R)) {
          pc <- tryCatch(stats::prcomp(R, center = TRUE, scale. = FALSE), error = function(e) NULL)
          if (!is.null(pc) && ncol(pc$x) >= 2) {
            tech_raw <- pc$x[, 1:2, drop = FALSE]
            colnames(tech_raw) <- paste0("PC", seq_len(ncol(tech_raw)))
            ## gate（避免吃到生物訊號）
            tech_ba <- tryCatch(gate_tech_pcs(tech_raw, DF_ba$predictor), error = function(e) NULL)
          }
        }
      }
    }
    
    ## NEW: SVA 回來後，強制補齊 rownames / colnames，確保能用樣本名對齊掛回 DF_ba
    if (!is.null(sv_ba)) {
      sv_ba <- as.matrix(sv_ba)
      
      ## 1) 行數應等於 SVA 使用的樣本數（samp_sva）
      if (nrow(sv_ba) == length(samp_sva)) {
        if (is.null(rownames(sv_ba)) || any(is.na(rownames(sv_ba))) || any(rownames(sv_ba) == "")) {
          rownames(sv_ba) <- samp_sva
        }
      } else if (ncol(sv_ba) == length(samp_sva) && nrow(sv_ba) < ncol(sv_ba)) {
        ## 有些 SVA/包裝會回傳 t(sv)；偵測到就轉置回來
        sv_ba <- t(sv_ba)
        rownames(sv_ba) <- samp_sva
      } else {
        logf("[SVA%s] WARNING: unexpected SV shape: %d x %d (samples=%d)",
             sprintf("|%s|BatchAdj", predictor_name),
             nrow(sv_ba), ncol(sv_ba), length(samp_sva))
      }
      
      ## 2) 確保有 SV1..SVk 的欄名
      if (is.null(colnames(sv_ba)) || any(colnames(sv_ba) == "")) {
        colnames(sv_ba) <- paste0("SV", seq_len(ncol(sv_ba)))
      }
    }
    
    ## NEW: gate SVA against biology (predictor & optional TP53)
    if (!is.null(sv_ba) && NROW(sv_ba) && NCOL(sv_ba) && exists("gate_tech_pcs", mode = "function")) {
      # 對 predictor（連續）做 gate
      sv_ba <- tryCatch(gate_tech_pcs(as.matrix(sv_ba)
                                      
                                      
                                      , DF_ba$predictor),
                        error = function(e){ if (exists("log_msg",mode="function")) log_msg("[gate_tech_pcs|pred] %s", conditionMessage(e)); sv_ba })
      
      ## [POLICY] TP53-only-at-SVA: gate SVs by TP53 in ALL (without using TP53 as covariate)
      if (is_ALL && exists("tp53_num_all") && !is.null(tp53_num_all)) {
        try({
          sv_ba <- gate_tech_pcs(as.matrix(sv_ba), suppressWarnings(as.numeric(tp53_num_all[colnames(sv_ba)])))
          logf("[gate-SVA] secondary TP53 gate applied")
        }, silent = TRUE)
      }
      
      # 若有 TP53（二元/因子）再做一次 gate（以 0/1 numeric 近似）
      tp53_col <- intersect(c("TP53_mutant","TP53_status","TP53"), colnames(DF_ba))
      if (length(tp53_col)) {
        tp53_vec <- tryCatch(suppressWarnings(as.numeric(DF_ba[[tp53_col[1]]])),
                             error = function(e) NULL)
        if (!is.null(tp53_vec)) {
          sv_ba <- tryCatch(gate_tech_pcs(as.matrix(sv_ba), tp53_vec),
                            error = function(e){ if (exists("log_msg",mode="function")) log_msg("[gate_tech_pcs|TP53] %s", conditionMessage(e)); sv_ba })
        }
      }
    }
    if (!is.null(sv_ba)) {
      kept <- colnames(as.matrix(sv_ba))
      if (is.null(kept) || !length(kept)) kept <- "none"
      if (exists("log_msg", mode="function")) {
        try(log_msg(sprintf("[gate-SVA] kept SV = {%s}", paste(kept, collapse=","))), silent = FALSE)
      }
    }
    
    ## === /SVA 稽核 ===
    
    ## 3) 把抓到的 SV 或 tech PCs 掛回 DF_ba
    ## 3) 把抓到的 SV 或 tech PCs 掛回 DF_ba（用樣本名對齊；缺者補 NA）
    attach_by_name <- function(target_df, add_mat) {
      if (is.null(add_mat)) return(target_df)
      add_mat <- as.matrix(add_mat)
      rn_tgt  <- rownames(target_df)
      ## OPTIONAL: 若 add_mat 沒有 rownames 但行數吻合，直接套用 target 的 rownames
      if (is.null(rownames(add_mat)) && nrow(add_mat) == length(rn_tgt)) {
        rownames(add_mat) <- rn_tgt
      }
      
      for (j in seq_len(ncol(add_mat))) {
        v <- rep(NA_real_, length(rn_tgt)); names(v) <- rn_tgt
        common <- intersect(rn_tgt, rownames(add_mat))
        if (length(common)) v[common] <- add_mat[common, j]
        target_df[[colnames(add_mat)[j]]] <- v
      }
      target_df
    }
    
    DF_ba <- attach_by_name(DF_ba, sv_ba)
    DF_ba <- attach_by_name(DF_ba, tech_ba)
  }
  ## NEW: 刪除全 NA 欄（例如被 coverage 丟掉但殘留的 purity 版本）
  all_na_col <- vapply(DF_ba, function(v) all(is.na(v)), logical(1))
  if (any(all_na_col)) {
    logf("[BatchAdj] drop all-NA columns: %s",
         paste(names(all_na_col)[all_na_col], collapse = ","))
    DF_ba <- DF_ba[, !all_na_col, drop = FALSE]
  }
  
  ## NEW: 取這一輪會真正用到的欄位（predictor + 已選 covars + SV/PC）
  use_cols <- unique(c(
    "predictor",
    intersect(c("purity","sex","age","batch","TP53_mutant"), colnames(DF_ba)),
    grep("^SV\\d+$|^PC\\d+$", colnames(DF_ba), value = TRUE)
  ))
  use_cols <- use_cols[use_cols %in% colnames(DF_ba)]  # 防守
  
  ## 以「要用到的欄」檢查完整列
  ok_rows <- complete.cases(DF_ba[, use_cols, drop = FALSE])
  n_ok    <- sum(ok_rows)
  
  if (n_ok < 16L) {
    logf("  [%s|BatchAdj] 可用樣本不足（%d < 16）→ skip", predictor_name, n_ok)
    return(invisible(NULL))   # ← 與你原本的控制流一致
  }
  
  ## 對齊到完整列的樣本（這一步很關鍵）
  DF_ba <- DF_ba[ok_rows, , drop = FALSE]
  
  ## ——中心化（保留你原有的設計）
  DF_ba$predictor <- DF_ba$predictor - mean(DF_ba$predictor, na.rm = TRUE)
  
  ## 用最終樣本子集出 M_ba（保留你原有的寫法）
  M_ba <- as.matrix(mat[, rownames(DF_ba), drop = FALSE])
  
  form_ba <- if (ncol(DF_ba) > 0) {
    as.formula(paste("~ 1 +", paste(colnames(DF_ba), collapse = " + ")))
  } else {
    ~ 1
  }
  des_ba <- stats::model.matrix(form_ba, data = DF_ba, na.action = stats::na.fail)
  stopifnot(nrow(des_ba) == ncol(M_ba))
  audit_covars_coverage(
    tag        = sprintf("%s|BatchAdj", predictor_name),
    ds_id      = ds_id,
    stratum    = basename(out_root),
    su         = predictor_name,
    sample_ids = rownames(DF_ba),
    batch      = if (!is.null(batch_all)) droplevels(as.factor(batch_all[rownames(DF_ba)])) else NULL,
    covars     = DF_ba[, intersect(c("purity","sex","age","batch"), colnames(DF_ba)), drop = FALSE],
    sv         = DF_ba[, grep("^SV\\d+$", colnames(DF_ba), perl = TRUE),  drop = FALSE],
    tech       = DF_ba[, grep("^PC\\d+$", colnames(DF_ba), perl = TRUE),  drop = FALSE]
  )
  
  ## --- NEW: 只在 RESIDUAL_*（BatchAdj）加上 dataset/stratum prefix；其餘維持原樣 ---
  tag_ba_align <- if (grepl("^RESIDUAL_", predictor_name)) {
    sprintf("%s_%s_%s|BatchAdj", ds_id, basename(out_root), predictor_name)
  } else {
    sprintf("%s|BatchAdj", predictor_name)
  }
  audit_design_alignment(
    tag         = tag_ba_align,
    samples     = colnames(M_ba),
    mod_interest= des_ba,
    mod_nuisance= NULL,
    out_dir     = file.path("run_info","covars_audit")
  )
  
  ## -------- Spearman（RAW/BatchAdj）：維持「未中心化」predictor --------
  if (.RUN_SPEARMAN) {
    audit_spearman_pairs(
      predictor_name       = predictor_name,
      predictor_vec        = pred,               # names(pred) == sample_order
      mat0                 = mat0,               # 原始表達矩陣（未插補；函式內會用 rownames(mat0)）
      out_dir              = file.path("run_info","covars_audit"),
      min_pairs_spearman   = max(16L, get0("min_pairs_spearman", ifnotfound = 16L))
    )  # → 會寫 spearman_pairs_audit.csv 與 summary.txt  :contentReference[oaicite:19]{index=19}
    
    ## -------- Spearman（RAW/BatchAdj）：維持「未中心化」predictor --------
    # 先做配對數稽核（不變）
    audit_spearman_pairs(
      predictor_name = predictor_name,
      predictor_vec  = pred,
      mat0           = mat0,
      out_dir        = file.path("run_info","covars_audit"),
      min_pairs_spearman = max(16L, get0("min_pairs_spearman", ifnotfound = 16L))
    )
    
    # 對 genesets_by_group 的每一個 group 都跑（H, C6, C5:BP, ...）
    for (grp in names(genesets_by_group)) {
      pw <- genesets_by_group[[grp]]
      if (is.null(pw) || !length(pw)) next
      
      ## RAW
      if (!("RAW" %in% getOption("csn.run_passes", c("BatchAdj")))) { log_msg("[spearman:RAW] skip: pass 未選") } else
        try({
          log_msg("[spearman:RAW] start; subunit=%s | group=%s", predictor_name, grp)
          expr_raw <- as.matrix(mat0[, colnames(M_raw), drop = FALSE])
          pred_raw <- pred[rownames(DF_raw)]; names(pred_raw) <- rownames(DF_raw)
          stats_raw_sp <- spearman_prerank(expr_raw, predictor = pred_raw,
                                           min_non_na = max(16L, get0("min_pairs_spearman", ifnotfound = 16L)))
          if (!is.null(stats_raw_sp) && length(stats_raw_sp) && any(is.finite(stats_raw_sp))) {
            run_and_write_gsea_spearman(
              stats_raw_sp, pw, grp,
              out_dir_root = file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp), "RAW", ds_id, basename(out_root)),
              subunit      = predictor_name,
              pass_label   = "RAW"
            )
          } else {
            log_msg("[spearman:RAW] skip: empty/NA stats after QC")
          }
        }, silent = TRUE)
      
      ## BatchAdj（partial）
      if (!("BatchAdj" %in% getOption("csn.run_passes", c("BatchAdj")))) { log_msg("[spearman:BatchAdj|partial] skip: pass 未選") } else
        try({
          log_msg("[spearman:BatchAdj|partial] start; subunit=%s | group=%s", predictor_name, grp)
          stats_ba_sp <- spearman_prerank_partial(
            M_ba, DF_ba,
            residualize_predictor = TRUE,
            min_non_na = max(16L, get0("min_pairs_spearman", ifnotfound = 16L))
          )
          run_and_write_gsea_spearman(
            stats_ba_sp, pw, grp,
            out_dir_root = file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp), "BatchAdj", ds_id, basename(out_root)),
            subunit      = predictor_name,
            pass_label   = "BatchAdj"
          )
        }, silent = TRUE)
    }
  }
  
  
  ## --------（可選）TP53：interaction 與差異相關（此處皆已用中心化後的 DF_raw/DF_ba$predictor）--------
  if (isTRUE(is_ALL) && !is.null(tp53_num_all) &&
      length(genesets_by_group) &&
      any(vapply(genesets_by_group, function(x) length(x) > 0L, logical(1)))) {
    ## ---------- RAW branch（interaction + diffcorr）----------
    if (!("RAW" %in% getOption("csn.run_passes", c("BatchAdj")))) { log_msg("[RAW|interaction] skip: pass 未選") } else
      try({
        tp53_raw_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_raw)]))
        tp53_raw_fac <- factor(ifelse(tp53_raw_num == 1, "MT", "WT"), levels = c("WT","MT"))
        covars_only_raw <- DF_raw[, setdiff(colnames(DF_raw),
                                            c("predictor","TP53_mutant","TP53_status","TP53")), drop = FALSE]
        df_int_raw <- data.frame(pred = DF_raw$predictor, tp53 = tp53_raw_fac,
                                 covars_only_raw, row.names = rownames(DF_raw), check.names = FALSE)
        ## NEW: 互動守門—TP53 至少要有 2 個 level，且設計不可對 TP53/交互作用「不可估」
        if (nlevels(droplevels(tp53_raw_fac)) < 2) {
          log_msg("[RAW|interaction] skip: TP53 has a single level among usable samples")
        } else {
          des_int_raw <- stats::model.matrix(~ pred * tp53 + . , data = df_int_raw)  # ← 你原本這行
          ne <- tryCatch(limma::nonEstimable(des_int_raw), error = function(e) NULL)
          if (!is.null(ne) && any(grepl("^(tp53|pred:tp53)", ne))) {
            log_msg("[RAW|interaction] skip: design has non-estimable TP53 terms (%s)",
                    paste(ne[grepl("^(tp53|pred:tp53)", ne)], collapse = ","))
          } else {
            ## ↓↓↓ 你原本從這裡開始的 RAW 互動分析全段（fit_int_raw, topTable, 排名與 GSEA 等）↓↓↓
            if (.RUN_LIMMA) {
              fit_int_raw <- limma::eBayes(limma::lmFit(M_raw, des_int_raw))
              coef_int <- "pred:tp53MT"
              if (coef_int %in% colnames(coef(fit_int_raw))) {
                tt   <- limma::topTable(fit_int_raw, coef = coef_int, number = nrow(M_raw), sort.by = "none")
                tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
                names(tvec) <- rownames(tt)
                # 與主線一致：補名稱、過濾非有限/重複並降冪排序
                ranks <- ._ensure_stats_names(tvec, rownames(M_raw))
                ranks <- ._finite_rank_stats(ranks, label = paste0("A1-", predictor_name, "-TP53_interaction"))
                for (grp in names(genesets_by_group)) {
                  pw <- genesets_by_group[[grp]]
                  if (is.null(pw) || !length(pw)) next
                  out_dir_int_raw <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp), "RAW", ds_id, basename(out_root), predictor_name)
                  dir.create(out_dir_int_raw, recursive = TRUE, showWarnings = FALSE)
                  
                  # ---- 互動：fgsea/camera/mroast ----
                  res_fg <- .gsea_from_ranks(
                    pathways = pw,
                    stats    = ranks,
                    minSize  = opt("minSize", 15L),
                    maxSize  = opt("maxSize", 500L),
                    gsea_eps = 1e-10,
                    label    = paste0("A1-", predictor_name, "-TP53_interaction")
                  )
                  data.table::fwrite(as.data.frame(res_fg), file.path(out_dir_int_raw, "GSEA_limma_interaction.csv"))
                  
                  idx <- limma::ids2indices(pw, names(ranks), remove.empty = TRUE)
                  if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) { res_cam <- limma::cameraPR(ranks, idx); res_cam$pathway <- rownames(res_cam) }
                  if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) data.table::fwrite(as.data.frame(res_cam), file.path(out_dir_int_raw, "GSEA_limma_interaction_camera.csv"))
                  
                  if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) res_mr <- limma::mroast(y = M_raw, index = idx, design = des_int_raw,
                                                                                                 contrast = which(colnames(des_int_raw) == coef_int),
                                                                                                 nrot = 9999, set.statistic = "mean")
                  if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) { res_mr <- as.data.frame(res_mr); res_mr$pathway <- rownames(res_mr) }
                  data.table::fwrite(res_mr, file.path(out_dir_int_raw, "GSEA_limma_interaction_mroast.csv"))
                }
              }
            }
            ## ↑↑↑ 原有邏輯到這裡結束 ↑↑↑
          }
        }
        
        ## Spearman 差異相關（用未中心化 predictor 計算相關，不影響結果）
        expr_raw_z <- as.matrix(mat0[, colnames(M_raw), drop = FALSE])
        sWT <- which(tp53_raw_fac == "WT"); sMT <- which(tp53_raw_fac == "MT")
        min_pairs <- max(10L, get0("min_pairs_spearman", ifnotfound = 10L))
        if (length(sWT) >= min_pairs && length(sMT) >= min_pairs) {
          orig_pred_raw <- pred[rownames(DF_raw)]
          rho_wt <- apply(expr_raw_z, 1, function(g) suppressWarnings(cor(g[sWT], orig_pred_raw[sWT], method = "spearman", use = "pairwise.complete.obs")))
          rho_mt <- apply(expr_raw_z, 1, function(g) suppressWarnings(cor(g[sMT], orig_pred_raw[sMT], method = "spearman", use = "pairwise.complete.obs")))
          z_wt <- atanh(pmax(pmin(rho_wt,  0.999999), -0.999999))
          z_mt <- atanh(pmax(pmin(rho_mt,  0.999999), -0.999999))
          n_wt <- sum(is.finite(orig_pred_raw[sWT])); n_mt <- sum(is.finite(orig_pred_raw[sMT]))
          if (n_wt >= 4 && n_mt >= 4) {
            if (.RUN_SPEARMAN) {
              z_diff <- (z_mt - z_wt) / sqrt(1/(n_mt-3) + 1/(n_wt-3))
              z_diff <- z_diff[is.finite(z_diff)]
              z_diff <- sort(z_diff, decreasing = TRUE)
              for (grp in names(genesets_by_group)) {
                pw <- genesets_by_group[[grp]]
                if (is.null(pw) || !length(pw)) next
                out_dir_int_raw <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp), "RAW", ds_id, basename(out_root), predictor_name)
                dir.create(out_dir_int_raw, recursive = TRUE, showWarnings = FALSE)
                
                idx2 <- limma::ids2indices(pw, names(z_diff), remove.empty = TRUE)
                res_fg2 <- gsea_from_diffcorr(pw, z_diff, label = "diff-corr")
                data.table::fwrite(as.data.frame(res_fg2), file.path(out_dir_int_raw, "GSEA_spearman_diffcorr.csv"))
                if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) { res_cam2 <- limma::cameraPR(z_diff, idx2); res_cam2$pathway <- rownames(res_cam2) }
                if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) data.table::fwrite(as.data.frame(res_cam2), file.path(out_dir_int_raw, "GSEA_spearman_diffcorr_camera.csv"))
              }
            }
          }
        }
      }, silent = TRUE)
    
    ## ---------- BatchAdj branch（interaction + diffcorr）----------
    if (!("BatchAdj" %in% getOption("csn.run_passes", c("BatchAdj")))) { log_msg("[BatchAdj|interaction] skip: pass 未選") } else
      try({
        tp53_ba_num <- suppressWarnings(as.numeric(tp53_num_all[rownames(DF_ba)]))
        tp53_ba_fac <- factor(ifelse(tp53_ba_num == 1, "MT", "WT"), levels = c("WT","MT"))
        covars_only_ba <- DF_ba[, setdiff(colnames(DF_ba),
                                          c("predictor","TP53_mutant","TP53_status","TP53")), drop = FALSE]
        
        df_int_ba <- data.frame(pred = DF_ba$predictor, tp53 = tp53_ba_fac,
                                covars_only_ba, row.names = rownames(DF_ba), check.names = FALSE)
        des_int_ba <- stats::model.matrix(~ pred * tp53 + ., data = df_int_ba)
        if (.RUN_LIMMA) {
          fit_int_ba <- limma::eBayes(limma::lmFit(M_ba, des_int_ba))
          coef_int <- "pred:tp53MT"
          if (coef_int %in% colnames(coef(fit_int_ba))) {
            tt   <- limma::topTable(fit_int_ba, coef = coef_int, number = nrow(M_ba), sort.by = "none")
            tvec <- if ("t" %in% names(tt)) tt$t else tt$logFC
            names(tvec) <- rownames(tt)
            # 與主線一致：補名稱、過濾非有限/重複並降冪排序
            ranks <- ._ensure_stats_names(tvec, rownames(M_ba))
            ranks <- ._finite_rank_stats(ranks, label = paste0("A2-", predictor_name, "-TP53_interaction"))
            for (grp in names(genesets_by_group)) {
              pw <- genesets_by_group[[grp]]
              if (is.null(pw) || !length(pw)) next
              out_dir_int_ba <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp), "BatchAdj", ds_id, basename(out_root), predictor_name)
              dir.create(out_dir_int_ba, recursive = TRUE, showWarnings = FALSE)
              
              # ---- 互動：fgsea/camera/mroast ----
              res_fg <- .gsea_from_ranks(
                pathways = pw,
                stats    = ranks,
                minSize  = opt("minSize", 15L),
                maxSize  = opt("maxSize", 500L),
                gsea_eps = 1e-10,
                label    = paste0("A2-", predictor_name, "-TP53_interaction")
              )
              data.table::fwrite(as.data.frame(res_fg), file.path(out_dir_int_ba, "GSEA_limma_interaction.csv"))
              
              idx <- limma::ids2indices(pw, names(ranks), remove.empty = TRUE)
              if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) { res_cam <- limma::cameraPR(ranks, idx); res_cam$pathway <- rownames(res_cam) }
              if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) data.table::fwrite(as.data.frame(res_cam), file.path(out_dir_int_ba, "GSEA_limma_interaction_camera.csv"))
              
              if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) res_mr <- limma::mroast(y = M_ba, index = idx, design = des_int_ba,
                                                                                             contrast = which(colnames(des_int_ba) == coef_int),
                                                                                             nrot = 9999, set.statistic = "mean")
              if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) { res_mr <- as.data.frame(res_mr); res_mr$pathway <- rownames(res_mr) }
              data.table::fwrite(res_mr, file.path(out_dir_int_ba, "GSEA_limma_interaction_mroast.csv"))
            }
          }
        }
        ## Spearman 差異相關（未中心化 predictor）
        sWT <- which(tp53_ba_fac == "WT"); sMT <- which(tp53_ba_fac == "MT")
        min_pairs <- max(10L, get0("min_pairs_spearman", ifnotfound = 10L))
        if (length(sWT) >= min_pairs && length(sMT) >= min_pairs) {
          orig_pred_ba <- pred[rownames(DF_ba)]
          rho_wt <- apply(M_ba, 1, function(g) suppressWarnings(cor(g[sWT], orig_pred_ba[sWT], method = "spearman", use = "pairwise.complete.obs")))
          rho_mt <- apply(M_ba, 1, function(g) suppressWarnings(cor(g[sMT], orig_pred_ba[sMT], method = "spearman", use = "pairwise.complete.obs")))
          z_wt <- atanh(pmax(pmin(rho_wt,  0.999999), -0.999999))
          z_mt <- atanh(pmax(pmin(rho_mt,  0.999999), -0.999999))
          n_wt <- sum(is.finite(orig_pred_ba[sWT])); n_mt <- sum(is.finite(orig_pred_ba[sMT]))
          if (n_wt >= 4 && n_mt >= 4) {
            if (.RUN_SPEARMAN) {
              z_diff <- (z_mt - z_wt) / sqrt(1/(n_mt-3) + 1/(n_wt-3))
              z_diff <- z_diff[is.finite(z_diff)]
              z_diff <- sort(z_diff, decreasing = TRUE)
              for (grp in names(genesets_by_group)) {
                pw <- genesets_by_group[[grp]]
                if (is.null(pw) || !length(pw)) next
                out_dir_int_ba <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp), "BatchAdj", ds_id, basename(out_root), predictor_name)
                dir.create(out_dir_int_ba, recursive = TRUE, showWarnings = FALSE)
                
                idx2 <- limma::ids2indices(pw, names(z_diff), remove.empty = TRUE)
                res_fg2 <- gsea_from_diffcorr(pw, z_diff, label = "diff-corr")
                data.table::fwrite(as.data.frame(res_fg2), file.path(out_dir_int_ba, "GSEA_spearman_diffcorr.csv"))
                if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) { res_cam2 <- limma::cameraPR(z_diff, idx2); res_cam2$pathway <- rownames(res_cam2) }
                if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) data.table::fwrite(as.data.frame(res_cam2), file.path(out_dir_int_ba, "GSEA_spearman_diffcorr_camera.csv"))
              }
            }
          }
        }
      }, silent = TRUE)
  }
  
  ## -------- 小工具 --------
  coef_name <- function(des) if ("predictor" %in% colnames(des)) "predictor" else colnames(des)[ncol(des)]
  rank_finite <- function(v, exclude = NULL, tag = NULL){
    keep <- is.finite(v); if (any(!keep)) logf("[gsea-%s] drop non-finite: %d", ifelse(is.null(tag),"NA",tag), sum(!keep))
    v <- v[keep]
    if (!is.null(exclude) && length(exclude)) v <- v[ setdiff(names(v), exclude) ]
    if (!length(v)) return(NULL)
    v
  }
  gsea_from <- function(pw, stats){
    if (is.null(pw) || !length(pw) || is.null(stats) || !length(stats)) return(NULL)
    # 讀全域門檻
    minSize <- opt("minSize", 15L); maxSize <- opt("maxSize", 500L)
    # 顯式交集 + size 篩
    orig_n <- length(pw)
    pw_use <- ._intersect_and_filter_pathways(pw, names(stats), minSize = minSize, maxSize = maxSize)
    message(sprintf("[gsea] (gsea_from) |H| 原=%d → 用=%d", orig_n, length(pw_use)))
    if (!length(pw_use)) return(NULL)
    
    set.seed(1L)
    res <- tryCatch(
      {
        suppressWarnings(fgsea::fgseaMultilevel(
          pathways = pw_use, stats = stats,
          minSize  = minSize, maxSize = maxSize,
          eps = 1e-10
        ))
      },
      error = function(e) {
        message(sprintf("[gsea] Multilevel 失敗：%s → 改用 fgseaSimple", conditionMessage(e)))
        suppressWarnings(fgsea::fgseaSimple(
          pathways = pw_use, stats = stats,
          nperm = 10000, minSize = minSize, maxSize = maxSize
        ))
      }
    )
    res
  }
  
  ## -------- 逐 group：A1（RAW）/ A2（BatchAdj）----------
  if (.RUN_LIMMA) {
    for (grp_name in names(genesets_by_group)) {
      pw <- genesets_by_group[[grp_name]]
      
      ## [NEW] collection-first root: <collection>/<dataset>/<stratum>/<RAW|BatchAdj>/<predictor>
      out_root_coll <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))
      if ("RAW" %in% getOption("csn.run_passes", c("BatchAdj"))) {
        out_dir_A1 <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "RAW", ds_id, basename(out_root), predictor_name)
        dir.create(out_dir_A1, recursive = TRUE, showWarnings = FALSE)
        fit1 <- limma::eBayes(limma::lmFit(M_raw, des_raw))
        t1 <- fit1$t[, coef_name(des_raw)]
        t1 <- ._ensure_stats_names(t1, rownames(M_raw))
        if (exists("exclude_genes") && length(exclude_genes)) {
          t1 <- t1[setdiff(names(t1), exclude_genes)]
        }
        t1 <- ._finite_rank_stats(t1, label = paste0("A1-", predictor_name))
        if (!is.null(t1)) {
          res1 <- gsea_from(pw, t1)
          if (!is.null(res1) && nrow(res1)) data.table::fwrite(res1, file.path(out_dir_A1, "GSEA_limma_t_cont.csv"))
        }
        # camera / mroast（敏感度）
        if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) try(run_camera_from_ranks_save(t1, pw, out_dir_A1, outfile = "GSEA_camera.csv"), silent = TRUE)
        if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) try(run_mroast_save(M_raw, des_raw, coef_name(des_raw), pw, out_dir_A1,
                                                                                   outfile = "GSEA_mroast.csv", universe = names(t1)), silent = TRUE)
      }
      
      ## A2
      if ("BatchAdj" %in% getOption("csn.run_passes", c("BatchAdj"))) {
        out_dir_A2 <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root), predictor_name)
        dir.create(out_dir_A2, recursive = TRUE, showWarnings = FALSE)
        fit2 <- limma::eBayes(limma::lmFit(M_ba, des_ba))
        t2 <- fit2$t[, coef_name(des_ba)]
        t2 <- ._ensure_stats_names(t2, rownames(M_ba))
        if (exists("exclude_genes") && length(exclude_genes)) {
          t2 <- t2[setdiff(names(t2), exclude_genes)]
        }
        t2 <- ._finite_rank_stats(t2, label = paste0("A2-", predictor_name))
        if (!is.null(t2)) {
          res2 <- gsea_from(pw, t2)
          if (!is.null(res2) && nrow(res2)) data.table::fwrite(res2, file.path(out_dir_A2, "GSEA_limma_t_cont.csv"))
        }
        if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) try(run_camera_from_ranks_save(t2, pw, out_dir_A2, outfile = "GSEA_camera.csv"), silent = TRUE)
        if (isTRUE(getOption("csn.run_camera_mroast", FALSE))) try(run_mroast_save(M_ba, des_ba, coef_name(des_ba), pw, out_dir_A2,
                                                                                   outfile = "GSEA_mroast.csv", universe = names(t2)), silent = TRUE)
      }
    }
  }
}





## ===== CSN subunits 覆蓋率 & CSN_SCORE（PC1）可行性審核（穩健版）=====
audit_csn_score_feasibility <- function(ds_id, stratum, mat0, present_sub,
                                        min_members = 5L,
                                        pca_min_samples = 10L,   # 與 build_csn_score 內門檻一致
                                        min_per_group = 8L,      # 你的 limma 門檻
                                        out_dir = file.path("run_info","csn_score_audit")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  n_samples <- ncol(mat0)
  if (!is.matrix(mat0) || is.null(rownames(mat0)) || is.null(colnames(mat0))) {
    log_msg("[CSN-audit] %s | %s：mat0 結構不完整，略過", ds_id, stratum)
    return(invisible(NULL))
  }
  present_sub <- intersect(present_sub, rownames(mat0))
  if (!length(present_sub) || n_samples == 0) {
    log_msg("[CSN-audit] %s | %s：沒有可用的 CSN subunit 或樣本為 0，略過", ds_id, stratum)
    return(invisible(NULL))
  }
  
  ## 1) 每個 subunit 的覆蓋率
  sub_cov <- vapply(present_sub, function(g) mean(is.finite(mat0[g, ])) * 100, numeric(1))
  sub_tbl <- data.frame(
    dataset  = ds_id, 
    stratum  = stratum,
    subunit  = present_sub,
    nonNA_pct = round(sub_cov, 1),
    nonNA_n   = vapply(present_sub, function(g) sum(is.finite(mat0[g, ])), integer(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  cov_min <- if (length(sub_cov)) round(min(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_med <- if (length(sub_cov)) round(stats::median(sub_cov, na.rm = TRUE), 1) else NA_real_
  cov_max <- if (length(sub_cov)) round(max(sub_cov, na.rm = TRUE), 1) else NA_real_
  
  ## 2) 每個樣本有多少 subunits 非 NA；滿足 >= min_members 的樣本數
  sample_counts <- colSums(is.finite(mat0[present_sub, , drop = FALSE]))
  enough <- sample_counts >= min_members
  n_enough <- sum(enough)
  
  ## 3) 嘗試計算 CSN_SCORE（PC1），並記錄可行性
  csn_score <- build_csn_score(mat0, subunits = present_sub,
                               combine_7AB = TRUE, min_members = min_members)
  csn_nonNA   <- sum(is.finite(csn_score))
  csn_can_pca <- (n_enough >= pca_min_samples) && (csn_nonNA >= pca_min_samples)
  
  ## 額外：PC1 解釋變異（可行時）
  pc1_var_pct <- NA_real_
  if (isTRUE(csn_can_pca)) {
    ok_sam <- names(sample_counts)[sample_counts >= min_members]
    get_z <- function(v) {
      v <- as.numeric(v)
      mu <- mean(v[is.finite(v)], na.rm = TRUE)
      sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[!is.finite(v)] <- mu
      (v - mu) / sdv
    }
    Z <- do.call(rbind, lapply(present_sub, function(g) get_z(mat0[g, ok_sam, drop = FALSE])))
    rownames(Z) <- present_sub
    if (all(c("COPS7A","COPS7B") %in% rownames(Z))) {
      Z7 <- colMeans(Z[c("COPS7A","COPS7B"), , drop = FALSE], na.rm = TRUE)
      Z  <- rbind(Z[setdiff(rownames(Z), c("COPS7A","COPS7B")), , drop = FALSE],
                  "COPS7*" = Z7)
    }
    pc <- tryCatch(stats::prcomp(t(Z), center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (!is.null(pc)) {
      ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
      pc1_var_pct <- round(ve[1], 1)
    }
  }
  
  ## 4) 寫出檔案
  # 4a) 總表（append）
  sum_row <- data.frame(
    dataset = ds_id, stratum = stratum,
    n_samples = n_samples,
    n_subunits_present = length(present_sub),
    subunit_cov_min = cov_min, subunit_cov_median = cov_med, subunit_cov_max = cov_max,
    min_members = min_members,
    samples_with_ge_min_members = n_enough,
    pca_min_samples = pca_min_samples,
    csn_score_nonNA = csn_nonNA,
    csn_pc1_feasible = csn_can_pca,
    csn_pc1_var_pct = pc1_var_pct,
    can_form_high_low_groups = (csn_nonNA >= 2 * min_per_group),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  fp_sum <- file.path(out_dir, "csn_score_feasibility_summary.csv")
  data.table::fwrite(sum_row, fp_sum, append = file.exists(fp_sum))
  
  # 4b) subunit 覆蓋率（每分層一個檔）
  tag <- paste(ds_id, stratum, sep = "_")
  fp_sub <- file.path(out_dir, sprintf("%s_subunit_coverage.csv", tag))
  data.table::fwrite(sub_tbl, fp_sub)
  
  # 4c) 每樣本非 NA subunit 數（每分層一個檔）
  sample_tbl <- data.frame(
    dataset = ds_id, stratum = stratum,
    sample_id = colnames(mat0),
    nonNA_subunits_n = as.integer(sample_counts),
    ge_min_members = as.logical(enough),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  fp_sam <- file.path(out_dir, sprintf("%s_sample_subunit_counts.csv", tag))
  data.table::fwrite(sample_tbl, fp_sam)
  
  log_msg("[CSN-audit] %s | %s：subunits=%d；min/median/max 覆蓋率=%.1f/%.1f/%.1f%%；eligible_samples=%d；CSN_SCORE nonNA=%d；PC1可行=%s（PC1%%=%.1f）",
          ds_id, stratum, length(present_sub), cov_min %||% NaN, cov_med %||% NaN, cov_max %||% NaN,
          n_enough, csn_nonNA, as.character(csn_can_pca), pc1_var_pct %||% NaN)
  
  invisible(list(summary = sum_row, per_subunit = sub_tbl, per_sample = sample_tbl))
}



## ===== 針對「指定樣本集合」執行一次完整 GSEA 的小函式 =====
## =========================================================
## 每個 stratum：同時產 RAW 與 BatchAdj（新增版本）
## =========================================================
## ===== run_one_stratum：統一 Spearman 檔名、寫 sample_counts、批次備註補強（2,3,7）=====
## ===== 針對「指定樣本集合」執行一次完整 GSEA 的小函式（修正版）=====
## ===== 針對「指定樣本集合」執行一次完整 GSEA 的小函式（修正版；BatchAdj 也可帶入 age_missing/age_z_imputed）=====
run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group){
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("  -- stratum：%s | N(sample_keep)=%d", basename(out_root), length(sample_keep))
  
  ## 1) 子集樣本 + 規模檢查 + log 標度
  keep <- intersect(colnames(mat0_full), sample_keep)
  if (length(keep) < 4) { log_msg("  [跳過] 樣本過少：%d", length(keep)); return(invisible(NULL)) }
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  
  ## 2) limma 用的插補+過濾矩陣（Spearman 仍用 raw mat0） + CSN_SCORE 可行性稽核
  mat  <- impute_and_filter(mat0, min_frac = min_frac_complete)
  present_sub <- intersect(csn_subunits, rownames(mat0))
  if (!length(present_sub)) { log_msg("  [跳過] 此分層沒有任何 CSN subunit"); return(invisible(NULL)) }
  audit_csn_score_feasibility_safe(
    ds_id    = ds_id,
    stratum  = basename(out_root),
    mat0     = mat0,
    present_sub = present_sub,
    min_members     = 5L,
    pca_min_samples = 10L,
    min_per_group   = min_per_group,
    out_dir = file.path("run_info","csn_score_audit")
  )
  
  ## === 在 ALL 分析：準備 TP53 狀態（0/1；mutant=1）===
  is_ALL <- identical(basename(out_root), "ALL")
  tp53_num_all <- NULL
  if (is_ALL) {
    tp53_status_all <- get_tp53_status(ds_dir, colnames(mat0))
    tp53_num_all <- as.numeric(tp53_status_all == "TP53_mutant")
    names(tp53_num_all) <- colnames(mat0)
  }
  
  ## 3) 共變項 / 批次（對齊樣本）—— 建立 limma 用共變項（可選帶入 missing-indicator）
  bi_all     <- get_batch_factor(ds_dir, colnames(mat0))
  batch_all  <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
  sa_all     <- get_sex_age_covariates(ds_dir, colnames(mat0))  # data.frame(sex, age [, age_missing, age_z_imputed])
  sa_all_limma <- sa_all
  if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing","age_z_imputed") %in% colnames(sa_all))) {
    keep_cols <- intersect(c("sex","age","age_missing","age_z_imputed"), colnames(sa_all))
    sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
  }
  
  ## 4) 逐「預測向量」執行分析：subunit 原值、CSN_SCORE、RESIDUAL_<SU>
  {
    present_sub <- intersect(csn_subunits, rownames(mat0))
    if (!length(present_sub)) { log_msg("  [跳過] 此分層沒有任何 CSN subunit"); return(invisible(NULL)) }
    
    # 再取一次（確保與上方一致）
    is_ALL <- identical(basename(out_root), "ALL")
    tp53_num_all <- if (is_ALL) {
      status_all <- get_tp53_status(ds_dir, colnames(mat0))
      setNames(as.numeric(status_all == "TP53_mutant"), colnames(mat0))
    } else NULL
    
    bi_all     <- get_batch_factor(ds_dir, colnames(mat0))
    batch_all  <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
    purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
    sa_all     <- get_sex_age_covariates(ds_dir, colnames(mat0))
    sa_all_limma <- sa_all
    if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing","age_z_imputed") %in% colnames(sa_all))) {
      keep_cols <- intersect(c("sex","age","age_missing","age_z_imputed"), colnames(sa_all))
      sa_all_limma <- sa_all[, keep_cols, drop = FALSE]
    }
    
    ## 4a) 原本：逐 subunit（排除自身 gene）
    for (su in present_sub) {
      run_predictor_analyses(
        predictor_name = su,
        predictor_vec  = mat0[su, ],
        exclude_genes  = su,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    }
    
    ## 4b) CSN complex score（PC1；方向校正）；排名時排除所有 CSN 成員
    csn_score <- build_csn_score_safe(
      mat0, subunits = present_sub, combine_7AB = TRUE,
      min_members = 5L, pca_min_samples = 10L
    )
    
    if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
      run_predictor_analyses(
        predictor_name = "CSN_SCORE",
        predictor_vec  = csn_score,
        exclude_genes  = present_sub,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    } else {
      log_msg("  [CSN_SCORE] 非 NA 樣本不足，略過")
    }
    
    ## 4c) RESIDUAL_<SU>（加護欄：CSN score 非 NA 樣本數不足 → 整批略過）
    min_n_resid <- min_per_group
    csn_nonNA   <- sum(is.finite(csn_score))
    if (csn_nonNA < min_n_resid) {
      log_msg("  [RESIDUAL] CSN score 非 NA 樣本不足（%d < %d），整批略過 residual_*", csn_nonNA, min_n_resid)
    } else {
      # 只在需要殘差時才建 base covars（省計算）
      base_covars_all <- data.frame(
        purity = as.numeric(purity_all[colnames(mat0)]),
        sex    = as.numeric(sa_all[colnames(mat0), "sex"]),
        age    = as.numeric(sa_all[colnames(mat0), "age"]),
        row.names = colnames(mat0), check.names = FALSE
      )
      # （可選）missing-indicator 一致化：在殘差化中帶入
      if (isTRUE(USE_AGE_MISSING_INDICATOR) && all(c("age_missing","age_z_imputed") %in% colnames(sa_all))) {
        base_covars_all$age_missing   <- as.numeric(sa_all[colnames(mat0), "age_missing"])
        base_covars_all$age_z_imputed <- as.numeric(sa_all[colnames(mat0), "age_z_imputed"])
      }
      if (is_ALL && !is.null(tp53_num_all)) {
        base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
      }
      
      for (su in present_sub) {
        res_sub <- residualize_vector(
          y = mat0[su, ],
          csn_score = csn_score,
          batch = batch_all,
          covars = base_covars_all,
          min_n = min_n_resid
        )
        if (sum(is.finite(res_sub)) < (2 * min_per_group)) {
          log_msg("  [RESIDUAL_%s] 非 NA 樣本不足，略過", su); next
        }
        run_predictor_analyses(
          predictor_name = paste0("RESIDUAL_", su),
          predictor_vec  = res_sub,
          exclude_genes  = su,
          ds_id = ds_id, ds_dir = ds_dir,
          mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
          out_root = out_root,
          genesets_by_group = genesets_by_group,
          batch_all = batch_all, purity_all = purity_all, sa_all = sa_all_limma,
          tp53_num_all = tp53_num_all, is_ALL = is_ALL
        )
      }
    }
  }
  
  ## 5) 各版本各自彙整
  present_sub <- intersect(csn_subunits, rownames(mat0))
  sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))
  
  # RAW：limma_t_cont + interaction（只彙整現有 csv，不重跑）
  for (grp_name in names(genesets_by_group)) {
    out_root_coll <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))
    ver_root <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "RAW", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = c("GSEA_limma_t_cont", "GSEA_limma_interaction")
    )
  }
  
  
  # BatchAdj：limma_t_cont + spearman + interaction（只彙整現有 csv，不重跑）
  for (grp_name in names(genesets_by_group)) {
    out_root_coll <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), ds_id, basename(out_root))
    ver_root <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), "BatchAdj", ds_id, basename(out_root))
    summarize_all_groups(
      out_root = ver_root,
      csn_subunits = sum_units,
      genesets_by_group = setNames(list(genesets_by_group[[grp_name]]), grp_name),
      stat_tags = {
        st <- c("GSEA_limma_t_cont", "GSEA_spearman", "GSEA_limma_interaction")
        if (isTRUE(RUN_HILO_SENSITIVITY)) st <- c(st, "GSEA_limma_t_hilo")
        st
      }
    )
  }
  
  invisible(NULL)
}


## 檢視目前取到哪些 PLEX（原始 vs 清理後）
inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy   = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
  # 讀你實際有用到的樣本（以 protein 矩陣為準）
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(mat0)
  
  # 讀臨床表並對齊樣本
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  stopifnot(file.exists(meta_fp))
  meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  id_cols <- intersect(c("SAMPLE_ID","sample_id","Sample_ID","Sample","sample"), names(meta))
  stopifnot(length(id_cols) > 0)
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  
  # 取原始與清理後的 PLEX
  raw <- as.character(meta[[col]])
  clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)
  
  # 摘要
  cat("\n==== ", basename(ds_dir), " | 欄位：", col, " ====\n", sep = "")
  cat("[原始 PLEX levels]：\n")
  print(sort(table(raw), decreasing = TRUE))
  cat("\n[含 '|' 的原始值（前幾個）]：\n")
  print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))
  
  cat("\n[清理後 PLEX levels]（pipe_policy = ", pipe_policy, 
      ", min_per_level = ", min_per_level, ")：\n", sep = "")
  print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))
  
  # 輸出對照表（前 10 列示範）
  df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
  cat("\n[樣本對照表（前 10 列）]\n")
  print(utils::head(df_map, 10))
  
  invisible(list(raw_counts = sort(table(raw), decreasing = TRUE),
                 clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
                 map = df_map))
}

## ---- 全域開關：永遠加 tech covariates（即使有 batch）----
ALWAYS_ADD_TECH_COVARS <- FALSE

## ===== Missingness sensitivity（AGE）=====
USE_AGE_MISSING_INDICATOR <- FALSE  # 主分析預設：不用 missing-indicator，僅保留 NA

## 安全 z-score：只以有限值估計 mean/sd，保留 NA，不做均值補
.z_no_impute <- function(v){
  v <- suppressWarnings(as.numeric(v))
  fin <- is.finite(v)
  mu <- mean(v[fin], na.rm = TRUE)
  sdv <- stats::sd(v[fin], na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
  out <- (v - mu) / sdv
  out[!fin] <- NA_real_
  out
}

## 工具：欄名正規化、z-score、0~1
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+","_", x))
.zscore <- function(v){
  v <- as.numeric(v); mu <- mean(v[is.finite(v)], na.rm=TRUE)
  sdv <- stats::sd(v[is.finite(v)], na.rm=TRUE); if (!is.finite(sdv) || sdv==0) sdv <- 1
  v[!is.finite(v)] <- mu; (v - mu)/sdv
}
.to01 <- function(v){
  v <- suppressWarnings(as.numeric(v))
  if (sum(is.finite(v) & v>1, na.rm=TRUE) > sum(is.finite(v) & v<=1, na.rm=TRUE)) v <- v/100
  pmin(pmax(v,0),1)
}

## PAAD 用：以分號切的值取中位數
.median_from_semicolon <- function(x_chr){
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(vv)]; if (!length(vv)) return(NA_real_)
  stats::median(vv)
}


## 3) 依資料集抓 purity（named numeric，names=sample_ids，0~1）
get_purity_covariate <- function(ds_id, ds_dir, sample_ids){
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  samp <- if (file.exists(samp_fp)) suppressMessages(readr::read_tsv(samp_fp, show_col_types=FALSE, comment="#")) else NULL
  if (!is.null(samp)) { samp <- as.data.frame(samp); names(samp) <- .norm_names(names(samp)) }
  purity <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  
  if (ds_id == "brca_cptac_2020") {
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    if (file.exists(pat_fp) && !is.null(samp)) {
      pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types=FALSE, comment="#")) |> as.data.frame()
      names(pat) <- .norm_names(names(pat))
      sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
      pid_in_samp <- intersect(c("PATIENT_ID","PATIENT"), names(samp))[1]
      pid_in_pat  <- intersect(c("PATIENT_ID","PATIENT"), names(pat))[1]
      if (!is.na(sid) && !is.na(pid_in_samp) && !is.na(pid_in_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat)) {
        map_pt <- setNames(as.character(samp[[pid_in_samp]]), samp[[sid]])
        pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_in_pat]]))
        purity[] <- unname(pt_pur[ map_pt[sample_ids] ])
      }
    }
  } else if (ds_id == "luad_cptac_2020") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
      if (!is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp))
        purity[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][ match(sample_ids, samp[[sid]]) ])
    }
  } else if (ds_id == "lusc_cptac_2021") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
      if (!is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp))
        purity[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][ match(sample_ids, samp[[sid]]) ])
    }
  } else if (ds_id == "paad_cptac_2021") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
      if (!is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
        medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
        purity[] <- .to01(medv[ match(sample_ids, samp[[sid]]) ])
      }
    }
  } else if (ds_id == "ucec_cptac_2020") {
    if (!is.null(samp)) {
      sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
      if (!is.na(sid)) {
        pc <- suppressWarnings(as.numeric(samp[["PURITY_CANCER"]]))
        pi <- suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]]))
        ps <- suppressWarnings(as.numeric(samp[["PURITY_STROMA"]]))
        idx <- match(sample_ids, samp[[sid]])
        purity_calc <- ifelse(is.finite(pc[idx]), pc[idx], 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0))
        purity[] <- pmin(pmax(purity_calc,0),1)
      }
    }
  }
  purity
}


## ========================
## 取 sex/age（sample 對齊；sex: 0/1；age: z-score）
## ========================

if (!exists(".norm_names", inherits = FALSE)) {
  .norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+","_", x))
}
if (!exists(".zscore", inherits = FALSE)) {
  .zscore <- function(v){
    v <- suppressWarnings(as.numeric(v))
    mu <- mean(v[is.finite(v)], na.rm = TRUE)
    sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
    v[!is.finite(v)] <- mu
    (v - mu)/sdv
  }
}

## 主流程用：嚴格以 patient 檔的 SEX 與 AGE 產生 covariates（sex: 0/1；age: z-score）
get_sex_age_covariates <- function(ds_dir, sample_ids) {
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp  <- file.path(ds_dir, "data_clinical_patient.txt")
  
  if (!file.exists(samp_fp) || !file.exists(pat_fp)) {
    log_msg("[covars] 缺少 sample/patient 檔，sex/age 以 NA 代替")
    out <- cbind(sex = rep(NA_real_, length(sample_ids)),
                 age = rep(NA_real_, length(sample_ids)))
    rownames(out) <- sample_ids
    return(out)
  }
  
  .NN <- function(x) toupper(gsub("[^A-Za-z0-9]+","_", x))
  .z_no_impute <- function(v){
    v <- suppressWarnings(as.numeric(v))
    fin <- is.finite(v)
    mu <- mean(v[fin], na.rm = TRUE)
    sdv <- stats::sd(v[fin], na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out <- (v - mu) / sdv
    out[!fin] <- NA_real_
    out
  }
  # 若外部未定義 .zscore（均值補後 z），提供後備版本以維持一致邏輯
  if (!exists(".zscore", mode = "function")) {
    .zscore <- function(v){
      v <- suppressWarnings(as.numeric(v))
      mu <- mean(v[is.finite(v)], na.rm = TRUE)
      sdv <- stats::sd(v[is.finite(v)], na.rm = TRUE); if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[!is.finite(v)] <- mu
      (v - mu)/sdv
    }
  }
  
  samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  pat  <- suppressMessages(readr::read_tsv(pat_fp,  show_col_types = FALSE, comment = "#")) |> as.data.frame()
  names(samp) <- .NN(names(samp)); names(pat) <- .NN(names(pat))
  
  sid      <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
  pid_samp <- intersect(c("PATIENT_ID","PATIENT"), names(samp))[1]
  pid_pat  <- intersect(c("PATIENT_ID","PATIENT"), names(pat))[1]
  if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) {
    stop("[get_sex_age_covariates] 無法建立 sample↔patient 對應（缺 SAMPLE_ID/PATIENT_ID）")
  }
  map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])
  
  ## SEX: Male=1, Female=0
  sex_col <- intersect(c("SEX","GENDER"), names(pat))[1]
  sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  if (!is.na(sex_col)) {
    raw <- toupper(as.character(pat[[sex_col]]))
    val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
    names(val) <- as.character(pat[[pid_pat]])
    sex[] <- unname(val[ map_pt[sample_ids] ])
  }
  
  ## AGE: 主分析不插補；敏感度可選 missing-indicator + 均值補後 z
  age_col <- intersect(c("AGE","AGE_AT_DIAGNOSIS","AGE_AT_INDEX","AGE_YEARS"), names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_missing <- age_z_imputed <- NULL
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[ map_pt[sample_ids] ])
    age[] <- .z_no_impute(age_raw)  # NA 保留（主分析）
    if (!exists("USE_AGE_MISSING_INDICATOR", inherits = FALSE)) USE_AGE_MISSING_INDICATOR <- FALSE
    if (isTRUE(USE_AGE_MISSING_INDICATOR)) {
      age_missing   <- as.numeric(is.na(age_raw))
      age_z_imputed <- .zscore(age_raw)  # 均值補後 z（敏感度）
    }
  }
  
  cov_sex <- mean(is.finite(sex)) * 100
  cov_age <- mean(is.finite(age)) * 100
  if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
    cov_age_imp <- mean(is.finite(age_z_imputed)) * 100
    cov_age_mis <- mean(is.finite(age_missing)) * 100
    log_msg("    covariates coverage：sex %.1f%%, age(NA-as-NA) %.1f%%, age_z_imputed %.1f%%, age_missing %.1f%%",
            cov_sex, cov_age, cov_age_imp, cov_age_mis)
  } else {
    log_msg("    covariates coverage：sex %.1f%%, age(NA-as-NA) %.1f%%", cov_sex, cov_age)
  }
  
  out <- data.frame(sex = as.numeric(sex), age = as.numeric(age),
                    row.names = sample_ids, check.names = FALSE)
  out$age_z <- out$age
  if (isTRUE(get0("USE_AGE_MISSING_INDICATOR", ifnotfound = FALSE))) {
    out$age_missing   <- as.numeric(age_missing)
    out$age_z_imputed <- as.numeric(age_z_imputed)
  }
  out
}


select_covars_safely <- function(
    df,                    # 候選共變數；rownames=sample_order
    sample_order,
    label = "covars",
    y = NULL,              # 連續 predictor；可為 NULL
    min_cov_named = c(purity=0.60, sex=0.80, age=0.80),
    max_abs_cor   = 0.30,
    min_pairs     = 20L
){
  logf <- function(...) if (exists("log_msg", mode="function")) try(log_msg(...), silent=TRUE)
  ## NEW: 逐欄 drop 原因收集器
  drop_log <- list()
  .append_drop <- function(step, col, reason, metric = NA_real_, threshold = NA_real_, extra = NA_character_) {
    drop_log[[length(drop_log) + 1L]] <<- data.frame(
      step       = step,
      covariate  = col,
      reason     = reason,
      metric     = metric,
      threshold  = threshold,
      extra      = extra,
      stringsAsFactors = FALSE
    )
  }
  if (is.null(df) || !nrow(df)) return(NULL)
  
  
  # 1) 對齊樣本順序
  if (is.null(rownames(df))) { logf("  [covars-%s] df 無 rownames → skip", label); return(NULL) }
  so <- as.character(sample_order)
  logf("  [covars-%s] before align: C_dim=%d x %d; has_rownames=%s", label, NROW(df), NCOL(df), !is.null(rownames(df)))
  df <- df[so, , drop = FALSE]
  
  ## [POLICY] Never include TP53 as a covariate (ALL / MT / WT strata)
  tp_cols_idx <- grep("^TP53($|_|)", colnames(df), ignore.case = FALSE)
  if (length(tp_cols_idx) > 0L) {
    tp_cols <- colnames(df)[tp_cols_idx]
    for (cc in tp_cols) {
      .append_drop(step = "policy", col = cc, reason = "exclude_TP53_as_covariate",
                   metric = NA_real_, threshold = NA_real_, extra = "global policy")
    }
    df <- df[, setdiff(colnames(df), tp_cols), drop = FALSE]
    logf(sprintf("  [covars-%s] policy: drop TP53 columns from covariates → %s", label, paste(tp_cols, collapse=",")))
  }
  
  logf("  [covars-%s] after  align: C_dim=%d x %d", label, NROW(df), NCOL(df))
  
  
  # 2) 一律用 data.frame；保留因子
  df <- as.data.frame(df, stringsAsFactors = TRUE, check.names = FALSE)
  
  # 3) coverage：未指名欄位使用全域最小門檻（安全索引，不用 [[cn]]）
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  global_min <- if (length(min_cov_named)) min(min_cov_named, na.rm = TRUE) else 0
  keep_cov <- rep(TRUE, ncol(df)); names(keep_cov) <- colnames(df)
  
  for (cn in colnames(df)) {
    v <- df[[cn]]
    # 統一 coverage 定義：數值用 is.finite；非數值用 !is.na
    cover <- if (is.numeric(v)) mean(is.finite(v)) else mean(!is.na(v))
    
    thr <- if (!is.null(names(min_cov_named)) && (cn %in% names(min_cov_named))) {
      min_cov_named[[cn]]            # ★ 安全：只有存在名稱時才 [[cn]]
    } else {
      global_min
    }
    
    if (is.na(cover) || cover < thr) {
      logf("  [covars-%s] drop %s (coverage=%s < %s)", label, cn,
           sprintf("%.0f%%", 100*ifelse(is.na(cover), 0, cover)),
           sprintf("%.0f%%", 100*thr))
      keep_cov[cn] <- FALSE
      if (is.na(cover) || cover < thr) {
        logf("  [covars-%s] drop %s (coverage=%s < %s)", label, cn,
             sprintf("%.0f%%", 100*ifelse(is.na(cover), 0, cover)),
             sprintf("%.0f%%", 100*thr))
        keep_cov[cn] <- FALSE
        ## NEW:
        .append_drop("coverage", cn, "low_coverage",
                     metric = ifelse(is.na(cover), 0, cover), threshold = thr)
      }
    }
  }
  df <- df[, keep_cov, drop = FALSE]
  if (!ncol(df)) return(NULL)
  
  ## NEW: 重新初始化 keep_cov，讓它與「已縮小的 df」同長度且同欄名
  keep_cov <- rep(TRUE, ncol(df))
  names(keep_cov) <- colnames(df)
  ## /NEW
  
  # 4) 與 biology 的相關性過濾（只對數值欄）
  if (!is.null(y)) {
    y <- suppressWarnings(as.numeric(y))
    ## 不對這些共變數做 rho gate
    skip_rho_gate <- c("purity", "age", "sex")
    
    for (cn in colnames(df)) {
      v <- df[[cn]]
      ## 只在「數值且不在 skip 名單」時才做 |rho| 篩選
      if (is.numeric(v) && !(cn %in% skip_rho_gate)) {
        fin <- is.finite(v) & is.finite(y)
        if (sum(fin) >= min_pairs) {
          r <- suppressWarnings(stats::cor(v[fin], y[fin], method = "spearman"))
          if (is.finite(r) && abs(r) >= max_abs_cor) {
            logf("  [covars-%s] drop %s (|rho|=%.3f ≥ %.2f vs biology)", label, cn, abs(r), max_abs_cor)
            keep_cov[cn] <- FALSE
          }
        }
      }
    }
    df <- df[, keep_cov, drop = FALSE]
    if (!ncol(df)) return(NULL)
  }
  
  
  # 5) 回傳：列序 = sample_order；給 model.matrix 展開
  stopifnot(nrow(df) == length(so), identical(rownames(df), so))
  ## NEW: 把每欄 drop 的原因掛在屬性上（由呼叫端決定是否寫 CSV）
  attr(df, "drop_log") <- if (length(drop_log)) do.call(rbind, drop_log) else NULL
  df
}



## ==============================
## 小稽核：sex/age + purity + batch 整體檢查
## ==============================

# 安全：工具函式（若尚未定義）
if (!exists(".norm_names", inherits = FALSE)) {
  .norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+","_", x))
}
if (!exists("%||%", inherits = FALSE)) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b
}
if (!exists(".to01", inherits = FALSE)) {
  .to01 <- function(v){
    v <- suppressWarnings(as.numeric(v))
    # 若 >1 的點比 <=1 多，視為 0~100 的百分比，轉為 0~1
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v/100
    pmin(pmax(v, 0), 1)
  }
}
.median_from_semicolon <- function(x_chr){
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(vv)]; if (!length(vv)) return(NA_real_)
  stats::median(vv)
}

# 小工具：把 batch 的每層樣本數整理成字串
.format_batch_sizes <- function(fac) {
  if (is.null(fac)) return(NA_character_)
  tb <- sort(table(fac), decreasing = TRUE)
  paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}

# --- 新增：依資料集規則審核 purity 的來源與覆蓋率 ---
.audit_purity_for_dataset <- function(ds_id, ds_dir, sample_ids){
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp  <- file.path(ds_dir, "data_clinical_patient.txt")
  has_samp <- file.exists(samp_fp); has_pat <- file.exists(pat_fp)
  
  purity_col <- "NONE"
  purity_vec <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  
  samp <- pat <- NULL
  sid <- pid_samp <- pid_pat <- NA
  
  if (has_samp) {
    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .norm_names(names(samp))
    sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID","PATIENT"), names(samp))[1]
  }
  if (has_pat) {
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(pat) <- .norm_names(names(pat))
    pid_pat  <- intersect(c("PATIENT_ID","PATIENT"), names(pat))[1]
  }
  
  # 建 sample→patient 對應（有需要時）
  map_pt <- NULL
  if (!is.na(sid) && !is.na(pid_samp)) map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])
  
  ds <- ds_id
  
  if (ds == "brca_cptac_2020" && has_pat && !is.na(pid_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat) && !is.null(map_pt)) {
    purity_col <- "PATIENT:ESTIMATE_TUMORPURITY"
    pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_pat]]))
    purity_vec[] <- unname(pt_pur[ map_pt[sample_ids] ])
    
  } else if (ds == "luad_cptac_2020" && has_samp && !is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp)) {
    purity_col <- "SAMPLE:TUMOR_PURITY_BYESTIMATE_RNASEQ"
    purity_vec[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][ match(sample_ids, samp[[sid]]) ])
    
  } else if (ds == "lusc_cptac_2021" && has_samp && !is.na(sid) && "ESTIMATE_TUMORPURITY" %in% names(samp)) {
    purity_col <- "SAMPLE:ESTIMATE_TUMORPURITY"
    purity_vec[] <- .to01(samp[["ESTIMATE_TUMORPURITY"]][ match(sample_ids, samp[[sid]]) ])
    
  } else if (ds == "paad_cptac_2021" && has_samp && !is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
    purity_col <- "SAMPLE:NEOPLASTIC_CELLULARITY(median;)"
    medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
    purity_vec[] <- .to01(medv[ match(sample_ids, samp[[sid]]) ])
    
  } else if (ds == "ucec_cptac_2020" && has_samp && !is.na(sid)) {
    if ("PURITY_CANCER" %in% names(samp)) {
      purity_col <- "SAMPLE:PURITY_CANCER"
      purity_vec[] <- .to01(samp[["PURITY_CANCER"]][ match(sample_ids, samp[[sid]]) ])
    } else {
      # 嘗試由 immune/stroma 推回（若欄位存在）
      pi <- if ("PURITY_IMMUNE" %in% names(samp)) suppressWarnings(as.numeric(samp[["PURITY_IMMUNE"]])) else NULL
      ps <- if ("PURITY_STROMA" %in% names(samp)) suppressWarnings(as.numeric(samp[["PURITY_STROMA"]])) else NULL
      if (!is.null(pi) || !is.null(ps)) {
        idx <- match(sample_ids, samp[[sid]])
        purity_est <- 1 - (pi[idx] %||% 0) - (ps[idx] %||% 0)
        purity_vec[] <- .to01(purity_est)
        purity_col <- "SAMPLE:1-(PURITY_IMMUNE+PURITY_STROMA)"
      }
    }
  } else {
    purity_col <- "NONE"
  }
  
  nonNA <- mean(is.finite(purity_vec)) * 100
  pur_min <- if (any(is.finite(purity_vec))) min(purity_vec, na.rm = TRUE) else NA_real_
  pur_med <- if (any(is.finite(purity_vec))) stats::median(purity_vec, na.rm = TRUE) else NA_real_
  pur_max <- if (any(is.finite(purity_vec))) max(purity_vec, na.rm = TRUE) else NA_real_
  
  list(column = purity_col,
       nonNA_pct = nonNA,
       min = pur_min, median = pur_med, max = pur_max,
       vec = purity_vec)
}

## ==============================
## 小稽核：sex/age + purity + batch 整體檢查（改良版）
## —— 追加：policy/batch_min/limma_min 於日誌與輸出欄位
## ==============================
audit_one_dataset_sa_batch <- function(ds_dir, pipe_policy = "NA", min_per_level = 2) {
  ds_id <- basename(ds_dir)
  # 1) 取蛋白矩陣樣本
  m <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(m)
  
  # 2) 只從 patient 檔抓 SEX/AGE
  samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  pat_fp  <- file.path(ds_dir, "data_clinical_patient.txt")
  has_samp <- file.exists(samp_fp); has_pat <- file.exists(pat_fp)
  
  sex_nonNA <- age_nonNA <- NA_real_; sex_M <- sex_F <- NA_integer_
  age_min <- age_med <- age_max <- NA_real_
  sex_col <- age_col <- "MISSING"; sex_src <- age_src <- if (has_pat) "patient" else "NONE"
  
  if (has_samp && has_pat) {
    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    pat  <- suppressMessages(readr::read_tsv(pat_fp,  show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .norm_names(names(samp))
    names(pat)  <- .norm_names(names(pat))
    
    sid <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID","PATIENT"), names(samp))[1]
    pid_pat  <- intersect(c("PATIENT_ID","PATIENT"), names(pat))[1]
    
    if (!is.na(sid) && !is.na(pid_samp) && !is.na(pid_pat)) {
      map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])
      
      # SEX
      if ("SEX" %in% names(pat)) {
        sex_col <- "SEX"
        raw <- toupper(as.character(pat$SEX))
        val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
        names(val) <- as.character(pat[[pid_pat]])
        sex_vec <- unname(val[ map_pt[sample_ids] ])
        sex_nonNA <- mean(is.finite(sex_vec)) * 100
        sex_M <- sum(sex_vec == 1, na.rm = TRUE)
        sex_F <- sum(sex_vec == 0, na.rm = TRUE)
      }
      
      # AGE（原值統計；主流程再做 z-score）
      if ("AGE" %in% names(pat)) {
        age_col <- "AGE"
        v <- suppressWarnings(as.numeric(pat$AGE))
        names(v) <- as.character(pat[[pid_pat]])
        age_vec <- unname(v[ map_pt[sample_ids] ])
        age_nonNA <- mean(is.finite(age_vec)) * 100
        if (any(is.finite(age_vec))) {
          age_min <- min(age_vec, na.rm = TRUE)
          age_med <- stats::median(age_vec, na.rm = TRUE)
          age_max <- max(age_vec, na.rm = TRUE)
        }
      }
    }
  }
  
  # 3) purity 審核
  pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)
  
  # 4) batch 檢查（用你現有的 get_batch_factor）
  bi <- get_batch_factor(ds_dir, sample_ids,
                         pipe_policy   = pipe_policy,
                         min_per_level = min_per_level)
  batch_col    <- if (!is.null(bi)) bi$name else "NONE"
  batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
  batch_nonNA  <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
  batch_sizes  <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_
  
  # 印一行摘要（補上 policy/batch_min/limma_min）
  line <- sprintf(
    "[Audit] %-24s | SEX:%-7s nonNA=%5.1f%% (M=%d,F=%d) | AGE:%-7s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | PURITY:%-35s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | BATCH:%-18s levels=%-2s nonNA=%5.1f%% [%s] | policy=%s | batch_min=%d | limma_min=%d",
    ds_id,
    sex_col, sex_nonNA %||% NaN, sex_M %||% NA_integer_, sex_F %||% NA_integer_,
    age_col, age_nonNA %||% NaN,
    ifelse(is.finite(age_min), round(age_min,1), "NA"),
    ifelse(is.finite(age_med), round(age_med,1), "NA"),
    ifelse(is.finite(age_max), round(age_max,1), "NA"),
    pur$column, pur$nonNA_pct %||% NaN,
    ifelse(is.finite(pur$min), round(pur$min,3), "NA"),
    ifelse(is.finite(pur$median), round(pur$median,3), "NA"),
    ifelse(is.finite(pur$max), round(pur$max,3), "NA"),
    batch_col, as.integer(batch_levels), batch_nonNA %||% NaN, batch_sizes %||% "NA",
    pipe_policy, min_per_level, min_per_group
  )
  if (exists("log_msg", mode = "function")) log_msg(line) else cat(line, "\n")
  
  data.frame(
    dataset = ds_id,
    sex_column = sex_col, sex_nonNA_pct = round(sex_nonNA,1), sex_M = sex_M, sex_F = sex_F,
    age_column = age_col, age_nonNA_pct = round(age_nonNA,1),
    age_min = ifelse(is.finite(age_min), round(age_min,1), NA),
    age_median = ifelse(is.finite(age_med), round(age_med,1), NA),
    age_max = ifelse(is.finite(age_max), round(age_max,1), NA),
    purity_column = pur$column, purity_nonNA_pct = round(pur$nonNA_pct,1),
    purity_min = ifelse(is.finite(pur$min), round(pur$min,3), NA),
    purity_median = ifelse(is.finite(pur$median), round(pur$median,3), NA),
    purity_max = ifelse(is.finite(pur$max), round(pur$max,3), NA),
    batch_column = batch_col, batch_levels = as.integer(batch_levels),
    batch_nonNA_pct = round(batch_nonNA,1), batch_sizes = batch_sizes,
    batch_policy = pipe_policy,
    batch_min_per_level = as.integer(min_per_level),
    limma_min_per_group = as.integer(min_per_group),
    stringsAsFactors = FALSE
  )
}


audit_all_datasets_sa_batch <- function(dataset_dirs, pipe_policy = "NA", min_per_level = 2) {
  out <- lapply(dataset_dirs, function(p) {
    tryCatch(audit_one_dataset_sa_batch(p, pipe_policy, min_per_level),
             error = function(e) {
               if (exists("log_msg", mode = "function")) log_msg("[Audit] %s ERROR: %s", basename(p), e$message)
               data.frame(
                 dataset = basename(p),
                 sex_column = "ERROR", sex_nonNA_pct = NA, sex_M = NA, sex_F = NA,
                 age_column = "ERROR", age_nonNA_pct = NA, age_min = NA, age_median = NA, age_max = NA,
                 purity_column = "ERROR", purity_nonNA_pct = NA, purity_min = NA, purity_median = NA, purity_max = NA,
                 batch_column = "ERROR", batch_levels = NA, batch_nonNA_pct = NA, batch_sizes = NA,
                 stringsAsFactors = FALSE
               )
             })
  })
  dplyr::bind_rows(out)
}


# ===== Helpers: cameraPR (from ranks) & mroast (from expression+design) =====
run_camera_from_ranks_save <- function(stats_named_vec, pathways, out_dir, outfile = "GSEA_camera.csv",
                                       min_size = 15L, max_size = 500L) {
  if (is.null(pathways) || !length(pathways) || is.null(stats_named_vec) || !length(stats_named_vec)) return(invisible(NULL))
  stats_named_vec <- stats_named_vec[is.finite(stats_named_vec)]
  if (!length(stats_named_vec)) return(invisible(NULL))
  genes <- names(stats_named_vec)
  idx <- try(limma::ids2indices(pathways, genes, remove.empty = TRUE), silent = TRUE)
  if (inherits(idx, "try-error")) return(invisible(NULL))
  # 篩掉過大/過小的 gene sets
  idx <- idx[vapply(idx, function(i) length(i) >= min_size && length(i) <= max_size, logical(1))]
  if (!length(idx)) return(invisible(NULL))
  cam <- try(limma::cameraPR(statistic = stats_named_vec, index = idx, sort = TRUE), silent = TRUE)
  if (inherits(cam, "try-error") || is.null(cam)) return(invisible(NULL))
  cam <- as.data.frame(cam); cam$pathway <- rownames(cam)
  out_csv <- file.path(out_dir, outfile)
  data.table::fwrite(cam, out_csv)
  message(sprintf("[cameraPR] -> %s (n_pathways=%d)", out_csv, nrow(cam)))
  invisible(out_csv)
}

run_mroast_save <- function(expr_mat, design, coef_or_contrast, pathways, out_dir,
                            outfile = "GSEA_mroast.csv",
                            nrot = 9999L, set_statistic = "mean",
                            min_size = 15L, max_size = 500L,
                            universe = NULL) {
  if (is.null(pathways) || !length(pathways) || is.null(expr_mat) || !nrow(expr_mat)) return(invisible(NULL))
  if (is.null(rownames(expr_mat))) stop("expr_mat 需要 rownames() 作為基因 ID")
  
  # 1) gene set 對齊 + 大小過濾（與 ranks 同步的 universe）
  row_genes <- rownames(expr_mat)
  idx <- try(limma::ids2indices(pathways, row_genes, remove.empty = TRUE), silent = TRUE)
  if (inherits(idx, "try-error") || is.null(idx) || !length(idx)) return(invisible(NULL))
  
  # 若指定 universe（例如 names(t-stat)），把 index 映射到 expr_mat 的 row 索引後再做交集
  if (!is.null(universe) && length(universe)) {
    keep_rows <- match(unique(universe), row_genes, nomatch = 0L)
    keep_rows <- unique(keep_rows[keep_rows > 0L])
    if (length(keep_rows)) {
      idx <- lapply(idx, function(i) intersect(i, keep_rows))
    }
  }
  
  # 大小篩選
  idx <- idx[vapply(idx, function(i) length(i) >= min_size && length(i) <= max_size, logical(1))]
  if (!length(idx)) return(invisible(NULL))
  
  # 2) 準備 contrast（可給欄名或索引）
  if (is.character(coef_or_contrast)) {
    ci <- match(coef_or_contrast, colnames(design))
    if (is.na(ci)) stop(sprintf("mroast: 找不到設計矩陣欄名：%s", coef_or_contrast))
    contrast_arg <- ci
  } else if (is.numeric(coef_or_contrast) && length(coef_or_contrast) == 1L) {
    contrast_arg <- as.integer(coef_or_contrast)
  } else if (is.numeric(coef_or_contrast) && length(coef_or_contrast) == ncol(design)) {
    # 直接給對比向量（如需要）
    contrast_arg <- coef_or_contrast
  } else {
    stop("mroast: coef_or_contrast 必須是 欄名/單一整數索引/或等長對比向量")
  }
  
  # 3) 執行 mroast（注意：是 contrast=，不是 coef=）
  mr <- try(limma::mroast(y = expr_mat, index = idx, design = design,
                          contrast = contrast_arg, nrot = nrot,
                          set.statistic = set_statistic),
            silent = TRUE)
  if (inherits(mr, "try-error") || is.null(mr)) return(invisible(NULL))
  
  mr <- as.data.frame(mr); mr$pathway <- rownames(mr)
  out_csv <- file.path(out_dir, outfile)
  data.table::fwrite(mr, out_csv)
  message(sprintf("[mroast] -> %s (n_pathways=%d; contrast=%s)",
                  out_csv, nrow(mr),
                  if (is.numeric(contrast_arg)) paste(contrast_arg, collapse=",") else "vector"))
  invisible(out_csv)
}

## ===== DO_AUDIT 開關包住 QC/探索（4）=====
DO_AUDIT <- TRUE  # <— 預設關閉；需要時改 TRUE

if (DO_AUDIT) {
  results <- list()
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    if (!dir.exists(ds_dir)) { log_msg("略過：找不到資料夾 %s", ds_dir); next }
    results[[ds]] <- screen_batch_need(ds_dir)
  }
  
  ds_dir <- file.path(getwd(), "brca_cptac_2020")
  if (dir.exists(ds_dir)) try(inspect_plex(ds_dir), silent = TRUE)
  
  sex_age_batch_audit <- audit_all_datasets_sa_batch(dataset_dirs)
  print(sex_age_batch_audit)
  readr::write_csv(sex_age_batch_audit, file.path(getwd(), "sex_age_batch_audit_summary.csv"))
}

msigdbr_version <- write_geneset_manifest(genesets_by_group)
yaml::write_yaml(list(
  seed = 1234,
  min_frac_complete   = min_frac_complete,
  hi_lo_quantile      = hi_lo_quantile,
  minSize = minSize, maxSize = maxSize, fgsea_eps = fgsea_eps,
  min_per_group      = min_per_group, 
  min_pairs_spearman = min_pairs_spearman,
  MAKE_PLOTS = MAKE_PLOTS,
  datasets_root = datasets_root,
  dataset_ids = dataset_ids,
  strata = strata,
  geneset_groups_selected = GENESET_GROUPS_TO_RUN,
  msigdbr_version = msigdbr_version
), file.path("run_info","run_manifest.yml"))


## 想從哪一個開始重跑（可改）
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix  <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' 不在 dataset_ids 裡", start_from))
dataset_dirs_run <- dataset_dirs[ ord[ix:length(ord)] ]

exists("coerce_covariates_safely")
getAnywhere("coerce_covariates_safely")
exists("opt")          # 應回 TRUE
getAnywhere("opt")     # 應顯示在 .GlobalEnv

options(csn.audit_sva_sanity = FALSE)  # 想關掉就設 FALSE

## ---- 可選：測試用，只跑特定 strata ----
#only_strata <- c("TP53_mutant")   # 想同時跑多個就寫 c("ALL","TP53_mutant") 等

## 小工具：若未設定 only_strata，則預設全跑
.should_run <- function(tag) {
  if (!exists("only_strata") || is.null(only_strata) || !length(only_strata)) return(TRUE)
  tag %in% only_strata
}


## =========================================================
## 依序處理本輪要跑的 datasets（補齊主迴圈＋定義各分層樣本集）
## =========================================================
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  log_msg("== 開始資料集：%s ==", ds)
  
  ## ---- [NEW | COMBO per-dataset RUN_PASSES] ----
  if (!is.null(COMBO_MODE)) {
    if (COMBO_MODE == "combo_1") {
      options(csn.run_passes = c("BatchAdj"))
    } else if (COMBO_MODE == "combo_2") {
      if (ds %in% .combo2_raw_ds) options(csn.run_passes = "RAW") else options(csn.run_passes = c("BatchAdj"))
    } else if (COMBO_MODE == "combo_3") {
      options(csn.run_passes = "RAW")
    }
    options(csn.run_camera_mroast = FALSE)
  }
  ## ---- [END NEW] ----
  
  ## 蛋白矩陣 & TP53 狀態
  mat0_full   <- load_matrix_from_dataset_dir(ds_dir)
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))
  
  ## 三個分層的樣本集
  samples_ALL <- colnames(mat0_full)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT  <- names(tp53_status)[tp53_status == "TP53_wild_type"]
  
  ## 每個分層的輸出根目錄
  base_tp53_root <- if (is.null(COMBO_PREFIX)) file.path(ds_dir, "csn_gsea_results_TP53") else file.path(COMBO_PREFIX, ds, "csn_gsea_results_TP53")
  
  
  ## 依序跑三個 strata（RAW & BatchAdj 都會在 run_one_stratum 裡寫好）
  if (.should_run("ALL")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_ALL,
      out_root = file.path(base_tp53_root, "ALL"),
      genesets_by_group = genesets_by_group
    )
  }
  
  if (.should_run("TP53_mutant")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_MUT,
      out_root = file.path(base_tp53_root, "TP53_mutant"),
      genesets_by_group = genesets_by_group
    )
  }
  
  if (.should_run("TP53_wild_type")) {
    run_one_stratum(
      ds_id = ds, ds_dir = ds_dir,
      mat0_full = mat0_full,
      sample_keep = samples_WT,
      out_root = file.path(base_tp53_root, "TP53_wild_type"),
      genesets_by_group = genesets_by_group
    )
  }
  
  log_msg("== 完成資料集：%s（TP53 分層輸出 → %s）==", ds, base_tp53_root)
}


## ---------- [NEW] WT vs MT ΔNES 聚合（從既有輸出彙整） ----------
tp53_delta_nes_aggregate <- function(datasets_root,
                                     dataset_ids = NULL,
                                     versions    = c("RAW","BatchAdj"),
                                     groups      = names(genesets_by_group),
                                     stat_tags   = c("GSEA_limma_t_cont","GSEA_spearman")) {
  if (is.null(dataset_ids)) {
    dataset_ids <- list.dirs(datasets_root, full.names = FALSE, recursive = FALSE)
  }
  sfn <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  for (ds in dataset_ids) {
    base_dir <- file.path(datasets_root, ds, "csn_gsea_results_TP53")
    for (ver in versions) {
      # 自動偵測 subunits：讀 WT/MT 資料夾下的子目錄名稱
      subunits <- unique(unlist(lapply(groups, function(g){
        b1 <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, "TP53_wild_type")
        b2 <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, "TP53_mutant")
        c(list.dirs(b1, full.names = FALSE, recursive = FALSE),
          list.dirs(b2, full.names = FALSE, recursive = FALSE))
      })))
      subunits <- subunits[nzchar(subunits)]
      subunits <- subunits[nzchar(subunits)]
      if (!length(subunits)) next
      for (su in subunits) {
        for (grp in groups) {
          grp_safe <- safe_fs_name(grp)
          for (st in stat_tags) {
            fp_wt <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "TP53_wild_type", su, paste0(st, ".csv"))
            fp_mt <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "TP53_mutant", su, paste0(st, ".csv"))
            if (!file.exists(fp_wt) || !file.exists(fp_mt)) next
            wt <- tryCatch(data.table::fread(fp_wt, na.strings=c("NA","NaN","")), error=function(e) NULL)
            mt <- tryCatch(data.table::fread(fp_mt, na.strings=c("NA","NaN","")), error=function(e) NULL)
            if (is.null(wt) || is.null(mt) || !nrow(wt) || !nrow(mt)) next
            wt <- setNames(as.data.frame(wt), tolower(names(wt)))
            mt <- setNames(as.data.frame(mt), tolower(names(mt)))
            # 標準化欄名
            if (!("pathway" %in% names(wt)) && ("term" %in% names(wt))) names(wt)[names(wt)=="term"] <- "pathway"
            if (!("pathway" %in% names(mt)) && ("term" %in% names(mt))) names(mt)[names(mt)=="term"] <- "pathway"
            if (!("nes" %in% names(wt)) && ("enrichment_score" %in% names(wt))) names(wt)[names(wt)=="enrichment_score"] <- "nes"
            if (!("nes" %in% names(mt)) && ("enrichment_score" %in% names(mt))) names(mt)[names(mt)=="enrichment_score"] <- "nes"
            for (nm in c("nes","padj","pval")) if (!nm %in% names(wt)) wt[[nm]] <- NA_real_
            for (nm in c("nes","padj","pval")) if (!nm %in% names(mt)) mt[[nm]] <- NA_real_
            df <- dplyr::full_join(wt[, c("pathway","nes","padj","pval")],
                                   mt[, c("pathway","nes","padj","pval")],
                                   by = "pathway", suffix = c("_WT","_MT"))
            if (!nrow(df)) next
            df <- df |>
              dplyr::mutate(
                NES_WT = as.numeric(nes_WT),
                NES_MT = as.numeric(nes_MT),
                padj_WT = as.numeric(padj_WT),
                padj_MT = as.numeric(padj_MT),
                pval_WT = as.numeric(pval_WT),
                pval_MT = as.numeric(pval_MT),
                delta_NES = NES_MT - NES_WT,
                direction_consistency = dplyr::case_when(
                  is.finite(NES_WT) & is.finite(NES_MT) & NES_MT > 0 & NES_WT > 0 ~ "BOTH_POS",
                  is.finite(NES_WT) & is.finite(NES_MT) & NES_MT < 0 & NES_WT < 0 ~ "BOTH_NEG",
                  is.finite(NES_WT) & is.finite(NES_MT)                           ~ "MIXED",
                  TRUE ~ "NA"
                )
              ) |>
              dplyr::select(pathway, NES_WT, NES_MT, delta_NES, padj_WT, padj_MT, pval_WT, pval_MT, direction_consistency) |>
              dplyr::arrange(dplyr::desc(abs(delta_NES)), pathway)
            out_dir <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", grp_safe, ver, ds, "DeltaWT_MT", su)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out_csv <- file.path(out_dir, paste0(st, "_deltaWTMT.csv"))
            data.table::fwrite(df, out_csv)
          }
        }
      }
      message(sprintf("[ΔNES] %s/%s 完成", ds, ver))
    }
  }
  invisible(TRUE)
}

## =========================================================
## [NEW-1] 產出跨資料集的 meta-FDR（Stouffer → BH）
## 條件：上方 for 迴圈已全部完成且各資料集/分層/統計已落地
## =========================================================
meta_fdr_stouffer(
  dataset_dirs = dataset_dirs_run,                 # ← 用本輪實際跑過且可用的集合
  strata       = strata,                           # c("ALL","TP53_mutant","TP53_wild_type")
  stat_tags    = c("GSEA_limma_t_cont",
                   "GSEA_spearman",
                   "GSEA_limma_interaction"),      # ← 新增這一個
  groups       = names(genesets_by_group),
  out_root = if (is.null(COMBO_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr"
  else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr")
)


## [NEW-1b] WT vs MT ΔNES 聚合（從已落地結果彙整）
try(tp53_delta_nes_aggregate(
  datasets_root = datasets_root,
  dataset_ids   = names(dataset_dirs_run),
  versions      = c("RAW","BatchAdj"),
  groups        = names(genesets_by_group),
  stat_tags     = c("GSEA_limma_t_cont","GSEA_spearman")
), silent = TRUE)

## =========================================================
## [NEW-2] 插補敏感度示範：MinProb vs QRILC（單一 dataset/stratum/predictor）
## 位置：for 迴圈外，避免重複示範；若要換隊列或 predictor，複製改參數即可
## =========================================================
.demo_ds <- "brca_cptac_2020"
if (!is.null(dataset_dirs_run[[.demo_ds]])) {
  run_imputation_sensitivity_demo(
    ds_id = .demo_ds,
    ds_dir = dataset_dirs_run[[.demo_ds]],
    stratum = "ALL",
    predictor_name = "COPS7B", # 也可用 "COPS5" 或 "RESIDUAL_COPS5"
    geneset_group = "H"           # 與目前預設一致（Hallmark）
  )
} else {
  log_msg("[impute-demo] 指定示範隊列 %s 不在 dataset_dirs_run 中，略過示範", .demo_ds)
}




## -----------------------------------------------
## 自動偵測本輪可用的次單元（含 CSN_SCORE 與所有 RESIDUAL_*）
## -----------------------------------------------
# 更穩健版本（建議）
detect_subunits_tp53 <- function(dataset_dirs, strata, versions, groups) {
  blacklist <- c("summary", "run_info", "tables")
  subs <- list()
  for (ds in names(dataset_dirs)) {
    for (st in strata) for (ver in versions) for (g in groups) {
      base <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(g), ver, ds, st)
      if (!dir.exists(base)) next
      cand <- list.dirs(base, full.names = FALSE, recursive = FALSE)
      cand <- cand[cand != "" & !cand %in% blacklist]
      if (length(cand)) subs[[length(subs) + 1L]] <- cand
    }
  }
  if (!length(subs)) return(character(0))
  sort(unique(unlist(subs)))
}


## ------- 你原本的設定 --------
stat_tags <- c("GSEA_limma_t_cont", "GSEA_spearman", "GSEA_limma_interaction")
versions  <- c("RAW","BatchAdj")
groups    <- names(genesets_by_group)
dataset_ids_have <- names(dataset_dirs_run)

## -----------------------------------------------
## I/O：讀單一 CSV（自動容錯 pval 缺欄）
## -----------------------------------------------
read_one_result_tp53 <- function(ds_id, stratum, version, subunit, grp_name, stat_tag){
  fp <- file.path(COMBO_PREFIX %||% "csn_gsea_results_TP53_by_collection", safe_fs_name(grp_name), version, ds_id, stratum, subunit, paste0(stat_tag, ".csv"))
  if (!file.exists(fp)) return(NULL)
  dt <- tryCatch(data.table::fread(fp, na.strings=c("NA","NaN","")), error=function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  need <- intersect(c("pathway","NES","padj","pval"), names(dt))
  if (!all(c("pathway","NES","padj") %in% need)) return(NULL)
  out <- as.data.frame(dt[, need, with = FALSE])
  if (!"pval" %in% names(out)) out$pval <- NA_real_
  transform(out,
            dataset = ds_id, stratum = stratum, version = version,
            group   = grp_name, stat = stat_tag, subunit = subunit,
            NES = as.numeric(NES), padj = as.numeric(padj), pval = as.numeric(pval)) |>
    dplyr::select(dataset, stratum, version, group, stat, pathway, subunit, NES, padj, pval)
}

## -----------------------------------------------
## 收集所有結果（自動偵測 subunits；含 CSN_SCORE 與所有 RESIDUAL_*）
## -----------------------------------------------
collect_all_long_tp53 <- function(dataset_ids, strata, versions, groups, stat_tags, subunits = NULL){
  if (is.null(subunits)) {
    subunits <- detect_subunits_tp53(dataset_dirs = setNames(file.path(datasets_root, dataset_ids), dataset_ids),
                                     strata = strata, versions = versions, groups = groups)
  }
  res <- list(); k <- 1L
  for (ds in dataset_ids)
    for (stt in strata)
      for (verx in versions)
        for (g in groups)
          for (st in stat_tags)
            for (su in subunits) {
              tmp <- read_one_result_tp53(ds, stt, verx, su, g, st)
              if (!is.null(tmp)) { res[[k]] <- tmp; k <- k + 1L }
            }
  if (!length(res)) return(NULL)
  dplyr::bind_rows(res)
}






## =========================================================
## PAN 匯總段落（擴充版）
## —— 自動納入 CSN_SCORE 與所有 RESIDUAL_*；
## —— 追加 subunit_class（BASE / RESIDUAL / SCORE）維度彙整；
## —— 仍輸出你原本的 per_dataset / pan_summary / 表1 / 表2。
## =========================================================
for (ver in versions) {
  log_msg("PAN(TP53,%s) 聚合開始", ver)
  long_df <- collect_all_long_tp53(dataset_ids_have, strata, versions = ver, groups, stat_tags, subunits = NULL)
  if (is.null(long_df) || !nrow(long_df)) { log_msg("  無可用結果（%s）", ver); next }
  
  ## 標記 subunit 類別：BASE / RESIDUAL / SCORE
  long_df <- long_df |>
    dplyr::mutate(
      subunit_class = dplyr::case_when(
        grepl("^RESIDUAL_", subunit) ~ "RESIDUAL",
        subunit == "CSN_SCORE"       ~ "SCORE",
        TRUE                         ~ "BASE"
      ),
      sig_005 = !is.na(padj) & padj < 0.05,
      sig_025 = !is.na(padj) & padj < 0.25,
      pos     = NES > 0,
      neg     = NES < 0
    )
  
  pan_out_root <- if (is.null(COMBO_PREFIX)) file.path(getwd(), "csn_gsea_pan_summary_TP53", ver) else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53", ver)
  dir.create(pan_out_root, recursive = TRUE, showWarnings = FALSE)
  # ---- groups 標籤與版本層級輸出檔名（避免互相覆蓋）----
  groups_tag <- paste(sort(unique(vapply(groups, safe_fs_name, FUN.VALUE = character(1)))), collapse = "-")
  f_per  <- file.path(pan_out_root, paste0("per_dataset_summary__groups_", groups_tag, ".csv"))
  f_pan  <- file.path(pan_out_root, paste0("pan_summary__groups_", groups_tag, ".csv"))
  f_cls  <- file.path(pan_out_root, paste0("pan_summary_by_subunit_class__groups_", groups_tag, ".csv"))
  f_xlsx <- file.path(pan_out_root, paste0("pan_summary__groups_", groups_tag, ".xlsx"))
  
  # 每癌別 × 路徑 摘要（per_dataset）
  per_ds_term <- long_df |>
    dplyr::group_by(dataset, stratum, version, group, stat, pathway) |>
    dplyr::summarise(
      sig_n_padj_0_05 = sum(sig_005, na.rm=TRUE),
      pos_n_padj_0_05 = sum(sig_005 & pos, na.rm=TRUE),
      neg_n_padj_0_05 = sum(sig_005 & neg, na.rm=TRUE),
      sig_n_padj_0_25 = sum(sig_025, na.rm=TRUE),
      pos_n_padj_0_25 = sum(sig_025 & pos, na.rm=TRUE),
      neg_n_padj_0_25 = sum(sig_025 & neg, na.rm=TRUE),
      median_NES_sig_0_05 = suppressWarnings(stats::median(NES[sig_005], na.rm=TRUE)),
      median_NES_sig_0_25 = suppressWarnings(stats::median(NES[sig_025], na.rm=TRUE)),
      median_NES_all      = suppressWarnings(stats::median(NES, na.rm=TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      direction_0_05 = dplyr::case_when(
        pos_n_padj_0_05 > 0 & neg_n_padj_0_05 == 0 ~ "POS",
        neg_n_padj_0_05 > 0 & pos_n_padj_0_05 == 0 ~ "NEG",
        pos_n_padj_0_05 > 0 & neg_n_padj_0_05 > 0  ~ "MIXED",
        TRUE                                        ~ "NONE"
      ),
      direction_0_25 = dplyr::case_when(
        pos_n_padj_0_25 > 0 & neg_n_padj_0_25 == 0 ~ "POS",
        neg_n_padj_0_25 > 0 & pos_n_padj_0_25 == 0 ~ "NEG",
        pos_n_padj_0_25 > 0 & neg_n_padj_0_25 > 0  ~ "MIXED",
        TRUE                                        ~ "NONE"
      )
    )
  
  # 輸出「每癌別」彙整（總表）
  data.table::fwrite(per_ds_term, f_per)
  wb_top <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_top, "per_dataset")
  openxlsx::writeData(wb_top, "per_dataset", per_ds_term)
  
  # pan-cancer 匯總（跨 dataset 聚合）
  pan_term <- per_ds_term |>
    dplyr::group_by(stratum, version, group, stat, pathway) |>
    dplyr::summarise(
      datasets_with_sig_0_05 = sum(sig_n_padj_0_05 > 0, na.rm = TRUE),
      total_sig_0_05         = sum(sig_n_padj_0_05,    na.rm = TRUE),
      total_pos_0_05         = sum(pos_n_padj_0_05,    na.rm = TRUE),
      total_neg_0_05         = sum(neg_n_padj_0_05,    na.rm = TRUE),
      median_NES_sig_0_05    = suppressWarnings(stats::median(median_NES_sig_0_05, na.rm = TRUE)),
      datasets_with_sig_0_25 = sum(sig_n_padj_0_25 > 0, na.rm = TRUE),
      total_sig_0_25         = sum(sig_n_padj_0_25,    na.rm = TRUE),
      total_pos_0_25         = sum(pos_n_padj_0_25,    na.rm = TRUE),
      total_neg_0_25         = sum(neg_n_padj_0_25,    na.rm = TRUE),
      median_NES_sig_0_25    = suppressWarnings(stats::median(median_NES_sig_0_25, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      prop_pos_all = ifelse((total_pos_0_05 + total_neg_0_05) > 0,
                            total_pos_0_05 / (total_pos_0_05 + total_neg_0_05), NA_real_),
      prop_neg_all = ifelse((total_pos_0_05 + total_neg_0_05) > 0,
                            total_neg_0_05 / (total_pos_0_05 + total_neg_0_05), NA_real_),
      direction_consistency = dplyr::case_when(
        total_pos_0_05 > 0 & total_neg_0_05 == 0 ~ "ALL_POS",
        total_neg_0_05 > 0 & total_pos_0_05 == 0 ~ "ALL_NEG",
        total_pos_0_05 > 0 & total_neg_0_05 > 0  ~ "MIXED",
        TRUE                                     ~ "NONE"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      dplyr::desc(datasets_with_sig_0_05),
      dplyr::desc(total_sig_0_05),
      pathway
    )
  
  data.table::fwrite(pan_term, f_pan)
  openxlsx::addWorksheet(wb_top, "pan_summary")
  openxlsx::writeData(wb_top, "pan_summary", pan_term)
  
  ## ------ 新增：按 subunit_class 的 pan 彙總（BASE / RESIDUAL / SCORE）
  pan_by_class <- long_df |>
    dplyr::group_by(stratum, version, group, stat, pathway, subunit_class) |>
    dplyr::summarise(
      datasets_with_sig_0_05 = dplyr::n_distinct(dataset[sig_005], na.rm = TRUE),
      total_sig_0_05         = sum(sig_005, na.rm = TRUE),
      total_pos_0_05         = sum(sig_005 & pos, na.rm = TRUE),
      total_neg_0_05         = sum(sig_005 & neg, na.rm = TRUE),
      median_NES_sig_0_05    = suppressWarnings(stats::median(NES[sig_005], na.rm = TRUE)),
      datasets_with_sig_0_25 = dplyr::n_distinct(dataset[sig_025], na.rm = TRUE),
      total_sig_0_25         = sum(sig_025, na.rm = TRUE),
      total_pos_0_25         = sum(sig_025 & pos, na.rm = TRUE),
      total_neg_0_25         = sum(sig_025 & neg, na.rm = TRUE),
      median_NES_sig_0_25    = suppressWarnings(stats::median(NES[sig_025], na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(stratum, group, stat, pathway, subunit_class)
  
  data.table::fwrite(pan_by_class, f_cls)
  openxlsx::addWorksheet(wb_top, "pan_by_class")
  openxlsx::writeData(wb_top, "pan_by_class", pan_by_class)
  
  openxlsx::saveWorkbook(wb_top, f_xlsx, overwrite = TRUE)
  
  ## 每個 stratum × gene-set group × stat 的輸出（保留原本表1/表2；新增把 class 結果也寫進去）
  for (stt in strata) {
    for (g in groups) {
      for (st in stat_tags) {
        sub_out <- file.path(pan_out_root, stt, safe_fs_name(g), st)
        dir.create(sub_out, recursive = TRUE, showWarnings = FALSE)
        
        long_sub   <- long_df    |> dplyr::filter(stratum == stt, group == g, stat == st)
        if (nrow(long_sub) == 0) next
        
        per_ds_sub <- per_ds_term |> dplyr::filter(stratum == stt, group == g, stat == st)
        pan_sub    <- pan_term    |> dplyr::filter(stratum == stt, group == g, stat == st)
        
        ## —— I²（等權近似）計算
        calc_I2_unweighted <- function(v){
          v <- v[is.finite(v)]
          k <- length(v); if (k <= 1) return(NA_real_)
          m <- mean(v); Q <- sum((v - m)^2)
          I2 <- ifelse(Q > 0, max(0, (Q - (k - 1)) / Q), 0)
          100 * I2
        }
        nes_by_term <- long_sub |>
          dplyr::group_by(pathway, dataset) |>
          dplyr::summarise(medNES = suppressWarnings(stats::median(NES, na.rm=TRUE)), .groups="drop")
        I2_tbl <- nes_by_term |>
          dplyr::group_by(pathway) |>
          dplyr::summarise(I2_unweighted = calc_I2_unweighted(medNES), .groups="drop")
        
        pan_sub <- pan_sub |>
          dplyr::left_join(I2_tbl, by = "pathway") |>
          dplyr::arrange(
            dplyr::desc(datasets_with_sig_0_05),
            dplyr::desc(total_pos_0_05 - total_neg_0_05),
            pathway
          )
        
        data.table::fwrite(long_sub,   file.path(sub_out, "all_pairs_long.csv"))
        data.table::fwrite(per_ds_sub, file.path(sub_out, "per_dataset_per_term_summary.csv"))
        data.table::fwrite(pan_sub,    file.path(sub_out, "pan_term_summary.csv"))
        data.table::fwrite(pan_sub |> dplyr::filter(datasets_with_sig_0_05 > 0),
                           file.path(sub_out, "pan_term_summary_padjLT0.05.csv"))
        data.table::fwrite(pan_sub |> dplyr::filter(datasets_with_sig_0_25 > 0),
                           file.path(sub_out, "pan_term_summary_padjLT0.25.csv"))
        
        ## 表1/表2（原樣保留）
        out_tables <- file.path(sub_out, "tables")
        dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
        
        df_gs <- long_sub |>
          dplyr::mutate(
            sig_005 = !is.na(padj) & padj < 0.05,
            sig_025 = !is.na(padj) & padj < 0.25,
            pos = NES > 0, neg = NES < 0
          )
        
        # 表1：各 pathway 在各 subunit 的 avg NES + 各癌別平均顯著計數
        tab1_avgNES <- df_gs |>
          dplyr::group_by(pathway, subunit) |>
          dplyr::summarise(avg_NES = mean(NES, na.rm=TRUE), .groups="drop") |>
          tidyr::pivot_wider(names_from = subunit, values_from = avg_NES, names_prefix = "avg_NES_")
        
        per_ds_counts <- df_gs |>
          dplyr::group_by(dataset, pathway) |>
          dplyr::summarise(
            sig_n_padj_0_05 = sum(sig_005, na.rm=TRUE),
            pos_n_padj_0_05 = sum(sig_005 & pos, na.rm=TRUE),
            neg_n_padj_0_05 = sum(sig_005 & neg, na.rm=TRUE),
            sig_n_padj_0_25 = sum(sig_025, na.rm=TRUE),
            pos_n_padj_0_25 = sum(sig_025 & pos, na.rm=TRUE),
            neg_n_padj_0_25 = sum(sig_025 & neg, na.rm=TRUE),
            .groups = "drop"
          )
        tab1_avgCounts <- per_ds_counts |>
          dplyr::group_by(pathway) |>
          dplyr::summarise(
            dplyr::across(
              c(sig_n_padj_0_05, pos_n_padj_0_05, neg_n_padj_0_05,
                sig_n_padj_0_25, pos_n_padj_0_25, neg_n_padj_0_25),
              ~ mean(.x, na.rm=TRUE),
              .names="avg_{.col}"
            ),
            .groups="drop"
          )
        table1 <- dplyr::left_join(tab1_avgNES, tab1_avgCounts, by="pathway")
        data.table::fwrite(table1, file.path(out_tables, "table1_pathway_avgNES_per_subunit_and_avgCounts_across_datasets.csv"))
        wb1 <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb1, "table1")
        openxlsx::writeData(wb1, "table1", table1)
        openxlsx::saveWorkbook(wb1, file.path(out_tables, "table1_pathway_avgNES_per_subunit_and_avgCounts_across_datasets.xlsx"), overwrite=TRUE)
        
        # 表2：各 pathway 在各 dataset 的 NES（每個 subunit 一張表）
        ds_order <- sort(unique(df_gs$dataset))
        for (su in unique(df_gs$subunit)) {
          df_su <- df_gs |> dplyr::filter(subunit == su)
          if (!nrow(df_su)) next
          tab2_nes_per_ds <- df_su |>
            dplyr::select(pathway, dataset, NES) |>
            tidyr::pivot_wider(names_from = dataset, values_from = NES, names_prefix = "NES_")
          nes_cols_order <- paste0("NES_", ds_order)
          nes_cols_order <- nes_cols_order[nes_cols_order %in% names(tab2_nes_per_ds)]
          tab2_nes_per_ds <- tab2_nes_per_ds |>
            dplyr::select(pathway, tidyselect::all_of(nes_cols_order))
          tab2_counts_int <- df_su |>
            dplyr::group_by(pathway) |>
            dplyr::summarise(
              sig_n_padj_0_05 = as.integer(sum(sig_005, na.rm = TRUE)),
              pos_n_padj_0_05 = as.integer(sum(sig_005 & pos, na.rm = TRUE)),
              neg_n_padj_0_05 = as.integer(sum(sig_005 & neg, na.rm = TRUE)),
              sig_n_padj_0_25 = as.integer(sum(sig_025, na.rm = TRUE)),
              pos_n_padj_0_25 = as.integer(sum(sig_025 & pos, na.rm = TRUE)),
              neg_n_padj_0_25 = as.integer(sum(sig_025 & neg, na.rm = TRUE)),
              .groups = "drop"
            )
          table2 <- dplyr::left_join(tab2_nes_per_ds, tab2_counts_int, by = "pathway")
          base <- paste0("table2_pathway_datasetNES_for_", su)
          data.table::fwrite(table2, file.path(out_tables, paste0(base, ".csv")))
          wb2 <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb2, "table2")
          openxlsx::writeData(wb2, "table2", table2)
          openxlsx::saveWorkbook(wb2, file.path(out_tables, paste0(base, ".xlsx")), overwrite = TRUE)
        }
        
        ## 額外把「subunit_class 匯總」也分檔輸出
        pan_sub_by_class <- pan_by_class |>
          dplyr::filter(stratum == stt, group == g, stat == st)
        data.table::fwrite(pan_sub_by_class, file.path(sub_out, "pan_term_summary_by_subunit_class.csv"))
      }
    }
  }
  
  if (exists("pan_out_root")) {
    message(sprintf("[PAN-TP53] 整合完成（%s）→ %s", ver, pan_out_root))
  }
}

## =========================================================
## 後處理：只用已輸出的 CSV 生成每個 dataset/stratum/version 下
##         的 summary/H/GSEA_limma_interaction（不重跑任何模型）
## =========================================================
posthoc_summary_interaction <- function(
    dataset_dirs = dataset_dirs_run,
    strata = c("ALL","TP53_mutant","TP53_wild_type"),
    versions = c("RAW","BatchAdj"),
    groups = names(genesets_by_group)
){
  # 這裡不讀矩陣、不跑模型；僅用 summarize_all_groups() 從現有 CSV 彙總
  # summarize_all_groups 會呼叫 read_gsea_table()，可容錯不同路徑排版
  # （包含 ./<group>/<subunit>/H/STAT.csv 或 ./<subunit>/<group>/STAT.csv）：
  # 見 read_gsea_table 的候選路徑設計。:contentReference[oaicite:1]{index=1}
  
  for (ds in names(dataset_dirs)) {
    ds_dir <- dataset_dirs[[ds]]
    base_root <- file.path(ds_dir, "csn_gsea_results_TP53")
    
    for (st in strata) for (ver in versions) {
      ver_root <- file.path(base_root, st, ver)
      if (!dir.exists(ver_root)) next
      
      # 不探勘現場子目錄，直接用全套 subunits + CSN_SCORE + RESIDUAL_*；
      # read_gsea_table 讀不到的會回傳 NULL，merge_subunit_tables 會自動略過。
      sum_units <- c(csn_subunits, "CSN_SCORE", paste0("RESIDUAL_", csn_subunits))
      
      if (exists("log_msg", mode="function")) {
        try(log_msg("後處理彙整：ds=%s | stratum=%s | version=%s | stat=GSEA_limma_interaction",
                    ds, st, ver), silent = TRUE)
      }
      
      summarize_all_groups(
        out_root          = ver_root,
        csn_subunits      = sum_units,
        genesets_by_group = genesets_by_group,
        stat_tags         = c("GSEA_limma_interaction")
      )
      # summarize_all_groups 會把輸出寫到：
      # <ver_root>/summary/<group>/GSEA_limma_interaction/...
      # 並自動生成 *_ALL.csv、padjLT0.05/0.25 與 .xlsx；詳見其內部寫檔：:contentReference[oaicite:2]{index=2} 和 :contentReference[oaicite:3]{index=3}
    }
  }
  invisible(TRUE)
}

## 直接執行一次，從既有 CSV 產出 interaction 的 summary
posthoc_summary_interaction()






















## =========================================================
## [NEW] CSN subunits pairwise correlations per dataset
##   - 三版本：NoCovariate / RAW_covars / BatchAdj_covars
##   - 產出：CSV + XLSX；（可選）高畫質散佈圖 tiff/png/jpg（預設不畫）
##   - 不依賴原 GSEA 迴圈
## =========================================================

if (!exists("opt", mode = "function")) {
  opt <- function(nm, default) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else default
}

.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# 簡易殘差化（不含 csn_score；僅對 batch + covars 殘差）
residualize_to_covars <- function(y, batch = NULL, covars = NULL, min_n = 8L) {
  stopifnot(!is.null(names(y)))
  sam <- names(y)
  DF <- data.frame(row.names = sam, check.names = FALSE)
  if (!is.null(batch)) DF[["batch"]] <- droplevels(batch[sam])
  if (!is.null(covars)) {
    C <- as.data.frame(covars, check.names = FALSE)
    rn <- rownames(C); if (is.null(rn)) stop("[residualize_to_covars] covars 必須 rownames=樣本")
    C <- C[sam, , drop = FALSE]
    # 與你現有邏輯一致的安全轉型
    if (exists("coerce_covariates_safely", mode = "function")) C <- coerce_covariates_safely(C)
    for (cn in colnames(C)) DF[[cn]] <- C[[cn]]
  }
  yv <- suppressWarnings(as.numeric(y[sam]))
  # 清掉全 NA/常數/單一水準欄，避免設計矩陣奇異
  if (ncol(DF) > 0) {
    all_na <- vapply(DF, function(v) all(is.na(v)), logical(1))
    if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]
    is_const <- vapply(DF, function(v) {
      vv <- v[!is.na(v)]
      if (!length(vv)) TRUE else if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
    }, logical(1))
    if (any(is_const)) DF <- DF[, !is_const, drop = FALSE]
  }
  ok <- is.finite(yv) & (if (ncol(DF) == 0) TRUE else stats::complete.cases(DF))
  if (sum(ok) < min_n) return(setNames(rep(NA_real_, length(y)), names(y)))
  des <- if (ncol(DF) == 0) model.matrix(~ 1) else model.matrix(~ 1 + ., data = DF[ok, , drop = FALSE])
  fit <- lm.fit(x = des, y = yv[ok])
  out <- setNames(rep(NA_real_, length(y)), names(y))
  out[ok] <- fit$residuals
  out
}

# 單一 dataset：產生三版本 pairwise correlation
csn_pairwise_correlation_one_ds <- function(
    ds_id, ds_dir,
    out_root = "CSN_subunits_correlation_coefficient",
    subunits = csn_subunits,
    min_pairs = 10L,
    dpi = 600,
    width = 4.5, height = 4.5,  # inches（期刊友善）
    point_size = 1.6
){
  
  # 準備輸出目錄
  dir.create(out_root, showWarnings = FALSE)
  ds_out <- file.path(out_root, ds_id)
  dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)
  
  # 讀矩陣
  mat0 <- load_matrix_from_dataset_dir(ds_dir)  # 你現成的讀檔器
  present <- intersect(subunits, rownames(mat0))
  if (length(present) < 2) {
    log_msg("[pairwise] %s: 可用 CSN subunits < 2，略過", ds_id)
    return(invisible(NULL))
  }
  sam_all <- colnames(mat0)
  
  # 共變項候選（purity/sex/age）
  build_covars_df <- function(ds_id, ds_dir, sample_ids){
    pur <- get_purity_covariate(ds_id, ds_dir, sample_ids)        # 既有
    sa  <- get_sex_age_covariates(ds_dir, sample_ids)             # 既有
    df  <- data.frame(purity = as.numeric(pur),
                      sex    = as.numeric(sa[, "sex"]),
                      age    = as.numeric(sa[, "age"]),
                      row.names = sample_ids,
                      check.names = FALSE)
    if (exists("coerce_covariates_safely", mode = "function"))
      df <- coerce_covariates_safely(df)
    df
  }
  cov0 <- build_covars_df(ds_id, ds_dir, sam_all)
  
  # 批次（若有），並對齊樣本
  bi <- get_batch_factor(ds_dir, sam_all)   # 既有
  batch <- NULL
  if (!is.null(bi) && !is.null(bi$fac)) {
    b <- bi$fac
    if (is.null(names(b))) {
      if (length(b) == length(sam_all)) {
        b <- stats::setNames(b, sam_all)
      } else {
        b <- NULL
      }
    }
    if (!is.null(b)) batch <- droplevels(b[sam_all])
  }
  
  
  # 建一份供 SVA / tech PC 用的「插補矩陣」
  keep_rows <- rowMeans(is.finite(mat0)) >= opt("min_frac_complete", 0.75)
  Mimp <- mat0[keep_rows, , drop = FALSE]
  set.seed(1234); Mimp <- imputeLCMD::impute.MinProb(Mimp, q = 0.01)
  
  # Dataset-level（與特定 predictor 無關）SVA：mod= cov0(+batch)，mod0= ~1
  sv_df <- NULL
  # 只有在沒有 batch 的情況下才嘗試 SVA（符合 batch > SV > tech PC 優先序）
  if (is.null(batch) && exists("estimate_svs", mode = "function")) {
    # 先確保 Mimp 沒有遺漏值（理論上 impute.MinProb 後不會有，但保險）
    if (any(!is.finite(Mimp))) {
      for (i in seq_len(nrow(Mimp))) {
        v <- Mimp[i, ]
        if (any(!is.finite(v))) {
          med <- stats::median(v[is.finite(v)], na.rm = TRUE)
          if (!is.finite(med)) med <- 0
          v[!is.finite(v)] <- med
          Mimp[i, ] <- v
        }
      }
    }
    
    # 建設計資料框並 rownames 對齊樣本；去掉全 NA 與常數欄
    DFmm <- data.frame(cov0[sam_all, , drop = FALSE], row.names = sam_all, check.names = FALSE)
    if (ncol(DFmm) > 0) {
      all_na <- vapply(DFmm, function(v) all(is.na(v)), logical(1))
      if (any(all_na)) DFmm <- DFmm[, !all_na, drop = FALSE]
      is_const <- vapply(DFmm, function(v) {
        vv <- v[!is.na(v)]
        if (!length(vv)) TRUE else if (is.factor(v)) nlevels(droplevels(vv)) <= 1 else stats::var(as.numeric(vv)) == 0
      }, logical(1))
      if (any(is_const)) DFmm <- DFmm[, !is_const, drop = FALSE]
    }
    
    mod  <- try(model.matrix(~ 1 + ., data = DFmm, na.action = stats::na.pass), silent = TRUE)
    if (inherits(mod, "try-error")) mod <- model.matrix(~ 1, data = DFmm, na.action = stats::na.pass)
    mod0 <- try(model.matrix(~ 1, data = DFmm, na.action = stats::na.pass), silent = TRUE)
    if (inherits(mod0, "try-error")) mod0 <- model.matrix(~ 1, data = DFmm, na.action = stats::na.pass)
    
    sv <- try(estimate_svs(M = Mimp, mod_interest = mod, mod_nuisance = mod0, max_k = 5,
                           label = sprintf("%s-pair", ds_id)), silent = TRUE)
    if (!inherits(sv, "try-error") && !is.null(sv)) {
      sv_df <- as.data.frame(sv)
      if (ncol(sv_df) > 0) {
        colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
        rownames(sv_df) <- colnames(Mimp)
      } else {
        sv_df <- NULL
      }
    } else {
      sv_df <- NULL
    }
  }
  
  
  # Tech PCs（前 2 個）；per-pair 會用 gate_tech_pcs() 依 predictor 篩
  tech_df <- NULL
  pc <- try(stats::prcomp(t(Mimp), center = TRUE, scale. = TRUE), silent = TRUE)
  if (!inherits(pc, "try-error") && !is.null(pc$x)) {
    k <- min(2L, ncol(pc$x))
    tech_df <- as.data.frame(pc$x[, seq_len(k), drop = FALSE])
    colnames(tech_df) <- paste0("PC", seq_len(k))
    rownames(tech_df) <- rownames(tech_df)  # sample names already in rownames(pc$x)
  }
  
  # 兩兩配對
  pairs <- utils::combn(present, 2, simplify = FALSE)
  
  # 收集三版本結果
  all_rows <- list()
  
  # 圖檔輸出的 helper（預設不執行；由 MAKE_PLOTS 控制）
  .save_scatter <- function(df_xy, title, out_base){
    gp <- ggplot2::ggplot(df_xy, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
      ggplot2::labs(title = title, x = df_xy$gx[1], y = df_xy$gy[1]) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
      )
    if (isTRUE(opt("MAKE_PLOTS", FALSE))) {
      ggplot2::ggsave(paste0(out_base, ".tiff"), gp, width = width, height = height, dpi = dpi, compression = "lzw")
      ggplot2::ggsave(paste0(out_base, ".png"),  gp, width = width, height = height, dpi = dpi)
      ggplot2::ggsave(paste0(out_base, ".jpg"),  gp, width = width, height = height, dpi = dpi)
    }
  }
  
  # 主迴圈：每一對 subunits
  for (p in pairs) {
    gx <- p[1]; gy <- p[2]
    x <- as.numeric(mat0[gx, sam_all]); names(x) <- sam_all
    y <- as.numeric(mat0[gy, sam_all]); names(y) <- sam_all
    ok <- stats::complete.cases(x, y)
    if (sum(ok) < min_pairs) next
    
    ## ---- 版本 1: NoCovariate ----
    rr <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "pearson"))
    rs <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    fit <- stats::lm(y ~ x, data = data.frame(x = x[ok], y = y[ok]))
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      dataset = ds_id, version = "NoCovariate",
      gene_x = gx, gene_y = gy, n = sum(ok),
      pearson_r = as.numeric(rr$estimate), pearson_p = rr$p.value,
      spearman_r = as.numeric(rs$estimate), spearman_p = rs$p.value,
      R2 = summary(fit)$r.squared,
      slope = unname(stats::coef(fit)[["x"]]),
      intercept = unname(stats::coef(fit)[["(Intercept)"]]),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    
    .save_scatter(
      df_xy = data.frame(x = x[ok], y = y[ok], gx = gx, gy = gy),
      title = sprintf("%s | %s vs %s (NoCovariate)", ds_id, gx, gy),
      out_base = file.path(ds_out, sprintf("%s_vs_%s_NoCovariate", .safe_fs(gx), .safe_fs(gy)))
    )
    
    ## ---- 版本 2: RAW_covars（依 predictor = gx 選 covariates）----
    # RAW_covars：依你的新政策，直接保留 sex/age/purity，不做 coverage/|rho| 篩選
    cov_raw <- cov0
    
    xr <- residualize_to_covars(x, batch = NULL, covars = cov_raw)
    yr <- residualize_to_covars(y, batch = NULL, covars = cov_raw)
    ok2 <- stats::complete.cases(xr, yr)
    if (sum(ok2) >= min_pairs) {
      rr2 <- suppressWarnings(stats::cor.test(xr[ok2], yr[ok2], method = "pearson"))
      rs2 <- suppressWarnings(stats::cor.test(xr[ok2], yr[ok2], method = "spearman", exact = FALSE))
      fit2 <- stats::lm(yr ~ xr, data = data.frame(xr = xr[ok2], yr = yr[ok2]))
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        dataset = ds_id, version = "RAW_covars",
        gene_x = gx, gene_y = gy, n = sum(ok2),
        pearson_r = as.numeric(rr2$estimate), pearson_p = rr2$p.value,
        spearman_r = as.numeric(rs2$estimate), spearman_p = rs2$p.value,
        R2 = summary(fit2)$r.squared,
        slope = unname(stats::coef(fit2)[["xr"]]),
        intercept = unname(stats::coef(fit2)[["(Intercept)"]]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      .save_scatter(
        df_xy = data.frame(x = xr[ok2], y = yr[ok2], gx = gx, gy = gy),
        title = sprintf("%s | %s vs %s (RAW_covars)", ds_id, gx, gy),
        out_base = file.path(ds_out, sprintf("%s_vs_%s_RAW_covars", .safe_fs(gx), .safe_fs(gy)))
      )
    }
    
    ## ---- 版本 3: BatchAdj_covars（優先序：batch > SV > tech PCs；三者互斥）----
    cov_ba <- cov_raw               # 先沿用 RAW 的共變項選擇
    batch_for_res <- NULL           # 只在「使用 batch」時帶入 residualize_to_covars 的 batch 參數
    
    if (!is.null(batch)) {
      # 有 batch：只用 batch（不再加 SV/tech）
      batch_for_res <- batch
      
    } else if (!is.null(sv_df)) {
      # 沒有 batch，有 SV：只用 SV（不再加 tech）
      sv_aln <- sv_df[sam_all, , drop = FALSE]
      cov_ba <- cbind(cov_ba, sv_aln)
      
    } else if (!is.null(tech_df)) {
      # 沒有 batch、也沒有 SV：退而用 tech PCs（以 predictor=x 做 gating；若無 gate_tech_pcs 就直接用前2PC）
      if (exists("gate_tech_pcs", mode = "function")) {
        tech_g <- gate_tech_pcs(tech_df[sam_all, , drop = FALSE], v = as.numeric(x))
      } else {
        tech_g <- tech_df[sam_all, , drop = FALSE]
      }
      cov_ba <- cbind(cov_ba, tech_g)
    }
    
    xb <- residualize_to_covars(x, batch = batch_for_res, covars = cov_ba)
    yb <- residualize_to_covars(y, batch = batch_for_res, covars = cov_ba)
    ok3 <- stats::complete.cases(xb, yb)
    
    if (sum(ok3) >= min_pairs) {
      rr3  <- suppressWarnings(stats::cor.test(xb[ok3], yb[ok3], method = "pearson"))
      rs3  <- suppressWarnings(stats::cor.test(xb[ok3], yb[ok3], method = "spearman", exact = FALSE))
      fit3 <- stats::lm(yb ~ xb, data = data.frame(xb = xb[ok3], yb = yb[ok3]))
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        dataset = ds_id, version = "BatchAdj_covars",
        gene_x = gx, gene_y = gy, n = sum(ok3),
        pearson_r = as.numeric(rr3$estimate), pearson_p = rr3$p.value,
        spearman_r = as.numeric(rs3$estimate), spearman_p = rs3$p.value,
        R2 = summary(fit3)$r.squared,
        slope = unname(stats::coef(fit3)[["xb"]]),
        intercept = unname(stats::coef(fit3)[["(Intercept)"]]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      .save_scatter(
        df_xy = data.frame(x = xb[ok3], y = yb[ok3], gx = gx, gy = gy),
        title = sprintf("%s | %s vs %s (BatchAdj_covars)", ds_id, gx, gy),
        out_base = file.path(ds_out, sprintf("%s_vs_%s_BatchAdj_covars", .safe_fs(gx), .safe_fs(gy)))
      )
    }
  } # end for pairs
  
  # 彙整/校正與輸出
  if (!length(all_rows)) {
    log_msg("[pairwise] %s: 無可用配對（樣本不足或 NA 過多）", ds_id)
    return(invisible(NULL))
  }
  RES <- do.call(rbind, all_rows)
  # 依 dataset + version 做 BH 校正
  RES$pearson_padj <- ave(RES$pearson_p, interaction(RES$dataset, RES$version, drop = TRUE),
                          FUN = function(p) stats::p.adjust(p, method = "BH"))
  RES$spearman_padj <- ave(RES$spearman_p,
                           interaction(RES$dataset, RES$version, drop = TRUE),
                           FUN = function(p) stats::p.adjust(p, method = "BH"))
  # 落地 CSV 與 XLSX
  out_csv  <- file.path(ds_out, "pairwise_correlations_all_versions.csv")
  data.table::fwrite(RES, out_csv)
  
  out_xlsx <- file.path(ds_out, "pairwise_correlations_all_versions.xlsx")
  wb <- openxlsx::createWorkbook()
  for (ver in unique(RES$version)) {
    openxlsx::addWorksheet(wb, ver)
    openxlsx::writeData(wb, ver, RES[RES$version == ver, ], withFilter = TRUE)
  }
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  
  log_msg("[pairwise] %s: 完成。CSV: %s；XLSX: %s", ds_id, basename(out_csv), basename(out_xlsx))
  invisible(RES)
}

# 一鍵跑全部 dataset（不需要也不會呼叫你的 GSEA 迴圈）
run_csn_subunits_pairwise_correlations <- function(
    dataset_dirs_map = NULL,
    out_root = "CSN_subunits_correlation_coefficient"
){
  if (is.null(dataset_dirs_map)) {
    # 盡量重用你現成的 dataset_dirs_run；若不存在，退回 dataset_dirs
    if (exists("dataset_dirs_run")) {
      dataset_dirs_map <- get("dataset_dirs_run", inherits = TRUE)
    } else if (exists("dataset_dirs")) {
      # 僅挑有檔的
      dd <- get("dataset_dirs", inherits = TRUE)
      dataset_dirs_map <- dd[dir.exists(dd) & file.exists(file.path(dd, "data_protein_quantification.txt"))]
    } else stop("找不到 dataset_dirs_run 或 dataset_dirs")
  }
  dir.create(out_root, showWarnings = FALSE)
  for (ds in names(dataset_dirs_map)) {
    try(csn_pairwise_correlation_one_ds(ds_id = ds, ds_dir = dataset_dirs_map[[ds]], out_root = out_root), silent = FALSE)
  }
  invisible(TRUE)
}

## ----（選擇性）執行：預設僅定義函式，不自動執行 ----
# run_csn_subunits_pairwise_correlations()


## =========================================================
## 直接執行：CSN subunits 兩兩相關（不依賴你的 GSEA 迴圈）
##  - 預設不輸出散佈圖；若要同時輸出 .tiff/.png/.jpg，改 MAKE_PLOTS <- TRUE
## =========================================================
MAKE_PLOTS <- FALSE

run_csn_subunits_pairwise_correlations(
  dataset_dirs_map = if (exists("dataset_dirs_run")) dataset_dirs_run else NULL,
  out_root = "CSN_subunits_correlation_coefficient"
)








## =========================================================
## [NEW] Correlation coefficient heatmaps from pairwise CSVs
##   - For each dataset under CSN_subunits_correlation_coefficient/<ds>/
##   - Read: pairwise_correlations_all_versions.csv
##   - Make 3 versions: NoCovariate / RAW_covars / BatchAdj_covars
##   - Save: .tiff (LZW), .pdf, .jpg, .png (600 dpi)
## =========================================================
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr", repos = "https://cloud.r-project.org")
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(stringr)
})

CSN_SUB_ORDER <- c("GPS1","COPS2","COPS3","COPS4","COPS5","COPS6","COPS7A","COPS7B","COPS8","COPS9")
CELL_BLUE <- "#3B4CC0"; CELL_WHITE <- "#F7F7F7"; CELL_RED <- "#B40426"

# 安全檔名
.safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# 從 pairwise CSV（單 dataset）整理出某版本的方陣資料（r 與 padj），補齊缺格、做對稱
# 從 pairwise CSV（單 dataset）整理出某版本的方陣資料（r 與 padj）
# 修正：
#  - 只保留在該 dataset 該版本實際出現的 subunits（依 CSN_SUB_ORDER 排序）
#  - 自行補上對角線 r=1（且不標黑點）
#  - 僅輸出「下三角（含對角線）」的格子
.build_corr_grid <- function(csv_path, version, sub_order = CSN_SUB_ORDER,
                             triangle = getOption("csn_heatmap_triangle", "full"),
                             method   = getOption("csn_heatmap_method", "pearson")) {
  stopifnot(file.exists(csv_path))
  df <- suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE))
  req <- c("version","gene_x","gene_y","pearson_r","pearson_padj","spearman_r","spearman_padj")
  if (!all(req %in% names(df))) stop("CSV 缺欄位：", paste(setdiff(req, names(df)), collapse = ", "))
  
  m <- match.arg(tolower(method), c("pearson","spearman"))
  rcol <- if (m == "pearson") "pearson_r"   else "spearman_r"
  pcol <- if (m == "pearson") "pearson_padj" else "spearman_padj"
  
  dfv <- df %>% dplyr::filter(.data$version == !!version) %>%
    dplyr::mutate(
      gene_x = stringr::str_trim(gene_x),
      gene_y = stringr::str_trim(gene_y)
    )
  tmp_r <- dfv[[rcol]]
  tmp_p <- dfv[[pcol]]
  dfv   <- dfv %>% dplyr::transmute(gene_x, gene_y, r = tmp_r, padj = tmp_p)
  
  present <- intersect(sub_order, unique(c(dfv$gene_x, dfv$gene_y)))
  if (length(present) == 0L) stop("此版本在該 dataset 無任何 CSN subunits")
  
  grid <- tidyr::expand_grid(gene_x = present, gene_y = present) %>%
    dplyr::mutate(
      key     = paste(gene_x, gene_y, sep = "|"),
      key_rev = paste(gene_y, gene_x, sep = "|")
    )
  
  dfv_key <- dfv %>% dplyr::mutate(key = paste(gene_x, gene_y, sep = "|")) %>%
    dplyr::select(key, r, padj)
  
  grid2 <- grid %>%
    dplyr::left_join(dfv_key, by = "key") %>%
    dplyr::left_join(dplyr::rename(dfv_key, key_rev = key, r2 = r, padj2 = padj), by = "key_rev") %>%
    dplyr::mutate(
      r    = dplyr::coalesce(r, r2, ifelse(gene_x == gene_y, 1, NA_real_)),
      padj = dplyr::coalesce(padj, padj2, ifelse(gene_x == gene_y, NA_real_, NA_real_))
    ) %>%
    dplyr::select(gene_x, gene_y, r, padj) %>%
    dplyr::mutate(
      gene_x = factor(gene_x, levels = present),
      gene_y = factor(gene_y, levels = present),
      signif = !is.na(padj) & padj < 0.05 & abs(r) < 0.999999
    )
  
  tri <- match.arg(tolower(triangle), c("full","lower","upper"))
  if (tri != "full") {
    ix <- as.integer(grid2$gene_x); iy <- as.integer(grid2$gene_y)
    if (tri == "lower") grid2 <- grid2[iy >= ix, , drop = FALSE]
    if (tri == "upper") grid2 <- grid2[iy <= ix, , drop = FALSE]
  }
  grid2
}





# 繪圖（不顯示，只存檔）
.plot_corr_heatmap_save <- function(df_grid, title, out_base,
                                    width = 6.5, height = 6.5, dpi = 600) {
  # 以藍-白-紅發散色，限制 -1~1；NA 淡灰
  p <- ggplot(df_grid, aes(x = gene_x, y = gene_y, fill = r)) +
    geom_tile(color = "white", size = 0.3, na.rm = FALSE) +
    # 顯著性黑點
    geom_point(data = subset(df_grid, signif), aes(x = gene_x, y = gene_y),
               inherit.aes = FALSE, shape = 16, size = 2.0, color = "black", alpha = 0.9) +
    scale_fill_gradient2(low = CELL_BLUE, mid = CELL_WHITE, high = CELL_RED,
                         midpoint = 0, na.value = "grey90", limits = c(-1, 1), oob = scales::squish,
                         name = "Pearson r") +
    scale_x_discrete(position = "top") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(angle = 45, hjust = 0),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(),
      legend.position = "right"
    )
  
  # 儲存四種格式
  ggsave(paste0(out_base, ".tiff"), p, width = width, height = height, dpi = dpi, compression = "lzw")
  ggsave(paste0(out_base, ".png"),  p, width = width, height = height, dpi = dpi)
  ggsave(paste0(out_base, ".jpg"),  p, width = width, height = height, dpi = dpi)
  ggsave(paste0(out_base, ".pdf"),  p, width = width, height = height)  # 向量化
}

# 找到所有 dataset 的 CSV 並繪圖
run_csn_corrcoef_plots <- function(out_root = "CSN_subunits_correlation_coefficient") {
  if (!dir.exists(out_root)) {
    message("[corrplot] 找不到根目錄：", out_root)
    return(invisible(FALSE))
  }
  csvs <- list.files(out_root, pattern = "^pairwise_correlations_all_versions\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(csvs)) {
    message("[corrplot] 找不到任何 pairwise_correlations_all_versions.csv")
    return(invisible(FALSE))
  }
  
  versions <- c("NoCovariate","RAW_covars","BatchAdj_covars")
  for (csv in csvs) {
    ds_dir <- dirname(csv)
    ds_id  <- basename(ds_dir)
    plots_dir <- file.path(ds_dir)  # 直接放在 dataset 資料夾內（符合你的要求）
    # 若你想集中在子資料夾，可改成：file.path(ds_dir, "corr_plots")
    
    for (ver in versions) {
      ## Pearson（保留原本檔名）
      grid_p <- try(.build_corr_grid(csv, ver, CSN_SUB_ORDER,
                                     triangle = getOption("csn_heatmap_triangle", "full"),
                                     method   = "pearson"), silent = TRUE)
      if (!inherits(grid_p, "try-error")) {
        title_p <- sprintf("%s — %s (Pearson r)", ds_id, ver)
        out_base_p <- file.path(plots_dir, sprintf("%s_corrcoef_%s", .safe_fs(ds_id), .safe_fs(ver)))
        try(.plot_corr_heatmap_save(grid_p, title_p, out_base_p), silent = FALSE)
      }
      
      ## Spearman（新增；檔名加 _spearman 後綴）
      grid_s <- try(.build_corr_grid(csv, ver, CSN_SUB_ORDER,
                                     triangle = getOption("csn_heatmap_triangle", "full"),
                                     method   = "spearman"), silent = TRUE)
      if (!inherits(grid_s, "try-error")) {
        title_s <- sprintf("%s — %s (Spearman \u03C1)", ds_id, ver)
        out_base_s <- file.path(plots_dir, sprintf("%s_corrcoef_%s_spearman", .safe_fs(ds_id), .safe_fs(ver)))
        try(.plot_corr_heatmap_save(grid_s, title_s, out_base_s), silent = FALSE)
      }
    }
    message("[corrplot] 完成：", ds_id)
  }
  invisible(TRUE)
}

## ---- 直接執行繪圖（可保留在腳本底部） ----
run_csn_corrcoef_plots(out_root = "CSN_subunits_correlation_coefficient")













## =========================================================
## 批次將單一 dataset/stratum 的 GSEA Summary CSV 繪製成 heatmap (heatmap_gsea_ggplot_V4.R)
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
  ## [NEW] 非 H collection：依 CSN_SCORE 的 NES 篩選前/後 N 條（需 padj<0.05）
  .coll_tok <- sub("^Summary_([^_]+)_GSEA.*$", "\\1", basename(csv_file))
  if (!identical(toupper(.coll_tok), "H") && "NES_CSN_SCORE" %in% present_preds) {
    .top_n <- get0("DATASET_HEATMAP_TOP_N", ifnotfound = 25)
    .bot_n <- get0("DATASET_HEATMAP_BOTTOM_N", ifnotfound = 25)
    csn_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "NES_CSN_SCORE",
                    is.finite(.data$padj), .data$padj < 0.05)
    keep_up <- csn_tbl %>%
      dplyr::arrange(dplyr::desc(.data$NES)) %>%
      dplyr::slice_head(n = .top_n) %>% dplyr::pull(.data$pathway)
    keep_dn <- csn_tbl %>%
      dplyr::arrange(.data$NES) %>%
      dplyr::slice_head(n = .bot_n) %>% dplyr::pull(.data$pathway)
    keep <- unique(c(keep_up, keep_dn))
    if (length(keep)) {
      df_long <- df_long %>% dplyr::filter(.data$pathway %in% keep)
    }
  }
  
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
    pattern = "^Summary_[^_]+_GSEA_(limma_t(_cont)?|limma_interaction|spearman)_(ALL|TP53_mutant|TP53_wild_type)\\.csv$",
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
    scan_dir <- if (is.null(COMBO_PREFIX)) ds_dir else COMBO_PREFIX
    csvs <- .find_target_csvs_for_dataset(scan_dir)
    ## [NEW] combo 層級下依 dataset 名稱過濾（也要求路徑含 /summary/）
    if (!is.null(COMBO_PREFIX)) {
      ds_name <- basename(ds_dir)
      sep_pat <- "[/\\\\]"
      pat_ds  <- paste0(sep_pat, ds_name, sep_pat)
      csvs <- csvs[grepl(pat_ds, csvs, ignore.case = TRUE)]
      csvs <- csvs[grepl(paste0(sep_pat, "summary", sep_pat), tolower(csvs))]
      # 若有指定 groups（例如只跑 H），路徑再限縮到這些群組
      .grps <- get0("GENESET_GROUPS_TO_RUN", ifnotfound = NULL)
      if (!is.null(.grps) && length(.grps)) {
        .gs  <- unique(vapply(.grps, safe_fs_name, FUN.VALUE = character(1)))
        .pat <- paste(paste0(sep_pat, tolower(.gs), sep_pat), collapse="|")
        csvs <- csvs[grepl(.pat, tolower(csvs))]
      }
    }
    ## [NEW] 篩選要畫的 collections（單一 dataset 熱圖）
    .cols_ds <- get0("PLOT_DATASET_COLLECTIONS", ifnotfound = NULL)
    if (!is.null(.cols_ds) && length(.cols_ds)) {
      .al <- unique(vapply(.cols_ds, safe_fs_name, FUN.VALUE = character(1)))
      .bas <- basename(csvs)
      .pref <- paste0("Summary_", .al, "_")
      .keep <- Reduce(`|`, lapply(.pref, function(p) startsWith(.bas, p)))
      csvs <- csvs[.keep]
    }
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
      
      ## === 修正：limma_interaction 一般版也套用綠–橘 Cell 配色 ===
      bn_lc <- tolower(basename(csv))
      is_interaction_now <- grepl("gsea_limma_interaction", bn_lc, fixed = TRUE)
      
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
          # 依「同 dataset + 同 variant + 同 collection（若有）」構造 ALL 對應路徑：
          # 直接把當前路徑中的 stratum 目錄（ALL/TP53_mutant/TP53_wild_type）改成 ALL
          parts    <- strsplit(normalizePath(csv, winslash = "/"), "/")[[1]]
          parts_lc <- tolower(parts)
          sidx <- which(parts_lc %in% c("all","tp53_mutant","tp53_wild_type"))
          if (length(sidx)) {
            parts2 <- parts
            parts2[sidx[length(sidx)]] <- "ALL"  # 將實際的 stratum 目錄改成 ALL
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
  collect_root = if (is.null(COMBO_PREFIX)) "single_dataset_GSEA_heatmap" else file.path(COMBO_PREFIX, "single_dataset_GSEA_heatmap")
)












## =========================================================
## [NEW] Summarize meta-FDR (Stouffer → BH) across subunits
##   - Input:  csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/<RAW|BatchAdj>/<SUBUNIT>/<GROUP>/GSEA_limma_t_cont_meta_fdr.csv
##   - Output: .../summary/<STRATUM>/<RAW|BatchAdj>/<GROUP>/GSEA_limma_t_cont_meta_fdr/
##              ├─ Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_ALL.csv
##              ├─ Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_padjLT0.05.csv
##              └─ Summary_<GROUP>_GSEA_limma_t_cont_meta_fdr_padjLT0.25.csv
##   - Columns merged per subunit: Z_<SUBUNIT>, padj_meta_<SUBUNIT>
##   - Counting uses padj_meta thresholds; direction by sign(Z).
## =========================================================

.read_meta_fdr_table <- function(meta_root, stratum, version, subunit,
                                 group_name, stat_tag = "GSEA_limma_t_cont_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("請先安裝 data.table")
  grp <- safe_fs_name(group_name)
  fp  <- file.path(meta_root, stratum, version, subunit, grp, paste0(stat_tag, ".csv"))
  if (!file.exists(fp)) return(NULL)
  dt <- tryCatch(data.table::fread(fp, na.strings = c("NA","NaN","")), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  need <- c("pathway","Z","padj_meta")
  # 寬鬆對應
  nms <- tolower(names(dt))
  if (!"pathway"   %in% names(dt) && "pathway"   %in% nms) names(dt)[match("pathway", nms)]   <- "pathway"
  if (!"Z"         %in% names(dt) && "z"         %in% nms) names(dt)[match("z", nms)]         <- "Z"
  if (!"padj_meta" %in% names(dt) && "padj_meta" %in% nms) names(dt)[match("padj_meta", nms)] <- "padj_meta"
  if (!all(need %in% names(dt))) return(NULL)
  out <- as.data.frame(dt[, ..need])
  names(out)[names(out) == "Z"]         <- paste0("Z_", subunit)
  names(out)[names(out) == "padj_meta"] <- paste0("padj_meta_", subunit)
  out
}

.merge_subunit_tables_meta <- function(tbl_list) {
  keep <- tbl_list[!vapply(tbl_list, is.null, logical(1))]
  if (!length(keep)) return(NULL)
  out <- Reduce(function(x, y) dplyr::full_join(x, y, by = "pathway"), keep)
  num_cols <- setdiff(names(out), "pathway")
  out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
  out
}

.add_sig_counts_meta <- function(df, alphas = c(0.05, 0.25)) {
  if (is.null(df) || !nrow(df)) return(df)
  padj_cols <- grep("^padj_meta_", names(df), value = TRUE)
  subs      <- sub("^padj_meta_", "", padj_cols)
  z_cols    <- paste0("Z_", subs)
  keep_idx  <- z_cols %in% names(df)
  padj_cols <- padj_cols[keep_idx]; z_cols <- z_cols[keep_idx]
  if (!length(padj_cols)) return(df)
  
  for (a in alphas) {
    a_tag <- gsub("\\.", "_", sprintf("%.2f", a))
    sig_mat <- sapply(padj_cols, function(p) { pv <- df[[p]]; as.integer(is.finite(pv) & pv < a) })
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    df[[sprintf("sig_n_padj_meta_%s", a_tag)]] <- rowSums(sig_mat, na.rm = TRUE)
    
    pos_mat <- mapply(function(p, zc) { pv <- df[[p]]; zv <- df[[zc]]
    as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv > 0) },
    padj_cols, z_cols, SIMPLIFY = TRUE)
    if (is.null(dim(pos_mat))) pos_mat <- matrix(pos_mat, ncol = 1)
    df[[sprintf("pos_n_padj_meta_%s", a_tag)]] <- rowSums(pos_mat, na.rm = TRUE)
    
    neg_mat <- mapply(function(p, zc) { pv <- df[[p]]; zv <- df[[zc]]
    as.integer(is.finite(pv) & pv < a & is.finite(zv) & zv < 0) },
    padj_cols, z_cols, SIMPLIFY = TRUE)
    if (is.null(dim(neg_mat))) neg_mat <- matrix(neg_mat, ncol = 1)
    df[[sprintf("neg_n_padj_meta_%s", a_tag)]] <- rowSums(neg_mat, na.rm = TRUE)
  }
  df
}

.write_summary_outputs_meta_csv <- function(df, out_dir, group_name,
                                            stat_tag = "GSEA_limma_t_cont_meta_fdr") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("請先安裝 data.table")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("請先安裝 openxlsx")
  if (is.null(df) || !nrow(df)) { 
    if (exists("log_msg", mode="function")) try(log_msg("  [略過輸出] {group_name} | {stat_tag} 無可用結果"), silent = TRUE)
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0("Summary_", safe_fs_name(group_name), "_", stat_tag))
  
  ## ---- 排序：優先顯示顯著數量（0.05）、再顯示正/負向顯著數 ----
  ord_keys <- c("sig_n_padj_meta_0_05","pos_n_padj_meta_0_05","neg_n_padj_meta_0_05")
  ord_keys <- intersect(ord_keys, names(df))
  if (length(ord_keys)) {
    df <- df |>
      dplyr::arrange(dplyr::desc(.data[[ord_keys[1]]]),
                     dplyr::desc(ifelse(length(ord_keys) > 1, .data[[ord_keys[2]]], 0)),
                     dplyr::desc(ifelse(length(ord_keys) > 2, .data[[ord_keys[3]]], 0)),
                     .data[["pathway"]])
  } else {
    df <- df |> dplyr::arrange(.data[["pathway"]])
  }
  
  ## ---- CSV 輸出（ALL / padj<0.05 / padj<0.25）----
  data.table::fwrite(df, paste0(base, "_ALL.csv"))
  
  df005 <- NULL; df025 <- NULL
  if ("sig_n_padj_meta_0_05" %in% names(df)) {
    df005 <- df |> dplyr::filter(.data[["sig_n_padj_meta_0_05"]] > 0)
    dir.create(dirname(paste0(base, "_padjLT0.05.csv")), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df005, paste0(base, "_padjLT0.05.csv"))
  }
  if ("sig_n_padj_meta_0_25" %in% names(df)) {
    df025 <- df |> dplyr::filter(.data[["sig_n_padj_meta_0_25"]] > 0)
    dir.create(dirname(paste0(base, "_padjLT0.25.csv")), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df025, paste0(base, "_padjLT0.25.csv"))
  }
  
  ## ---- XLSX 輸出（ALL / padjLT0.05 / padjLT0.25 三工作表）----
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "ALL")
  openxlsx::writeData(wb, "ALL", df)
  openxlsx::freezePane(wb, "ALL", firstRow = TRUE)
  openxlsx::setColWidths(wb, "ALL", cols = 1:ncol(df), widths = "auto")
  
  openxlsx::addWorksheet(wb, "padjLT0.05")
  if (!is.null(df005)) {
    openxlsx::writeData(wb, "padjLT0.05", df005)
    openxlsx::freezePane(wb, "padjLT0.05", firstRow = TRUE)
    openxlsx::setColWidths(wb, "padjLT0.05", cols = 1:ncol(df), widths = "auto")
  }
  
  openxlsx::addWorksheet(wb, "padjLT0.25")
  if (!is.null(df025)) {
    openxlsx::writeData(wb, "padjLT0.25", df025)
    openxlsx::freezePane(wb, "padjLT0.25", firstRow = TRUE)
    openxlsx::setColWidths(wb, "padjLT0.25", cols = 1:ncol(df), widths = "auto")
  }
  
  openxlsx::saveWorkbook(wb, paste0(base, ".xlsx"), overwrite = TRUE)
  
  if (exists("log_msg", mode="function")) {
    try(log_msg("  [完成輸出] {group_name} | {stat_tag} -> {dirname(base)} (.csv + .xlsx)"), silent = TRUE)
  }
  invisible(NULL)
}


summarize_meta_fdr_across_subunits <- function(
    meta_root = if (is.null(COMBO_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata    = c("ALL","TP53_mutant","TP53_wild_type"),
    versions  = c("RAW","BatchAdj"),
    genesets_by_group = genesets_by_group,
    stat_tag  = "GSEA_limma_t_cont_meta_fdr"
){
  if (!requireNamespace("data.table", quietly = TRUE)) stop("請先安裝 data.table")
  for (st in strata) for (ver in versions) {
    base_dir <- file.path(meta_root, st, ver)
    if (!dir.exists(base_dir)) next
    # 自動偵測 subunit 清單：抓該層的第一層子目錄名稱
    subs <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
    subs <- subs[nzchar(subs)]
    if (!length(subs)) next
    
    for (grp in names(genesets_by_group)) {
      if (exists("log_msg", mode="function")) try(log_msg("== [meta-summary] stratum=%s | version=%s | group=%s ==", st, ver, grp), silent = TRUE)
      lst <- setNames(vector("list", length(subs)), subs)
      for (su in subs) lst[[su]] <- .read_meta_fdr_table(meta_root, st, ver, su, grp, stat_tag)
      wide <- .merge_subunit_tables_meta(lst)
      wide <- .add_sig_counts_meta(wide, alphas = c(0.05, 0.25))
      out_dir <- file.path(meta_root, "summary", st, ver, safe_fs_name(grp), stat_tag)
      .write_summary_outputs_meta_csv(wide, out_dir, grp, stat_tag)
    }
  }
  invisible(TRUE)
}

## ---- 一鍵執行（不重跑 GSEA；僅讀既有 *_meta_fdr.csv 後彙整）----
posthoc_summary_meta_fdr <- function(){
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(COMBO_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata           = c("ALL","TP53_mutant","TP53_wild_type"),
    versions         = c("RAW","BatchAdj"),
    genesets_by_group = genesets_by_group,
    stat_tag         = "GSEA_limma_t_cont_meta_fdr"
  )
  invisible(TRUE)
}

## 直接跑一次（可重複執行；會覆蓋同名 csv）
posthoc_summary_meta_fdr()




## =========================================================
## [NEW] Summarize meta-FDR for limma interaction across subunits
##   - Input:
##     csn_gsea_pan_summary_TP53/meta_fdr/<STRATUM>/<RAW|BatchAdj>/<SUBUNIT>/<GROUP>/GSEA_limma_interaction_meta_fdr.csv
##   - Output:
##     csn_gsea_pan_summary_TP53/meta_fdr/summary/<STRATUM>/<RAW|BatchAdj>/<GROUP>/GSEA_limma_interaction_meta_fdr/
##       ├─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr_ALL.csv
##       ├─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr_padjLT0.05.csv
##       ├─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr_padjLT0.25.csv
##       └─ Summary_<GROUP>_GSEA_limma_interaction_meta_fdr.xlsx
##   - 不重跑 GSEA；沿用你先前新增的 summarize_meta_fdr_across_subunits() 等工具。
## =========================================================

posthoc_summary_meta_fdr_interaction <- function(){
  summarize_meta_fdr_across_subunits(
    meta_root = if (is.null(COMBO_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr" else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr"),
    strata            = c("ALL","TP53_mutant","TP53_wild_type"),
    versions          = c("RAW","BatchAdj"),
    genesets_by_group = genesets_by_group,
    stat_tag          = "GSEA_limma_interaction_meta_fdr"
  )
  invisible(TRUE)
}

## ---- 直接執行一次（可重複執行；只會覆蓋輸出檔，不會重跑GSEA）----
posthoc_summary_meta_fdr_interaction()
















## =========================================================
## [NEW] Heatmaps for meta-FDR summary CSVs (Z + padj_meta)
##   - Inputs (四種):
##     1) ALL   : Summary_H_GSEA_limma_t_cont_meta_fdr_ALL.csv
##     2) ALL   : Summary_H_GSEA_limma_interaction_meta_fdr_ALL.csv
##     3) TP53_mutant    : Summary_H_GSEA_limma_t_cont_meta_fdr_ALL.csv
##     4) TP53_wild_type : Summary_H_GSEA_limma_t_cont_meta_fdr_ALL.csv
##   - Fill: Z_<subunit>; Dot: padj_meta_<subunit> < 0.05
##   - 输出：和輸入同資料夾；檔名加 prefix = heatmap_<STRATUM>_<VERSION>_
##   - ordered 版本：
##       a) (ALL) interaction  → 依「同版本」 ALL 的 limma_t_cont_meta_fdr 的 y-order
##       b) (TP53_mutant / TP53_wild_type) t_cont  → 依「同版本」 ALL 的 limma_t_cont_meta_fdr 的 y-order
## =========================================================

## ===== 安全載入舊版工具函式（若未定義才建立）=====

if (!exists(".read_csv_safe", mode = "function")) {
  .read_csv_safe <- function(path) {
    p <- path
    # Windows: 若原字串不存在，嘗試 / -> \
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
    
    # 3) 讀成全文字串再解析
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
}

if (!exists(".safe_fs", mode = "function")) {
  .safe_fs <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

if (!exists(".pred_order_all", mode = "any")) {
  # 目標 X 軸順序（沿用你原腳本）
  .pred_order_all <- c(
    "NES_CSN_SCORE", "NES_GPS1",
    "NES_COPS2", "NES_COPS3", "NES_COPS4", "NES_COPS5", "NES_COPS6",
    "NES_COPS7A", "NES_COPS7B", "NES_COPS8", "NES_COPS9",
    "NES_RESIDUAL_GPS1",
    "NES_RESIDUAL_COPS2", "NES_RESIDUAL_COPS3", "NES_RESIDUAL_COPS4", "NES_RESIDUAL_COPS5", "NES_RESIDUAL_COPS6",
    "NES_RESIDUAL_COPS7A", "NES_RESIDUAL_COPS7B", "NES_RESIDUAL_COPS8", "NES_RESIDUAL_COPS9"
  )
}


# Z / padj_meta 的 predictor 順序（沿用你原本 NES_* 的順序）
.pred_order_meta <- gsub("^NES_", "Z_", .pred_order_all)

# 把 Summary_*_meta_fdr_ALL.csv 讀成「長表」(pathway, predictor, Z, padj_meta)
.meta_read_summary_long <- function(csv_file) {
  df <- suppressMessages(.read_csv_safe(csv_file))
  # 自動找 pathway 欄位
  path_candidates <- c("pathway","Pathway","term","Term","gs_name","NAME","set","Set")
  path_col <- intersect(path_candidates, names(df))[1]
  if (is.na(path_col)) stop("找不到 pathway 欄位於: ", csv_file)
  
  z_cols     <- grep("^Z_",         names(df), value = TRUE)
  padj_cols  <- grep("^padj_meta_", names(df), value = TRUE)
  if (!length(z_cols)) stop("找不到 Z_* 欄位於: ", csv_file)
  
  df <- dplyr::rename(df, pathway = dplyr::all_of(path_col))
  z_map     <- tibble::tibble(z_col = z_cols, key = sub("^Z_", "", z_cols))
  padj_map  <- tibble::tibble(p_col = padj_cols, key = sub("^padj_meta_", "", padj_cols))
  pair_map  <- dplyr::left_join(z_map, padj_map, by = "key")
  
  lst <- lapply(seq_len(nrow(pair_map)), function(i){
    zc <- pair_map$z_col[i]
    pc <- pair_map$p_col[i]
    tibble::tibble(
      pathway   = df$pathway,
      predictor = paste0("Z_", pair_map$key[i]),
      Z         = suppressWarnings(as.numeric(df[[zc]])),
      padj_meta = if (!is.na(pc)) suppressWarnings(as.numeric(df[[pc]])) else NA_real_
    )
  })
  dplyr::bind_rows(lst)
}

# 以 meta 檔（Z）計算 y 軸順序：優先用 Z_CSN_SCORE 降冪，否則用行平均 Z
.meta_compute_y_order <- function(csv_file) {
  df_long <- .meta_read_summary_long(csv_file)
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("無可用 predictors 於: ", csv_file)
  
  if ("Z_CSN_SCORE" %in% present) {
    ord <- df_long %>%
      dplyr::filter(.data$predictor == "Z_CSN_SCORE") %>%
      dplyr::select(.data$pathway, .data$Z) %>% dplyr::distinct() %>%
      dplyr::arrange(dplyr::desc(.data$Z)) %>% dplyr::pull(.data$pathway)
  } else {
    ord <- df_long %>% dplyr::group_by(.data$pathway) %>%
      dplyr::summarise(m = mean(.data$Z, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$m)) %>% dplyr::pull(.data$pathway)
  }
  ord
}

# 建 meta-FDR heatmap（Z 做色階；padj_meta<0.05 打黑點）
.meta_make_heatmap_plot <- function(csv_file, y_order = NULL,
                                    palette = c(low = "#053061", mid = "#FFFFFF", high = "#67001F")) {
  df_long <- .meta_read_summary_long(csv_file)
  
  present <- intersect(.pred_order_meta, unique(df_long$predictor))
  if (!length(present)) stop("資料中沒有任何指定的 predictors：", csv_file)
  df_long <- dplyr::filter(df_long, .data$predictor %in% present)
  ## [NEW] 非 H collection：依 CSN_SCORE 的 Z 篩選前/後 N 條（需 padj_meta<0.05）
  .coll_tok <- sub("^Summary_([^_]+)_GSEA.*$", "\\1", basename(csv_file))
  if (!identical(toupper(.coll_tok), "H") && "Z_CSN_SCORE" %in% present) {
    .top_n <- get0("PAN_HEATMAP_TOP_N", ifnotfound = 25)
    .bot_n <- get0("PAN_HEATMAP_BOTTOM_N", ifnotfound = 25)
    csn_tbl <- df_long %>%
      dplyr::filter(.data$predictor == "Z_CSN_SCORE",
                    is.finite(.data$padj_meta), .data$padj_meta < 0.05)
    keep_up <- csn_tbl %>%
      dplyr::arrange(dplyr::desc(.data$Z)) %>%
      dplyr::slice_head(n = .top_n) %>% dplyr::pull(.data$pathway)
    keep_dn <- csn_tbl %>%
      dplyr::arrange(.data$Z) %>%
      dplyr::slice_head(n = .bot_n) %>% dplyr::pull(.data$pathway)
    keep <- unique(c(keep_up, keep_dn))
    if (length(keep)) {
      df_long <- dplyr::filter(df_long, .data$pathway %in% keep)
    }
  }
  
  # y-order：若外部有給就用，否則依 Z_CSN_SCORE 或平均 Z
  if (is.null(y_order)) {
    y_order <- .meta_compute_y_order(csv_file)
  }
  y_use <- intersect(y_order, unique(df_long$pathway))
  if (!length(y_use)) stop("y-order 與資料無交集：", csv_file)
  df_long <- dplyr::mutate(df_long, pathway = factor(.data$pathway, levels = rev(y_use)))
  
  # X 軸位置與 gap 規則（與你原規則一致；改用 Z_* 名稱）
  gap <- 0.4
  needs_gap1 <- all(c("Z_CSN_SCORE","Z_GPS1") %in% present)
  needs_gap2 <- "Z_RESIDUAL_GPS1" %in% present
  
  pos_map <- list(); pos <- 0
  for (p in present) {
    if (p == "Z_GPS1" && needs_gap1) pos <- pos + gap
    if (p == "Z_RESIDUAL_GPS1" && needs_gap2) pos <- pos + gap
    pos <- pos + 1
    pos_map[[p]] <- pos
  }
  pos_map <- unlist(pos_map)
  df_long <- dplyr::mutate(df_long, xpos = unname(pos_map[.data$predictor]))
  
  L <- max(abs(df_long$Z), na.rm = TRUE); if (!is.finite(L) || L == 0) L <- 1
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = xpos, y = pathway, fill = Z)) +
    ggplot2::geom_tile(width = 1, height = 0.9, color = NA) +
    ggplot2::geom_point(data = dplyr::filter(df_long, is.finite(.data$padj_meta), .data$padj_meta < 0.05),
                        ggplot2::aes(x = xpos, y = pathway),
                        shape = 16, size = 1.6, color = "black", inherit.aes = FALSE) +
    ggplot2::scale_fill_gradient2(low = palette[["low"]], mid = palette[["mid"]], high = palette[["high"]],
                                  limits = c(-L, L), midpoint = 0, oob = scales::squish, name = "Z") +
    ggplot2::scale_x_continuous(
      breaks = unname(pos_map[present]),
      labels = present,
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
  W <- max(8, length(present) * 0.45)
  H <- max(6, n_path * 0.22)
  list(plot = p, width = W, height = H)
}

# 解析 meta_fdr/summary 路徑，取得 version(=RAW/BatchAdj) 與 stratum
.meta_parse_summary_path <- function(csv_file) {
  parts    <- strsplit(normalizePath(csv_file, winslash = "/"), "/")[[1]]
  parts_lc <- tolower(parts)
  ver_idx  <- which(parts_lc %in% c("raw","batchadj"))
  version  <- if (length(ver_idx)) parts[ver_idx[1]] else NA_character_
  st_idx   <- which(parts_lc %in% c("all","tp53_mutant","tp53_wild_type"))
  stratum  <- if (length(st_idx)) parts[st_idx[length(st_idx)]] else NA_character_
  list(version = version, stratum = stratum)
}

# 儲存到「與輸入相同的資料夾」，檔名前加 prefix：heatmap_<STRATUM>_<VERSION>_
.meta_save_near <- function(p, width, height, csv_file, meta, suffix = "", dpi = 600) {
  out_dir <- dirname(csv_file)
  bn      <- basename(csv_file)            # 原檔名
  pre     <- paste0("heatmap_", .safe_fs(meta$stratum), "_", .safe_fs(meta$version), "_")
  base    <- file.path(out_dir, paste0(pre, bn, suffix))
  ggplot2::ggsave(paste0(base, ".tiff"), p, width = width, height = height, units = "in",
                  dpi = dpi, bg = "white", compression = "lzw")
  # ggplot2::ggsave(paste0(base, ".png"),  p, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  # ggplot2::ggsave(paste0(base, ".jpg"),  p, width = width, height = height, units = "in", dpi = dpi, bg = "white", quality = 100)
  # ggplot2::ggsave(paste0(base, ".pdf"),  p, width = width, height = height, units = "in", bg = "white")
}

# 找到四種目標 CSV（只抓 group=H 的 ALL / MUT / WT，且檔名完全匹配）
.meta_find_target_csvs <- function(root = "csn_gsea_pan_summary_TP53/meta_fdr/summary") {
  if (!dir.exists(root)) return(character(0))
  pat_all_t <- "^Summary_[^_]+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  pat_all_i <- "^Summary_[^_]+_GSEA_limma_interaction_meta_fdr_ALL\\.csv$"
  pat_mut_t   <- "^Summary_[^_]+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  pat_wt_t    <- "^Summary_[^_]+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$"
  
  files <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # 仍要求路徑含 RAW/BatchAdj
  keep <- grepl("(/|\\\\)(RAW)(/|\\\\)", files, ignore.case = TRUE) |
    grepl("(/|\\\\)(BatchAdj)(/|\\\\)", files, ignore.case = TRUE)
  files <- files[keep]
  if (!length(files)) return(files)
  
  bas <- basename(files)
  is_all   <- grepl("(/|\\\\)ALL(/|\\\\)", files, ignore.case = TRUE)
  is_mut   <- grepl("(/|\\\\)TP53_mutant(/|\\\\)", files, ignore.case = TRUE)
  is_wt    <- grepl("(/|\\\\)TP53_wild_type(/|\\\\)", files, ignore.case = TRUE)
  
  tgt <- c(
    files[ is_all & grepl(pat_all_t, bas) ],
    files[ is_all & grepl(pat_all_i, bas) ],
    files[ is_mut & grepl(pat_mut_t, bas) ],
    files[ is_wt  & grepl(pat_wt_t,  bas) ]
  )
  
  # 剔除 size=0
  if (length(tgt)) {
    finfo <- suppressWarnings(file.info(tgt))
    tgt <- tgt[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
  }
  unique(tgt)
}

# 若尚未定義，使用舊版的掃描函式
if (!exists(".find_target_csvs_for_dataset", mode = "function")) {
  .find_target_csvs_for_dataset <- function(ds_dir) {
    # 搜尋名稱匹配的 Summary 檔
    files <- list.files(
      ds_dir,
      pattern = "^Summary_[^_]+_GSEA_(limma_t(_cont)?|limma_interaction|spearman)_(ALL|TP53_mutant|TP53_wild_type)\\.csv$",
      recursive = TRUE, full.names = TRUE
    )
    if (!length(files)) return(files)
    
    # 仍要求路徑中要有 RAW 或 BatchAdj（同時接受 / 或 \）
    sep_pat <- "[/\\\\]"  # 同時匹配 / 與 \
    keep <- grepl(paste0(sep_pat, "raw", sep_pat), files, ignore.case = TRUE) |
      grepl(paste0(sep_pat, "batchadj", sep_pat), files, ignore.case = TRUE)
    files <- files[keep]
    
    # 剔除 size=0（尚未寫完或空檔）
    if (length(files)) {
      finfo <- suppressWarnings(file.info(files))
      files <- files[is.finite(finfo$size) & finfo$size > 0 & !finfo$isdir]
    }
    files
  }
}

# 若尚未定義，使用舊版的 .parse_meta_from_path()
if (!exists(".parse_meta_from_path", mode = "function")) {
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
}



# 主程式：對找到的 meta_fdr Summary 檔畫圖，並產生 ordered 版本（依 ALL/t_cont）
run_meta_fdr_heatmaps <- function(root = if (is.null(COMBO_PREFIX)) "csn_gsea_pan_summary_TP53/meta_fdr/summary" else file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr/summary")) {
  csvs <- .meta_find_target_csvs(root)
  ## [NEW] 篩選要畫的 collections（PAN 熱圖）
  .cols_pan <- get0("PLOT_PAN_COLLECTIONS", ifnotfound = NULL)
  if (!is.null(.cols_pan) && length(.cols_pan)) {
    .al <- unique(vapply(.cols_pan, safe_fs_name, FUN.VALUE = character(1)))
    .bas <- basename(csvs)
    .pref <- paste0("Summary_", .al, "_")
    .keep <- Reduce(`|`, lapply(.pref, function(p) startsWith(.bas, p)))
    csvs <- csvs[.keep]
  }
  if (!length(csvs)) {
    message("[meta-heatmap] 找不到目標 Summary_*_meta_fdr_ALL.csv，略過。")
    return(invisible(TRUE))
  }
  for (csv in csvs) {
    meta  <- .meta_parse_summary_path(csv)
    bn_lc <- tolower(basename(csv))
    is_interaction <- grepl("gsea_limma_interaction_meta_fdr", bn_lc, fixed = TRUE)
    
    ## 若是 interaction，改用綠–橘 Cell 調色
    pal_cell <- NULL
    if (is_interaction) {
      pal0 <- .get_interaction_palette()  # c(neg=..., mid=..., pos=...)
      pal_cell <- c(low = pal0[["neg"]], mid = pal0[["mid"]], high = pal0[["pos"]])
    }
    
    # 先畫「一般版」
    h <- try({
      if (is_interaction) {
        .meta_make_heatmap_plot(csv, y_order = NULL, palette = pal_cell)
      } else {
        .meta_make_heatmap_plot(csv)  # 其它維持原藍紅
      }
    }, silent = TRUE)
    
    if (inherits(h, "try-error")) {
      message("[meta-heatmap] 建圖失敗：", csv, " | ", as.character(h))
      next
    }
    .meta_save_near(h$plot, h$width, h$height, csv, meta, suffix = "")
    
    # ----- 檢查是否要畫 ordered 版 -----
    # Case A：ALL 的 interaction → 依同版本 ALL 的 limma_t_cont_meta_fdr 的順序
    is_interaction_all <- grepl("gsea_limma_interaction_meta_fdr_all\\.csv$", bn_lc)
    if (is_interaction_all && identical(tolower(meta$stratum), "all")) {
      dir_now    <- normalizePath(dirname(csv), winslash = "/", mustWork = FALSE)
      parent_dir <- dirname(dir_now)  # e.g. .../meta_fdr/summary/ALL/RAW/H
      cont_dir   <- file.path(parent_dir, "GSEA_limma_t_cont_meta_fdr")
      cont_csvs <- list.files(cont_dir, pattern = "^Summary_[^_]+_GSEA_limma_t_cont_meta_fdr_ALL\\.csv$", full.names = TRUE)
      cont_csv  <- if (length(cont_csvs)) cont_csvs[1L] else ""
      if (file.exists(cont_csv)) {
        yo <- try(.meta_compute_y_order(cont_csv), silent = TRUE)
        if (!inherits(yo, "try-error")) {
          h2 <- try(.meta_make_heatmap_plot(csv, y_order = yo, palette = pal_cell), silent = TRUE)
          if (!inherits(h2, "try-error")) {
            .meta_save_near(h2$plot, h2$width, h2$height, csv, meta, suffix = "_ordered")
          } else {
            message("[meta-heatmap][ordered] 建圖失敗：", csv, " | ", as.character(h2))
          }
        } else {
          message("[meta-heatmap][ordered] 取 y-order 失敗（ALL/t_cont）：", cont_csv, " | ", as.character(yo))
        }
      } else {
        message("[meta-heatmap][ordered] 找不到 ALL/t_cont 檔案：", cont_csv)
      }
    }
    
    # Case B：TP53_mutant / TP53_wild_type 的 t_cont → 依同版本 ALL 的 t_cont
    is_tcont_all <- grepl("gsea_limma_t_cont_meta_fdr_all\\.csv$", bn_lc)
    if (is_tcont_all && tolower(meta$stratum) %in% c("tp53_mutant","tp53_wild_type")) {
      # 把路徑中的 stratum 改成 ALL
      parts    <- strsplit(normalizePath(csv, winslash = "/"), "/")[[1]]
      parts_lc <- tolower(parts)
      st_idx   <- which(parts_lc %in% c("all","tp53_mutant","tp53_wild_type"))
      if (length(st_idx)) {
        parts2 <- parts; parts2[st_idx[length(st_idx)]] <- "ALL"
        all_csv <- paste(parts2, collapse = "/")
        if (file.exists(all_csv)) {
          yo <- try(.meta_compute_y_order(all_csv), silent = TRUE)
          if (!inherits(yo, "try-error")) {
            h3 <- try(.meta_make_heatmap_plot(csv, y_order = yo), silent = TRUE)  # t_cont 照舊藍紅
            if (!inherits(h3, "try-error")) {
              .meta_save_near(h3$plot, h3$width, h3$height, csv, meta, suffix = "_ordered")
            } else {
              message("[meta-heatmap][ordered] 建圖失敗：", csv, " | ", as.character(h3))
            }
          } else {
            message("[meta-heatmap][ordered] 取 ALL/t_cont 的 y-order 失敗：", all_csv, " | ", as.character(yo))
          }
        } else {
          message("[meta-heatmap][ordered] 找不到對應 ALL/t_cont 檔案：", all_csv)
        }
      }
    }
    
    try(closeAllConnections(), silent = TRUE)
    invisible(gc(FALSE))
    message("[meta-heatmap] 完成：", csv)
  }
  invisible(TRUE)
}



## ---- 執行（例）----
run_meta_fdr_heatmaps(
  root = if (is.null(COMBO_PREFIX)) 
    "csn_gsea_pan_summary_TP53/meta_fdr/summary" 
  else 
    file.path(COMBO_PREFIX, "csn_gsea_pan_summary_TP53/meta_fdr/summary")
)


