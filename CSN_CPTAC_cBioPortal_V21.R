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

## 建立 run_info 目錄（寫任何 run_info/* 檔前）
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)


## 附上多重假設層級說明（6）
cat(paste(
  "Multiple-testing policy:",
  " - Within each gene-set group (H/C2-KEGG-legacy/C3-TFT-legacy/C6) and per statistic",
  "   (limma t-rank, Spearman rho-rank), we control FDR via fgsea padj.",
  " - Pan-cancer summaries aggregate directions and counts across datasets/strata",
  "   (descriptive; no extra global FDR).",
  sep = "\n"),
  file = file.path("run_info","analysis_notes.txt"))


## ===== 參數 =====
set.seed(1234)

csn_subunits <- c("GPS1","COPS2","COPS3","COPS4","COPS5","COPS6","COPS7A","COPS7B","COPS8","COPS9")

min_frac_complete   <- 0.75     # 每基因至少 75% 非 NA
hi_lo_quantile      <- 0.25     # 高低各 25%
minSize <- 15; maxSize <- 500
fgsea_eps <- 1e-10              # fgseaMultilevel() 精度

min_per_group      <- 8         # High/Low 各至少多少樣本（不足就跳過 limma）
min_pairs_spearman <- 10        # Spearman 至少配對樣本數
MAKE_PLOTS <- FALSE

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


## ===== Gene set 選擇（預設只跑 Hallmark）=====
## 可選值： "H", "C2_KEGG", "C3_TFT", "C6"
GENESET_GROUPS_TO_RUN <- c("H")

## 也可以像這樣同時跑多組（示例）：
## GENESET_GROUPS_TO_RUN <- c("H", "C6")
## GENESET_GROUPS_TO_RUN <- c("H", "C2_KEGG", "C3_TFT", "C6")


## ===== Gene sets（只建置被選到的群組）=====
log_msg("準備 MSigDB gene sets；將依 GENESET_GROUPS_TO_RUN 建置：%s",
        paste(GENESET_GROUPS_TO_RUN, collapse = ", "))

genesets_by_group <- list()
.get_subcat_col <- function(df){
  cand <- c("gs_subcat","subcategory","sub_category","gs_subcategory")
  cand[cand %in% names(df)][1]
}

# HALLMARK
if ("H" %in% GENESET_GROUPS_TO_RUN) {
  df_H <- msigdbr(species = "Homo sapiens", category = "H")
  if (nrow(df_H) > 0) {
    genesets_by_group[["H"]] <- lapply(split(df_H$gene_symbol, df_H$gs_name), unique)
    log_msg("  取得 Hallmark：%d 個集合", length(genesets_by_group[["H"]]))
  } else {
    log_msg("  注意：msigdbr 版本中抓不到 Hallmark（H）")
  }
}

# C2 KEGG legacy（有些 msigdbr 版本可能無法提供；已做保護）
if ("C2_KEGG" %in% GENESET_GROUPS_TO_RUN) {
  df_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
  if (nrow(df_C2) > 0) {
    subcol <- .get_subcat_col(df_C2); subtxt <- if (!is.null(subcol)) df_C2[[subcol]] else NA_character_
    keep <- grepl("KEGG", paste(df_C2$gs_name, subtxt), ignore.case = TRUE)
    df_C2_ke <- df_C2[keep, , drop = FALSE]
    if (nrow(df_C2_ke) > 0) {
      grp <- if (!is.null(subcol) && all(!is.na(df_C2_ke[[subcol]]))) paste0("C2__", unique(df_C2_ke[[subcol]])[1]) else "C2__KEGG"
      genesets_by_group[[grp]] <- lapply(split(df_C2_ke$gene_symbol, df_C2_ke$gs_name), unique)
      log_msg("  取得 %s：%d 個集合", grp, length(genesets_by_group[[grp]]))
    } else {
      log_msg("  注意：此 msigdbr 版本找不到 C2-KEGG（可能因授權變動），略過")
    }
  } else {
    log_msg("  注意：此 msigdbr 版本無 C2 類別，略過 C2-KEGG")
  }
}

# C3 TFT legacy（TRANSFAC/JASPAR legacy；若無 legacy 就退而求其次用 TFT）
if ("C3_TFT" %in% GENESET_GROUPS_TO_RUN) {
  df_C3 <- msigdbr(species = "Homo sapiens", category = "C3")
  if (nrow(df_C3) > 0) {
    subcol <- .get_subcat_col(df_C3); subtxt <- if (!is.null(subcol)) df_C3[[subcol]] else NA_character_
    keep_tft <- grepl("TFT", subtxt, ignore.case = TRUE)
    keep_leg <- grepl("TRANSFAC|JASPAR|LEGACY", paste(subtxt, df_C3$gs_name), ignore.case = TRUE)
    df_C3_tft_legacy <- df_C3[keep_tft & keep_leg, , drop = FALSE]
    if (nrow(df_C3_tft_legacy) == 0) df_C3_tft_legacy <- df_C3[keep_tft, , drop = FALSE]
    if (nrow(df_C3_tft_legacy) > 0) {
      grp <- if (!is.null(subcol) && all(!is.na(df_C3_tft_legacy[[subcol]]))) paste0("C3__", unique(df_C3_tft_legacy[[subcol]])[1]) else "C3__TFT"
      genesets_by_group[[grp]] <- lapply(split(df_C3_tft_legacy$gene_symbol, df_C3_tft_legacy$gs_name), unique)
      log_msg("  取得 %s：%d 個集合", grp, length(genesets_by_group[[grp]]))
    } else {
      log_msg("  注意：此 msigdbr 版本抓不到 C3-TFT（legacy 或一般），略過")
    }
  } else {
    log_msg("  注意：此 msigdbr 版本無 C3 類別，略過 C3-TFT")
  }
}

# C6
if ("C6" %in% GENESET_GROUPS_TO_RUN) {
  df_C6 <- msigdbr(species = "Homo sapiens", category = "C6")
  if (nrow(df_C6) > 0) {
    genesets_by_group[["C6"]] <- lapply(split(df_C6$gene_symbol, df_C6$gs_name), unique)
    log_msg("  取得 C6：%d 個集合", length(genesets_by_group[["C6"]]))
  } else {
    log_msg("  注意：msigdbr 版本中抓不到 C6")
  }
}

if (!length(genesets_by_group)) stop("沒有可用的 gene-set（請設定 GENESET_GROUPS_TO_RUN，例如 GENESET_GROUPS_TO_RUN <- c('H')）")
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
audit_sva_sanity <- function(M, sv, batch_factor = NULL,
                             purity = NULL, sex = NULL, age = NULL,
                             y_bio = NULL, sample_order = NULL, label = "") {
  tag <- if (nzchar(label)) paste0(": ", label) else ""
  if (is.null(sv)) { log_msg("[SVA sanity%s] SV=NULL, skip", tag); return(invisible(NULL)) }
  
  M <- as.matrix(M)
  nm <- if (is.null(sample_order)) colnames(M) else sample_order
  
  # ---------- 對齊工具 ----------
  align_vec <- function(v, nm) {
    if (is.null(v)) return(NULL)
    v <- if (is.matrix(v) && ncol(v) == 1) v[,1] else v
    if (is.factor(v)) v <- as.character(v)
    if (is.null(names(v))) {
      if (length(v) == length(nm)) {
        names(v) <- nm
        log_msg("[SVA sanity%s] unnamed vector -> assume same order as M", tag)
      } else {
        log_msg("[SVA sanity%s] WARN: length mismatch for unnamed vector; drop", tag)
        return(NULL)
      }
    }
    out <- suppressWarnings(v[nm])
    out
  }
  align_fac <- function(f, nm) {
    if (is.null(f)) return(NULL)
    if (!is.factor(f)) f <- factor(f)
    if (is.null(names(f))) {
      if (length(f) == length(nm)) {
        names(f) <- nm
        log_msg("[SVA sanity%s] unnamed factor -> assume same order as M", tag)
      } else {
        log_msg("[SVA sanity%s] WARN: length mismatch for factor; drop", tag)
        return(NULL)
      }
    }
    droplevels(f[nm])
  }
  
  sv <- as.matrix(sv); rownames(sv) <- nm
  batch_factor <- align_fac(batch_factor, nm)
  purity <- suppressWarnings(as.numeric(align_vec(purity, nm)))
  sex    <- suppressWarnings(as.numeric(align_vec(sex, nm)))
  age    <- suppressWarnings(as.numeric(align_vec(age, nm)))
  y_bio  <- suppressWarnings(as.numeric(align_vec(y_bio, nm)))
  
  # ---------- (1) SV 與 batch 的關聯 ----------
  if (!is.null(batch_factor)) {
    mm_b <- model.matrix(~ 0 + batch_factor)
    cb   <- suppressWarnings(cor(sv, mm_b, use = "pairwise.complete.obs"))
    if (is.matrix(cb)) {
      mx <- round(apply(abs(cb), 1, max), 3)
      log_msg("[SVA sanity%s] max|cor(SV, batch dummies)| = %s", tag, paste(mx, collapse=", "))
    }
  }
  
  # ---------- (2) SV 與 purity/age/sex/biology 的 |rho| ----------
  sp <- function(a,b) suppressWarnings(cor(a, b, method="spearman", use="pairwise.complete.obs"))
  show_vec <- function(xx) paste(round(abs(xx), 3), collapse=", ")
  
  if (!all(is.na(purity))) {
    r <- apply(sv, 2, sp, purity); log_msg("[SVA sanity%s] |rho(SV, purity)| = %s", tag, show_vec(r))
  }
  if (!all(is.na(age))) {
    r <- apply(sv, 2, sp, age);    log_msg("[SVA sanity%s] |rho(SV, age)|    = %s", tag, show_vec(r))
  }
  if (!all(is.na(sex))) {
    r <- apply(sv, 2, sp, sex);    log_msg("[SVA sanity%s] |rho(SV, sex)|    = %s", tag, show_vec(r))
  }
  if (!all(is.na(y_bio))) {
    r <- apply(sv, 2, sp, y_bio);  log_msg("[SVA sanity%s] |rho(SV, biology)|= %s", tag, show_vec(r))
  }
  
  # ---------- (3) 殘差 PCA（加不加 SV 的對照） ----------
  build_design <- function(with_sv = FALSE) {
    parts <- list("(Intercept)" = rep(1, length(nm)))
    if (!all(is.na(y_bio))) parts$bio <- y_bio
    if (!is.null(batch_factor)) parts$batch <- model.matrix(~ 0 + batch_factor)
    if (!all(is.na(purity))) parts$purity <- purity
    if (!all(is.na(age)))    parts$age    <- age
    if (!all(is.na(sex)))    parts$sex    <- sex
    if (with_sv)             parts$sv     <- sv
    X <- do.call(cbind, parts)
    # 丟掉常數/全 NA 欄
    keep <- apply(X, 2, function(z) sd(z[is.finite(z)]) > 0)
    as.matrix(X[, keep, drop = FALSE])
  }
  
  safe_pca_top5 <- function(res) {
    pc <- tryCatch(stats::prcomp(res, center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (is.null(pc)) return("NA")
    ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
    paste(head(ve, 5), collapse = ", ")
  }
  
  X0 <- build_design(FALSE)
  X1 <- build_design(TRUE)
  Y  <- t(M[, nm, drop = FALSE])
  fit0 <- tryCatch(stats::lm.fit(x = X0, y = Y), error = function(e) NULL)
  fit1 <- tryCatch(stats::lm.fit(x = X1, y = Y), error = function(e) NULL)
  if (!is.null(fit0) && !is.null(fit1)) {
    r0 <- fit0$residuals; r1 <- fit1$residuals
    log_msg("[SVA sanity%s] residual var%% top5 (noSV) = %s",  tag, safe_pca_top5(r0))
    log_msg("[SVA sanity%s] residual var%% top5 (withSV)= %s", tag, safe_pca_top5(r1))
  }
  invisible(NULL)
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

run_fgsea_save <- function(stats, genesets, out_prefix, top_plot_n = 0, plot_title = "") {
  set.seed(1234)  # <— 確保 fgseaMultilevel/fgsea 的穩定性
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  
  # 先強制序列化執行，避免 BiocParallel 在 Windows 的 SOCK cluster 造成 serialize 連線錯誤
  bp <- tryCatch({
    if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::SerialParam()
    else NULL
  }, error = function(e) NULL)
  
  # 先試 fgseaMultilevel（序列化）；若失敗再回退到經典 fgsea()
  res <- tryCatch({
    suppressWarnings(
      fgseaMultilevel(
        pathways = genesets,
        stats    = stats,
        minSize  = minSize,
        maxSize  = maxSize,
        eps      = fgsea_eps,
        BPPARAM  = bp
      )
    )
  }, error = function(e) {
    log_msg("    [fgseaMultilevel] 失敗：%s → 改用 fgsea() 回退", e$message)
    suppressWarnings(
      fgsea(
        pathways = genesets,
        stats    = stats,
        minSize  = minSize,
        maxSize  = maxSize,
        nperm    = 10000,
        BPPARAM  = bp
      )
    )
  })
  
  # 之後流程照舊
  res <- res |> dplyr::arrange(padj, dplyr::desc(NES))
  data.table::fwrite(as.data.frame(res), paste0(out_prefix, ".csv"))
  wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "GSEA")
  openxlsx::writeData(wb, 1, as.data.frame(res))
  openxlsx::saveWorkbook(wb, paste0(out_prefix, ".xlsx"), overwrite = TRUE)
  
  if (isTRUE(MAKE_PLOTS) && top_plot_n > 0 && nrow(res) > 0) {
    topH <- res |> dplyr::filter(grepl("^HALLMARK_", pathway)) |> dplyr::arrange(padj) |> head(top_plot_n)
    for (pw in topH$pathway) {
      fn <- paste0(out_prefix, "_", gsub("[^A-Za-z0-9_]+","_", pw), ".pdf")
      try({
        pdf(fn, width=6, height=4)
        print(plotEnrichment(genesets[[pw]], stats) + ggplot2::ggtitle(glue::glue("{plot_title}\n{pw}")))
        dev.off()
      }, silent=TRUE)
    }
  }
  res
}


read_gsea_table <- function(out_root, subunit, group_name, stat_tag) {
  grp_safe <- safe_fs_name(group_name)
  fp <- file.path(out_root, subunit, grp_safe, paste0(stat_tag, ".csv"))
  if (!file.exists(fp)) { log_msg("    檔案不存在：{fp}"); return(NULL) }
  dt <- tryCatch(data.table::fread(fp, na.strings=c("NA","NaN","")), error=function(e) NULL)
  if (is.null(dt)) return(NULL)
  need_cols <- c("pathway","NES","padj")
  if (!all(need_cols %in% names(dt))) {
    if (!"pathway" %in% names(dt) && "pathway" %in% tolower(names(dt))) names(dt)[match("pathway", tolower(names(dt)))] <- "pathway"
    if (!"NES" %in% names(dt) && "nes" %in% tolower(names(dt))) names(dt)[match("nes", tolower(names(dt)))] <- "NES"
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
      out_dir <- file.path(sum_root, grp_safe, stat_tag)
      write_summary_outputs(wide, out_dir, grp_name, stat_tag)
    }
  }
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
residualize_vector <- function(y, csn_score, batch = NULL, covars = NULL, min_n = 8) {
  # 只留 y、CSN 都是有限值的樣本
  sam_all <- intersect(names(y), names(csn_score))
  yv  <- as.numeric(y[sam_all])
  sc  <- as.numeric(csn_score[sam_all])
  keep <- is.finite(yv) & is.finite(sc)
  if (sum(keep) < min_n) {
    return(setNames(rep(NA_real_, length(y)), names(y)))
  }
  sam <- sam_all[keep]; yv <- yv[keep]; sc <- sc[keep]
  
  # batch：只有 >=2 levels 才建設計；否則設為 NULL（避免 contrasts<- 炸掉）
  Xb <- NULL
  if (!is.null(batch)) {
    b <- droplevels(as.factor(batch[sam]))
    if (nlevels(b) >= 2) {
      Xb <- stats::model.matrix(~ 0 + b)
      colnames(Xb) <- paste0("b_", colnames(Xb))
      rownames(Xb) <- sam
    }
  }
  
  # covars：移除缺值過多/零變異欄
  Xc <- NULL
  if (!is.null(covars)) {
    cv <- as.matrix(covars[sam, , drop = FALSE])
    if (ncol(cv) > 0) {
      cv[!is.finite(cv)] <- NA
      var_ok <- apply(cv, 2, function(z) {
        z <- as.numeric(z); sum(is.finite(z)) >= 3 && stats::var(z[is.finite(z)]) > 0
      })
      if (any(var_ok)) {
        Xc <- cv[, var_ok, drop = FALSE]
        rownames(Xc) <- sam
      }
    }
  }
  
  # 組合設計矩陣；再做一次欄位篩選（防護）
  des <- cbind("(Intercept)" = 1, CSN = sc,
               if (!is.null(Xb)) Xb else NULL,
               if (!is.null(Xc)) Xc else NULL)
  keep_col <- apply(des, 2, function(z) {
    z <- as.numeric(z); sum(is.finite(z)) >= 3 && stats::var(z[is.finite(z)]) > 0
  })
  des <- as.matrix(des[, keep_col, drop = FALSE])
  
  ok <- is.finite(yv) & rowSums(is.finite(des)) == ncol(des)
  res <- setNames(rep(NA_real_, length(y)), names(y))
  if (sum(ok) < min_n) return(res)
  
  fit <- stats::lm.fit(x = des[ok, , drop = FALSE], y = yv[ok])
  res[sam[ok]] <- fit$residuals
  res
}


audit_covars_coverage <- function(tag, ds_id, stratum, su,
                                  sample_ids,
                                  batch = NULL,
                                  covars = NULL,
                                  sv = NULL,
                                  tech = NULL) {
  dir.create(file.path("run_info","covars_audit"), recursive = TRUE, showWarnings = FALSE)
  
  cov_m <- .ensure_mat_or_null(covars)
  cov_nms <- if (!is.null(cov_m)) colnames(cov_m) else character(0)
  
  get_cov <- function(k) {
    if (!is.null(cov_m) && k %in% colnames(cov_m)) mean(is.finite(cov_m[,k]))*100 else NA_real_
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
  X <- cbind("(Intercept)" = 1, X)
  
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
dataset_dirs <- setNames(file.path(datasets_root, dataset_ids), dataset_ids)

strata <- c("ALL","TP53_mutant","TP53_wild_type")

message("datasets_root = ", datasets_root)
stopifnot(all(dir.exists(dataset_dirs)))



## ==== PATCH: SVA (full-rank + nested mods + robust alignment v3; non-count safe) ====
estimate_svs <- function(M, mod_interest, mod_nuisance = NULL, max_k = 5, label = "") {
  tag <- if (nzchar(label)) paste0(":", label) else ""
  if (!requireNamespace("sva", quietly = TRUE)) {
    log_msg("[SVA%s] package 'sva' not available -> SV=NULL", tag)
    return(NULL)
  }
  
  tryCatch({
    M <- as.matrix(M)
    storage.mode(M) <- "double"
    samp <- colnames(M)  # 以 M 的樣本名為準
    
    ## ---- 更魯棒的設計矩陣對齊 ----
    fix_design <- function(X, nm, role = "") {
      if (is.null(X)) return(NULL)
      X <- as.matrix(X)
      rn <- rownames(X)
      
      if (is.null(rn)) {
        rownames(X) <- nm
        log_msg("[SVA%s] fix_design(%s): no rownames -> assign sample names (safe)", tag, role)
        return(X)
      }
      
      # 先嘗試用名稱對齊並重排
      idx <- match(nm, rn)
      if (!anyNA(idx)) {
        if (!all(idx == seq_along(idx))) {
          log_msg("[SVA%s] fix_design(%s): matched by names and permuted to M", tag, role)
        } else {
          log_msg("[SVA%s] fix_design(%s): names already aligned", tag, role)
        }
        out <- X[idx, , drop = FALSE]
        rownames(out) <- nm
        out[!is.finite(out)] <- 0
        return(out)
      }
      
      # 名稱對不上的 fallback（高風險）：只有在列數吻合時才覆寫
      if (nrow(X) == length(nm)) {
        nbad <- sum(is.na(idx))
        log_msg("[SVA%s] fix_design(%s): %d/%d names not found; ASSUME current row order==M and overwrite rownames (RISKY)",
                tag, role, nbad, length(nm))
        rownames(X) <- nm
        X[!is.finite(X)] <- 0
        return(X)
      }
      
      # 最後防呆：維度不合、直接放棄
      log_msg("[SVA%s] fix_design(%s): DIM MISMATCH (nrow=%d vs samples=%d) -> drop design",
              tag, role, nrow(X), length(nm))
      return(NULL)
    }
    
    mod_interest <- fix_design(mod_interest, samp, "interest")
    mod_nuisance <- fix_design(mod_nuisance, samp, "nuisance")
    
    ## ---- 先建 mod0（攔截 + 干擾），再擴充為滿秩的 mod（保證 mod0 ⊂ mod）----
    mod0 <- cbind("(Intercept)" = 1, mod_nuisance)
    qr0  <- qr(mod0)
    mod0 <- mod0[, qr0$pivot[seq_len(qr0$rank)], drop = FALSE]
    
    add_full_rank <- function(B, Xadd) {
      if (is.null(Xadd)) return(B)
      Xadd <- as.matrix(Xadd)
      if (ncol(Xadd) == 0) return(B)
      keep_nc <- apply(Xadd, 2, function(z) sd(z) > 0)
      if (any(!keep_nc)) Xadd <- Xadd[, keep_nc, drop = FALSE]
      if (ncol(Xadd) == 0) return(B)
      rB <- qr(B)$rank
      for (j in seq_len(ncol(Xadd))) {
        Tst <- cbind(B, Xadd[, j])
        if (qr(Tst)$rank > rB) { B <- Tst; rB <- rB + 1 }
      }
      B
    }
    mod <- add_full_rank(mod0, mod_interest)
    
    if (ncol(M) != nrow(mod)) {
      log_msg("[SVA%s] DIM MISMATCH: ncol(M)=%d vs nrow(mod)=%d -> SV=NULL", tag, ncol(M), nrow(mod))
      return(NULL)
    }
    
    r_mod  <- qr(mod)$rank
    r_mod0 <- qr(mod0)$rank
    res_df <- ncol(M) - r_mod
    log_msg("[SVA%s] dims: genes=%d, samples=%d | rank(mod)=%s, rank(mod0)=%s | residual_df=%s",
            tag, nrow(M), ncol(M), r_mod, r_mod0, res_df)
    
    if (!is.finite(res_df) || res_df < 5) {
      log_msg("[SVA%s] residual df too small -> k=0 (skip SVA)", tag)
      return(NULL)
    }
    
    ## ---- k（只用 be / leek；be 失敗就退 leek）----
    k_be <- tryCatch(sva::num.sv(dat = M, mod = mod, method = "be"),
                     error = function(e){ log_msg("[SVA%s] num.sv(be) ERROR: %s", tag, e$message); NA_integer_ })
    k_le <- tryCatch(sva::num.sv(dat = M, mod = mod, method = "leek"),
                     error = function(e){ log_msg("[SVA%s] num.sv(leek) ERROR: %s", tag, e$message); NA_integer_ })
    log_msg("[SVA%s] num.sv -> be=%s, leek=%s", tag, k_be, k_le)
    
    k <- if (is.finite(k_be) && k_be > 0) k_be else if (is.finite(k_le) && k_le > 0) k_le else 0L
    k <- as.integer(max(0, min(max_k, k, res_df - 1L)))
    if (k == 0) {
      R <- tryCatch({ fit <- stats::lm.fit(x = mod, y = t(M)); fit$residuals }, error = function(e) NULL)
      if (!is.null(R)) {
        pc <- tryCatch(stats::prcomp(R, center = TRUE, scale. = FALSE), error = function(e) NULL)
        if (!is.null(pc)) {
          ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
          log_msg("[SVA%s] residual PCA var%% (top5) = %s", tag, paste(ve[seq_len(min(5, length(ve)))], collapse = ", "))
        }
      }
      return(NULL)
    }
    
    ## ---- 選擇引擎：非 counts → sva；count-like → 優先 svaseq，失敗回退 sva ----
    is_counts <- all(is.finite(M)) && min(M) >= 0 && mean(abs(M - round(M))) < 0.05
    engine <- if (is_counts) "svaseq" else "sva"
    
    if (engine == "svaseq") {
      sv <- tryCatch(
        sva::svaseq(dat = as.matrix(M), mod = as.matrix(mod), mod0 = as.matrix(mod0), n.sv = k)$sv,
        error = function(e){
          log_msg("[SVA%s] svaseq error: %s -> fallback to sva", tag, e$message)
          NULL
        }
      )
      if (is.null(sv)) {
        sv <- tryCatch(sva::sva(dat = as.matrix(M), mod = as.matrix(mod), mod0 = as.matrix(mod0), n.sv = k)$sv,
                       error = function(e){ log_msg("[SVA%s] sva fallback error: %s", tag, e$message); NULL })
        engine <- "sva(fallback)"
      }
    } else {
      sv <- tryCatch(sva::sva(dat = as.matrix(M), mod = as.matrix(mod), mod0 = as.matrix(mod0), n.sv = k)$sv,
                     error = function(e){ log_msg("[SVA%s] sva error: %s", tag, e$message); NULL })
    }
    
    if (is.null(sv)) { log_msg("[SVA%s] returned NULL with k=%d", tag, k); return(NULL) }
    sv <- as.matrix(sv); colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
    log_msg("[SVA%s] %s produced %d SV(s)", tag, engine, ncol(sv))
    sv
  }, error = function(e) {
    log_msg("[SVA%s] ERROR: %s -> SV=NULL", tag, e$message)
    NULL
  })
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

## === NEW: generalized runner for any predictor vector 
run_predictor_analyses <- function(predictor_name,
                                   predictor_vec,          # named numeric by sample
                                   exclude_genes,          # character genes to drop from ranking target
                                   ds_id, ds_dir,
                                   mat0,                   # raw (log2 可能已處理) matrix used by Spearman
                                   mat,                    # imputed+filtered matrix used by limma
                                   out_root,
                                   genesets_by_group,
                                   batch_all, purity_all, sa_all,
                                   tp53_num_all = NULL,
                                   is_ALL = FALSE) {
  
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("      >> predictor：%s", predictor_name)
  
  # ---------- 初始化 ----------
  sub_raw <- predictor_vec
  stats_t_RAW <- NULL
  stats_t_ADJ <- NULL
  stats_rho_ADJ <- NULL
  nH <- NA_integer_
  nL <- NA_integer_
  limma_design_cols <- NA_character_
  spearman_covar_cols <- NA_character_
  
  # ---------- LIMMA 子流程（早退風格；無 else） ----------
  do_limma <- function() {
    out <- list(stats_raw=NULL, stats_adj=NULL, nH=NA_integer_, nL=NA_integer_, design_cols=NA_character_)
    
    v <- sub_raw[is.finite(sub_raw)]
    if (length(v) < (2 * min_per_group)) {
      log_msg("        可用樣本過少（非NA=%d），略過 limma", length(v))
      return(out)
    }
    
    qs <- stats::quantile(v, probs = c(hi_lo_quantile, 1 - hi_lo_quantile), names = FALSE)
    grp <- setNames(rep("Mid", length(sub_raw)), names(sub_raw))
    grp[names(v)[v <= qs[1]]] <- "Low"
    grp[names(v)[v >= qs[2]]] <- "High"
    
    keep_s <- names(sub_raw)[grp != "Mid" & is.finite(sub_raw)]
    nH_loc <- sum(grp[keep_s] == "High")
    nL_loc <- sum(grp[keep_s] == "Low")
    
    if (nH_loc < min_per_group || nL_loc < min_per_group) {
      k <- max(min_per_group, ceiling(hi_lo_quantile * length(v)))
      ord_asc  <- names(sort(v, decreasing = FALSE))
      ord_desc <- names(sort(v, decreasing = TRUE))
      idx_low  <- head(ord_asc,  k)
      idx_high <- head(ord_desc, k)
      grp <- setNames(rep("Mid", length(sub_raw)), names(sub_raw))
      grp[idx_low]  <- "Low"
      grp[idx_high] <- "High"
      keep_s <- union(idx_low, idx_high)
      nH_loc <- sum(grp[keep_s] == "High")
      nL_loc <- sum(grp[keep_s] == "Low")
      log_msg("        ties → 改用排名：High=%d, Low=%d", nH_loc, nL_loc)
    }
    
    if (!(nH_loc >= min_per_group && nL_loc >= min_per_group)) {
      log_msg("        limma 條件不足，略過（High=%d, Low=%d）", nH_loc, nL_loc)
      return(out)
    }
    
    grp2 <- factor(ifelse(grp[keep_s] == "High","High","Low"), levels = c("Low","High"))
    mat_limma <- mat[setdiff(rownames(mat), exclude_genes), keep_s, drop = FALSE]
    
    # RAW
    stats_raw <- tryCatch(
      limma_t_with_covars(mat_limma, grp2),
      error = function(e){ log_msg("        limma RAW 失敗：%s", e$message); NULL }
    )
    
    # 共變項（base）
    batch_subset <- if (!is.null(batch_all)) droplevels(batch_all[keep_s]) else NULL
    base_covars_raw <- data.frame(
      purity = as.numeric(purity_all[keep_s]),
      sex    = as.numeric(sa_all[keep_s, "sex"]),
      age    = as.numeric(sa_all[keep_s, "age"]),
      row.names = keep_s, check.names = FALSE
    )
    if (is_ALL && !is.null(tp53_num_all)) base_covars_raw$TP53_mutant <- as.numeric(tp53_num_all[keep_s])
    
    audit_covars_coverage("limma-base:before-select", ds_id, basename(out_root), predictor_name,
                          keep_s, batch_subset, base_covars_raw, sv = NULL, tech = NULL)
    
    covars <- select_covars_safely(
      base_covars_raw, keep_s, label = "limma-base",
      y = as.numeric(grp2 == "High"),
      min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
      max_abs_cor = 0.30
    )
    
    # 技術/批次（SVA 或 PC1-2）
    svs <- NULL
    tech <- NULL
    include_tech <- is.null(batch_subset) || nlevels(batch_subset) < 2 || length(keep_s) >= 30
    if (include_tech) {
      nuis_mm <- cbind(
        if (!is.null(batch_subset)) model.matrix(~ 0 + batch_subset) else NULL,
        covars
      )
      
      mod_interest_l <- as.matrix(model.matrix(~ 0 + grp2))
      if (nrow(mod_interest_l) == length(keep_s)) rownames(mod_interest_l) <- keep_s
      if (is_ALL && !is.null(tp53_num_all)) {
        tp53_l <- matrix(as.numeric(tp53_num_all[keep_s]), ncol = 1,
                         dimnames = list(keep_s, "TP53_mutant"))
        mod_interest_l <- cbind(mod_interest_l, tp53_l)
      }
      
      if (!is.null(nuis_mm)) {
        nuis_mm <- as.matrix(nuis_mm)
        if (nrow(nuis_mm) == length(keep_s)) rownames(nuis_mm) <- keep_s
      }
      
      audit_design_alignment(
        tag      = sprintf("%s|%s|%s|limma:preSVA", ds_id, basename(out_root), predictor_name),
        samples  = colnames(mat_limma),
        mod_interest = mod_interest_l,
        mod_nuisance = nuis_mm,
        out_dir  = out_root
      )
      
      svs <- estimate_svs(
        M            = mat_limma,
        mod_interest = mod_interest_l,
        mod_nuisance = nuis_mm,
        label        = sprintf("%s|%s|%s|limma", ds_id, basename(out_root), predictor_name)
      )
      
      audit_sva_sanity(
        M = mat_limma, sv = svs,
        batch_factor = batch_subset,
        purity = purity_all, sex = sa_all$sex, age = sa_all$age,
        y_bio = as.numeric(grp2 == "High"),
        label = sprintf("%s|%s|%s|limma", ds_id, basename(out_root), predictor_name)
      )
      
      if (is.null(svs)) {
        tmp <- build_tech_covars(mat_limma)[, c("PC1","PC2"), drop = FALSE]
        tech <- gate_tech_pcs(tmp, as.numeric(grp2 == "High"))
      }
      if (!is.null(svs)) tech <- svs
      
      tech <- orthogonalize_to(tech, nuis_mm)
      tech <- .ensure_mat_or_null(tech)
      
      if (!is.null(tech)) {
        if (is.null(colnames(tech))) colnames(tech) <- paste0("tech", seq_len(ncol(tech)))
        covars <- if (is.null(covars)) tech else cbind(covars, tech)
      }
      
      audit_covars_coverage("limma-tech", ds_id, basename(out_root), predictor_name,
                            keep_s, batch_subset, covars, sv = svs, tech = tech)
    }
    
    # 最後保守篩選 + TP53（ALL）
    covars <- select_covars_safely(
      covars, keep_s, label = "limma-final",
      y = as.numeric(grp2 == "High"),
      min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
      max_abs_cor = 0.30
    )
    if (is_ALL && !is.null(tp53_num_all)) {
      tp53_c <- matrix(as.numeric(tp53_num_all[keep_s]), ncol = 1,
                       dimnames = list(keep_s, "TP53_mutant"))
      if (is.null(covars)) covars <- tp53_c
      if (!is.null(covars) && !("TP53_mutant" %in% colnames(covars))) covars <- cbind(covars, tp53_c)
    }
    
    audit_covars_coverage("limma-final:design-ready", ds_id, basename(out_root), predictor_name,
                          keep_s, batch_subset, covars, sv = svs, tech = tech)
    
    # ADJ
    stats_adj <- tryCatch(
      limma_t_with_covars(mat_limma, grp2, batch = batch_subset, covars = covars),
      error = function(e){ log_msg("        limma BatchAdj 失敗：%s", e$message); NULL }
    )
    
    design_cols <- attr(stats_adj, "design_cols", exact = TRUE)
    if (is.null(design_cols)) design_cols <- NA_character_
    
    out$stats_raw   <- stats_raw
    out$stats_adj   <- stats_adj
    out$nH          <- nH_loc
    out$nL          <- nL_loc
    out$design_cols <- design_cols
    out
  }
  
  lim <- do_limma()
  stats_t_RAW <- lim$stats_raw
  stats_t_ADJ <- lim$stats_adj
  nH <- lim$nH
  nL <- lim$nL
  limma_design_cols <- lim$design_cols
  
  # ---------- SPEARMAN 子流程（早退風格；無 else） ----------
  do_spearman <- function() {
    out <- list(stats=NULL, covar_cols=NA_character_, svs=NULL, tech=NULL)
    
    covars_all_base <- data.frame(
      purity = as.numeric(purity_all[colnames(mat0)]),
      sex    = as.numeric(sa_all[colnames(mat0), "sex"]),
      age    = as.numeric(sa_all[colnames(mat0), "age"]),
      row.names = colnames(mat0), check.names = FALSE
    )
    if (is_ALL && !is.null(tp53_num_all)) covars_all_base$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
    
    audit_covars_coverage("spearman-base:before-select", ds_id, basename(out_root), predictor_name,
                          colnames(mat0), batch_all, covars_all_base, sv = NULL, tech = NULL)
    
    covars_all <- select_covars_safely(
      covars_all_base, colnames(mat0), label = "spearman-base",
      y = as.numeric(sub_raw),
      min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
      max_abs_cor = 0.30
    )
    
    svs_all  <- NULL
    tech_all <- NULL
    include_tech_all <- is.null(batch_all) || nlevels(batch_all) < 2 || ncol(mat) >= 30
    if (include_tech_all) {
      nuis_mm_all <- cbind(
        if (!is.null(batch_all)) model.matrix(~ 0 + batch_all) else NULL,
        covars_all
      )
      
      x_bio <- as.numeric(sub_raw[colnames(mat)])
      mod_interest_s <- matrix(x_bio, ncol = 1, dimnames = list(colnames(mat), "sub_raw"))
      if (is_ALL && !is.null(tp53_num_all)) {
        tp53_s <- matrix(as.numeric(tp53_num_all[colnames(mat)]), ncol = 1,
                         dimnames = list(colnames(mat), "TP53_mutant"))
        mod_interest_s <- cbind(mod_interest_s, tp53_s)
      }
      
      if (!is.null(nuis_mm_all)) {
        nuis_mm_all <- as.matrix(nuis_mm_all)
        if (nrow(nuis_mm_all) == ncol(mat)) rownames(nuis_mm_all) <- colnames(mat)
      }
      
      audit_design_alignment(
        tag      = sprintf("%s|%s|%s|spearman:preSVA", ds_id, basename(out_root), predictor_name),
        samples  = colnames(mat),
        mod_interest = mod_interest_s,
        mod_nuisance = nuis_mm_all,
        out_dir  = out_root
      )
      
      svs_all <- estimate_svs(
        M            = mat,
        mod_interest = mod_interest_s,
        mod_nuisance = nuis_mm_all,
        label        = sprintf("%s|%s|%s|spearman", ds_id, basename(out_root), predictor_name)
      )
      
      audit_sva_sanity(
        M = mat, sv = svs_all,
        batch_factor = batch_all,
        purity = purity_all, sex = sa_all$sex, age = sa_all$age,
        y_bio = as.numeric(sub_raw),
        label = sprintf("%s|%s|%s|spearman", ds_id, basename(out_root), predictor_name)
      )
      
      tech_all <- if (is.null(svs_all)) {
        tmp <- build_tech_covars(mat)[, c("PC1","PC2"), drop = FALSE]
        gate_tech_pcs(tmp, sub_raw[colnames(mat)])
      } else {
        svs_all
      }
      
      tech_all <- orthogonalize_to(tech_all, nuis_mm_all)
      tech_all <- .ensure_mat_or_null(tech_all)
      
      if (!is.null(tech_all)) {
        if (is.null(colnames(tech_all))) colnames(tech_all) <- paste0("tech", seq_len(ncol(tech_all)))
        covars_all <- if (is.null(covars_all)) tech_all else cbind(covars_all, tech_all)
        
        audit_covars_coverage("spearman-tech", ds_id, basename(out_root), predictor_name,
                              colnames(mat0), batch_all, covars_all,
                              sv = svs_all, tech = tech_all)
      }
      
      if (is.null(tech_all)) {
        audit_covars_coverage("spearman-tech", ds_id, basename(out_root), predictor_name,
                              colnames(mat0), batch_all, covars_all,
                              sv = svs_all, tech = NULL)
      }
    }
    
    # 最後保守篩選 + TP53（ALL）
    covars_all <- select_covars_safely(
      covars_all, colnames(mat0), label = "spearman-final",
      y = as.numeric(sub_raw),
      min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
      max_abs_cor = 0.30
    )
    if (is_ALL && !is.null(tp53_num_all)) {
      tp53_all_c <- matrix(as.numeric(tp53_num_all[colnames(mat0)]), ncol = 1,
                           dimnames = list(colnames(mat0), "TP53_mutant"))
      if (is.null(covars_all)) covars_all <- tp53_all_c
      if (!is.null(covars_all) && !("TP53_mutant" %in% colnames(covars_all))) covars_all <- cbind(covars_all, tp53_all_c)
    }
    
    audit_covars_coverage("spearman-final:design-ready", ds_id, basename(out_root), predictor_name,
                          colnames(mat0), batch_all, covars_all,
                          sv = svs_all, tech = tech_all)
    
    mat0_no_ex <- mat0[setdiff(rownames(mat0), exclude_genes), , drop = FALSE]
    stats_rho <- tryCatch(
      spearman_rank_with_covars(
        mat_raw = mat0_no_ex,
        sub_raw = sub_raw,
        batch = batch_all, covars = covars_all,
        min_pairs = min_pairs_spearman, impute_for_rbe = TRUE
      ),
      error = function(e){ log_msg("        Spearman BatchAdj 無法進行：%s", e$message); NULL }
    )
    
    covar_cols <- if (!is.null(stats_rho)) attr(stats_rho, "covar_cols", exact = TRUE) else NULL
    if (is.null(covar_cols)) covar_cols <- NA_character_
    
    out$stats <- stats_rho
    out$covar_cols <- covar_cols
    out$svs <- svs_all
    out$tech <- tech_all
    out
  }
  
  sp <- do_spearman()
  stats_rho_ADJ <- sp$stats
  spearman_covar_cols <- sp$covar_cols
  svs_all <- sp$svs
  tech_all <- sp$tech
  
  # ---------- 輸出 ----------
  out_roots <- list(RAW = file.path(out_root, "RAW"),
                    BatchAdj = file.path(out_root, "BatchAdj"))
  for (ver in names(out_roots)) {
    dir.create(file.path(out_roots[[ver]], predictor_name), recursive = TRUE, showWarnings = FALSE)
  }
  
  write_one <- function(stats_list, ver_tag){
    if (is.null(stats_list$limma) && is.null(stats_list$spearman)) return(invisible(NULL))
    out_dir <- file.path(out_roots[[ver_tag]], predictor_name)
    for (grp_name in names(genesets_by_group)) {
      subdir <- file.path(out_dir, safe_fs_name(grp_name))
      dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      if (!is.null(stats_list$limma)) {
        run_fgsea_save(
          stats = stats_list$limma,
          genesets = genesets_by_group[[grp_name]],
          out_prefix = file.path(subdir, "GSEA_limma_t"),
          top_plot_n = 0,
          plot_title = glue::glue("{ds_id} | {predictor_name} | High vs Low | {grp_name} | {ver_tag}")
        )
      }
      if (!is.null(stats_list$spearman)) {
        run_fgsea_save(
          stats = stats_list$spearman,
          genesets = genesets_by_group[[grp_name]],
          out_prefix = file.path(subdir, "GSEA_spearman"),
          top_plot_n = 0,
          plot_title = glue::glue("{ds_id} | {predictor_name} | Spearman | {grp_name} | {ver_tag}")
        )
      }
    }
  }
  write_one(list(limma = stats_t_RAW, spearman = NULL), "RAW")
  write_one(list(limma = stats_t_ADJ,  spearman = stats_rho_ADJ), "BatchAdj")
  
  # ---------- 記錄 ----------
  counts_row <- data.frame(
    dataset = ds_id,
    stratum = basename(out_root),
    subunit = predictor_name,
    limma_High = nH, limma_Low = nL,
    spearman_nonNA = sum(is.finite(sub_raw)),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(counts_row, file.path(out_root, "sample_counts.csv"),
                     append = file.exists(file.path(out_root, "sample_counts.csv")))
  
  note_fp <- file.path(out_roots$BatchAdj, predictor_name, "batch_correction_notes.txt")
  info <- c(
    sprintf("dataset=%s; predictor=%s", ds_id, predictor_name),
    sprintf("exclude_genes_from_ranking=%s", paste(exclude_genes, collapse=",")),
    sprintf("batch_levels=%s", if (!is.null(batch_all)) nlevels(batch_all) else NA),
    "limma: design = ~ 0 + group (+ batch/SV/PC1-2 + optional purity/sex/age + TP53 in ALL)",
    sprintf("limma_design_cols=%s", if (length(limma_design_cols) && all(!is.na(limma_design_cols))) paste(limma_design_cols, collapse=",") else "NA"),
    "spearman: removeBatchEffect(batch, covariates=SV or PC1-2 + optional purity/sex/age + TP53 in ALL; design=~predictor_raw (+ TP53 protected in SVA))",
    sprintf("spearman_covar_cols=%s", if (length(spearman_covar_cols) && all(!is.na(spearman_covar_cols))) paste(spearman_covar_cols, collapse=",") else "NA")
  )
  cat(paste(info, collapse = "\n"), file = note_fp)
  
  invisible(NULL)
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
run_one_stratum <- function(ds_id, ds_dir, mat0_full, sample_keep, out_root, genesets_by_group){
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  log_msg("  -- stratum：%s | N(sample_keep)=%d", basename(out_root), length(sample_keep))
  
  ## 1) 子集樣本 + 規模檢查 + log 標度
  keep <- intersect(colnames(mat0_full), sample_keep)
  if (length(keep) < 4) { log_msg("  [跳過] 樣本過少：%d", length(keep)); return(invisible(NULL)) }
  mat0 <- mat0_full[, keep, drop = FALSE]
  mx <- suppressWarnings(max(mat0, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) mat0 <- log2(mat0 + 1)
  
  ## 2) limma 用的插補+過濾矩陣（Spearman 仍用 raw mat0）
  mat  <- impute_and_filter(mat0, min_frac = min_frac_complete)
  present_sub <- intersect(csn_subunits, rownames(mat0))
  if (!length(present_sub)) { log_msg("  [跳過] 此分層沒有任何 CSN subunit"); return(invisible(NULL)) }
  ## —— 加入 CSN subunits 覆蓋率 & CSN_SCORE 可行性審核 —— ##
  audit_csn_score_feasibility(
    ds_id    = ds_id,
    stratum  = basename(out_root),
    mat0     = mat0,
    present_sub = present_sub,
    min_members     = 5L,              # 與 build_csn_score 同步
    pca_min_samples = 10L,             # 與 build_csn_score 同步
    min_per_group   = min_per_group,   # 你的 limma 門檻（高/低各至少 8）
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
  
  ## 3) 共變項 / 批次（對齊樣本）
  bi_all     <- get_batch_factor(ds_dir, colnames(mat0))
  batch_all  <- if (!is.null(bi_all)) droplevels(bi_all$fac[colnames(mat0)]) else NULL
  purity_all <- get_purity_covariate(ds_id, ds_dir, colnames(mat0))
  sa_all     <- get_sex_age_covariates(ds_dir, colnames(mat0))  # data.frame(sex, age)
  
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
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    }
    
    ## 4b) CSN complex score（PC1；方向校正）；排名時排除所有 CSN 成員
    csn_score <- build_csn_score(mat0, subunits = present_sub, combine_7AB = TRUE, min_members = 5L)
    if (sum(is.finite(csn_score)) >= (2 * min_per_group)) {
      run_predictor_analyses(
        predictor_name = "CSN_SCORE",
        predictor_vec  = csn_score,
        exclude_genes  = present_sub,
        ds_id = ds_id, ds_dir = ds_dir,
        mat0 = mat0, mat = impute_and_filter(mat0, min_frac = min_frac_complete),
        out_root = out_root,
        genesets_by_group = genesets_by_group,
        batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
        tp53_num_all = tp53_num_all, is_ALL = is_ALL
      )
    } else {
      log_msg("  [CSN_SCORE] 非 NA 樣本不足，略過")
    }
    
    ## 4c) RESIDUAL_<SU>（加護欄：CSN score 非 NA 樣本數不足 → 整批略過）
    min_n_resid <- min_per_group  # 與你的 limma 門檻對齊（目前為 8）
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
      if (is_ALL && !is.null(tp53_num_all)) {
        base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat0)])
      }
      
      for (su in present_sub) {
        res_sub <- residualize_vector(
          y = mat0[su, ],
          csn_score = csn_score,
          batch = batch_all,
          covars = base_covars_all,
          min_n = min_n_resid      # 明確傳入門檻
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
          batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
          tp53_num_all = tp53_num_all, is_ALL = is_ALL
        )
      }
    }
  }
  
  ## 5) 各版本各自彙整
  present_sub <- intersect(csn_subunits, rownames(mat0))
  sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))
  
  # RAW：只有 limma_t
  ver_root <- file.path(out_root, "RAW")
  summarize_all_groups(
    out_root = ver_root,
    csn_subunits = sum_units,
    genesets_by_group = genesets_by_group,
    stat_tags = c("GSEA_limma_t")
  )
  
  # BatchAdj：limma_t + spearman
  ver_root <- file.path(out_root, "BatchAdj")
  summarize_all_groups(
    out_root = ver_root,
    csn_subunits = sum_units,
    genesets_by_group = genesets_by_group,
    stat_tags = c("GSEA_limma_t", "GSEA_spearman")
  )
  
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
  .z  <- function(v){ v <- suppressWarnings(as.numeric(v))
  m <- mean(v[is.finite(v)], na.rm=TRUE)
  s <- stats::sd(v[is.finite(v)], na.rm=TRUE); if (!is.finite(s) || s==0) s <- 1
  v[!is.finite(v)] <- m; (v - m)/s
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
  
  # --- SEX（Male=1, Female=0；其餘 NA）---
  sex_col_candidates <- c("SEX","GENDER")
  sex_col <- intersect(sex_col_candidates, names(pat))[1]
  sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  if (!is.na(sex_col)) {
    raw <- toupper(as.character(pat[[sex_col]]))
    val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
    names(val) <- as.character(pat[[pid_pat]])
    sex[] <- unname(val[ map_pt[sample_ids] ])
  }
  
  # --- AGE（優先 AGE，其次常見替代欄位）---
  age_col_candidates <- c("AGE","AGE_AT_DIAGNOSIS","AGE_AT_INDEX","AGE_YEARS")
  age_col <- intersect(age_col_candidates, names(pat))[1]
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[ map_pt[sample_ids] ])
    age[] <- .z(age_raw)
  }
  
  # 記錄 coverage（全體樣本）
  cov_sex <- mean(is.finite(sex)) * 100
  cov_age <- mean(is.finite(age)) * 100
  log_msg("    covariates（patient→sample 對齊） coverage：sex %.1f%%, age %.1f%%", cov_sex, cov_age)
  
  out <- data.frame(sex = as.numeric(sex), age = as.numeric(age),
                    row.names = sample_ids, check.names = FALSE)
  return(out)
}

select_covars_safely <- function(df, sample_order, label = "",
                                 y = NULL,
                                 min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
                                 max_abs_cor = 0.30,
                                 min_pairs = 20) {
  if (is.null(df) || ncol(as.data.frame(df)) == 0) return(NULL)
  
  # 轉成 data.frame 並保留欄名
  mat <- as.data.frame(df, check.names = FALSE)
  
  # 對齊樣本（用 match，能抓到不存在的樣本 → NA）
  if (is.null(rownames(mat))) rownames(mat) <- sample_order
  idx <- match(sample_order, rownames(mat))
  if (anyNA(idx)) {
    log_msg("  [covars-%s] 警告：%d/%d 個樣本在共變數表中找不到 → 以 NA 補",
            label, sum(is.na(idx)), length(idx))
  }
  mat <- mat[idx, , drop = FALSE]
  rownames(mat) <- sample_order
  
  # 每欄轉為數值；把 NaN 改成 NA（as.numeric 可能產生 NaN）
  for (cn in colnames(mat)) {
    v <- mat[[cn]]
    if (is.matrix(v) && ncol(v) == 1) v <- v[, 1]
    if (!is.numeric(v)) v <- suppressWarnings(as.numeric(v))
    v[is.nan(v)] <- NA_real_
    mat[[cn]] <- v
  }
  
  keep_cols <- list()
  for (cn in colnames(mat)) {
    x <- mat[[cn]]
    
    # 覆蓋率：用 is.finite，±Inf 也當缺值
    cover <- sum(is.finite(x)) / length(x)
    thr   <- if (cn %in% names(min_cov_named)) min_cov_named[[cn]] else min(min_cov_named)
    if (!is.finite(cover) || cover < thr) {
      log_msg("  [covars-%s] drop %s (coverage=%.0f%% < %.0f%%)", label, cn, 100*cover, 100*thr)
      next
    }
    
    fin <- is.finite(x)
    if (sum(fin) < 3 || sd(x[fin]) == 0) {
      log_msg("  [covars-%s] drop %s (near-constant)", label, cn)
      next
    }
    
    if (!is.null(y)) {
      ok <- fin & is.finite(y)
      if (sum(ok) >= min_pairs) {
        rshow <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
        if (is.finite(rshow) && abs(rshow) >= max_abs_cor) {
          log_msg("  [covars-%s] drop %s (|rho|=%.3f ≥ %.2f vs biology)", label, cn, abs(rshow), max_abs_cor)
          next
        }
      }
    }
    
    keep_cols[[cn]] <- x
  }
  
  if (!length(keep_cols)) return(NULL)
  out <- do.call(cbind, keep_cols)
  out <- as.matrix(out)
  rownames(out) <- sample_order
  colnames(out) <- names(keep_cols)
  out
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
  
  # 印一行摘要（加上 PURITY）
  line <- sprintf(
    "[Audit] %-24s | SEX:%-7s nonNA=%5.1f%% (M=%d,F=%d) | AGE:%-7s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | PURITY:%-35s nonNA=%5.1f%% (min/med/max=%s/%s/%s) | BATCH:%-18s levels=%-2s nonNA=%5.1f%% [%s]",
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
    batch_col, as.integer(batch_levels), batch_nonNA %||% NaN, batch_sizes %||% "NA"
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



## ===== DO_AUDIT 開關包住 QC/探索（4）=====
DO_AUDIT <- FALSE  # <— 預設關閉；需要時改 TRUE

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
  geneset_groups_selected = GENESET_GROUPS_TO_RUN
), file.path("run_info","run_manifest.yml"))


## 想從哪一個開始重跑（可改）
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix  <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' 不在 dataset_ids 裡", start_from))
dataset_dirs_run <- dataset_dirs[ ord[ix:length(ord)] ]

## ===== 主要流程：每個 data set × strata（ALL / TP53_mutant / TP53_wild_type）=====
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  if (!dir.exists(ds_dir)) { log_msg("略過：找不到資料夾 {ds_dir}"); next }
  if (!file.exists(file.path(ds_dir, "data_protein_quantification.txt"))) {
    log_msg("略過：{ds} 缺少 data_protein_quantification.txt"); next
  }
  log_msg("== 開始資料集：{ds}（含 TP53 分層）==")
  
  ## 先讀入完整矩陣
  mat0_full <- load_matrix_from_dataset_dir(ds_dir)
  
  ## 取得 TP53 狀態 → 三個 strata 的樣本清單
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))
  samples_ALL <- colnames(mat0_full)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT  <- names(tp53_status)[tp53_status == "TP53_wild_type"]
  
  ## 每個分層的輸出根目錄
  base_tp53_root <- file.path(ds_dir, "csn_gsea_results_TP53")
  
  ## 依序跑三個 strata（RAW & BatchAdj 都會在 run_one_stratum 裡寫好）
  run_one_stratum(ds_id = ds, ds_dir = ds_dir,
                  mat0_full = mat0_full,
                  sample_keep = samples_ALL,
                  out_root = file.path(base_tp53_root, "ALL"),
                  genesets_by_group = genesets_by_group)
  
  run_one_stratum(ds_id = ds, ds_dir = ds_dir,
                  mat0_full = mat0_full,
                  sample_keep = samples_MUT,
                  out_root = file.path(base_tp53_root, "TP53_mutant"),
                  genesets_by_group = genesets_by_group)
  
  run_one_stratum(ds_id = ds, ds_dir = ds_dir,
                  mat0_full = mat0_full,
                  sample_keep = samples_WT,
                  out_root = file.path(base_tp53_root, "TP53_wild_type"),
                  genesets_by_group = genesets_by_group)
  
  log_msg("== 完成資料集：{ds}（TP53 分層輸出 → {base_tp53_root}）==")
}


## ===== PAN：跨癌別 × 跨 CSN 次單元（TP53 分層版；RAW/BatchAdj 各自做）=====
versions <- c("RAW","BatchAdj")
groups    <- names(genesets_by_group)
stat_tags <- c("GSEA_limma_t","GSEA_spearman")
cand <- list.dirs(datasets_root, full.names=FALSE, recursive=FALSE)
dataset_ids_have <- intersect(cand, dataset_ids)

read_one_result_tp53 <- function(ds_id, stratum, version, subunit, grp_name, stat_tag){
  ds_dir <- file.path(datasets_root, ds_id)
  fp <- file.path(ds_dir, "csn_gsea_results_TP53", stratum, version,
                  subunit, safe_fs_name(grp_name), paste0(stat_tag, ".csv"))
  if (!file.exists(fp)) return(NULL)
  dt <- tryCatch(data.table::fread(fp, na.strings=c("NA","NaN","")), error=function(e) NULL)
  if (is.null(dt)) return(NULL)
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

collect_all_long_tp53 <- function(dataset_ids, strata, versions, groups, stat_tags, subunits){
  res <- list(); k <- 1L
  for (ds in dataset_ids) for (stt in strata) for (ver in versions) for (g in groups) for (st in stat_tags) for (su in subunits) {
    tmp <- read_one_result_tp53(ds, stt, ver, su, g, st); if (!is.null(tmp)) { res[[k]] <- tmp; k <- k + 1L }
  }
  if (!length(res)) return(NULL)
  dplyr::bind_rows(res)
}

for (ver in versions) {
  log_msg("PAN(TP53,%s) 聚合開始", ver)
  long_df <- collect_all_long_tp53(dataset_ids_have, strata, versions = ver, groups, stat_tags, csn_subunits)
  if (is.null(long_df) || !nrow(long_df)) { log_msg("  無可用結果（%s）", ver); next }
  
  long_df <- long_df |>
    dplyr::mutate(
      sig_005 = !is.na(padj) & padj < 0.05,
      sig_025 = !is.na(padj) & padj < 0.25,
      pos     = NES > 0,
      neg     = NES < 0
    )
  
  pan_out_root <- file.path(getwd(), "csn_gsea_pan_summary_TP53", ver)
  dir.create(pan_out_root, recursive = TRUE, showWarnings = FALSE)
  
  # 每癌別 × 路徑 摘要
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
  data.table::fwrite(per_ds_term, file.path(pan_out_root, "per_dataset_summary.csv"))
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
    dplyr::arrange(
      dplyr::desc(datasets_with_sig_0_05),
      dplyr::desc(total_sig_0_05),
      pathway
    )
  
  data.table::fwrite(pan_term, file.path(pan_out_root, "pan_summary.csv"))
  openxlsx::addWorksheet(wb_top, "pan_summary")
  openxlsx::writeData(wb_top, "pan_summary", pan_term)
  openxlsx::saveWorkbook(wb_top, file.path(pan_out_root, "pan_summary.xlsx"), overwrite = TRUE)
  
  ## 每個 stratum × gene-set group × stat 的輸出（沿用你原本的表1/表2）
  for (stt in strata) {
    for (g in groups) {
      for (st in stat_tags) {
        sub_out <- file.path(pan_out_root, stt, safe_fs_name(g), st)
        dir.create(sub_out, recursive = TRUE, showWarnings = FALSE)
        
        long_sub   <- long_df    |> dplyr::filter(stratum == stt, group == g, stat == st)
        if (nrow(long_sub) == 0) next
        
        per_ds_sub <- per_ds_term |> dplyr::filter(stratum == stt, group == g, stat == st)
        pan_sub    <- pan_term    |> dplyr::filter(stratum == stt, group == g, stat == st) |>
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
        
        ## 表1/表2
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
        for (su in csn_subunits) {
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
      }
    }
  }
  # 迴圈最末尾（該 ver 的輸出都完成後）
  if (exists("pan_out_root")) {
    message(sprintf("[PAN-TP53] 整合完成（%s）→ %s", ver, pan_out_root))
  }
}
