#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(openxlsx)
})

min_frac_complete   <- 0.75     # 每基因至少 75% 非 NA
min_per_group      <- 8         # High/Low 各至少多少樣本（不足就跳過 limma）
csn_subunits <- c("GPS1","COPS2","COPS3","COPS4","COPS5","COPS6","COPS7A","COPS7B","COPS8","COPS9")
RUN_HILO_SENSITIVITY <- FALSE

## =========================
## 共同小工具（單一定義）
## =========================
`%||%` <- function(a,b) if (!is.null(a)) a else b
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+","_", x))
.zscore <- function(v){
  v <- suppressWarnings(as.numeric(v))
  mu <- mean(v[is.finite(v)], na.rm=TRUE)
  sdv <- stats::sd(v[is.finite(v)], na.rm=TRUE); if (!is.finite(sdv) || sdv==0) sdv <- 1
  v[is.nan(v)] <- NA_real_
  v[!is.finite(v)] <- mu
  (v - mu)/sdv
}
.to01 <- function(v){
  v <- suppressWarnings(as.numeric(v))
  if (sum(is.finite(v) & v>1, na.rm=TRUE) > sum(is.finite(v) & v<=1, na.rm=TRUE)) v <- v/100
  pmin(pmax(v,0),1)
}
.median_from_semicolon <- function(x_chr){
  vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
  vv <- vv[is.finite(v)]; if (!length(vv)) return(NA_real_)
  stats::median(vv)
}

## 全域開關：永遠加 tech covariates（即使有 batch）
ALWAYS_ADD_TECH_COVARS <- FALSE

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

impute_and_filter <- function(mat, min_frac=0.75) {
  keep <- rowMeans(!is.na(mat)) >= min_frac
  m <- mat[keep, , drop=FALSE]
  if (any(is.na(m))) { set.seed(1234); m <- imputeLCMD::impute.MinProb(m, q=0.01) }
  m
}

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
  
  # ---------- LIMMA 子流程（主：連續；副：High/Low 敏感度可選） ----------
  do_limma <- function() {
    out <- list(stats_cont_raw=NULL, stats_cont_adj=NULL,
                stats_hilo_raw=NULL, stats_hilo_adj=NULL,
                nH=NA_integer_, nL=NA_integer_, design_cols=NA_character_)
    
    # ---- (A1) 連續 predictor：RAW（只 ~ x）----
    mat_limma_all <- mat[setdiff(rownames(mat), exclude_genes), , drop = FALSE]
    x_all <- sub_raw[colnames(mat_limma_all)]
    if (sum(is.finite(x_all)) >= 4) {
      out$stats_cont_raw <- tryCatch(
        limma_cont_with_covars(mat_limma_all, x_all, batch = NULL, covars = NULL),
        error = function(e){ log_msg("        limma-cont RAW 失敗：%s", e$message); NULL }
      )
    } else {
      log_msg("        連續 predictor 可用樣本過少，跳過 limma-cont RAW")
    }
    
    # ---- (A2) 連續 predictor：BatchAdj（~ x + batch + covars [+ tech/SV]）----
    batch_subset <- batch_all
    base_covars_all <- data.frame(
      purity = as.numeric(purity_all[colnames(mat)]),
      sex    = as.numeric(sa_all[colnames(mat), "sex"]),
      age    = as.numeric(sa_all[colnames(mat), "age"]),
      row.names = colnames(mat), check.names = FALSE
    )
    if (is_ALL && !is.null(tp53_num_all)) base_covars_all$TP53_mutant <- as.numeric(tp53_num_all[colnames(mat)])
    
    audit_covars_coverage("limma-cont:base-before-select", ds_id, basename(out_root), predictor_name,
                          colnames(mat), batch_subset, base_covars_all, sv = NULL, tech = NULL)
    
    covars_all <- select_covars_safely(
      base_covars_all, colnames(mat), label = "limma-cont:base",
      y = as.numeric(sub_raw[colnames(mat)]),
      min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
      max_abs_cor = 0.30
    )
    
    svs <- NULL; tech <- NULL
    include_tech <- is.null(batch_subset) || nlevels(batch_subset) < 2 || ncol(mat) >= 30
    if (include_tech) {
      nuis_mm <- cbind(
        if (!is.null(batch_subset)) model.matrix(~ 0 + batch_subset) else NULL,
        covars_all
      )
      x_bio <- as.numeric(sub_raw[colnames(mat)])
      mod_interest_l <- matrix(x_bio, ncol = 1, dimnames = list(colnames(mat), "x"))
      if (is_ALL && !is.null(tp53_num_all)) {
        tp53_l <- matrix(as.numeric(tp53_num_all[colnames(mat)]), ncol = 1,
                         dimnames = list(colnames(mat), "TP53_mutant"))
        mod_interest_l <- cbind(mod_interest_l, tp53_l)
      }
      
      audit_design_alignment(
        tag      = sprintf("%s|%s|%s|limma-cont:preSVA", ds_id, basename(out_root), predictor_name),
        samples  = colnames(mat),
        mod_interest = mod_interest_l,
        mod_nuisance = nuis_mm,
        out_dir  = out_root
      )
      
      svs <- estimate_svs(
        M            = mat,
        mod_interest = mod_interest_l,
        mod_nuisance = nuis_mm,
        label        = sprintf("%s|%s|%s|limma-cont", ds_id, basename(out_root), predictor_name)
      )
      
      audit_sva_sanity(
        M = mat, sv = svs,
        batch_factor = batch_subset,
        purity = purity_all, sex = sa_all$sex, age = sa_all$age,
        y_bio = as.numeric(sub_raw[colnames(mat)]),
        label = sprintf("%s|%s|%s|limma-cont", ds_id, basename(out_root), predictor_name)
      )
      
      tech <- if (is.null(svs)) {
        tmp <- build_tech_covars(mat)[, c("PC1","PC2"), drop = FALSE]
        gate_tech_pcs(tmp, sub_raw[colnames(mat)])
      } else svs
      
      tech <- orthogonalize_to(tech, nuis_mm)
      tech <- .ensure_mat_or_null(tech)
      if (!is.null(tech)) {
        if (is.null(colnames(tech))) colnames(tech) <- paste0("tech", seq_len(ncol(tech)))
        covars_all <- if (is.null(covars_all)) tech else cbind(covars_all, tech)
        audit_covars_coverage("limma-cont:tech", ds_id, basename(out_root), predictor_name,
                              colnames(mat), batch_subset, covars_all, sv = svs, tech = tech)
      } else {
        audit_covars_coverage("limma-cont:tech", ds_id, basename(out_root), predictor_name,
                              colnames(mat), batch_subset, covars_all, sv = svs, tech = NULL)
      }
    }
    
    covars_all <- select_covars_safely(
      covars_all, colnames(mat), label = "limma-cont:final",
      y = as.numeric(sub_raw[colnames(mat)]),
      min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
      max_abs_cor = 0.30
    )
    if (is_ALL && !is.null(tp53_num_all)) {
      tp53_c <- matrix(as.numeric(tp53_num_all[colnames(mat)]), ncol = 1,
                       dimnames = list(colnames(mat), "TP53_mutant"))
      if (is.null(covars_all)) covars_all <- tp53_c
      if (!is.null(covars_all) && !("TP53_mutant" %in% colnames(covars_all))) covars_all <- cbind(covars_all, tp53_c)
    }
    
    audit_covars_coverage("limma-cont:final-design", ds_id, basename(out_root), predictor_name,
                          colnames(mat), batch_subset, covars_all, sv = svs, tech = tech)
    
    mat_limma_all <- mat[setdiff(rownames(mat), exclude_genes), , drop = FALSE]
    x_all <- sub_raw[colnames(mat_limma_all)]
    out$stats_cont_adj <- tryCatch(
      limma_cont_with_covars(mat_limma_all, x_all, batch = batch_subset, covars = covars_all),
      error = function(e){ log_msg("        limma-cont BatchAdj 失敗：%s", e$message); NULL }
    )
    design_cols <- attr(out$stats_cont_adj, "design_cols", exact = TRUE)
    if (is.null(design_cols)) design_cols <- NA_character_
    out$design_cols <- design_cols
    
    # ---- (A3) 可選：High/Low 敏感度（沿用舊函式）----
    if (isTRUE(RUN_HILO_SENSITIVITY)) {
      v <- sub_raw[is.finite(sub_raw)]
      if (length(v) >= (2 * min_per_group)) {
        qs <- stats::quantile(v, probs = c(hi_lo_quantile, 1 - hi_lo_quantile), names = FALSE)
        grp <- setNames(rep("Mid", length(sub_raw)), names(sub_raw))
        grp[names(v)[v <= qs[1]]] <- "Low"
        grp[names(v)[v >= qs[2]]] <- "High"
        keep_s <- names(sub_raw)[grp != "Mid" & is.finite(sub_raw)]
        nH_loc <- sum(grp[keep_s] == "High"); nL_loc <- sum(grp[keep_s] == "Low")
        
        if (nH_loc >= min_per_group && nL_loc >= min_per_group) {
          grp2 <- factor(ifelse(grp[keep_s] == "High","High","Low"), levels = c("Low","High"))
          mat_hilo <- mat[setdiff(rownames(mat), exclude_genes), keep_s, drop = FALSE]
          
          out$stats_hilo_raw <- tryCatch(
            limma_t_with_covars(mat_hilo, grp2),
            error = function(e){ log_msg("        limma Hi/Lo RAW 失敗：%s", e$message); NULL }
          )
          batch_sub2 <- if (!is.null(batch_subset)) droplevels(batch_subset[keep_s]) else NULL
          covars_sub2 <- if (!is.null(covars_all)) as.matrix(covars_all[keep_s, , drop = FALSE]) else NULL
          out$stats_hilo_adj <- tryCatch(
            limma_t_with_covars(mat_hilo, grp2, batch = batch_sub2, covars = covars_sub2),
            error = function(e){ log_msg("        limma Hi/Lo BatchAdj 失敗：%s", e$message); NULL }
          )
          out$nH <- nH_loc; out$nL <- nL_loc
        } else {
          log_msg("        Hi/Lo 不足（High=%d, Low=%d），略過敏感度", nH_loc, nL_loc)
        }
      }
    }
    
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
    if (all(vapply(stats_list, is.null, logical(1)))) return(invisible(NULL))
    out_dir <- file.path(out_roots[[ver_tag]], predictor_name)
    for (grp_name in names(genesets_by_group)) {
      subdir <- file.path(out_dir, safe_fs_name(grp_name))
      dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      
      if (!is.null(stats_list$limma_cont)) {
        run_fgsea_save(
          stats = stats_list$limma_cont,
          genesets = genesets_by_group[[grp_name]],
          out_prefix = file.path(subdir, "GSEA_limma_t_cont"),
          top_plot_n = 0,
          plot_title = glue::glue("{ds_id} | {predictor_name} | limma (continuous) | {grp_name} | {ver_tag}")
        )
      }
      if (!is.null(stats_list$limma_hilo)) {
        run_fgsea_save(
          stats = stats_list$limma_hilo,
          genesets = genesets_by_group[[grp_name]],
          out_prefix = file.path(subdir, "GSEA_limma_t_hilo"),
          top_plot_n = 0,
          plot_title = glue::glue("{ds_id} | {predictor_name} | limma (High vs Low) | {grp_name} | {ver_tag}")
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
  
  write_one(
    list(limma_cont = lim$stats_cont_raw,
         limma_hilo = if (RUN_HILO_SENSITIVITY) lim$stats_hilo_raw else NULL,
         spearman   = NULL),
    "RAW"
  )
  write_one(
    list(limma_cont = lim$stats_cont_adj,
         limma_hilo = if (RUN_HILO_SENSITIVITY) lim$stats_hilo_adj else NULL,
         spearman   = stats_rho_ADJ),
    "BatchAdj"
  )
  
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
    "limma-cont: design = ~ x (+ batch/SV or PC1-2 + optional purity/sex/age + TP53 in ALL); ranking = moderated t for 'x'",
    sprintf("limma_design_cols=%s", if (length(limma_design_cols) && all(!is.na(limma_design_cols))) paste(limma_design_cols, collapse=",") else "NA"),
    if (isTRUE(RUN_HILO_SENSITIVITY)) "limma-hilo: sensitivity; design = ~ 0 + group (+ covars)" else "limma-hilo: NOT RUN",
    "spearman: removeBatchEffect(...) then rho with continuous predictor; TP53 protected in SVA when ALL",
    sprintf("spearman_covar_cols=%s", if (length(spearman_covar_cols) && all(!is.na(spearman_covar_cols))) paste(spearman_covar_cols, collapse=",") else "NA")
  )
  
  cat(paste(info, collapse = "\n"), file = note_fp)
  
  invisible(NULL)
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

genesets_by_group <- list()
.get_subcat_col <- function(df){
  cand <- c("gs_subcat","subcategory","sub_category","gs_subcategory")
  cand[cand %in% names(df)][1]
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


## =========================================================
## 小函式：針對「指定樣本集合」執行一次完整 GSEA（修正版）
## 依賴的外部函式/物件（請先定義或載入）：
## 函式：
##   log_msg, impute_and_filter, audit_csn_score_feasibility, build_csn_score,
##   get_tp53_status, get_batch_factor, get_purity_covariate, get_sex_age_covariates,
##   residualize_vector, run_predictor_analyses, summarize_all_groups
## 全域參數/物件：
##   csn_subunits, min_frac_complete, min_per_group, RUN_HILO_SENSITIVITY
## 傳入參數：
##   genesets_by_group
## =========================================================

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
          batch_all = batch_all, purity_all = purity_all, sa_all = sa_all,
          tp53_num_all = tp53_num_all, is_ALL = is_ALL
        )
      }
    }
  }
  
  ## 5) 各版本各自彙整
  present_sub <- intersect(csn_subunits, rownames(mat0))
  sum_units <- c(present_sub, "CSN_SCORE", paste0("RESIDUAL_", present_sub))
  
  # RAW：只有 limma_t_cont
  ver_root <- file.path(out_root, "RAW")
  summarize_all_groups(
    out_root = ver_root,
    csn_subunits = sum_units,
    genesets_by_group = genesets_by_group,
    stat_tags = c("GSEA_limma_t_cont")
  )
  
  # BatchAdj：limma_t_cont + spearman (+ optional hilo)
  ver_root <- file.path(out_root, "BatchAdj")
  summarize_all_groups(
    out_root = ver_root,
    csn_subunits = sum_units,
    genesets_by_group = genesets_by_group,
    stat_tags = {
      st <- c("GSEA_limma_t_cont", "GSEA_spearman")
      if (isTRUE(RUN_HILO_SENSITIVITY)) st <- c(st, "GSEA_limma_t_hilo")
      st
    }
  )
  
  invisible(NULL)
}

## =========================================================
## 檢視目前取到哪些 PLEX（原始 vs 清理後）
## 依賴：load_matrix_from_dataset_dir, sanitize_batch_levels
## =========================================================
inspect_plex <- function(ds_dir, col = "TMT_PLEX",
                         pipe_policy   = BATCH_PIPE_POLICY,
                         min_per_level = BATCH_MIN_PER_LEVEL) {
  mat0 <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(mat0)
  meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
  stopifnot(file.exists(meta_fp))
  meta <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  id_cols <- intersect(c("SAMPLE_ID","sample_id","Sample_ID","Sample","sample"), names(meta))
  stopifnot(length(id_cols) > 0)
  meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
  meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
  raw <- as.character(meta[[col]])
  clean <- sanitize_batch_levels(raw, pipe_policy = pipe_policy, min_per_level = min_per_level)
  cat("\n==== ", basename(ds_dir), " | 欄位：", col, " ====\n", sep = "")
  cat("[原始 PLEX levels]：\n"); print(sort(table(raw), decreasing = TRUE))
  cat("\n[含 '|' 的原始值（前幾個）]：\n"); print(utils::head(unique(raw[grepl("\\|", raw %||% "")]), 10))
  cat("\n[清理後 PLEX levels]（pipe_policy = ", pipe_policy, 
      ", min_per_level = ", min_per_level, ")：\n", sep = "")
  print(sort(table(clean, useNA = "ifany"), decreasing = TRUE))
  df_map <- data.frame(SAMPLE_ID = sample_ids, raw_plex = raw, clean_plex = as.character(clean))
  cat("\n[樣本對照表（前 10 列）]\n"); print(utils::head(df_map, 10))
  invisible(list(raw_counts = sort(table(raw), decreasing = TRUE),
                 clean_counts = sort(table(clean, useNA = "ifany"), decreasing = TRUE),
                 map = df_map))
}

## =========================================================
## 取 purity（依資料集規則；named numeric, 0~1）
## 依賴：.to01, .median_from_semicolon
## =========================================================
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

## =========================================================
## 取 sex/age（sample 對齊；sex: 0/1；age: z-score；以 patient 檔為準）
## =========================================================
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
  
  samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
  pat  <- suppressMessages(readr::read_tsv(pat_fp,  show_col_types = FALSE, comment = "#")) |> as.data.frame()
  names(samp) <- .norm_names(names(samp)); names(pat) <- .norm_names(names(pat))
  
  sid      <- intersect(c("SAMPLE_ID","SAMPLE"), names(samp))[1]
  pid_samp <- intersect(c("PATIENT_ID","PATIENT"), names(samp))[1]
  pid_pat  <- intersect(c("PATIENT_ID","PATIENT"), names(pat))[1]
  if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) {
    stop("[get_sex_age_covariates] 無法建立 sample↔patient 對應（缺 SAMPLE_ID/PATIENT_ID）")
  }
  map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])
  
  # SEX（Male=1, Female=0；其餘 NA）
  sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  if ("SEX" %in% names(pat) || "GENDER" %in% names(pat)) {
    sex_col <- intersect(c("SEX","GENDER"), names(pat))[1]
    raw <- toupper(as.character(pat[[sex_col]]))
    val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
    names(val) <- as.character(pat[[pid_pat]])
    sex[] <- unname(val[ map_pt[sample_ids] ])
  }
  
  # AGE（優先 AGE，其次常見替代欄位）
  age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  age_col <- intersect(c("AGE","AGE_AT_DIAGNOSIS","AGE_AT_INDEX","AGE_YEARS"), names(pat))[1]
  if (!is.na(age_col)) {
    v <- suppressWarnings(as.numeric(pat[[age_col]]))
    names(v) <- as.character(pat[[pid_pat]])
    age_raw <- unname(v[ map_pt[sample_ids] ])
    age[] <- .zscore(age_raw)
  }
  
  cov_sex <- mean(is.finite(sex)) * 100
  cov_age <- mean(is.finite(age)) * 100
  log_msg("    covariates（patient→sample 對齊） coverage：sex %.1f%%, age %.1f%%", cov_sex, cov_age)
  
  data.frame(sex = as.numeric(sex), age = as.numeric(age),
             row.names = sample_ids, check.names = FALSE)
}

## =========================================================
## 共變量選擇（覆蓋率/與生物預測向量相關性過大則剔除）
## =========================================================
select_covars_safely <- function(df, sample_order, label = "",
                                 y = NULL,
                                 min_cov_named = c(purity = 0.60, sex = 0.80, age = 0.80),
                                 max_abs_cor = 0.30,
                                 min_pairs = 20) {
  if (is.null(df) || ncol(as.data.frame(df)) == 0) return(NULL)
  mat <- as.data.frame(df, check.names = FALSE)
  if (is.null(rownames(mat))) rownames(mat) <- sample_order
  idx <- match(sample_order, rownames(mat))
  if (anyNA(idx)) {
    log_msg("  [covars-%s] 警告：%d/%d 個樣本在共變數表中找不到 → 以 NA 補",
            label, sum(is.na(idx)), length(idx))
  }
  mat <- mat[idx, , drop = FALSE]
  rownames(mat) <- sample_order
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
    cover <- sum(is.finite(x)) / length(x)
    thr   <- if (cn %in% names(min_cov_named)) min_cov_named[[cn]] else min(min_cov_named)
    if (!is.finite(cover) || cover < thr) { log_msg("  [covars-%s] drop %s (coverage=%.0f%% < %.0f%%)", label, cn, 100*cover, 100*thr); next }
    fin <- is.finite(x)
    if (sum(fin) < 3 || sd(x[fin]) == 0) { log_msg("  [covars-%s] drop %s (near-constant)", label, cn); next }
    if (!is.null(y)) {
      ok <- fin & is.finite(y)
      if (sum(ok) >= min_pairs) {
        rshow <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
        if (is.finite(rshow) && abs(rshow) >= max_abs_cor) { log_msg("  [covars-%s] drop %s (|rho|=%.3f ≥ %.2f vs biology)", label, cn, abs(rshow), max_abs_cor); next }
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

## =========================================================
## 小稽核：sex/age + purity + batch 整體檢查
## 依賴：get_batch_factor, screen_batch_need（可選）
## =========================================================
.format_batch_sizes <- function(fac) {
  if (is.null(fac)) return(NA_character_)
  tb <- sort(table(fac), decreasing = TRUE)
  paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
}

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
  m <- load_matrix_from_dataset_dir(ds_dir)
  sample_ids <- colnames(m)
  
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
  
  pur <- .audit_purity_for_dataset(ds_id, ds_dir, sample_ids)
  
  bi <- get_batch_factor(ds_dir, sample_ids,
                         pipe_policy   = pipe_policy,
                         min_per_level = min_per_level)
  batch_col    <- if (!is.null(bi)) bi$name else "NONE"
  batch_levels <- if (!is.null(bi)) nlevels(bi$fac) else 0L
  batch_nonNA  <- if (!is.null(bi)) mean(!is.na(bi$fac)) * 100 else NA_real_
  batch_sizes  <- if (!is.null(bi)) .format_batch_sizes(bi$fac) else NA_character_
  
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

## =========================================================
## DO_AUDIT：需要時開啟（預設 FALSE）
## 依賴：screen_batch_need, dataset_dirs
## =========================================================
DO_AUDIT <- FALSE
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

## =========================================================
## run manifest（依賴：多個全域變數）
## =========================================================
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

## =========================================================
## 想從哪一個開始重跑（可改）
## =========================================================
start_from <- "brca_cptac_2020"
ord <- dataset_ids
ix  <- match(start_from, ord)
if (is.na(ix)) stop(sprintf("'%s' 不在 dataset_ids 裡", start_from))
dataset_dirs_run <- dataset_dirs[ ord[ix:length(ord)] ]

## =========================================================
## 依序處理本輪要跑的 datasets（TP53 三個分層）
## 依賴：load_matrix_from_dataset_dir, get_tp53_status, genesets_by_group
## =========================================================
for (ds in names(dataset_dirs_run)) {
  ds_dir <- dataset_dirs_run[[ds]]
  log_msg("== 開始資料集：%s ==", ds)
  
  mat0_full   <- load_matrix_from_dataset_dir(ds_dir)
  tp53_status <- get_tp53_status(ds_dir, colnames(mat0_full))
  
  samples_ALL <- colnames(mat0_full)
  samples_MUT <- names(tp53_status)[tp53_status == "TP53_mutant"]
  samples_WT  <- names(tp53_status)[tp53_status == "TP53_wild_type"]
  
  base_tp53_root <- file.path(ds_dir, "csn_gsea_results_TP53")
  
  run_one_stratum(
    ds_id = ds, ds_dir = ds_dir,
    mat0_full = mat0_full,
    sample_keep = samples_ALL,
    out_root = file.path(base_tp53_root, "ALL"),
    genesets_by_group = genesets_by_group
  )
  run_one_stratum(
    ds_id = ds, ds_dir = ds_dir,
    mat0_full = mat0_full,
    sample_keep = samples_MUT,
    out_root = file.path(base_tp53_root, "TP53_mutant"),
    genesets_by_group = genesets_by_group
  )
  run_one_stratum(
    ds_id = ds, ds_dir = ds_dir,
    mat0_full = mat0_full,
    sample_keep = samples_WT,
    out_root = file.path(base_tp53_root, "TP53_wild_type"),
    genesets_by_group = genesets_by_group
  )
  
  log_msg("== 完成資料集：%s（TP53 分層輸出 → %s）==", ds, base_tp53_root)
}

## =========================================================
## [NEW-1] 產出跨資料集的 meta-FDR（Stouffer → BH）
## 依賴：meta_fdr_stouffer
## =========================================================
meta_fdr_stouffer(
  dataset_dirs = dataset_dirs_run,
  strata = strata,                 # c("ALL","TP53_mutant","TP53_wild_type")
  stat_tags = c("GSEA_limma_t_cont","GSEA_spearman"),
  groups = names(genesets_by_group),
  out_root = "csn_gsea_pan_summary_TP53/meta_fdr"
)

## =========================================================
## [NEW-2] 插補敏感度示範：MinProb vs QRILC（單一 dataset/stratum/predictor）
## 依賴：run_imputation_sensitivity_demo
## =========================================================
.demo_ds <- "luad_cptac_2020"
if (!is.null(dataset_dirs_run[[.demo_ds]])) {
  run_imputation_sensitivity_demo(
    ds_id = .demo_ds,
    ds_dir = dataset_dirs_run[[.demo_ds]],
    stratum = "ALL",
    predictor_name = "CSN_SCORE", # 也可用 "COPS5" 或 "RESIDUAL_COPS5"
    geneset_group = "H"           # Hallmark
  )
} else {
  log_msg("[impute-demo] 指定示範隊列 %s 不在 dataset_dirs_run 中，略過示範", .demo_ds)
}

## =========================================================
## PAN：跨癌別 × 跨 CSN 次單元（TP53 分層版；RAW/BatchAdj 各自做）
## 依賴：safe_fs_name, csn_subunits
## =========================================================
versions <- c("RAW","BatchAdj")
groups    <- names(genesets_by_group)
dataset_ids_have <- names(dataset_dirs_run)

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

for (ver in versions) {
  log_msg("PAN(TP53,%s) 聚合開始", ver)
  
  ## 依版本設定正確的 stat_tags
  stat_tags_ver <- if (ver == "RAW") {
    c("GSEA_limma_t_cont")
  } else {
    c("GSEA_limma_t_cont", "GSEA_spearman",
      if (isTRUE(RUN_HILO_SENSITIVITY)) "GSEA_limma_t_hilo" else NULL)
  }
  
  long_df <- collect_all_long_tp53(
    dataset_ids = dataset_ids_have,
    strata      = strata,
    versions    = ver,
    groups      = groups,
    stat_tags   = stat_tags_ver,
    subunits    = csn_subunits
  )
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
  
  data.table::fwrite(per_ds_term, file.path(pan_out_root, "per_dataset_summary.csv"))
  wb_top <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_top, "per_dataset")
  openxlsx::writeData(wb_top, "per_dataset", per_ds_term)
  
  # pan-cancer 匯總
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
  
  ## 每個 stratum × gene-set group × stat 的輸出（表1/表2）
  for (stt in strata) {
    for (g in groups) {
      for (st in stat_tags_ver) {
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
  
  if (exists("pan_out_root")) {
    message(sprintf("[PAN-TP53] 整合完成（%s）→ %s", ver, pan_out_root))
  }
}
