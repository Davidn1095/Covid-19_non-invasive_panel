# ===============================================================
# External-only heatmap, value ± SD
# Split boundary = 2020-04-15, tuning + threshold metric = MCC
# WHO_score_admission_mod: controls' NA -> 0
# Calibration removed, raw probabilities are used
# C4.5 fixed at J48 defaults (no fine tuning)
# Seed fixed to 123
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(kknn)
})

# ---------------- config ----------------
cfg <- list(
  file_path     = "biomarkers_acuteCOVID_meta.xlsx",
  sheet         = "meta",
  outcome       = "Hospital_ID",
  pos_label     = "Yes",
  neg_label     = "No",
  boundary_date = as.Date("2020-04-15"),
  inner_k       = 3,
  tune_len      = 3,
  seed_cv       = 123,
  B_boot        = 1000
)

# ---------------- features ----------------
feat_triage <- c("Diagnosis","WHO_score_admission_mod","Age","Gender","SpO2_admission")
feat_full <- c(
  feat_triage,
  "CRP","D_Dimer","albumin",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc"
)

# ---------------- labels, order, palette ----------------
metric_order <- c("MCC","AUC-ROC","F1-score","Accuracy","Precision","Sensitivity","Specificity")
algo_order   <- c("C4.5","k-NN","SVM","RF","LR")
pal_red_green <- c("#F4CCCC","#D9EAD3")

# ---------------- helpers ----------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))
sanitize_name <- function(x) gsub("[^A-Za-z0-9]", "", tolower(x))

to_YN <- function(x, pos = cfg$pos_label, neg = cfg$neg_label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(factor(ifelse(x, pos, neg), levels = c(pos, neg)))
  if (is.numeric(x)) {
    stopifnot(all(x %in% c(0,1), na.rm = TRUE))
    return(factor(ifelse(x == 1, pos, neg), levels = c(pos, neg)))
  }
  if (is.character(x)) {
    s <- trimws(tolower(x))
    map_yes <- c("1","yes","y","true","pos","positive")
    map_no  <- c("0","no","n","false","neg","negative")
    out <- ifelse(s %in% map_yes, pos, ifelse(s %in% map_no, neg, NA_character_))
    if (any(is.na(out))) stop("Outcome not coercible to Yes/No")
    return(factor(out, levels = c(pos, neg)))
  }
  stop("Unsupported outcome type")
}

parse_excel_date <- function(v) {
  if (inherits(v, "Date")) return(v)
  if (inherits(v, "POSIXct") || inherits(v, "POSIXt")) return(as.Date(v))
  if (is.numeric(v)) return(as.Date(v, origin = "1899-12-30"))
  if (is.character(v)) {
    fmts <- c("%d/%m/%Y","%Y-%m-%d","%m/%d/%Y","%d-%m-%Y","%d.%m.%Y",
              "%d/%m/%Y %H:%M","%Y-%m-%d %H:%M","%m/%d/%Y %H:%M")
    for (f in fmts) {
      d <- suppressWarnings(as.Date(v, format = f))
      if (!all(is.na(d))) return(d)
    }
  }
  as.Date(NA)
}

resolve_date_col <- function(df) {
  sn <- sanitize_name(names(df))
  prefs <- c("samplingdate","admissiondate","date","sampling_date","acqdate","acq_date")
  for (p in prefs) {
    hit <- which(sn == p)
    if (length(hit)) return(names(df)[hit[1]])
  }
  contains <- which(grepl("sampling|admission|date", sn))
  if (length(contains)) return(names(df)[contains[1]])
  stop("A date column was not found")
}

# ---------------- data loader + cohort filter ----------------
load_data <- function(path, sheet, outcome, keep_predictors) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet = sheet) %>% as.data.frame()
  
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df) == grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  for (v in intersect(c("group","Gender","Diagnosis","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' missing", outcome))
  df[[outcome]] <- to_YN(df[[outcome]])
  
  if ("WHO_score_admission_mod" %in% names(df)) {
    if (!is.numeric(df$WHO_score_admission_mod)) df$WHO_score_admission_mod <- as_num(df$WHO_score_admission_mod)
    idx_ctrl_na <- with(df, is.na(WHO_score_admission_mod) & group == "CTRL_noCOVID")
    df$WHO_score_admission_mod[idx_ctrl_na] <- 0
  }
  
  num_candidates <- c("Age","SpO2_admission","CRP","D_Dimer","albumin",
                      "monocyte_abs_number","monocytes_perc",
                      "lymphocyte_abs_number","lymphocytes_perc",
                      "neutrophil_abs_number","neutrophils_perc",
                      "WHO_score_admission_mod")
  for (v in intersect(num_candidates, names(df))) {
    if (!is.numeric(df[[v]])) df[[v]] <- as_num(df[[v]])
  }
  
  dcol <- resolve_date_col(df)
  message(sprintf("[date] Using '%s' as the admission or sampling date.", dcol))
  df[[dcol]] <- parse_excel_date(df[[dcol]])
  
  keepx <- unique(c(keep_predictors, outcome, "group", dcol, "data_split"))
  keepx <- intersect(keepx, names(df))
  list(df = df[, keepx, drop = FALSE], admit_col = dcol)
}

# ---------------- temporal split ----------------
enforce_temporal_split <- function(df, admit_col, boundary_date = cfg$boundary_date) {
  stopifnot(admit_col %in% names(df))
  ad <- parse_excel_date(df[[admit_col]])
  df$data_split <- factor(ifelse(ad <= boundary_date, "train", "external"),
                          levels = c("train","external"))
  message("\n[Temporal split]\n  Boundary date: ", format(boundary_date), "\n")
  print(table(df$data_split, useNA = "ifany"))
  df
}

# ---------------- recipes (method-specific) ----------------
make_recipe <- function(dat, yvar, method) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c("data_split"), names(dat))
  if (length(ign)) rec <- rec %>% update_role(all_of(ign), new_role = "ignore")
  
  if (identical(method, "J48")) {
    rec %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_impute_mode(all_nominal_predictors())  %>%
      step_novel(all_nominal_predictors())        %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_zv(all_predictors())
  } else {
    rec %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_impute_mode(all_nominal_predictors())  %>%
      step_novel(all_nominal_predictors())        %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_dummy(all_nominal_predictors())        %>%
      step_zv(all_predictors())                   %>%
      step_normalize(all_numeric_predictors())
  }
}

# ---------------- metrics ----------------
mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(NA_real_)
  num/den
}

binary_auc <- function(obs, p_pos, pos_label = cfg$pos_label) {
  y <- as.integer(obs == pos_label)
  if (length(unique(y)) < 2) return(NA_real_)
  suppressMessages(as.numeric(pROC::auc(y, p_pos, quiet = TRUE)))
}

compute_metrics_binary <- function(obs, pred, p_pos, pos_label = cfg$pos_label, neg_label = cfg$neg_label) {
  y    <- factor(obs,  levels = c(neg_label, pos_label))
  pred <- factor(pred, levels = c(neg_label, pos_label))
  tab  <- table(y, pred)
  TP <- as_num(tab[pos_label, pos_label] %||% 0)
  TN <- as_num(tab[neg_label, neg_label] %||% 0)
  FP <- as_num(tab[neg_label, pos_label] %||% 0)
  FN <- as_num(tab[pos_label, neg_label] %||% 0)
  preci <- if ((TP+FP)==0) NA_real_ else TP/(TP+FP)
  sens  <- if ((TP+FN)==0) NA_real_ else TP/(TP+FN)
  spec  <- if ((TN+FP)==0) NA_real_ else TN/(TN+FP)
  acc   <- (TP+TN)/sum(tab)
  f1    <- if (is.na(preci) || is.na(sens) || (preci+sens)==0) NA_real_ else 2*preci*sens/(preci+sens)
  mcc   <- mcc_from_counts(TP, FP, FN, TN)
  aucv  <- binary_auc(y, as_num(p_pos), pos_label)
  c(MCC = mcc, AUC = aucv, `F1` = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec)
}

# ---------------- MCC threshold on a dense grid ----------------
best_threshold_mcc <- function(obs, p_pos, pos = cfg$pos_label, neg = cfg$neg_label, grid = seq(0.01,0.99,by=0.001)) {
  y <- factor(obs, levels = c(neg, pos))
  best_t <- 0.5; best_m <- -Inf
  for (t in grid) {
    pred <- factor(ifelse(p_pos >= t, pos, neg), levels = c(neg,pos))
    tab <- table(y,pred)
    TP <- as_num(tab[pos,pos] %||% 0); TN <- as_num(tab[neg,neg] %||% 0)
    FP <- as_num(tab[neg,pos] %||% 0); FN <- as_num(tab[pos,neg] %||% 0)
    m  <- mcc_from_counts(TP, FP, FN, TN)
    if (is.finite(m) && m > best_m) { best_m <- m; best_t <- t }
  }
  list(t = best_t, mcc = best_m)
}

# ---------------- model registry ----------------
algo_requires <- function(method) {
  switch(method,
         "glmnet"    = "glmnet",
         "rf"        = "randomForest",
         "svmRadial" = "kernlab",
         "kknn"      = "kknn",
         "J48"       = "RWeka",
         NULL)
}
method_available <- function(method) {
  pkg <- algo_requires(method)
  if (is.null(pkg)) return(TRUE)
  requireNamespace(pkg, quietly = TRUE)
}

# C4.5 fixed, no tuning, others keep tune_len
model_specs <- function(tune_len, algos = c("C4.5","k-NN","SVM","RF","LR")) {
  all <- list(
    "LR"   = list(name = "LR",   method = "glmnet",    tuneLength = tune_len),
    "RF"   = list(name = "RF",   method = "rf",        tuneLength = tune_len),
    "SVM"  = list(name = "SVM",  method = "svmRadial", tuneLength = tune_len),
    "k-NN" = list(name = "k-NN", method = "kknn",      tuneLength = tune_len),
    "C4.5" = list(
      name   = "C4.5",
      method = "J48",
      grid   = data.frame(C = 0.25, M = 2)  # Weka defaults
    )
  )
  avail <- names(caret::getModelInfo())
  keep  <- names(all)[vapply(all, \(s) s$method %in% avail && method_available(s$method), logical(1))]
  out   <- all[keep]
  out[names(out)[order(match(names(out), algo_order))]]
}

# ---------------- custom MCC summary for caret tuning ----------------
mccSummary <- function(data, lev = NULL, model = NULL) {
  pos <- cfg$pos_label; neg <- cfg$neg_label
  obs  <- factor(data$obs,  levels = c(neg, pos))
  pred <- factor(data$pred, levels = c(neg, pos))
  tab <- table(obs, pred)
  TP <- as_num(tab[pos, pos] %||% 0)
  TN <- as_num(tab[neg, neg] %||% 0)
  FP <- as_num(tab[neg, pos] %||% 0)
  FN <- as_num(tab[pos, neg] %||% 0)
  mcc <- mcc_from_counts(TP, FP, FN, TN)
  c(MCC = mcc)
}

# ---------------- shared CV folds across algorithms ----------------
make_cv_indices <- function(y, k = cfg$inner_k, seed = cfg$seed_cv) {
  set.seed(seed)
  idx_tr <- caret::createFolds(y, k = k, returnTrain = TRUE)
  idx_te <- lapply(idx_tr, function(tr) setdiff(seq_along(y), tr))
  list(index = idx_tr, indexOut = idx_te)
}

# ---------------- training and holdout, MCC everywhere ----------------
train_inner <- function(dat_train, yvar, spec, inner_k, cv_idx) {
  rec <- make_recipe(dat_train, yvar, method = spec$method)
  tr_ctrl <- caret::trainControl(
    method = "cv", number = inner_k,
    classProbs = TRUE,
    summaryFunction = mccSummary,
    savePredictions = "final", allowParallel = TRUE,
    index = cv_idx$index, indexOut = cv_idx$indexOut
  )
  if (!is.null(spec$grid)) {
    caret::train(rec, data = dat_train, method = spec$method,
                 trControl = tr_ctrl, metric = "MCC", tuneGrid = spec$grid)
  } else {
    caret::train(rec, data = dat_train, method = spec$method,
                 trControl = tr_ctrl, metric = "MCC", tuneLength = spec$tuneLength)
  }
}

run_holdout <- function(df, yvar, features, algo_name) {
  specs <- model_specs(cfg$tune_len, algos = algo_name)
  if (!length(specs)) return(NULL)
  spec  <- specs[[algo_name]]
  
  use_cols <- unique(c(features, yvar, "group", "data_split"))
  use_cols <- intersect(use_cols, names(df))
  dtrain <- df[df$data_split == "train", use_cols, drop = FALSE]
  dext   <- df[df$data_split == "external", use_cols, drop = FALSE]
  if (!nrow(dtrain) || !nrow(dext)) return(NULL)
  
  cv_idx <- make_cv_indices(dtrain[[yvar]], k = cfg$inner_k, seed = cfg$seed_cv)
  fit <- train_inner(dtrain, yvar, spec, cfg$inner_k, cv_idx)
  
  p_tr <- as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label])
  thr_info <- best_threshold_mcc(dtrain[[yvar]], p_tr,
                                 pos = cfg$pos_label, neg = cfg$neg_label)
  thr <- thr_info$t
  
  p_ex <- as.numeric(predict(fit, newdata = dext, type = "prob")[, cfg$pos_label])
  neg  <- cfg$neg_label
  pred <- factor(ifelse(p_ex >= thr, cfg$pos_label, neg), levels = c(neg, cfg$pos_label))
  
  list(
    obs = factor(dext[[yvar]], levels = c(neg, cfg$pos_label)),
    p = p_ex, pred = pred, threshold = thr
  )
}

# ---------------- bootstrap SD on external ----------------
bootstrap_sd_all_metrics <- function(obs, pred, p, B = cfg$B_boot, seed = cfg$seed_cv + 99) {
  set.seed(seed)
  n <- length(obs)
  mets <- c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")
  M <- matrix(NA_real_, nrow = B, ncol = length(mets), dimnames = list(NULL, mets))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    m <- compute_metrics_binary(obs[idx], pred[idx], p[idx],
                                pos_label = cfg$pos_label, neg_label = cfg$neg_label)
    M[b, "MCC"]         <- as_num(m["MCC"])
    M[b, "AUC"]         <- as_num(m["AUC"])
    M[b, "F1"]          <- as_num(m["F1"])
    M[b, "Accuracy"]    <- as_num(m["Accuracy"])
    M[b, "Precision"]   <- as_num(m["Precision"])
    M[b, "Sensitivity"] <- as_num(m["Sensitivity"])
    M[b, "Specificity"] <- as_num(m["Specificity"])
  }
  apply(M, 2, sd, na.rm = TRUE)
}

# ---------------- External summary for heatmap (value + sd) ----------------
external_summary_for_heatmap <- function(df, yvar, feat_full, feat_triage) {
  specs <- model_specs(cfg$tune_len)
  
  build_set <- function(set_feats, set_name){
    out <- list()
    for (algo in names(specs)) {
      o <- run_holdout(df, yvar, set_feats, algo)
      if (is.null(o)) next
      
      mets_point <- compute_metrics_binary(o$obs, o$pred, o$p,
                                           pos_label = cfg$pos_label,
                                           neg_label = cfg$neg_label)
      sds <- bootstrap_sd_all_metrics(o$obs, o$pred, o$p)
      
      vals <- c(
        MCC         = as.numeric(mets_point["MCC"]),
        `AUC-ROC`   = as.numeric(mets_point["AUC"]),
        `F1-score`  = as.numeric(mets_point["F1"]),
        Accuracy    = as.numeric(mets_point["Accuracy"]),
        Precision   = as.numeric(mets_point["Precision"]),
        Sensitivity = as.numeric(mets_point["Sensitivity"]),
        Specificity = as.numeric(mets_point["Specificity"])
      )
      sds2 <- c(
        MCC         = sds["MCC"],
        `AUC-ROC`   = sds["AUC"],
        `F1-score`  = sds["F1"],
        Accuracy    = sds["Accuracy"],
        Precision   = sds["Precision"],
        Sensitivity = sds["Sensitivity"],
        Specificity = sds["Specificity"]
      )
      
      out[[algo]] <- tibble::tibble(
        Algorithm = algo,
        Metric    = names(vals),
        Mean      = as.numeric(vals),
        SD        = as.numeric(sds2),
        Set       = set_name,
        Split     = "External"
      )
    }
    dplyr::bind_rows(out)
  }
  
  dplyr::bind_rows(
    build_set(feat_full,   "Complete feature set"),
    build_set(feat_triage, "Triage feature set")
  )
}

# ---------------- External-only heatmap ----------------
plot_external_heatmap <- function(ext_df) {
  df <- ext_df %>%
    dplyr::mutate(
      Metric    = factor(Metric, levels = rev(metric_order), ordered = TRUE),
      Algorithm = factor(Algorithm, levels = algo_order, ordered = TRUE),
      Set       = factor(Set, levels = c("Complete feature set","Triage feature set"), ordered = TRUE),
      Mean_clip = pmin(pmax(Mean, 0), 1)
    ) %>%
    dplyr::group_by(Set, Metric) %>%
    dplyr::mutate(
      is_max = if (all(is.na(Mean_clip))) FALSE else Mean_clip == max(Mean_clip, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      is_max = dplyr::coalesce(is_max, FALSE),
      LabelText = sprintf("%.4f \u00B1 %.4f", Mean_clip, dplyr::coalesce(SD, 0))
    )
  
  ggplot2::ggplot(df, ggplot2::aes(Algorithm, Metric, fill = Mean_clip)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = LabelText, fontface = ifelse(is_max, "bold", "plain")),
      size = 3
    ) +
    ggplot2::scale_fill_gradient(limits = c(0,1),
                                 low = pal_red_green[1],
                                 high = pal_red_green[2],
                                 name = "Score") +
    ggplot2::facet_wrap(~Set, ncol = 1) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing.y = grid::unit(1.0, "lines")
    )
}

# ========================= DRIVER =========================
set.seed(123)

ld <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
df <- enforce_temporal_split(ld$df, admit_col = ld$admit_col, boundary_date = cfg$boundary_date)
yvar <- cfg$outcome

ext_df <- external_summary_for_heatmap(df, yvar, feat_full, feat_triage)
p_ext  <- plot_external_heatmap(ext_df)
print(p_ext)

# ggsave("external_heatmap_mcc_everywhere_boundary_2020-04-15.png", p_ext, width = 10, height = 8, dpi = 300, bg = "transparent")
