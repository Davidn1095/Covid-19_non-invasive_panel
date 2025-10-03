# ============================================================
# PREDICTIONS: full vs triage, binary metrics for Yes/No, no disk writes
# ============================================================

suppressPackageStartupMessages({
  library(readxl);  library(dplyr);  library(tidyr);   library(tibble)
  library(recipes); library(caret);  library(pROC);    library(purrr)
  library(stringr); library(RWeka);  library(kernlab); library(randomForest)
  library(glmnet);  library(kknn)
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflict_prefer("select","dplyr", quiet = TRUE)
    conflicted::conflict_prefer("filter","dplyr", quiet = TRUE)
    conflicted::conflict_prefer("mutate","dplyr", quiet = TRUE)
    conflicted::conflict_prefer("lag",   "dplyr", quiet = TRUE)
  }
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
set.seed(123)

# ---------------- config ----------------
cfg <- list(
  fast_mode = TRUE,
  cv_k      = 10,
  cv_R      = 1,
  inner_k   = 5,
  tune_len  = 5,
  print_fold_diag = TRUE,
  algos     = c("LR","RF","SVM","k-NN","C4.5"),
  file_path = "biomarkers_acuteCOVID_meta.xlsx",
  sheet     = "meta",
  outcome   = "Hospital_ID",
  pos_label = "Yes",
  seed_cv   = 444
)
if (cfg$fast_mode) {
  cfg$cv_k    <- 3
  cfg$inner_k <- 3
  cfg$tune_len<- 3
}

# exact predictors
feat_full <- c(
  "Diagnosis","severity_admission","Age","Gender",
  "SpO2_admission","albumin","CRP","D_Dimer",
  "monocyte_abs_number","neutrophil_abs_number","lymphocyte_abs_number",
  "monocytes_perc","neutrophils_perc","lymphocytes_perc"
)
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- data ----------------
load_data <- function(path, sheet, outcome, keep) {
  df <- readxl::read_excel(path, sheet = sheet) |>
    dplyr::select(dplyr::any_of(c(outcome, keep))) |>
    dplyr::mutate(
      Hosp_Bin = factor(dplyr::case_when(
        .data[[outcome]] %in% c("Yes","No") ~ as.character(.data[[outcome]]),
        TRUE ~ NA_character_
      ), levels = c("No","Yes"))
    ) |>
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) |>
    dplyr::relocate(Hosp_Bin, .before = dplyr::everything()) |>
    droplevels()
  df$.rid <- seq_len(nrow(df))
  df
}

# ---------------- recipe ----------------
make_recipe <- function(dat, yvar) {
  recipe(stats::as.formula(paste(yvar, "~ .")), data = dat) |>
    update_role(.rid, new_role = "id") |>
    step_impute_median(all_numeric_predictors()) |>
    step_impute_mode(all_nominal_predictors()) |>
    step_novel(all_nominal_predictors(), new_level = "NOVEL") |>
    step_other(all_nominal_predictors(), threshold = 0.01, other = "OTHER") |>
    step_zv(all_predictors()) |>
    step_YeoJohnson(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_dummy(all_nominal_predictors(), one_hot = TRUE)
}

# ---------------- binary metrics (pos = "Yes", neg = "No") ----------------
binary_auc <- function(obs, p_pos, pos_label = "Yes") {
  y <- as.integer(obs == pos_label)
  if (length(unique(y)) < 2 || is.null(p_pos) || all(is.na(p_pos))) return(NA_real_)
  r <- try(pROC::roc(y, as.numeric(p_pos), quiet = TRUE), silent = TRUE)
  if (inherits(r, "try-error")) NA_real_ else as.numeric(pROC::auc(r))
}
compute_metrics_binary <- function(obs, pred, p_pos, pos_label = "Yes", neg_label = "No") {
  y    <- factor(obs,  levels = c(neg_label, pos_label))
  yhat <- factor(pred, levels = c(neg_label, pos_label))
  cm <- table(y, yhat)
  TP <- cm[pos_label, pos_label]; TN <- cm[neg_label, neg_label]
  FP <- cm[neg_label, pos_label]; FN <- cm[pos_label, neg_label]
  eps <- 1e-9
  Sensitivity <- as.numeric(TP / pmax(TP + FN, eps))   # recall for Yes
  Specificity <- as.numeric(TN / pmax(TN + FP, eps))   # TNR for No
  Precision   <- as.numeric(TP / pmax(TP + FP, eps))
  F1          <- as.numeric(2 * TP / pmax(2 * TP + FP + FN, eps))
  Accuracy    <- as.numeric((TP + TN) / pmax(sum(cm), eps))
  MCC         <- as.numeric((TP*TN - FP*FN) / sqrt(pmax((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN), eps)))
  AUC         <- binary_auc(y, p_pos, pos_label)
  c(MCC = MCC, F1 = F1, Accuracy = Accuracy, Precision = Precision,
    Sensitivity = Sensitivity, Specificity = Specificity, AUC = AUC)
}
metricSummaryInner_binary <- function(data, lev = NULL, model = NULL, pos_label = "Yes") {
  lev <- lev %||% levels(data$obs)
  pcol <- if (pos_label %in% colnames(data)) pos_label else if (length(lev) == 2) lev[2] else NA_character_
  p_pos <- if (!is.na(pcol) && pcol %in% colnames(data)) data[[pcol]] else NULL
  auc <- binary_auc(data$obs, p_pos, pos_label = pcol %||% pos_label)
  acc <- mean(data$obs == data$pred, na.rm = TRUE)
  c(AUC = auc, Accuracy = acc)
}

# ---------------- model registry ----------------
model_specs <- function(tune_len, algos) {
  all <- list(
    "C4.5" = list(method = "J48",       tuneLength = NULL),
    "k-NN" = list(method = "kknn",      tuneLength = tune_len),
    "SVM"  = list(method = "svmRadial", tuneLength = tune_len),
    "RF"   = list(method = "rf",        tuneLength = tune_len),
    "LR"   = list(method = "glmnet",    tuneLength = tune_len)
  )
  all[names(all) %in% algos]
}

# ---------------- CV utilities ----------------
choose_k_for_task <- function(y, k_desired) {
  tab <- table(y); max_k <- max(2, min(tab))
  kk <- min(k_desired, max_k)
  if (kk < k_desired) message(sprintf("Folds reduced to %d for class balance", kk))
  kk
}
build_cv_splits <- function(y, R, k_desired, seed = 999) {
  set.seed(seed)
  splits <- list()
  for (r in seq_len(R)) {
    kk <- choose_k_for_task(y, k_desired)
    fold_list <- caret::createFolds(y, k = kk, returnTrain = FALSE)
    for (j in seq_len(kk)) {
      test_idx  <- fold_list[[j]]
      train_idx <- setdiff(seq_along(y), test_idx)
      splits[[paste0("r", r, "_f", j)]] <- list(train_idx = train_idx, test_idx = test_idx)
    }
  }
  splits
}

# ---------------- training wrappers ----------------
train_inner <- function(dat_train, yvar, spec, inner_k) {
  rec <- make_recipe(dat_train, yvar)
  tr_ctrl <- caret::trainControl(
    method = "cv", number = inner_k,
    classProbs = TRUE, savePredictions = "final",
    summaryFunction = metricSummaryInner_binary
  )
  fit <- suppressWarnings(tryCatch({
    if (!is.null(spec$tuneLength)) {
      caret::train(rec, data = dat_train, method = spec$method,
                   tuneLength = spec$tuneLength, trControl = tr_ctrl, metric = "AUC")
    } else {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = tr_ctrl, metric = "AUC")
    }
  }, error = function(e) NULL))
  if (is.null(fit)) return(NULL)
  pred <- fit$pred
  if (!is.null(fit$bestTune) && ncol(fit$bestTune) > 0) {
    pred <- dplyr::semi_join(pred, fit$bestTune, by = names(fit$bestTune))
  }
  list(fit = fit, oof = pred)
}
train_final <- function(dat_train, yvar, spec, bestTune, inner_k) {
  rec <- make_recipe(dat_train, yvar)
  tr_none <- caret::trainControl(method = "none", classProbs = TRUE,
                                 summaryFunction = metricSummaryInner_binary)
  suppressWarnings(tryCatch({
    if (!is.null(bestTune) && ncol(bestTune) > 0) {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = tr_none, metric = "AUC", tuneGrid = bestTune)
    } else if (!is.null(spec$tuneLength)) {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = caret::trainControl(method="cv", number=inner_k,
                                                   classProbs=TRUE, savePredictions="final",
                                                   summaryFunction = metricSummaryInner_binary),
                   tuneLength = spec$tuneLength, metric = "AUC")
    } else {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = tr_none, metric = "AUC")
    }
  }, error = function(e) NULL))
}

# ---------------- single fold eval (uses binary metrics) ----------------
eval_fold_binary <- function(dat_tr, dat_va, spec, yvar, pos_label, algo_name, fold_tag, inner_k, verbose) {
  lev <- levels(dat_tr[[yvar]])
  neg_label <- setdiff(lev, pos_label)[1]
  
  inner <- train_inner(dat_tr, yvar, spec, inner_k)
  if (is.null(inner) || !all(lev %in% colnames(inner$oof))) return(NULL)
  
  t_res <- best_threshold_mcc(
    inner$oof$obs,
    as.numeric(inner$oof[, pos_label]),
    pos = pos_label,
    neg = neg_label
  )
  bestTune <- inner$fit$bestTune %||% NULL
  fit <- train_final(dat_tr, yvar, spec, bestTune, inner_k)
  if (is.null(fit)) return(NULL)
  
  pr <- try(predict(fit, newdata = dat_va, type = "prob"), silent = TRUE)
  if (inherits(pr, "try-error") || is.null(pr) || !all(lev %in% colnames(pr))) return(NULL)
  
  pred <- factor(ifelse(pr[[pos_label]] >= t_res$t, pos_label, neg_label), levels = lev)
  
  m <- compute_metrics_binary(
    obs   = dat_va[[yvar]],
    pred  = pred,
    p_pos = pr[[pos_label]],
    pos_label = pos_label,
    neg_label = neg_label
  )
  
  if (isTRUE(verbose)) {
    cat(sprintf("[Fold %s] %s | %s | tuned.t=%.3f | pos=%s, neg=%s\n",
                fold_tag, yvar, algo_name, t_res$t, pos_label, neg_label))
  }
  m
}

# ---------------- unified CV runner ----------------
run_cv <- function(df, yvar, features, algos, type = "binary",
                   pos_label = NULL, cv_R = 1, cv_k = 10, inner_k = 5,
                   tune_len = 5, seed = 444, verbose = TRUE,
                   return = c("folds","agg")) {
  return <- match.arg(return)
  dat <- df |>
    dplyr::filter(!is.na(.data[[yvar]])) |>
    dplyr::select(dplyr::all_of(c(yvar, features, ".rid")))
  splits <- build_cv_splits(dat[[yvar]], R = cv_R, k_desired = cv_k, seed = seed)
  specs  <- model_specs(tune_len, algos)
  
  rows <- list()
  for (algo in names(specs)) {
    spec <- specs[[algo]]
    for (nm in names(splits)) {
      tr_idx <- splits[[nm]]$train_idx
      va_idx <- splits[[nm]]$test_idx
      dat_tr <- dat[tr_idx, , drop = FALSE]
      dat_va <- dat[va_idx, , drop = FALSE]
      m <- switch(type,
                  binary = eval_fold_binary(dat_tr, dat_va, spec, yvar, pos_label, algo, nm, inner_k, verbose),
                  stop("Only binary implemented"))
      if (is.null(m)) next
      rows <- c(rows, list(tibble::tibble(Algorithm = algo, Fold = nm,
                                          !!!as.list(m))))
    }
  }
  out <- if (length(rows)) dplyr::bind_rows(rows) else tibble::tibble()
  if (!nrow(out)) return(out)
  out_long <- out |>
    tidyr::pivot_longer(cols = -c(Algorithm, Fold), names_to = "Metric", values_to = "Value")
  
  if (return == "folds") return(out_long)
  
  out_long |>
    dplyr::group_by(Algorithm, Metric) |>
    dplyr::summarise(Mean = mean(Value, na.rm = TRUE),
                     SD   = stats::sd(Value, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(factor(Metric, levels = c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")))
}

# ---------------- reporting helper ----------------
pretty_task_table <- function(agg_tbl, title) {
  cat(sprintf("\n=== %s METRICS (mean ± SD) ===\n", toupper(title)))
  if (!nrow(agg_tbl)) return(cat("(no results)\n"))
  pretty <- agg_tbl |>
    tidyr::pivot_wider(names_from = Algorithm,
                       values_from = c(Mean, SD)) |>
    tidyr::pivot_longer(cols = -Metric, names_to = c("Stat","Algorithm"),
                        names_sep = "_", values_to = "Val") |>
    tidyr::pivot_wider(names_from = Stat, values_from = Val) |>
    dplyr::mutate(Value = sprintf("%.4f ± %.4f", Mean, SD)) |>
    dplyr::select(Metric, Algorithm, Value) |>
    tidyr::pivot_wider(names_from = Algorithm, values_from = Value)
  print(pretty, n = nrow(pretty))
}

# ---------------- main ----------------
df <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
yvar <- "Hosp_Bin"

full_folds <- run_cv(df, yvar, features = feat_full, algos = cfg$algos,
                     type = "binary", pos_label = cfg$pos_label,
                     cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
                     tune_len = cfg$tune_len, seed = cfg$seed_cv,
                     verbose = cfg$print_fold_diag, return = "folds")
full_agg <- run_cv(df, yvar, features = feat_full, algos = cfg$algos,
                   type = "binary", pos_label = cfg$pos_label,
                   cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
                   tune_len = cfg$tune_len, seed = cfg$seed_cv,
                   verbose = FALSE, return = "agg")

triage_folds <- run_cv(df, yvar, features = feat_triage, algos = cfg$algos,
                       type = "binary", pos_label = cfg$pos_label,
                       cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
                       tune_len = cfg$tune_len, seed = cfg$seed_cv,
                       verbose = cfg$print_fold_diag, return = "folds")
triage_agg <- run_cv(df, yvar, features = feat_triage, algos = cfg$algos,
                     type = "binary", pos_label = cfg$pos_label,
                     cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
                     tune_len = cfg$tune_len, seed = cfg$seed_cv,
                     verbose = FALSE, return = "agg")

pretty_task_table(full_agg,   "Hospitalization baseline, no resampling")
pretty_task_table(triage_agg, "Hospitalization triage, bedside only")

# Objects in memory: full_folds, triage_folds, full_agg, triage_agg

