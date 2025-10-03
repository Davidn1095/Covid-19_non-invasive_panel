# ============================================================
# Clinical prediction on acute COVID biomarkers
# Tasks:
# 1) Diagnosis, binary, MCC-thresholded on inner OOF (pos = COVID)
# 2) Severity at admission, 3-class, argmax
# 3) Hospitalization, binary, MCC-thresholded on inner OOF (pos = Yes)
# 4) Worst severity, 3-class, argmax
# Same algorithms, same CV, same metrics
# ============================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(tibble)
  library(recipes); library(caret); library(pROC); library(purrr); library(stringr)
  library(RWeka);   library(kernlab); library(randomForest); library(glmnet); library(kknn)
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("mutate", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("lag",    "dplyr", quiet = TRUE)
  }
})

set.seed(123)

# ---------------- switches ----------------
fast_mode <- TRUE
cv_k      <- if (fast_mode) 3 else 10
cv_R      <- 1
inner_k   <- if (fast_mode) 3 else 5
tune_len  <- if (fast_mode) 3 else 5
print_fold_diag <- TRUE

algos_to_run <- c("LR","RF","SVM","k-NN","C4.5")

# ---------------- I/O ----------------
file_path   <- "biomarkers_acuteCOVID_meta.xlsx"
sheet_name  <- "meta"
sev_adm_var <- "severity_admission"
worst_var   <- "worst_severity"

keep_cols <- c(
  "Group", sev_adm_var, worst_var, "Hospital_ID",
  "Age","Gender","SpO2_admission","albumin","CRP","D_Dimer",
  "monocytes_perc","neutrophils_perc","lymphocytes_perc",
  "monocyte_abs_number","neutrophil_abs_number","lymphocyte_abs_number"
)

# ---------------- load + label engineering ----------------
df <- readxl::read_excel(file_path, sheet = sheet_name) %>%
  dplyr::select(dplyr::any_of(keep_cols)) %>%
  dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
  dplyr::mutate(
    NLR_abs     = ifelse(is.na(neutrophil_abs_number) | is.na(lymphocyte_abs_number), NA_real_,
                         neutrophil_abs_number / ifelse(lymphocyte_abs_number == 0, NA, lymphocyte_abs_number)),
    CRP_albumin = ifelse(is.na(CRP) | is.na(albumin), NA_real_, CRP / ifelse(albumin == 0, NA, albumin)),
    Diag_Bin = factor(dplyr::case_when(
      Group %in% c("CTRL_noCOVID","COVID") ~ as.character(Group),
      TRUE ~ NA_character_
    ), levels = c("CTRL_noCOVID","COVID")),
    Sev3_Admission = factor(dplyr::case_when(
      Group == "CTRL_noCOVID" ~ "CTRL_noCOVID",
      Group == "COVID" & !is.na(.data[[sev_adm_var]]) ~ as.character(.data[[sev_adm_var]]),
      TRUE ~ NA_character_
    ), levels = c("CTRL_noCOVID","Mild","Severe")),
    Hosp_Bin = factor(dplyr::case_when(
      Hospital_ID %in% c("Yes","No") ~ as.character(Hospital_ID),
      TRUE ~ NA_character_
    ), levels = c("No","Yes")),
    Worst3 = factor(dplyr::case_when(
      Group == "CTRL_noCOVID" ~ "CTRL_noCOVID",
      Group == "COVID" & !is.na(.data[[worst_var]]) ~ as.character(.data[[worst_var]]),
      TRUE ~ NA_character_
    ), levels = c("CTRL_noCOVID","Mild","Severe"))
  ) %>%
  droplevels()

# prefer absolute counts over percentages
df <- df %>% dplyr::select(-dplyr::any_of(c("monocytes_perc","neutrophils_perc","lymphocytes_perc")))

# stable row ids
df$.rid <- seq_len(nrow(df))

# feature set
feature_vars <- setdiff(colnames(df), c("Group","Diag_Bin","Sev3_Admission","Hosp_Bin","Worst3",".rid", sev_adm_var, worst_var, "Hospital_ID"))

# ---------------- shared helpers ----------------
`%||%` <- function(a,b) if (!is.null(a)) a else b

choose_k_for_task <- function(y) {
  tab <- table(y); max_k <- max(2, min(tab))
  kk <- min(cv_k, max_k)
  if (kk < cv_k) message(sprintf("Folds reduced to %d for class balance", kk))
  kk
}

build_cv_splits <- function(y, R = cv_R, seed = 999) {
  set.seed(seed)
  splits <- list()
  for (r in seq_len(R)) {
    kk <- choose_k_for_task(y)
    fold_list <- caret::createFolds(y, k = kk, returnTrain = FALSE)
    for (j in seq_len(kk)) {
      test_idx  <- fold_list[[j]]
      train_idx <- setdiff(seq_along(y), test_idx)
      splits[[paste0("r", r, "_f", j)]] <- list(train_idx = train_idx, test_idx = test_idx)
    }
  }
  splits
}

make_recipe <- function(dat, yvar) {
  recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat) %>%
    update_role(.rid, new_role = "id") %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_novel(all_nominal_predictors(), new_level = "NOVEL") %>%
    step_other(all_nominal_predictors(), threshold = 0.01, other = "OTHER") %>%
    step_zv(all_predictors()) %>%
    step_YeoJohnson(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE)
}

# metrics
calc_conf_metrics <- function(obs, pred, lev) {
  cm <- table(factor(obs, levels = lev), factor(pred, levels = lev))
  TP <- diag(cm); FN <- rowSums(cm) - TP; FP <- colSums(cm) - TP
  TN <- sum(cm) - TP - FN - FP
  eps <- 1e-9
  sens <- TP / pmax(TP + FN, eps)
  spec <- TN / pmax(TN + FP, eps)
  prec <- TP / pmax(TP + FP, eps)
  f1   <- 2*TP / pmax(2*TP + FP + FN, eps)
  acc  <- sum(TP) / pmax(sum(cm), eps)
  mcc  <- (TP*TN - FP*FN) / sqrt(pmax((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN), eps))
  c(
    MCC = mean(mcc, na.rm = TRUE),
    F1 = mean(f1, na.rm = TRUE),
    Accuracy = acc,
    Precision = mean(prec, na.rm = TRUE),
    Sensitivity = mean(sens, na.rm = TRUE),
    Specificity = mean(spec, na.rm = TRUE)
  )
}

macro_auc_ovr <- function(obs, prob_mat, lev) {
  if (is.null(prob_mat)) return(NA_real_)
  keep <- intersect(colnames(prob_mat), lev)
  if (length(keep) < 1) return(NA_real_)
  aucs <- sapply(keep, function(cls) {
    y <- as.integer(obs == cls)
    p <- as.numeric(prob_mat[[cls]])
    if (length(unique(y)) < 2 || all(is.na(p))) return(NA_real_)
    r <- try(pROC::roc(y, p, quiet = TRUE), silent = TRUE)
    if (inherits(r, "try-error")) NA_real_ else as.numeric(pROC::auc(r))
  })
  if (all(is.na(aucs))) NA_real_ else mean(aucs, na.rm = TRUE)
}

compute_metrics <- function(obs, pred, prob_mat, lev) {
  base <- calc_conf_metrics(obs, pred, lev)
  auc  <- macro_auc_ovr(obs, prob_mat, lev)
  c(base, AUC = auc)
}

metricSummaryInner <- function(data, lev = NULL, model = NULL) {
  lev <- lev %||% levels(data$obs)
  auc <- macro_auc_ovr(data$obs, data[, intersect(colnames(data), lev), drop=FALSE], lev)
  acc <- mean(data$obs == data$pred, na.rm = TRUE)
  c(AUC = auc, Accuracy = acc)
}

safe_train <- function(expr) suppressWarnings(tryCatch(expr, error = function(e) NULL))

fit_fixed <- function(dat_train, yvar, spec, bestTune) {
  rec <- make_recipe(dat_train, yvar)
  tr_none <- caret::trainControl(method = "none", classProbs = TRUE,
                                 summaryFunction = metricSummaryInner)
  safe_train({
    if (!is.null(bestTune) && ncol(bestTune) > 0) {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = tr_none, metric = "AUC", tuneGrid = bestTune)
    } else if (!is.null(spec$tuneLength)) {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = caret::trainControl(method="cv", number=inner_k,
                                                   classProbs=TRUE, savePredictions="final",
                                                   summaryFunction = metricSummaryInner),
                   tuneLength = spec$tuneLength, metric = "AUC")
    } else {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = tr_none, metric = "AUC")
    }
  })
}

# ---------------- model specs ----------------
model_specs_full <- list(
  "C4.5" = list(method = "J48",       tuneLength = NULL),
  "k-NN" = list(method = "kknn",      tuneLength = tune_len),
  "SVM"  = list(method = "svmRadial", tuneLength = tune_len),
  "RF"   = list(method = "rf",        tuneLength = tune_len),
  "LR"   = list(method = "glmnet",    tuneLength = tune_len)
)
model_specs <- model_specs_full[names(model_specs_full) %in% algos_to_run]

# ---------------- thresholding helpers for binary ----------------
mcc_binary_vec <- function(y_true, y_pred, pos, neg) {
  y_true <- factor(y_true, levels = c(neg, pos))
  y_pred <- factor(y_pred, levels = c(neg, pos))
  cm <- table(y_true, y_pred)
  TP <- cm[pos,pos]; TN <- cm[neg,neg]
  FP <- cm[neg,pos]; FN <- cm[pos,neg]
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(0)
  as.numeric(num/den)
}

best_threshold_mcc <- function(y_true, p_pos, pos, neg, grid = seq(0.05, 0.95, by = 0.01)) {
  scores <- sapply(grid, function(t) {
    pred <- ifelse(p_pos >= t, pos, neg)
    mcc_binary_vec(y_true, pred, pos, neg)
  })
  t_star <- grid[which.max(scores)]
  list(t = t_star, mcc = max(scores))
}

# ---------------- generic inner training ----------------
train_with_inner <- function(dat_train, yvar, spec) {
  rec <- make_recipe(dat_train, yvar)
  tr_ctrl <- caret::trainControl(
    method = "cv", number = inner_k,
    classProbs = TRUE, savePredictions = "final",
    summaryFunction = metricSummaryInner
  )
  fit <- safe_train({
    if (!is.null(spec$tuneLength)) {
      caret::train(rec, data = dat_train, method = spec$method,
                   tuneLength = spec$tuneLength, trControl = tr_ctrl, metric = "AUC")
    } else {
      caret::train(rec, data = dat_train, method = spec$method,
                   trControl = tr_ctrl, metric = "AUC")
    }
  })
  if (is.null(fit)) return(NULL)
  pred <- fit$pred
  if (!is.null(fit$bestTune) && ncol(fit$bestTune) > 0) {
    pred <- dplyr::semi_join(pred, fit$bestTune, by = names(fit$bestTune))
  }
  list(fit = fit, oof = pred)
}

# ---------------- fold evaluators ----------------
fold_eval_binary <- function(dat_tr, dat_va, spec, yvar, pos_label, algo_name, fold_tag) {
  lev <- levels(dat_tr[[yvar]])
  neg_label <- setdiff(lev, pos_label)[1]
  inner <- train_with_inner(dat_tr, yvar, spec)
  if (is.null(inner)) return(NULL)
  if (!all(lev %in% colnames(inner$oof))) return(NULL)
  
  t_res <- best_threshold_mcc(inner$oof$obs, as.numeric(inner$oof[, pos_label]),
                              pos = pos_label, neg = neg_label)
  bestTune <- inner$fit$bestTune %||% NULL
  fit <- fit_fixed(dat_tr, yvar, spec, bestTune)
  if (is.null(fit)) return(NULL)
  
  pr <- try(predict(fit, newdata = dat_va, type = "prob"), silent = TRUE)
  if (inherits(pr, "try-error") || is.null(pr) || !all(lev %in% colnames(pr))) return(NULL)
  
  pred <- factor(ifelse(pr[[pos_label]] >= t_res$t, pos_label, neg_label), levels = lev)
  m    <- compute_metrics(dat_va[[yvar]], pred, pr[, lev, drop = FALSE], lev)
  
  if (print_fold_diag) {
    cat(sprintf("[Fold %s] %s | %s | tuned.t=%.3f | prob.cols={%s}\n",
                fold_tag, yvar, algo_name, t_res$t, paste(intersect(colnames(pr), lev), collapse=",")))
  }
  m
}

fold_eval_multiclass <- function(dat_tr, dat_va, spec, yvar, algo_name, fold_tag) {
  lev <- levels(dat_tr[[yvar]])
  inner <- train_with_inner(dat_tr, yvar, spec)
  if (is.null(inner)) return(NULL)
  bestTune <- inner$fit$bestTune %||% NULL
  fit <- fit_fixed(dat_tr, yvar, spec, bestTune)
  if (is.null(fit)) return(NULL)
  
  pr <- try(predict(fit, newdata = dat_va, type = "prob"), silent = TRUE)
  if (inherits(pr, "try-error") || is.null(pr) || !all(lev %in% colnames(pr))) return(NULL)
  
  pr_use <- pr[, lev, drop = FALSE]
  pred <- factor(colnames(pr_use)[max.col(as.matrix(pr_use), ties.method = "first")], levels = lev)
  m    <- compute_metrics(dat_va[[yvar]], pred, pr_use, lev)
  
  if (print_fold_diag) {
    cat(sprintf("[Fold %s] %s | %s | prob.cols={%s}\n",
                fold_tag, yvar, algo_name, paste(intersect(colnames(pr), lev), collapse=",")))
  }
  m
}

# ---------------- generic runner ----------------
run_task <- function(task_name, yvar, type, pos_label = NULL, seed = 100) {
  dat <- df %>%
    dplyr::filter(!is.na(.data[[yvar]])) %>%
    dplyr::select(dplyr::all_of(c(yvar, feature_vars, ".rid")))
  lev <- levels(dat[[yvar]])
  splits <- build_cv_splits(dat[[yvar]], R = cv_R, seed = seed)
  
  out_rows <- list()
  for (algo in names(model_specs)) {
    spec <- model_specs[[algo]]
    metrics <- list()
    for (nm in names(splits)) {
      tr_idx <- splits[[nm]]$train_idx
      va_idx <- splits[[nm]]$test_idx
      dat_tr <- dat[tr_idx, , drop = FALSE]
      dat_va <- dat[va_idx, , drop = FALSE]
      
      m <- if (type == "binary") {
        fold_eval_binary(dat_tr, dat_va, spec, yvar, pos_label, algo, nm)
      } else {
        fold_eval_multiclass(dat_tr, dat_va, spec, yvar, algo, nm)
      }
      if (!is.null(m)) metrics[[nm]] <- m
    }
    if (!length(metrics)) next
    M <- do.call(rbind, metrics)
    row <- tibble::tibble(Task = task_name, Algorithm = algo)
    for (mm in colnames(M)) {
      row[[paste0(mm,"_Mean")]] <- round(mean(M[, mm], na.rm = TRUE), 4)
      row[[paste0(mm,"_SD")]]   <- round(sd(M[, mm],   na.rm = TRUE), 4)
    }
    out_rows <- c(out_rows, list(row))
  }
  if (length(out_rows)) dplyr::bind_rows(out_rows) else tibble::tibble()
}

# ---------------- pretty printer ----------------
print_task_table <- function(tbl, task_label) {
  cat(sprintf("\n=== %s METRICS (mean ± SD) ===\n", toupper(task_label)))
  if (!nrow(tbl)) {
    cat("(no results)\n")
    return(invisible(NULL))
  }
  pretty <- tbl %>%
    tidyr::pivot_longer(cols = -c(Task,Algorithm),
                        names_to = c("Metric",".value"),
                        names_pattern = "(.*)_(Mean|SD)") %>%
    dplyr::mutate(Value = sprintf("%.4f ± %.4f", Mean, SD)) %>%
    dplyr::select(Task, Algorithm, Metric, Value) %>%
    tidyr::pivot_wider(names_from = Algorithm, values_from = Value) %>%
    dplyr::arrange(factor(Metric,
                          levels=c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")))
  print(pretty, n = nrow(pretty))
}

# ============================================================
# Run all tasks with one-liners
# ============================================================

# Task 1 Diagnosis, pos = COVID
diag_tbl <- run_task(task_name = "Diagnosis", yvar = "Diag_Bin", type = "binary", pos_label = "COVID", seed = 777)
print_task_table(diag_tbl, "Diagnosis")

# Task 2 Severity at admission, 3-class
sev_tbl  <- run_task(task_name = "SeverityAdmission3", yvar = "Sev3_Admission", type = "multiclass", seed = 999)
print_task_table(sev_tbl, "Severity-at-admission")

# Task 3 Hospitalization, pos = Yes
hosp_tbl <- run_task(task_name = "Hospitalization", yvar = "Hosp_Bin", type = "binary", pos_label = "Yes", seed = 444)
print_task_table(hosp_tbl, "Hospitalization")

# Task 4 Worst severity, 3-class
worst_tbl <- run_task(task_name = "WorstSeverity3", yvar = "Worst3", type = "multiclass", seed = 333)
print_task_table(worst_tbl, "Worst severity")
