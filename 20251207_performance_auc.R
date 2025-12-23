suppressPackageStartupMessages({
  library(tidyverse)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(kknn)
  library(randomForest)
  library(kernlab)
  library(RWeka)
})

options(warn = -1)

# ======================================================================
# CONFIG
# ======================================================================

max_n_per_dataset <- 5000L
use_class_weights <- TRUE

t_low_fixed  <- 0.33
t_high_fixed <- 0.66

cfg <- list(
  pos_label      = "Yes",
  neg_label      = "No",
  seed_cv        = 123,
  t_low_default  = t_low_fixed,
  t_high_default = t_high_fixed
)

bands <- expand.grid(
  t_low_band  = c(0.20, 0.25, 0.30, 0.33, 0.35),
  t_high_band = c(0.65, 0.66, 0.70, 0.75, 0.80)
) %>%
  dplyr::filter(t_low_band < t_high_band)

message(
  ">>> FINAL BENCHMARK (",
  max_n_per_dataset, " per dataset, ",
  ifelse(use_class_weights, "class weights, ", "no class weights, "),
  "engineered features, repeated CV, Platt scaling fitted on development out of fold predictions and applied to external, ",
  "thresholds selected on development only, cascade with fixed deferral thresholds ",
  t_low_fixed, " and ", t_high_fixed,
  ", models optimised for area under the receiver operating characteristic curve) <<<"
)

# ======================================================================
# 1. LOAD DATA
# ======================================================================

dev_path <- "mimic_final_22_features_60pct.csv"
ext_path <- "mcmed_final_22_features_60pct.csv"

if (!file.exists(dev_path) || !file.exists(ext_path)) {
  stop("Input CSV files were not found: ", dev_path, " and or ", ext_path)
}

dev_raw <- readr::read_csv(dev_path, show_col_types = FALSE)
ext_raw <- readr::read_csv(ext_path, show_col_types = FALSE)

prep_df <- function(df) {
  if ("outcome" %in% names(df)) {
    df %>%
      dplyr::mutate(outcome = factor(outcome, levels = c("No", "Yes"))) %>%
      dplyr::select(outcome, dplyr::everything())
  } else if ("group" %in% names(df)) {
    df %>%
      dplyr::filter(.data$group %in% c("CTRL_noCOVID", "COVID")) %>%
      dplyr::mutate(
        outcome = dplyr::case_when(
          .data$group == "COVID" ~ "Yes",
          .data$group == "CTRL_noCOVID" ~ "No",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(.data$outcome)) %>%
      dplyr::mutate(outcome = factor(.data$outcome, levels = c("No", "Yes"))) %>%
      dplyr::select(outcome, dplyr::everything())
  } else {
    stop("A required column named outcome or group was not found in the data.")
  }
}

set.seed(cfg$seed_cv)

df_dev_full <- prep_df(dev_raw)
df_ext_full <- prep_df(ext_raw)

n_dev  <- min(max_n_per_dataset, nrow(df_dev_full))
n_ext  <- min(max_n_per_dataset, nrow(df_ext_full))

df_train <- df_dev_full %>% dplyr::slice_sample(n = n_dev)
df_test  <- df_ext_full %>% dplyr::slice_sample(n = n_ext)

y_test <- df_test$outcome

# ======================================================================
# 2. ALIGN COLUMNS AND IMPUTE NUMERICS (TRAIN MEDIANS APPLIED TO TEST)
# ======================================================================

align_test_to_train <- function(train, test) {
  missing <- setdiff(names(train), names(test))
  if (length(missing) > 0L) {
    for (c in missing) {
      x <- train[[c]]
      if (is.factor(x)) {
        test[[c]] <- factor(NA_character_, levels = levels(x))
      } else if (is.character(x)) {
        test[[c]] <- NA_character_
      } else if (inherits(x, "Date")) {
        test[[c]] <- as.Date(NA)
      } else if (inherits(x, "POSIXct")) {
        test[[c]] <- as.POSIXct(NA)
      } else if (is.logical(x)) {
        test[[c]] <- NA
      } else {
        test[[c]] <- NA_real_
      }
    }
  }
  
  extra <- setdiff(names(test), names(train))
  if (length(extra) > 0L) {
    test <- test %>% dplyr::select(-dplyr::all_of(extra))
  }
  
  test %>% dplyr::select(dplyr::all_of(names(train)))
}

median_impute_numeric <- function(train, test) {
  med_tbl <- train %>%
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        ~ stats::median(., na.rm = TRUE)
      )
    )
  med_list <- as.list(med_tbl)
  
  train2 <- train %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ ifelse(is.na(.), med_list[[dplyr::cur_column()]], .)
      )
    )
  
  test2 <- test %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ ifelse(is.na(.), med_list[[dplyr::cur_column()]], .)
      )
    )
  
  list(train = train2, test = test2)
}

df_test <- align_test_to_train(df_train, df_test)
imp1 <- median_impute_numeric(df_train, df_test)
df_train <- imp1$train
df_test  <- imp1$test

# ======================================================================
# 3. FEATURE ENGINEERING, THEN IMPUTE AGAIN (FOR NEW FEATURES)
# ======================================================================

add_features <- function(df) {
  df2 <- df
  
  safe_ratio <- function(a, b) {
    out <- a / b
    out[!is.finite(out)] <- NA_real_
    out
  }
  
  if (all(c("neutrophils_raw", "lymphocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(NLR = safe_ratio(.data$neutrophils_raw, .data$lymphocytes_raw))
  }
  
  if (all(c("monocytes_raw", "lymphocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(MLR = safe_ratio(.data$monocytes_raw, .data$lymphocytes_raw))
  }
  
  if (all(c("neutrophils_raw", "monocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(NMR = safe_ratio(.data$neutrophils_raw, .data$monocytes_raw))
  }
  
  if (all(c("neutrophils_raw", "monocytes_raw", "lymphocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(SIRI = safe_ratio(.data$neutrophils_raw * .data$monocytes_raw, .data$lymphocytes_raw))
  }
  
  if (all(c("crp", "lymphocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(CLR = safe_ratio(.data$crp, .data$lymphocytes_raw))
  }
  
  if (all(c("crp", "albumin") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(CAR = safe_ratio(.data$crp, .data$albumin))
  }
  
  if (all(c("albumin", "lymphocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(PNI = .data$albumin + 5 * .data$lymphocytes_raw)
  }
  
  if (all(c("d_dimer", "lymphocytes_raw") %in% names(df2))) {
    df2 <- df2 %>% dplyr::mutate(DLR = safe_ratio(.data$d_dimer, .data$lymphocytes_raw))
  }
  
  df2
}

df_train <- add_features(df_train)
df_test  <- add_features(df_test)

df_test <- align_test_to_train(df_train, df_test)
imp2 <- median_impute_numeric(df_train, df_test)
df_train <- imp2$train
df_test  <- imp2$test

all_cols <- setdiff(names(df_train), "outcome")

triage_candidates <- c(
  "age",
  "sex",
  "spo2",
  "heart_rate",
  "resp_rate",
  "temp",
  "sbp",
  "dbp",
  "map"
)

feat_triage <- intersect(triage_candidates, all_cols)
if (length(feat_triage) == 0L) {
  feat_triage <- grep("age|sex|spo2|heart|resp|temp|sbp|dbp|map", all_cols, value = TRUE)
}

feat_full <- all_cols

# ======================================================================
# 4. METRICS
# ======================================================================

confusion_counts <- function(obs, pred, pos, neg) {
  tab <- table(
    factor(obs,  levels = c(neg, pos)),
    factor(pred, levels = c(neg, pos))
  )
  
  TN <- tab[neg, neg]
  FP <- tab[neg, pos]
  FN <- tab[pos, neg]
  TP <- tab[pos, pos]
  
  c(
    TN = as.numeric(TN),
    FP = as.numeric(FP),
    FN = as.numeric(FN),
    TP = as.numeric(TP)
  )
}

compute_metrics <- function(obs, pred, prob, pos, neg) {
  cc <- confusion_counts(obs, pred, pos, neg)
  TN <- cc["TN"]; FP <- cc["FP"]; FN <- cc["FN"]; TP <- cc["TP"]
  
  acc <- (TP + TN) / (TP + TN + FP + FN)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (denom > 0) ((TP * TN) - (FP * FN)) / denom else NA_real_
  
  f1 <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) {
    2 * prec * sens / (prec + sens)
  } else {
    NA_real_
  }
  
  auc_val <- NA_real_
  if (!all(is.na(prob))) {
    roc_obj <- try(
      pROC::roc(
        response  = obs,
        predictor = as.numeric(prob),
        levels    = c(neg, pos),
        direction = "<"
      ),
      silent = TRUE
    )
    if (!inherits(roc_obj, "try-error")) {
      auc_val <- as.numeric(pROC::auc(roc_obj))
    }
  }
  
  list(
    AUC         = auc_val,
    Accuracy    = acc,
    MCC         = mcc,
    F1          = f1,
    Sensitivity = sens,
    Specificity = spec,
    Precision   = prec,
    TN          = TN,
    FP          = FP,
    FN          = FN,
    TP          = TP
  )
}

best_thresh_roc <- function(obs, prob, pos = "Yes", neg = "No") {
  mask <- !is.na(prob)
  if (sum(mask) < 2L) return(0.5)
  
  roc_obj <- try(
    pROC::roc(
      response  = obs[mask],
      predictor = as.numeric(prob[mask]),
      levels    = c(neg, pos),
      direction = "<"
    ),
    silent = TRUE
  )
  if (inherits(roc_obj, "try-error")) return(0.5)
  
  thr <- try(
    pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold"),
    silent = TRUE
  )
  if (inherits(thr, "try-error") || length(thr) == 0L || is.na(thr[1])) return(0.5)
  
  as.numeric(thr[1])
}

fit_calibrator <- function(obs, prob, pos_label) {
  y  <- ifelse(obs == pos_label, 1, 0)
  p  <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  df <- data.frame(y = y, p = p)
  
  suppressWarnings(
    stats::glm(y ~ qlogis(p), data = df, family = stats::binomial(link = "logit"))
  )
}

apply_calibration <- function(prob, calib_model) {
  if (is.null(calib_model)) return(prob)
  p <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  lp <- predict(calib_model, newdata = data.frame(p = p), type = "link")
  1 / (1 + exp(-lp))
}

# ======================================================================
# 4.1 caret summary function to optimise area under the receiver operating characteristic curve
# ======================================================================

aucSummary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) lev <- levels(data$obs)
  if (length(lev) != 2L) stop("aucSummary requires a binary outcome.")
  
  pos <- lev[2]
  prob_pos <- data[[pos]]
  
  roc_obj <- try(
    pROC::roc(
      response  = data$obs,
      predictor = prob_pos,
      levels    = lev,
      direction = "<"
    ),
    silent = TRUE
  )
  
  roc_auc <- if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj))
  c(ROC = roc_auc)
}

ctrl <- caret::trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = aucSummary,
  savePredictions = "final",
  verboseIter     = FALSE,
  allowParallel   = TRUE
)

models <- c("LR", "RF", "SVM", "k-NN", "C4.5")

# ======================================================================
# 5. TRAIN MODEL, FIT PLATT SCALING ON OUT OF FOLD PREDICTIONS, SELECT THRESHOLD ON DEVELOPMENT ONLY
# ======================================================================

extract_oof_best <- function(fit) {
  pred <- fit$pred
  if (is.null(pred) || nrow(pred) == 0L) stop("Out of fold predictions were not available from caret.")
  
  bt <- fit$bestTune
  for (nm in names(bt)) {
    pred <- pred %>% dplyr::filter(.data[[nm]] == bt[[nm]])
  }
  
  pred
}

run_algo <- function(algo, feat_names, df_train, df_test = NULL) {
  dat_train <- df_train[, c("outcome", feat_names), drop = FALSE]
  dat_train$outcome <- droplevels(dat_train$outcome)
  
  if (!is.null(df_test)) {
    dat_test <- df_test[, c("outcome", feat_names), drop = FALSE]
    dat_test$outcome <- droplevels(dat_test$outcome)
  } else {
    dat_test <- NULL
  }
  
  y <- dat_train$outcome
  
  if (use_class_weights) {
    tab   <- table(y)
    n_pos <- as.numeric(tab[cfg$pos_label])
    n_neg <- as.numeric(tab[cfg$neg_label])
    w_pos <- if (!is.na(n_pos) && n_pos > 0) n_neg / n_pos else 1
    w <- ifelse(y == cfg$pos_label, w_pos, 1)
  } else {
    w <- rep(1, length(y))
  }
  
  method <- switch(
    algo,
    "LR"   = "glmnet",
    "RF"   = "rf",
    "SVM"  = "svmRadial",
    "k-NN" = "kknn",
    "C4.5" = "J48",
    stop("Unknown algorithm: ", algo)
  )
  
  fit <- switch(
    algo,
    "LR" = caret::train(
      outcome ~ .,
      data       = dat_train,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl,
      tuneLength = 20,
      weights    = w,
      family     = "binomial"
    ),
    "RF" = caret::train(
      outcome ~ .,
      data       = dat_train,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl,
      tuneLength = 5,
      weights    = w,
      ntree      = 500
    ),
    "SVM" = caret::train(
      outcome ~ .,
      data       = dat_train,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl,
      tuneLength = 10,
      weights    = w
    ),
    "k-NN" = caret::train(
      outcome ~ .,
      data       = dat_train,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl,
      tuneLength = 10,
      weights    = w
    ),
    "C4.5" = caret::train(
      outcome ~ .,
      data       = dat_train,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl,
      tuneLength = 5,
      weights    = w
    )
  )
  
  oof <- extract_oof_best(fit)
  
  prob_oof_raw <- oof[[cfg$pos_label]]
  obs_oof      <- oof$obs
  
  calib <- fit_calibrator(obs_oof, prob_oof_raw, cfg$pos_label)
  prob_oof_cal <- apply_calibration(prob_oof_raw, calib)
  
  thr_train <- best_thresh_roc(
    obs  = obs_oof,
    prob = prob_oof_cal,
    pos  = cfg$pos_label,
    neg  = cfg$neg_label
  )
  
  if (is.null(dat_test)) {
    return(list(
      fit       = fit,
      calib     = calib,
      thr_train = thr_train
    ))
  }
  
  prob_test_raw <- stats::predict(fit, newdata = dat_test, type = "prob")[, cfg$pos_label]
  prob_test_cal <- apply_calibration(prob_test_raw, calib)
  
  pred_test <- ifelse(prob_test_cal >= thr_train, cfg$pos_label, cfg$neg_label)
  pred_test <- factor(pred_test, levels = c(cfg$neg_label, cfg$pos_label))
  
  list(
    fit           = fit,
    calib         = calib,
    thr_train     = thr_train,
    prob_test_raw = prob_test_raw,
    prob_test_cal = prob_test_cal,
    pred_test     = pred_test
  )
}

# ======================================================================
# 6. FIT TRIAGE AND FULL MODELS ONCE
# ======================================================================

triage_res   <- vector("list", length(models)); names(triage_res)   <- models
complete_res <- vector("list", length(models)); names(complete_res) <- models

for (algo in models) {
  message(sprintf("Fitting %s (triage)...", algo))
  triage_res[[algo]] <- run_algo(algo, feat_triage, df_train, df_test)
  
  message(sprintf("Fitting %s (full)...", algo))
  complete_res[[algo]] <- run_algo(algo, feat_full, df_train, df_test)
}

# ======================================================================
# 7. CASCADE FOR A GIVEN DEFERRAL BAND
# ======================================================================

run_cascade_for_band <- function(t_low, t_high) {
  out_list <- vector("list", length(models) * 12L)
  idx <- 1L
  
  for (algo in models) {
    prob_T <- triage_res[[algo]]$prob_test_cal
    pred_T <- triage_res[[algo]]$pred_test
    
    prob_C <- complete_res[[algo]]$prob_test_cal
    pred_C <- complete_res[[algo]]$pred_test
    thr_C  <- complete_res[[algo]]$thr_train
    
    decision <- ifelse(
      prob_T < t_low,  cfg$neg_label,
      ifelse(prob_T > t_high, cfg$pos_label, "Defer")
    )
    
    cascade_pred <- character(length(decision))
    mask_nodefer <- decision != "Defer"
    cascade_pred[mask_nodefer] <- decision[mask_nodefer]
    
    mask_defer <- !mask_nodefer
    cascade_pred[mask_defer] <- ifelse(
      prob_C[mask_defer] >= thr_C,
      cfg$pos_label,
      cfg$neg_label
    )
    
    cascade_pred <- factor(cascade_pred, levels = c(cfg$neg_label, cfg$pos_label))
    defer_rate   <- mean(decision == "Defer")
    
    prob_S <- ifelse(decision == "Defer", prob_C, prob_T)
    
    mets_T <- compute_metrics(y_test, pred_T, prob_T, pos = cfg$pos_label, neg = cfg$neg_label)
    mets_C <- compute_metrics(y_test, pred_C, prob_C, pos = cfg$pos_label, neg = cfg$neg_label)
    mets_S <- compute_metrics(y_test, cascade_pred, prob_S, pos = cfg$pos_label, neg = cfg$neg_label)
    
    metrics_to_keep <- c(
      "AUC", "MCC", "F1", "Accuracy",
      "Sensitivity", "Specificity", "Precision"
    )
    
    for (mname in metrics_to_keep) {
      out_list[[idx]] <- tibble::tibble(
        Algorithm   = algo,
        t_low_band  = t_low,
        t_high_band = t_high,
        Metric_raw  = mname,
        Mean_T      = as.numeric(mets_T[[mname]]),
        Mean_C      = as.numeric(mets_C[[mname]]),
        Mean_S      = as.numeric(mets_S[[mname]])
      )
      idx <- idx + 1L
    }
    
    cc_S <- c(TN = mets_S$TN, FP = mets_S$FP, FN = mets_S$FN, TP = mets_S$TP)
    
    out_list[[idx]] <- tibble::tibble(
      Algorithm   = algo,
      t_low_band  = t_low,
      t_high_band = t_high,
      Metric_raw  = "DeferRate",
      Mean_T      = NA_real_,
      Mean_C      = NA_real_,
      Mean_S      = defer_rate
    )
    idx <- idx + 1L
    
    for (nm in names(cc_S)) {
      out_list[[idx]] <- tibble::tibble(
        Algorithm   = algo,
        t_low_band  = t_low,
        t_high_band = t_high,
        Metric_raw  = nm,
        Mean_T      = NA_real_,
        Mean_C      = NA_real_,
        Mean_S      = as.numeric(cc_S[[nm]])
      )
      idx <- idx + 1L
    }
    
    out_list[[idx]] <- tibble::tibble(
      Algorithm   = algo,
      t_low_band  = t_low,
      t_high_band = t_high,
      Metric_raw  = "t_low",
      Mean_T      = NA_real_,
      Mean_C      = NA_real_,
      Mean_S      = t_low
    )
    idx <- idx + 1L
    
    out_list[[idx]] <- tibble::tibble(
      Algorithm   = algo,
      t_low_band  = t_low,
      t_high_band = t_high,
      Metric_raw  = "t_high",
      Mean_T      = NA_real_,
      Mean_C      = NA_real_,
      Mean_S      = t_high
    )
    idx <- idx + 1L
  }
  
  dplyr::bind_rows(out_list)
}

options(warn = 0)

# ======================================================================
# 8. RUN CASCADE FOR ALL BANDS AND SUMMARISE
# ======================================================================

res_grid_list <- vector("list", nrow(bands))

for (i in seq_len(nrow(bands))) {
  message(
    ">>> Cascade band t_low = ",
    bands$t_low_band[i],
    ", t_high = ",
    bands$t_high_band[i],
    " <<<"
  )
  res_grid_list[[i]] <- run_cascade_for_band(bands$t_low_band[i], bands$t_high_band[i])
}

res_grid <- dplyr::bind_rows(res_grid_list)

res_all_num <- res_grid %>%
  dplyr::filter(
    Metric_raw %in% c(
      "AUC",
      "MCC",
      "F1",
      "Accuracy",
      "Precision",
      "Sensitivity",
      "Specificity",
      "DeferRate"
    )
  )

cat(
  "\n--- Grid results over deferral thresholds (",
  max_n_per_dataset, " per dataset, ",
  ifelse(use_class_weights, "class weights, ", "no class weights, "),
  "engineered features, repeated CV, Platt scaling and thresholds fitted on development only) ---\n",
  sep = ""
)

res_auc_best <- res_all_num %>%
  dplyr::filter(Metric_raw == "AUC") %>%
  dplyr::group_by(Algorithm) %>%
  dplyr::arrange(dplyr::desc(Mean_S), .by_group = TRUE) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::ungroup() %>%
  dplyr::select(Algorithm, t_low_band, t_high_band, Mean_S_AUC = Mean_S)

print(res_auc_best, n = Inf)

baseline_all_num <- res_all_num %>%
  dplyr::filter(
    abs(t_low_band  - t_low_fixed)  < 1e-9,
    abs(t_high_band - t_high_fixed) < 1e-9
  )

res_print <- baseline_all_num %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("Mean"),
      ~ ifelse(is.na(.), NA_character_, sprintf("%.4f", .))
    )
  )

cat(
  "\n--- Final results for baseline band (",
  max_n_per_dataset, " per dataset, ",
  ifelse(use_class_weights, "class weights, ", "no class weights, "),
  "deferral thresholds ",
  t_low_fixed, " and ", t_high_fixed,
  ", thresholds fitted on development only, with TP TN FP FN) ---\n",
  sep = ""
)

print(res_print, n = Inf)

band_metrics <- baseline_all_num %>%
  dplyr::filter(
    Metric_raw %in% c(
      "AUC",
      "MCC",
      "F1",
      "Accuracy",
      "Precision",
      "Sensitivity",
      "Specificity"
    )
  )

cat(
  "\n--- Detailed metrics for chosen band t_low = ",
  t_low_fixed,
  ", t_high = ",
  t_high_fixed,
  " ---\n",
  sep = ""
)

print(band_metrics, n = Inf)

auc_table <- band_metrics %>%
  dplyr::filter(Metric_raw == "AUC") %>%
  dplyr::select(Algorithm, Mean_T, Mean_C, Mean_S)

cat("\n--- AUCs per algorithm and model for chosen band ---\n")
print(auc_table, n = Inf)
