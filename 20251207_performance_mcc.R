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

max_n_per_dataset <- 10000L          # cap per dataset (set to Inf to use all)
use_class_weights <- TRUE            # use class weights

t_low_fixed  <- 0.33                 # chosen deferral band
t_high_fixed <- 0.66

cfg <- list(
  pos_label      = "Yes",
  neg_label      = "No",
  seed_cv        = 123,
  t_low_default  = t_low_fixed,
  t_high_default = t_high_fixed
)

# deferral band grid (INCLUDES 0.33 / 0.66)
bands <- expand.grid(
  t_low_band  = c(0.20, 0.25, 0.30, 0.33, 0.35),
  t_high_band = c(0.65, 0.66, 0.70, 0.75, 0.80)
) %>%
  dplyr::as_tibble() %>%
  dplyr::filter(t_low_band < t_high_band)

message(
  ">>> FINAL BENCHMARK (",
  max_n_per_dataset, " per dataset, ",
  ifelse(use_class_weights, "class weights, ", "no class weights, "),
  "engineered features, repeated CV, cascade with fixed deferral thresholds ",
  t_low_fixed, " and ", t_high_fixed,
  ", triage/full MCC thresholds, cascade threshold MCC-optimised) <<<"
)

# ======================================================================
# 1. LOAD DATA (CAP AT max_n_per_dataset INSTANCES PER DATASET)
# ======================================================================

if (file.exists("mimic_final_22_features_60pct.csv")) {
  dev_raw <- readr::read_csv("mimic_final_22_features_60pct.csv", show_col_types = FALSE)
  val_raw <- readr::read_csv("mcmed_final_22_features_60pct.csv",  show_col_types = FALSE)
} else {
  stop("Missing data files. Run data generation script first.")
}

prep_df <- function(df) {
  df %>%
    dplyr::mutate(outcome = factor(outcome, levels = c("No", "Yes"))) %>%
    dplyr::select(outcome, dplyr::everything())
}

set.seed(cfg$seed_cv)

df_dev_full  <- prep_df(dev_raw)
df_test_full <- prep_df(val_raw)

n_dev  <- min(max_n_per_dataset, nrow(df_dev_full))
n_test <- min(max_n_per_dataset, nrow(df_test_full))

df_train <- df_dev_full  %>% dplyr::sample_n(size = n_dev)
df_test  <- df_test_full %>% dplyr::sample_n(size = n_test)

# ======================================================================
# 2. IMPUTATION (MEDIANS ON TRAIN APPLIED TO TEST)
# ======================================================================

tr_medians <- df_train %>%
  dplyr::summarise(dplyr::across(where(is.numeric), median, na.rm = TRUE))

df_train <- df_train %>%
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ ifelse(is.na(.), tr_medians[[dplyr::cur_column()]], .)
    )
  )

missing <- setdiff(names(df_train), names(df_test))
if (length(missing) > 0) {
  for (c in missing) df_test[[c]] <- NA_real_
}

df_test <- df_test %>%
  dplyr::select(dplyr::all_of(names(df_train))) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ ifelse(is.na(.), tr_medians[[dplyr::cur_column()]], .)
    )
  )

# ======================================================================
# 3. FEATURE ENGINEERING (ENGINEERED LAB AND VITAL FEATURES)
# ======================================================================

add_features <- function(df) {
  nm <- names(df)
  
  has_age  <- "age"        %in% nm
  has_spo2 <- "spo2"       %in% nm
  has_hr   <- "heart_rate" %in% nm
  has_rr   <- "resp_rate"  %in% nm
  has_sbp  <- "sbp"        %in% nm
  has_dbp  <- "dbp"        %in% nm
  has_map  <- "map"        %in% nm
  
  idx_neut  <- grep("neut",     nm, ignore.case = TRUE)[1]
  idx_lymph <- grep("lymph",    nm, ignore.case = TRUE)[1]
  idx_mono  <- grep("mono",     nm, ignore.case = TRUE)[1]
  idx_crp   <- grep("crp",      nm, ignore.case = TRUE)[1]
  idx_dd    <- grep("d.?dimer", nm, ignore.case = TRUE)[1]
  idx_alb   <- grep("albumin",  nm, ignore.case = TRUE)[1]
  
  has_nlr <- !is.na(idx_neut)  && !is.na(idx_lymph)
  has_mlr <- !is.na(idx_mono)  && !is.na(idx_lymph)
  has_car <- !is.na(idx_crp)   && !is.na(idx_alb)
  has_clr <- !is.na(idx_crp)   && !is.na(idx_lymph)
  has_crp <- !is.na(idx_crp)
  has_dd  <- !is.na(idx_dd)
  has_alb <- !is.na(idx_alb)
  
  df %>%
    dplyr::mutate(
      age_sq    = if (has_age)  age^2  else NULL,
      spo2_sq   = if (has_spo2) spo2^2 else NULL,
      age_spo2  = if (has_age  && has_spo2) age * spo2        else NULL,
      age_hr    = if (has_age  && has_hr)   age * heart_rate  else NULL,
      age_rr    = if (has_age  && has_rr)   age * resp_rate   else NULL,
      spo2_rr   = if (has_spo2 && has_rr)   spo2 * resp_rate  else NULL,
      spo2_hr   = if (has_spo2 && has_hr)   spo2 * heart_rate else NULL,
      hr_rr     = if (has_hr   && has_rr)   heart_rate * resp_rate else NULL,
      sbp_dbp   = if (has_sbp  && has_dbp)  sbp * dbp         else NULL,
      map_rr    = if (has_map  && has_rr)   map * resp_rate   else NULL,
      
      log_crp   = if (has_crp) log1p(df[[nm[idx_crp]]]) else NULL,
      log_dd    = if (has_dd)  log1p(df[[nm[idx_dd]]])  else NULL,
      log_alb   = if (has_alb) log1p(df[[nm[idx_alb]]]) else NULL,
      
      nlr       = if (has_nlr)
        df[[nm[idx_neut]]] / pmax(df[[nm[idx_lymph]]], 1e-3)
      else NULL,
      mlr       = if (has_mlr)
        df[[nm[idx_mono]]] / pmax(df[[nm[idx_lymph]]], 1e-3)
      else NULL,
      car       = if (has_car)
        df[[nm[idx_crp]]] / pmax(df[[nm[idx_alb]]], 1e-3)
      else NULL,
      clr       = if (has_clr)
        df[[nm[idx_crp]]] / pmax(df[[nm[idx_lymph]]], 1e-3)
      else NULL,
      
      spo2_crp  = if (has_spo2 && has_crp)
        spo2 * df[[nm[idx_crp]]]
      else NULL,
      age_crp   = if (has_age  && has_crp)
        age * df[[nm[idx_crp]]]
      else NULL
    )
}

df_train <- add_features(df_train)
df_test  <- add_features(df_test)

# ======================================================================
# 4. METRICS UTILITIES
# ======================================================================

all_cols <- names(df_train)[names(df_train) != "outcome"]

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
if (length(feat_triage) == 0) {
  feat_triage <- grep("age|sex|spo2|heart|resp|temp|sbp|dbp|map",
                      all_cols, value = TRUE)
}
feat_full <- all_cols

compute_metrics <- function(obs, pred, prob, pos, neg) {
  tab <- table(
    factor(obs,  levels = c(neg, pos)),
    factor(pred, levels = c(neg, pos))
  )
  
  TP <- as.numeric(tab[pos, pos])
  TN <- as.numeric(tab[neg, neg])
  FP <- as.numeric(tab[neg, pos])
  FN <- as.numeric(tab[pos, neg])
  
  sens <- if ((TP + FN) == 0) 0 else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) 0 else TN / (TN + FP)
  ppv  <- if ((TP + FP) == 0) 0 else TP / (TP + FP)
  f1   <- if ((ppv + sens) == 0) 0 else 2 * ppv * sens / (ppv + sens)
  acc  <- if (sum(tab) == 0) 0 else sum(diag(tab)) / sum(tab)
  
  mcc_num <- (TP * TN - FP * FN)
  mcc_den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc     <- if (is.na(mcc_den) || mcc_den == 0) 0 else mcc_num / mcc_den
  
  auc_v <- tryCatch(
    as.numeric(pROC::auc(pROC::roc(obs, prob, quiet = TRUE))),
    error = function(e) 0.5
  )
  
  c(
    AUC         = auc_v,
    MCC         = mcc,
    F1          = f1,
    Accuracy    = acc,
    Precision   = ppv,
    Sensitivity = sens,
    Specificity = spec
  )
}

confusion_counts <- function(obs, pred, pos, neg) {
  tab <- table(
    factor(obs,  levels = c(neg, pos)),
    factor(pred, levels = c(neg, pos))
  )
  TP <- as.numeric(tab[pos, pos])
  TN <- as.numeric(tab[neg, neg])
  FP <- as.numeric(tab[neg, pos])
  FN <- as.numeric(tab[pos, neg])
  c(TP = TP, TN = TN, FP = FP, FN = FN)
}

best_thresh_mcc <- function(obs, prob, pos = "Yes", neg = "No") {
  mask  <- !is.na(prob)
  probs <- as.numeric(prob[mask])
  if (length(probs) == 0) return(0.5)
  
  obs_fac <- factor(obs[mask], levels = c(neg, pos))
  
  thr <- unique(probs)
  if (length(thr) > 200) {
    thr <- stats::quantile(
      thr,
      probs = seq(0.01, 0.99, length.out = 200),
      na.rm = TRUE
    )
    thr <- unique(as.numeric(thr))
  }
  
  best_m <- -Inf
  best_t <- 0.5
  
  for (t in thr) {
    pred <- factor(ifelse(probs >= t, pos, neg), levels = c(neg, pos))
    tab  <- table(obs_fac, pred)
    
    TP <- as.numeric(tab[pos, pos])
    TN <- as.numeric(tab[neg, neg])
    FP <- as.numeric(tab[neg, pos])
    FN <- as.numeric(tab[pos, neg])
    
    num <- TP * TN - FP * FN
    den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- if (is.na(den) || den == 0) NA_real_ else num / den
    
    if (!is.na(mcc) && mcc > best_m) {
      best_m <- mcc
      best_t <- t
    }
  }
  
  best_t
}

fit_calibrator <- function(obs, prob, pos_label) {
  y  <- ifelse(obs == pos_label, 1, 0)
  df <- data.frame(y = y, p = prob)
  df$p <- pmin(pmax(df$p, 1e-6), 1 - 1e-6)
  
  tryCatch(
    stats::glm(y ~ p, data = df, family = stats::binomial),
    error = function(e) NULL
  )
}

apply_calibration <- function(prob, calib_model) {
  if (is.null(calib_model)) return(prob)
  p <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  stats::predict(calib_model, newdata = data.frame(p = p), type = "response")
}

mccSummary <- function(data, lev = NULL, model = NULL) {
  prob_pos <- data[[cfg$pos_label]]
  
  mets <- compute_metrics(
    obs  = data$obs,
    pred = data$pred,
    prob = prob_pos,
    pos  = cfg$pos_label,
    neg  = cfg$neg_label
  )
  
  c(
    MCC = as.numeric(mets["MCC"]),
    ROC = as.numeric(mets["AUC"])
  )
}

# ======================================================================
# 4.2 MODELS AND TRAIN CONTROL
# ======================================================================

models <- list(
  "LR"   = list(method = "glmnet",    tuneLength = 5),
  "RF"   = list(method = "rf",        tuneLength = 5),
  "SVM"  = list(method = "svmRadial", tuneLength = 5),
  "k-NN" = list(method = "kknn",      tuneLength = 5),
  "C4.5" = list(method = "J48",       tuneLength = 5)
)

ctrl <- caret::trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = mccSummary,
  savePredictions = "final"
)

# ======================================================================
# 5. RUN ALGORITHM (CLASS WEIGHTS + MCC THRESHOLD FOR T AND C)
# ======================================================================

run_algo <- function(algo, feats, df_train_cur, df_test_cur) {
  set.seed(cfg$seed_cv)
  
  tr_x <- df_train_cur %>% dplyr::select(dplyr::all_of(feats))
  tr_y <- df_train_cur$outcome
  
  preProc <- caret::preProcess(tr_x, method = c("center", "scale"))
  tr_x_pp <- stats::predict(preProc, tr_x)
  
  n_pos <- sum(tr_y == cfg$pos_label)
  n_neg <- sum(tr_y == cfg$neg_label)
  n_tot <- n_pos + n_neg
  
  w_pos <- n_tot / (2 * n_pos)
  w_neg <- n_tot / (2 * n_neg)
  
  w_vec   <- ifelse(tr_y == cfg$pos_label, w_pos, w_neg)
  classwt <- c(w_neg, w_pos)
  names(classwt) <- c(cfg$neg_label, cfg$pos_label)
  
  if (algo == "LR") {
    fit <- tryCatch(
      caret::train(
        x          = tr_x_pp,
        y          = tr_y,
        method     = models[[algo]]$method,
        trControl  = ctrl,
        metric     = "MCC",
        tuneLength = models[[algo]]$tuneLength,
        family     = "binomial",
        weights    = if (use_class_weights) w_vec else NULL
      ),
      error = function(e) NULL
    )
  } else if (algo == "RF") {
    fit <- tryCatch(
      caret::train(
        x          = tr_x_pp,
        y          = tr_y,
        method     = models[[algo]]$method,
        trControl  = ctrl,
        metric     = "MCC",
        tuneLength = models[[algo]]$tuneLength,
        ntree      = 500,
        classwt    = if (use_class_weights) classwt else NULL
      ),
      error = function(e) NULL
    )
  } else if (algo == "SVM") {
    fit <- tryCatch(
      caret::train(
        x             = tr_x_pp,
        y             = tr_y,
        method        = models[[algo]]$method,
        trControl     = ctrl,
        metric        = "MCC",
        tuneLength    = models[[algo]]$tuneLength,
        class.weights = if (use_class_weights) classwt else NULL
      ),
      error = function(e) NULL
    )
  } else {
    fit <- tryCatch(
      caret::train(
        x          = tr_x_pp,
        y          = tr_y,
        method     = models[[algo]]$method,
        trControl  = ctrl,
        metric     = "MCC",
        tuneLength = models[[algo]]$tuneLength
      ),
      error = function(e) NULL
    )
  }
  
  if (is.null(fit)) return(NULL)
  
  pred_cv <- fit$pred
  if (!is.null(fit$bestTune)) {
    for (nm_bt in names(fit$bestTune)) {
      pred_cv <- pred_cv[pred_cv[[nm_bt]] == fit$bestTune[[nm_bt]], , drop = FALSE]
    }
  }
  
  prob_cv_raw <- pred_cv[[cfg$pos_label]]
  obs_cv      <- pred_cv$obs
  
  use_calibration <- algo %in% c("RF", "SVM", "k-NN", "C4.5")
  calib_model <- NULL
  
  if (use_calibration) {
    calib_model <- fit_calibrator(obs_cv, prob_cv_raw, cfg$pos_label)
    prob_cv     <- apply_calibration(prob_cv_raw, calib_model)
  } else {
    prob_cv <- prob_cv_raw
  }
  
  thr_mcc <- best_thresh_mcc(
    obs  = obs_cv,
    prob = prob_cv,
    pos  = cfg$pos_label,
    neg  = cfg$neg_label
  )
  
  tr_prob_raw <- stats::predict(fit, newdata = tr_x_pp, type = "prob")[, cfg$pos_label]
  tr_prob     <- if (use_calibration) apply_calibration(tr_prob_raw, calib_model) else tr_prob_raw
  
  te_x    <- df_test_cur %>% dplyr::select(dplyr::all_of(feats))
  te_x_pp <- stats::predict(preProc, te_x)
  
  prob_test_raw <- stats::predict(fit, newdata = te_x_pp, type = "prob")[, cfg$pos_label]
  prob_test     <- if (use_calibration) apply_calibration(prob_test_raw, calib_model) else prob_test_raw
  
  pred_test <- factor(
    ifelse(prob_test >= thr_mcc, cfg$pos_label, cfg$neg_label),
    levels = c(cfg$neg_label, cfg$pos_label)
  )
  
  list(
    p_train   = tr_prob,
    p_test    = prob_test,
    pred_test = pred_test
  )
}

# ======================================================================
# 6. FIT TRIAGE AND COMPLETE MODELS ONCE
# ======================================================================

df_train_use <- df_train
y_train_vec  <- df_train_use$outcome

triage_res   <- vector("list", length(models))
complete_res <- vector("list", length(models))
names(triage_res)   <- names(models)
names(complete_res) <- names(models)

for (algo in names(models)) {
  message(sprintf("Fitting %s (triage)...", algo))
  res_T <- run_algo(algo, feat_triage, df_train_use, df_test)
  message(sprintf("Fitting %s (full)...", algo))
  res_C <- run_algo(algo, feat_full,   df_train_use, df_test)
  
  triage_res[[algo]]   <- res_T
  complete_res[[algo]] <- res_C
}

# ======================================================================
# 7. CASCADE AS FUNCTION FOR A GIVEN DEFERRAL BAND
# ======================================================================

run_cascade_for_band <- function(t_low, t_high) {
  res_list <- vector("list", length(models))
  k <- 1L
  
  for (algo in names(models)) {
    res_T <- triage_res[[algo]]
    res_C <- complete_res[[algo]]
    
    if (is.null(res_T) || is.null(res_C)) next
    
    # triage metrics
    mT <- compute_metrics(
      df_test$outcome,
      res_T$pred_test,
      res_T$p_test,
      cfg$pos_label,
      cfg$neg_label
    )
    conf_T <- confusion_counts(
      df_test$outcome,
      res_T$pred_test,
      cfg$pos_label,
      cfg$neg_label
    )
    
    # complete metrics
    mC <- compute_metrics(
      df_test$outcome,
      res_C$pred_test,
      res_C$p_test,
      cfg$pos_label,
      cfg$neg_label
    )
    conf_C <- confusion_counts(
      df_test$outcome,
      res_C$pred_test,
      cfg$pos_label,
      cfg$neg_label
    )
    
    # cascade on train (MCC-optimised threshold)
    defer_tr  <- (res_T$p_train > t_low & res_T$p_train < t_high)
    p_train_S <- ifelse(defer_tr, res_C$p_train, res_T$p_train)
    
    best_thr_S <- best_thresh_mcc(
      obs  = y_train_vec,
      prob = p_train_S,
      pos  = cfg$pos_label,
      neg  = cfg$neg_label
    )
    
    # cascade on test
    defer_te <- (res_T$p_test > t_low & res_T$p_test < t_high)
    p_S      <- ifelse(defer_te, res_C$p_test, res_T$p_test)
    pred_S   <- factor(
      ifelse(p_S >= best_thr_S, cfg$pos_label, cfg$neg_label),
      levels = c(cfg$neg_label, cfg$pos_label)
    )
    
    mS <- compute_metrics(
      df_test$outcome,
      pred_S,
      p_S,
      cfg$pos_label,
      cfg$neg_label
    )
    conf_S <- confusion_counts(
      df_test$outcome,
      pred_S,
      cfg$pos_label,
      cfg$neg_label
    )
    
    defer_rate_test <- mean(defer_te)
    
    metrics_names <- names(mT)
    
    res_list[[k]] <- tibble::tibble(
      Algorithm   = algo,
      t_low_band  = t_low,
      t_high_band = t_high,
      Metric_raw  = c(
        metrics_names,
        "TP", "TN", "FP", "FN",
        "DeferRate", "t_low", "t_high"
      ),
      Mean_T = c(
        mT,
        conf_T["TP"], conf_T["TN"], conf_T["FP"], conf_T["FN"],
        NA_real_, NA_real_, NA_real_
      ),
      Mean_C = c(
        mC,
        conf_C["TP"], conf_C["TN"], conf_C["FP"], conf_C["FN"],
        NA_real_, NA_real_, NA_real_
      ),
      Mean_S = c(
        mS,
        conf_S["TP"], conf_S["TN"], conf_S["FP"], conf_S["FN"],
        defer_rate_test, t_low, t_high
      )
    )
    
    k <- k + 1L
  }
  
  dplyr::bind_rows(res_list)
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
  res_grid_list[[i]] <- run_cascade_for_band(
    bands$t_low_band[i],
    bands$t_high_band[i]
  )
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
      "TP",
      "TN",
      "FP",
      "FN",
      "DeferRate",
      "t_low",
      "t_high"
    )
  )

res_all <- res_all_num %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("Mean"),
      ~ ifelse(is.na(.), NA_character_, sprintf("%.4f", .))
    )
  )

cat(
  "\n--- Grid results over deferral thresholds (",
  max_n_per_dataset, " per dataset, ",
  ifelse(use_class_weights, "class weights, ", "no class weights, "),
  "engineered features, repeated CV, cascade with MCC-optimised cascade threshold) ---\n",
  sep = ""
)

# top 3 bands per algorithm according to cascade MCC (Mean_S)
res_mcc_best <- res_all_num %>%
  dplyr::filter(Metric_raw == "MCC") %>%
  dplyr::group_by(Algorithm) %>%
  dplyr::slice_max(order_by = Mean_S, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(Algorithm, t_low_band, t_high_band, Mean_S_MCC = Mean_S)

print(res_mcc_best, n = Inf)

# detailed table for the chosen baseline band t_low_fixed and t_high_fixed
baseline_all_num <- res_all_num %>%
  dplyr::filter(
    abs(t_low_band  - t_low_fixed)  < 1e-9,
    abs(t_high_band - t_high_fixed) < 1e-9
  ) %>%
  dplyr::arrange(Algorithm, Metric_raw)

res_print <- baseline_all_num %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::starts_with("Mean"),
      ~ ifelse(is.na(.), NA_character_, sprintf("%.4f", .))
    )
  )

cat(
  "\n--- Final Results for baseline band (",
  max_n_per_dataset, " per dataset, ",
  ifelse(use_class_weights, "class weights, ", "no class weights, "),
  "engineered features, repeated CV, cascade with deferral thresholds ",
  t_low_fixed, " and ", t_high_fixed,
  ", triage/full MCC thresholds, cascade threshold MCC-optimised, with TP TN FP FN) ---\n",
  sep = ""
)

print(res_print, n = Inf)

# optional explicit metrics block for the chosen band (same as baseline)
chosen_t_low  <- t_low_fixed
chosen_t_high <- t_high_fixed

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
  chosen_t_low,
  ", t_high = ",
  chosen_t_high,
  " ---\n",
  sep = ""
)

print(band_metrics, n = Inf)

# AUCs per algorithm and model for the chosen band
auc_table <- band_metrics %>%
  dplyr::filter(Metric_raw == "AUC") %>%
  dplyr::select(Algorithm, Mean_T, Mean_C, Mean_S) %>%
  dplyr::mutate(
    dplyr::across(dplyr::starts_with("Mean"), ~ round(., 4))
  )

cat("\n--- AUCs per algorithm and model for chosen band ---\n")
print(auc_table, n = Inf)
