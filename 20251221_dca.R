# ===============================================================
# COMPLETE SCRIPT — DCA curves (Complete, Triage, Cascade) per algorithm
# Uses YOUR current datasets:
#   dev_path = "mimic_final_22_features_60pct.csv"
#   ext_path = "mcmed_final_22_features_60pct.csv"
#
# - Development = df_train (from dev CSV)
# - External    = df_test  (from ext CSV)
# - tuneLength  = 5 for ALL algorithms
# - Class weights optional
# - Platt scaling fit on development OOF, applied to external
# - Cascade uses FIXED deferral band by default (t_low_fixed/t_high_fixed)
#   (optionally tune the band on triage OOF via tune_band = TRUE)
#
# NOTE: Does NOT write files unless save_png <- TRUE
# ===============================================================

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
  library(patchwork)
  library(cowplot)
  library(scales)
  library(grid)
})

options(warn = -1)
options(na.action = "na.pass")
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ===============================================================
# 0) CONFIG
# ===============================================================
max_n_per_dataset  <- 5000L
use_class_weights  <- TRUE
tune_len_all       <- 5L

t_low_fixed  <- 0.33
t_high_fixed <- 0.66

cfg <- list(
  pos_label      = "Yes",
  neg_label      = "No",
  seed_cv        = 123,
  t_low_default  = t_low_fixed,
  t_high_default = t_high_fixed,
  cascade = list(
    enabled   = TRUE,
    tune_band = FALSE,         # <- set TRUE if you want to tune band on TRIAGE OOF
    t_low     = t_low_fixed,
    t_high    = t_high_fixed,
    grid_tl   = seq(0.05, 0.45, by = 0.05),
    grid_th   = seq(0.55, 0.95, by = 0.05),
    max_defer = 1.0
  )
)

models_plot_order <- c("C4.5", "k-NN", "SVM", "RF", "LR")  # for facet/grid order
models_fit        <- c("LR", "RF", "SVM", "k-NN", "C4.5")  # for fitting loop

# DCA bootstrap ribbons
B_boot    <- 200L     # increase if you want smoother CI (e.g., 1000)
save_png  <- FALSE
png_path  <- "dca_mimic_mcmed.png"

set.seed(cfg$seed_cv)

# ===============================================================
# 1) LOAD DATA (YOUR DATASETS)
# ===============================================================
dev_path <- "mimic_final_22_features_60pct.csv"
ext_path <- "mcmed_final_22_features_60pct.csv"

if (!file.exists(dev_path) || !file.exists(ext_path)) {
  stop("Input CSV files were not found: ", dev_path, " and/or ", ext_path)
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
          .data$group == "COVID"        ~ "Yes",
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

df_dev_full <- prep_df(dev_raw)
df_ext_full <- prep_df(ext_raw)

n_dev <- min(max_n_per_dataset, nrow(df_dev_full))
n_ext <- min(max_n_per_dataset, nrow(df_ext_full))

df_train <- df_dev_full %>% dplyr::slice_sample(n = n_dev)
df_test  <- df_ext_full %>% dplyr::slice_sample(n = n_ext)

# enforce outcome levels
df_train$outcome <- factor(df_train$outcome, levels = c(cfg$neg_label, cfg$pos_label))
df_test$outcome  <- factor(df_test$outcome,  levels = c(cfg$neg_label, cfg$pos_label))

# ===============================================================
# 2) ALIGN COLUMNS + MEDIAN IMPUTE NUMERICS (TRAIN MEDIANS -> TEST)
# ===============================================================
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
  if (length(extra) > 0L) test <- test %>% dplyr::select(-dplyr::all_of(extra))
  
  test %>% dplyr::select(dplyr::all_of(names(train)))
}

median_impute_numeric <- function(train, test) {
  med_tbl <- train %>%
    dplyr::summarise(dplyr::across(where(is.numeric), ~ stats::median(., na.rm = TRUE)))
  med_list <- as.list(med_tbl)
  
  train2 <- train %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.), med_list[[dplyr::cur_column()]], .)))
  
  test2 <- test %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.), med_list[[dplyr::cur_column()]], .)))
  
  list(train = train2, test = test2)
}

df_test <- align_test_to_train(df_train, df_test)
imp1 <- median_impute_numeric(df_train, df_test)
df_train <- imp1$train
df_test  <- imp1$test

# ===============================================================
# 3) FEATURE ENGINEERING (then align + impute again)
# ===============================================================
add_features <- function(df) {
  df2 <- df
  
  safe_ratio <- function(a, b) {
    out <- a / b
    out[!is.finite(out)] <- NA_real_
    out
  }
  
  if (all(c("neutrophils_raw", "lymphocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(NLR = safe_ratio(.data$neutrophils_raw, .data$lymphocytes_raw))
  
  if (all(c("monocytes_raw", "lymphocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(MLR = safe_ratio(.data$monocytes_raw, .data$lymphocytes_raw))
  
  if (all(c("neutrophils_raw", "monocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(NMR = safe_ratio(.data$neutrophils_raw, .data$monocytes_raw))
  
  if (all(c("neutrophils_raw", "monocytes_raw", "lymphocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(SIRI = safe_ratio(.data$neutrophils_raw * .data$monocytes_raw, .data$lymphocytes_raw))
  
  if (all(c("crp", "lymphocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(CLR = safe_ratio(.data$crp, .data$lymphocytes_raw))
  
  if (all(c("crp", "albumin") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(CAR = safe_ratio(.data$crp, .data$albumin))
  
  if (all(c("albumin", "lymphocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(PNI = .data$albumin + 5 * .data$lymphocytes_raw)
  
  if (all(c("d_dimer", "lymphocytes_raw") %in% names(df2)))
    df2 <- df2 %>% dplyr::mutate(DLR = safe_ratio(.data$d_dimer, .data$lymphocytes_raw))
  
  df2
}

df_train <- add_features(df_train)
df_test  <- add_features(df_test)

df_test <- align_test_to_train(df_train, df_test)
imp2 <- median_impute_numeric(df_train, df_test)
df_train <- imp2$train
df_test  <- imp2$test

all_cols <- setdiff(names(df_train), "outcome")

triage_candidates <- c("age","sex","spo2","heart_rate","resp_rate","temp","sbp","dbp","map")
feat_triage <- intersect(triage_candidates, all_cols)
if (!length(feat_triage)) {
  feat_triage <- grep("age|sex|spo2|heart|resp|temp|sbp|dbp|map", all_cols, value = TRUE)
}
feat_full <- all_cols

# ===============================================================
# 4) caret summary (AUC), weights, Platt scaling, ROC threshold
# ===============================================================
to_num <- function(x) as.numeric(unlist(x, use.names = FALSE))

aucSummary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) lev <- levels(data$obs)
  pos <- lev[2]
  prob_pos <- to_num(data[[pos]])
  mask <- is.finite(prob_pos) & !is.na(data$obs)
  if (sum(mask) < 2L) return(c(ROC = NA_real_))
  roc_obj <- try(pROC::roc(response = data$obs[mask], predictor = prob_pos[mask],
                           levels = lev, direction = "<"), silent = TRUE)
  if (inherits(roc_obj, "try-error")) c(ROC = NA_real_) else c(ROC = as.numeric(pROC::auc(roc_obj)))
}

ctrl_auc <- caret::trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = aucSummary,
  savePredictions = "final",
  allowParallel   = TRUE
)

method_from_algo <- function(algo) {
  switch(algo,
         "LR"   = "glmnet",
         "RF"   = "rf",
         "SVM"  = "svmRadial",
         "k-NN" = "kknn",
         "C4.5" = "J48",
         stop("Unknown algo: ", algo))
}

compute_case_weights <- function(y) {
  if (!isTRUE(use_class_weights)) return(rep(1, length(y)))
  tab <- table(y)
  n_pos <- as.numeric(tab[cfg$pos_label])
  n_neg <- as.numeric(tab[cfg$neg_label])
  w_pos <- if (!is.na(n_pos) && n_pos > 0) n_neg / n_pos else 1
  ifelse(y == cfg$pos_label, w_pos, 1)
}

extract_oof_best <- function(fit) {
  oof <- fit$pred
  if (is.null(oof) || !nrow(oof)) stop("Out-of-fold predictions were not available from caret.")
  bt <- fit$bestTune
  for (nm in names(bt)) oof <- dplyr::filter(oof, .data[[nm]] == bt[[nm]])
  oof
}

fit_calibrator <- function(obs, prob, pos_label) {
  y <- ifelse(obs == pos_label, 1, 0)
  p <- pmin(pmax(to_num(prob), 1e-6), 1 - 1e-6)
  suppressWarnings(stats::glm(y ~ qlogis(p), data = data.frame(y = y, p = p),
                              family = stats::binomial(link = "logit")))
}

apply_calibration <- function(prob, calib_model) {
  p <- pmin(pmax(to_num(prob), 1e-6), 1 - 1e-6)
  if (is.null(calib_model)) return(p)
  lp <- suppressWarnings(predict(calib_model, newdata = data.frame(p = p), type = "link"))
  to_num(1 / (1 + exp(-lp)))
}

best_thresh_roc <- function(obs, prob, pos = cfg$pos_label, neg = cfg$neg_label) {
  prob <- to_num(prob)
  mask <- is.finite(prob) & !is.na(obs)
  if (sum(mask) < 2L) return(0.5)
  roc_obj <- try(pROC::roc(response = obs[mask], predictor = prob[mask],
                           levels = c(neg, pos), direction = "<"), silent = TRUE)
  if (inherits(roc_obj, "try-error")) return(0.5)
  thr <- try(pROC::coords(roc_obj, x = "best", best.method = "youden",
                          ret = "threshold", transpose = FALSE), silent = TRUE)
  if (inherits(thr, "try-error")) return(0.5)
  thr_num <- suppressWarnings(as.numeric(unlist(thr))[1])
  if (!is.finite(thr_num)) 0.5 else thr_num
}

# Align factor levels between train/test for predictors (avoid predict() failures)
align_train_test_levels <- function(dtr, dte, feats) {
  feats <- intersect(feats, intersect(names(dtr), names(dte)))
  for (nm in feats) {
    xtr <- dtr[[nm]]
    xte <- dte[[nm]]
    
    is_cat <- is.factor(xtr) || is.character(xtr) || is.factor(xte) || is.character(xte)
    if (!is_cat) next
    
    xtr <- as.character(xtr)
    xte <- as.character(xte)
    
    xtr[is.na(xtr)] <- "Unknown"
    xte[is.na(xte)] <- "Unknown"
    
    levs <- union(unique(xtr), unique(xte))
    levs <- unique(c(levs, "Unknown"))
    
    # map unseen in test to Unknown
    xte[!(xte %in% levs)] <- "Unknown"
    xtr[!(xtr %in% levs)] <- "Unknown"
    
    dtr[[nm]] <- factor(xtr, levels = levs)
    dte[[nm]] <- factor(xte, levels = levs)
  }
  list(train = dtr, test = dte)
}

# Numeric NA safety (train medians -> test)
impute_numeric_train_median <- function(dtr, dte, feats) {
  feats <- intersect(feats, intersect(names(dtr), names(dte)))
  num_feats <- feats[vapply(dtr[feats], is.numeric, logical(1))]
  if (!length(num_feats)) return(list(train = dtr, test = dte))
  meds <- vapply(num_feats, function(nm) stats::median(dtr[[nm]], na.rm = TRUE), numeric(1))
  for (nm in num_feats) {
    med <- meds[[nm]]
    if (!is.finite(med)) med <- 0
    dtr[[nm]][is.na(dtr[[nm]])] <- med
    dte[[nm]][is.na(dte[[nm]])] <- med
  }
  list(train = dtr, test = dte)
}

# LR needs >=2 columns for glmnet
add_lr_shadow_if_needed <- function(d, feats, algo) {
  if (identical(algo, "LR") && length(feats) == 1L) {
    nm <- "shadow_LR"
    while (nm %in% names(d)) nm <- paste0("shadow_LR_", sample(1000:9999, 1))
    set.seed(cfg$seed_cv + 17)
    d[[nm]] <- rnorm(nrow(d), sd = 1e-8)
    feats <- c(feats, nm)
  }
  list(df = d, feats = feats)
}

# ===============================================================
# 5) Holdout runner: fit on dev, calibrate on OOF, score external
# ===============================================================
run_holdout_dca <- function(df_train, df_test, yvar, features, algo, set_label) {
  pos <- cfg$pos_label; neg <- cfg$neg_label
  feats <- intersect(features, setdiff(intersect(names(df_train), names(df_test)), yvar))
  if (!length(feats)) return(NULL)
  
  dtr <- df_train[, c(yvar, feats), drop = FALSE]
  dte <- df_test[,  c(yvar, feats), drop = FALSE]
  
  # outcome levels
  dtr[[yvar]] <- factor(dtr[[yvar]], levels = c(neg, pos))
  dte[[yvar]] <- factor(dte[[yvar]], levels = c(neg, pos))
  
  # align categoricals + numeric NA safety
  al1 <- align_train_test_levels(dtr, dte, feats); dtr <- al1$train; dte <- al1$test
  al2 <- impute_numeric_train_median(dtr, dte, feats); dtr <- al2$train; dte <- al2$test
  
  # LR shadow if needed
  sh <- add_lr_shadow_if_needed(dtr, feats, algo); dtr <- sh$df; feats <- sh$feats
  # keep test aligned with shadow if added
  if (any(grepl("^shadow_LR", feats)) && !("shadow_LR" %in% names(dte))) {
    new_cols <- setdiff(names(dtr), names(dte))
    for (nm in new_cols) dte[[nm]] <- dtr[[nm]] * 0
    dte <- dte[, names(dtr), drop = FALSE]
  }
  
  w <- compute_case_weights(dtr[[yvar]])
  method <- method_from_algo(algo)
  
  set.seed(cfg$seed_cv + sum(utf8ToInt(paste0(algo, set_label))))
  fit <- switch(
    algo,
    "LR" = caret::train(
      outcome ~ .,
      data       = dtr,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tune_len_all,
      weights    = w,
      family     = "binomial"
    ),
    "RF" = caret::train(
      outcome ~ .,
      data       = dtr,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tune_len_all,
      weights    = w,
      ntree      = 500
    ),
    "SVM" = caret::train(
      outcome ~ .,
      data       = dtr,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tune_len_all,
      weights    = w
    ),
    "k-NN" = caret::train(
      outcome ~ .,
      data       = dtr,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tune_len_all,
      weights    = w
    ),
    "C4.5" = caret::train(
      outcome ~ .,
      data       = dtr,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tune_len_all,
      weights    = w
    )
  )
  
  oof <- extract_oof_best(fit)
  prob_oof_raw <- to_num(oof[[pos]])
  obs_oof      <- factor(oof$obs, levels = c(neg, pos))
  
  calib <- fit_calibrator(obs_oof, prob_oof_raw, pos)
  prob_oof_cal <- apply_calibration(prob_oof_raw, calib)
  thr_train    <- best_thresh_roc(obs_oof, prob_oof_cal, pos = pos, neg = neg)
  
  prob_test_raw <- to_num(predict(fit, newdata = dte, type = "prob")[, pos])
  prob_test_cal <- apply_calibration(prob_test_raw, calib)
  
  list(
    algo      = algo,
    set       = set_label,
    obs       = factor(dte[[yvar]], levels = c(neg, pos)),
    p         = prob_test_cal,
    oof       = tibble::tibble(obs = obs_oof, p = prob_oof_cal),
    threshold = thr_train
  )
}

# ===============================================================
# 6) DCA tables + stratified bootstrap ribbons
# ===============================================================
decision_curve_table <- function(obs, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  pos <- cfg$pos_label
  y_raw <- as.integer(obs == pos)
  p <- to_num(p)
  keep <- is.finite(p) & !is.na(y_raw)
  y <- y_raw[keep]; p <- p[keep]
  
  if (!length(y)) {
    return(tibble::tibble(
      threshold = thresholds,
      Net_Benefit = NA_real_,
      NB_TreatAll = NA_real_,
      NB_TreatNone = 0,
      prevalence = NA_real_
    ))
  }
  
  N <- length(y)
  prev <- mean(y)
  kvec <- thresholds / (1 - thresholds)
  
  purrr::map_dfr(seq_along(thresholds), function(j) {
    pt <- thresholds[j]
    pred <- as.integer(p >= pt)
    TP <- sum(pred == 1 & y == 1)
    FP <- sum(pred == 1 & y == 0)
    NB <- TP / N - FP / N * kvec[j]
    NB_all <- prev - (1 - prev) * kvec[j]
    tibble::tibble(
      threshold  = pt,
      Net_Benefit = NB,
      NB_TreatAll = NB_all,
      NB_TreatNone = 0,
      prevalence = prev
    )
  })
}

dca_boot_band_stratified <- function(obj,
                                     B = B_boot,
                                     thresholds = seq(0.01, 0.99, by = 0.01)) {
  if (is.null(obj)) return(NULL)
  
  y <- as.integer(obj$obs == cfg$pos_label)
  p <- to_num(obj$p)
  
  keep <- is.finite(p) & !is.na(y)
  y <- y[keep]; p <- p[keep]
  
  idx_pos <- which(y == 1L)
  idx_neg <- which(y == 0L)
  P  <- length(idx_pos)
  Nn <- length(idx_neg)
  N  <- P + Nn
  if (P == 0L || Nn == 0L) return(NULL)
  
  kvec <- thresholds / (1 - thresholds)
  nb_mat <- matrix(NA_real_, nrow = B, ncol = length(thresholds))
  
  for (b in seq_len(B)) {
    ip   <- sample(idx_pos, P,  replace = TRUE)
    ineg <- sample(idx_neg, Nn, replace = TRUE)
    idx  <- c(ip, ineg)
    yb   <- y[idx]
    pb   <- p[idx]
    
    for (j in seq_along(thresholds)) {
      pt <- thresholds[j]
      pred1 <- (pb >= pt)
      TP <- sum(pred1 & (yb == 1L))
      FP <- sum(pred1 & (yb == 0L))
      nb_mat[b, j] <- TP / N - FP / N * kvec[j]
    }
  }
  
  qs <- apply(nb_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  tibble::tibble(
    threshold = thresholds,
    NB_L = qs[1, ],
    NB_U = qs[2, ],
    set  = obj$set,
    algo = obj$algo
  )
}

# Optional: tune cascade band on TRIAGE OOF (net benefit at triage threshold)
tune_band_from_triage_oof <- function(oof_df, thr_tri,
                                      grid_tl   = cfg$cascade$grid_tl,
                                      grid_th   = cfg$cascade$grid_th,
                                      max_defer = cfg$cascade$max_defer) {
  stopifnot(all(c("obs", "p") %in% names(oof_df)))
  pos <- cfg$pos_label; neg <- cfg$neg_label
  y <- factor(oof_df$obs, levels = c(neg, pos))
  p <- to_num(oof_df$p)
  
  keep0 <- is.finite(p) & !is.na(y)
  y <- y[keep0]; p <- p[keep0]
  if (length(y) < 2L || length(unique(y)) < 2L) return(list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high))
  
  best <- NULL
  for (tl in grid_tl) for (th in grid_th) if (tl < th) {
    defer <- (p > tl & p < th)
    defer_rate <- mean(defer)
    if (defer_rate > max_defer) next
    keep <- !defer
    if (!any(keep)) next
    
    pred_keep <- factor(ifelse(p[keep] >= thr_tri, pos, neg), levels = c(neg, pos))
    obs_keep  <- factor(y[keep], levels = c(neg, pos))
    N <- length(obs_keep)
    yy <- as.integer(obs_keep == pos)
    pr <- as.integer(pred_keep == pos)
    TP <- sum(pr == 1 & yy == 1)
    FP <- sum(pr == 1 & yy == 0)
    NB <- TP / N - FP / N * (thr_tri / (1 - thr_tri))
    
    cand <- list(t_low = tl, t_high = th, score = NB, defer = defer_rate, n_keep = sum(keep))
    if (is.null(best) ||
        cand$score > best$score + 1e-12 ||
        (abs(cand$score - best$score) < 1e-12 &&
         (cand$defer < best$defer - 1e-12 ||
          (abs(cand$defer - best$defer) < 1e-12 && cand$n_keep > best$n_keep)))) {
      best <- cand
    }
  }
  
  if (is.null(best)) list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high) else best
}

# ===============================================================
# 7) FIT MODELS (Complete + Triage) AND BUILD CASCADE
# ===============================================================
holds_base <- purrr::map(models_fit, function(a) {
  full <- run_holdout_dca(df_train, df_test, "outcome", feat_full,   a, "Complete feature set")
  tri  <- run_holdout_dca(df_train, df_test, "outcome", feat_triage, a, "Triage feature set")
  list(full = full, tri = tri)
})

# flatten + drop NULL
holds_external <- purrr::compact(unlist(holds_base, recursive = FALSE))

# build cascade holds
by_algo <- split(holds_external, vapply(holds_external, function(x) x$algo, character(1)))
cascade_holds <- list()

for (a in names(by_algo)) {
  hx   <- by_algo[[a]]
  full <- purrr::detect(hx, ~ .x$set == "Complete feature set")
  tri  <- purrr::detect(hx, ~ .x$set == "Triage feature set")
  if (is.null(full) || is.null(tri)) next
  
  # choose band
  t_low  <- cfg$cascade$t_low
  t_high <- cfg$cascade$t_high
  if (isTRUE(cfg$cascade$enabled) && isTRUE(cfg$cascade$tune_band) &&
      !is.null(tri$oof) && nrow(tri$oof) > 0 && is.finite(tri$threshold)) {
    tuned <- tune_band_from_triage_oof(tri$oof, thr_tri = tri$threshold)
    t_low  <- tuned$t_low  %||% t_low
    t_high <- tuned$t_high %||% t_high
  }
  
  p_tri  <- to_num(tri$p)
  p_full <- to_num(full$p)
  
  # strict safety
  n <- min(length(p_tri), length(p_full), length(tri$obs))
  p_tri  <- p_tri[seq_len(n)]
  p_full <- p_full[seq_len(n)]
  obs    <- tri$obs[seq_len(n)]
  
  use_full <- (p_tri > t_low & p_tri < t_high)
  p_cas <- ifelse(use_full, p_full, p_tri)
  
  cascade_holds[[length(cascade_holds) + 1L]] <- list(
    obs  = obs,
    p    = p_cas,
    algo = a,
    set  = "Cascade"
  )
}

holds_external <- c(holds_external, cascade_holds)

# ===============================================================
# 8) BUILD DCA TABLES + RIBBONS
# ===============================================================
by_algo <- split(holds_external, vapply(holds_external, function(x) x$algo, character(1)))

combo <- lapply(names(by_algo), function(a) {
  hx <- by_algo[[a]]
  full <- purrr::detect(hx, ~ .x$set == "Complete feature set")
  tri  <- purrr::detect(hx, ~ .x$set == "Triage feature set")
  cas  <- purrr::detect(hx, ~ .x$set == "Cascade")
  
  d_full <- if (!is.null(full)) decision_curve_table(full$obs, full$p) %>% mutate(set = full$set, algo = a) else NULL
  d_tri  <- if (!is.null(tri))  decision_curve_table(tri$obs,  tri$p)  %>% mutate(set = tri$set,  algo = a) else NULL
  d_cas  <- if (!is.null(cas))  decision_curve_table(cas$obs,  cas$p)  %>% mutate(set = cas$set,  algo = a) else NULL
  
  ref <- full %||% tri %||% cas
  t_df <- decision_curve_table(ref$obs, ref$p) %>% transmute(algo = a, threshold, NB_TreatAll, NB_TreatNone)
  
  list(dca = bind_rows(d_full, d_tri, d_cas), treat = t_df)
})

dca_df   <- bind_rows(purrr::map(combo, "dca"))
treat_df <- bind_rows(purrr::map(combo, "treat"))

bands_df <- holds_external %>%
  purrr::map(dca_boot_band_stratified) %>%
  purrr::compact() %>%
  { if (!length(.)) NULL else bind_rows(.) }

# factor levels
dca_df$algo   <- factor(as.character(dca_df$algo), levels = models_plot_order, ordered = TRUE)
treat_df$algo <- factor(as.character(treat_df$algo), levels = models_plot_order, ordered = TRUE)
if (!is.null(bands_df)) bands_df$algo <- factor(as.character(bands_df$algo), levels = models_plot_order, ordered = TRUE)

dca_df$set <- factor(dca_df$set, levels = c("Complete feature set","Triage feature set","Cascade"))
if (!is.null(bands_df)) bands_df$set <- factor(bands_df$set, levels = levels(dca_df$set))

# ===============================================================
# 9) PLOT
# ===============================================================
cb_orange <- "#D55E00"  # Complete
cb_blue   <- "#0072B2"  # Triage
cb_green  <- "#009E73"  # Cascade

y_cap <- max(c(dca_df$Net_Benefit, treat_df$NB_TreatAll, 0), na.rm = TRUE)
if (!is.finite(y_cap) || y_cap <= 0) y_cap <- 0.05
y0 <- -0.02

axis_breaks_y <- function(y_cap) {
  if (y_cap <= 0.06) seq(0, y_cap, by = 0.01) else scales::breaks_extended()(c(0, y_cap))
}

theme_ref <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.line          = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "plain", color = "grey20"),
      legend.title       = ggplot2::element_blank()
    )
}

safe_ribbon_layers <- function(df_b, algo_name) {
  if (is.null(df_b) || !nrow(df_b)) return(list())
  df_b <- df_b %>% filter(algo == algo_name, is.finite(NB_L), is.finite(NB_U), !is.na(set))
  if (!nrow(df_b)) return(list())
  
  out <- list()
  for (s in levels(dca_df$set)) {
    d <- df_b %>% filter(set == s)
    if (!nrow(d)) next
    fill_col <- if (s == "Complete feature set") cb_orange else if (s == "Triage feature set") cb_blue else cb_green
    out[[length(out) + 1]] <- ggplot2::geom_ribbon(
      data = d, aes(threshold, ymin = NB_L, ymax = NB_U),
      fill = fill_col, alpha = 0.18, inherit.aes = FALSE
    )
  }
  out
}

dca_panel <- function(algo_name, show_y = FALSE, show_x = FALSE) {
  df_d <- dca_df %>% filter(algo == algo_name, is.finite(Net_Benefit))
  df_t <- treat_df %>% filter(algo == algo_name, is.finite(NB_TreatAll), is.finite(NB_TreatNone))
  
  ggplot() +
    safe_ribbon_layers(bands_df, algo_name) +
    geom_line(data = df_t, aes(threshold, NB_TreatAll, linetype = "Treat-all"),
              linewidth = 0.6, color = "black", na.rm = TRUE) +
    geom_line(data = df_t, aes(threshold, NB_TreatNone, linetype = "Treat-none"),
              linewidth = 0.6, color = "grey40", na.rm = TRUE) +
    geom_line(data = df_d, aes(threshold, Net_Benefit, color = set),
              linewidth = 0.9, na.rm = TRUE) +
    geom_segment(aes(x = 0, xend = 0, y = y0, yend = y_cap),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_segment(aes(x = 0, xend = 1, y = y0, yend = y0),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    scale_color_manual(values = c("Complete feature set" = cb_orange,
                                  "Triage feature set"   = cb_blue,
                                  "Cascade"              = cb_green),
                       name = NULL) +
    scale_linetype_manual(values = c("Treat-all" = "solid", "Treat-none" = "dashed"),
                          name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    scale_y_continuous(breaks = axis_breaks_y(y_cap),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(y0, y_cap)) +
    labs(
      x = if (show_x) "Threshold probability" else NULL,
      y = if (show_y) "Net benefit" else NULL,
      title = algo_name
    ) +
    theme_ref(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11),
      axis.text  = element_text(size = 8),
      axis.title = element_text(size = 9)
    )
}

# panels in a 2x3 grid + legend at bottom-right
algos_in_plot <- models_plot_order[models_plot_order %in% unique(as.character(dca_df$algo))]

p_list <- lapply(seq_along(algos_in_plot), function(i) {
  a <- algos_in_plot[i]
  show_y <- (i %in% c(1, 4))
  show_x <- (i %in% c(4, 5))
  dca_panel(a, show_y = show_y, show_x = show_x)
})

blank <- ggplot() + theme_void()
fill_panel <- function(i) if (i <= length(p_list)) p_list[[i]] else blank

p1 <- fill_panel(1)
p2 <- fill_panel(2)
p3 <- fill_panel(3)
p4 <- fill_panel(4)
p5 <- fill_panel(5)

legend_src <- ggplot(dca_df %>% filter(algo == algos_in_plot[1]),
                     aes(threshold, Net_Benefit, color = set)) +
  geom_line() +
  scale_color_manual(values = c("Complete feature set" = cb_orange,
                                "Triage feature set"   = cb_blue,
                                "Cascade"              = cb_green),
                     name = NULL) +
  geom_line(data = treat_df %>% filter(algo == algos_in_plot[1]),
            aes(threshold, NB_TreatAll, linetype = "Treat-all"),
            color = "black", inherit.aes = FALSE) +
  geom_line(data = treat_df %>% filter(algo == algos_in_plot[1]),
            aes(threshold, NB_TreatNone, linetype = "Treat-none"),
            color = "grey40", inherit.aes = FALSE) +
  scale_linetype_manual(values = c("Treat-all" = "solid", "Treat-none" = "dashed"), name = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right")

leg_plot <- patchwork::wrap_elements(cowplot::get_legend(legend_src))

p_dca_grid <- (p1 | p2 | p3) / (p4 | p5 | leg_plot)
print(p_dca_grid)

if (isTRUE(save_png)) {
  ggsave(filename = png_path, plot = p_dca_grid,
         width = 250, height = 120, units = "mm", dpi = 300, bg = "white")
}

# Data behind curves (optional inspect)
dca_df
