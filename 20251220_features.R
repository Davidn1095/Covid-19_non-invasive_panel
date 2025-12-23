suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(grid)
  
  library(caret)
  library(pROC)
  library(glmnet)
  library(kknn)
  library(randomForest)
  library(kernlab)
  library(RWeka)
})

options(warn = -1)

`%||%` <- function(a, b) if (!is.null(a)) a else b

seed_from_name <- function(name, offset = 0L) {
  as.integer(offset + sum(utf8ToInt(as.character(name))))
}

# ======================================================================
# 0. REQUIRED OBJECTS
# Development was taken from df_train, external was taken from df_test
# ======================================================================

stopifnot(exists("df_train", inherits = TRUE))
stopifnot(exists("df_test", inherits = TRUE))
stopifnot(exists("feat_triage", inherits = TRUE))
stopifnot(exists("feat_full", inherits = TRUE))
stopifnot(exists("cfg", inherits = TRUE))

df_train <- get("df_train", inherits = TRUE)
df_test  <- get("df_test",  inherits = TRUE)
feat_triage <- get("feat_triage", inherits = TRUE)
feat_full   <- get("feat_full",   inherits = TRUE)
cfg <- get("cfg", inherits = TRUE)

if (is.null(cfg$pos_label) || is.null(cfg$neg_label)) {
  stop("cfg$pos_label and cfg$neg_label are required.")
}

use_class_weights <- if (exists("use_class_weights", inherits = TRUE)) {
  isTRUE(get("use_class_weights", inherits = TRUE))
} else {
  TRUE
}

yvar <- "outcome"
if (!(yvar %in% names(df_train)) || !(yvar %in% names(df_test))) {
  stop("outcome was not found in df_train or df_test.")
}

df_train[[yvar]] <- factor(df_train[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
df_test[[yvar]]  <- factor(df_test[[yvar]],  levels = c(cfg$neg_label, cfg$pos_label))

feat_triage <- intersect(feat_triage, setdiff(names(df_train), yvar))
feat_triage <- intersect(feat_triage, setdiff(names(df_test),  yvar))

feat_full <- intersect(feat_full, setdiff(names(df_train), yvar))
feat_full <- intersect(feat_full, setdiff(names(df_test),  yvar))

if (!length(feat_triage)) stop("feat_triage had no columns present in both df_train and df_test.")
if (!length(feat_full))   stop("feat_full had no columns present in both df_train and df_test.")

algo_order <- c("LR", "RF", "SVM", "k-NN", "C4.5")

max_n_train <- 5000L

cascade_cfg <- list(
  grid_tl   = seq(0.05, 0.45, 0.05),
  grid_th   = seq(0.55, 0.95, 0.05),
  max_defer = 1
)

# ======================================================================
# 1. caret utilities, AUC optimisation, tuneLength fixed to 5 for all
# ======================================================================

aucSummary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) lev <- levels(data$obs)
  if (length(lev) != 2L) stop("aucSummary required a binary outcome.")
  pos <- lev[2]
  roc_obj <- try(
    pROC::roc(
      response  = data$obs,
      predictor = data[[pos]],
      levels    = lev,
      direction = "<"
    ),
    silent = TRUE
  )
  roc_auc <- if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj))
  c(ROC = roc_auc)
}

ctrl_auc <- caret::trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = aucSummary,
  savePredictions = "final",
  verboseIter     = FALSE,
  allowParallel   = TRUE
)

method_from_algo <- function(algo) {
  switch(
    algo,
    "LR"   = "glmnet",
    "RF"   = "rf",
    "SVM"  = "svmRadial",
    "k-NN" = "kknn",
    "C4.5" = "J48",
    stop("Unknown algorithm: ", algo)
  )
}

tune_len_from_algo <- function(algo) {
  5L
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
  if (is.null(oof) || !nrow(oof)) stop("Out of fold predictions were not available from caret.")
  bt <- fit$bestTune
  if (!is.null(bt) && nrow(bt)) {
    for (nm in names(bt)) oof <- oof[oof[[nm]] == bt[[nm]], , drop = FALSE]
  }
  oof
}

# ======================================================================
# 2. Platt scaling by GLM, MCC threshold selection, NA safe
# ======================================================================

fit_platt_glm_from_oof <- function(obs, p_raw, pos_label = cfg$pos_label) {
  y <- ifelse(obs == pos_label, 1, 0)
  p <- pmin(pmax(p_raw, 1e-6), 1 - 1e-6)
  df_cal <- data.frame(y = y, p = p)
  
  fit <- try(
    suppressWarnings(stats::glm(y ~ qlogis(p), data = df_cal, family = stats::binomial(link = "logit"))),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(function(p_new) {
      p2 <- pmin(pmax(p_new, 1e-6), 1 - 1e-6)
      p2
    })
  }
  
  function(p_new) {
    p2 <- pmin(pmax(p_new, 1e-6), 1 - 1e-6)
    lp <- try(stats::predict(fit, newdata = data.frame(p = p2), type = "link"), silent = TRUE)
    if (inherits(lp, "try-error")) return(p2)
    out <- 1 / (1 + exp(-lp))
    out[!is.finite(out)] <- p2[!is.finite(out)]
    out
  }
}

mcc_from_preds <- function(obs, pred, pos_label, neg_label) {
  obs  <- factor(obs,  levels = c(neg_label, pos_label))
  pred <- factor(pred, levels = c(neg_label, pos_label))
  tab <- table(obs, pred)
  
  TN <- as.numeric(tab[neg_label, neg_label])
  FP <- as.numeric(tab[neg_label, pos_label])
  FN <- as.numeric(tab[pos_label, neg_label])
  TP <- as.numeric(tab[pos_label, pos_label])
  
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  ((TP * TN) - (FP * FN)) / denom
}

best_threshold_mcc_med <- function(obs, p_pos,
                                   pos_label = cfg$pos_label,
                                   neg_label = cfg$neg_label,
                                   grid = seq(0.01, 0.99, 0.01)) {
  obs <- factor(obs, levels = c(neg_label, pos_label))
  p_pos[is.na(p_pos)] <- 0.5
  
  mcc <- vapply(grid, function(thr) {
    pred <- ifelse(p_pos >= thr, pos_label, neg_label)
    mcc_from_preds(obs, pred, pos_label, neg_label)
  }, numeric(1))
  
  if (all(is.na(mcc))) return(0.5)
  
  best <- suppressWarnings(max(mcc, na.rm = TRUE))
  if (!is.finite(best)) return(0.5)
  
  keep <- which(abs(mcc - best) < 1e-12)
  if (!length(keep)) return(0.5)
  
  out <- stats::median(grid[keep])
  if (!is.finite(out)) 0.5 else out
}

# ======================================================================
# 3. optional SMOTE and ENN, used if installed
# ======================================================================

maybe_smote_enn <- function(dtrain, yvar) {
  if (!requireNamespace("UBL", quietly = TRUE) || !requireNamespace("NoiseFiltersR", quietly = TRUE)) {
    return(dtrain)
  }
  
  dtrain[[yvar]] <- factor(dtrain[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  sm_formula <- stats::as.formula(paste(yvar, "~ ."))
  
  d_sm <- try(
    UBL::SmoteClassif(sm_formula, dtrain, C.perc = "balance", k = 5),
    silent = TRUE
  )
  if (!inherits(d_sm, "try-error") && nrow(d_sm) > 0) {
    dtrain <- d_sm
    dtrain[[yvar]] <- factor(dtrain[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  }
  
  enn_res <- try(
    NoiseFiltersR::ENN(sm_formula, dtrain, k = 3),
    silent = TRUE
  )
  if (!inherits(enn_res, "try-error") &&
      !is.null(enn_res$cleanData) &&
      nrow(enn_res$cleanData) > 0) {
    dtrain <- enn_res$cleanData
    dtrain[[yvar]] <- factor(dtrain[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  }
  
  dtrain
}

# ======================================================================
# 4. fit holdout model, calibrate from OOF, score external probabilities
# tuneLength = 5 was used for all algorithms
# ======================================================================

fit_holdout_for_importance <- function(df_train, df_test, yvar, features, algo_name, set_label) {
  features <- features[!is.na(features)]
  model_cols <- unique(c(yvar, features))
  model_cols <- intersect(model_cols, names(df_train))
  model_cols <- intersect(model_cols, names(df_test))
  
  dtrain <- df_train %>% dplyr::select(dplyr::all_of(model_cols))
  dext   <- df_test  %>% dplyr::select(dplyr::all_of(model_cols))
  
  if (!nrow(dtrain) || !nrow(dext)) return(NULL)
  
  if (nrow(dtrain) > max_n_train) {
    set.seed(seed_from_name(paste(algo_name, set_label, "subsample", sep = "|")))
    dtrain <- dplyr::sample_n(dtrain, max_n_train)
  }
  
  dtrain[[yvar]] <- factor(dtrain[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  dext[[yvar]]   <- factor(dext[[yvar]],   levels = c(cfg$neg_label, cfg$pos_label))
  
  dtrain <- maybe_smote_enn(dtrain, yvar)
  
  method <- method_from_algo(algo_name)
  tl     <- tune_len_from_algo(algo_name)
  
  y <- dtrain[[yvar]]
  w <- compute_case_weights(y)
  
  set.seed(seed_from_name(paste(algo_name, set_label, "fit", sep = "|")))
  
  fit <- switch(
    algo_name,
    "LR" = caret::train(
      stats::as.formula(paste(yvar, "~ .")),
      data       = dtrain,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tl,
      weights    = w,
      family     = "binomial"
    ),
    "RF" = caret::train(
      stats::as.formula(paste(yvar, "~ .")),
      data       = dtrain,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tl,
      weights    = w,
      ntree      = 500
    ),
    caret::train(
      stats::as.formula(paste(yvar, "~ .")),
      data       = dtrain,
      method     = method,
      metric     = "ROC",
      trControl  = ctrl_auc,
      tuneLength = tl,
      weights    = w
    )
  )
  
  oof <- extract_oof_best(fit)
  pos_col <- cfg$pos_label
  if (!(pos_col %in% names(oof))) stop("Positive class column was not found in out of fold predictions.")
  
  p_oof_raw <- oof[[pos_col]]
  p_oof_raw[is.na(p_oof_raw)] <- 0.5
  
  cal_fun <- fit_platt_glm_from_oof(oof$obs, p_oof_raw, cfg$pos_label)
  p_cal   <- cal_fun(p_oof_raw)
  p_cal[is.na(p_cal)] <- 0.5
  
  thr <- best_threshold_mcc_med(oof$obs, p_cal)
  if (!is.finite(thr)) thr <- 0.5
  
  p_ex_mat <- predict(fit, newdata = dext, type = "prob", na.action = na.pass)
  p_ex_raw <- if (is.vector(p_ex_mat)) as.numeric(p_ex_mat) else p_ex_mat[, cfg$pos_label]
  p_ex_raw[is.na(p_ex_raw)] <- 0.5
  p_ex <- cal_fun(p_ex_raw)
  p_ex[is.na(p_ex)] <- 0.5
  
  list(
    algo      = algo_name,
    set       = set_label,
    fit       = fit,
    threshold = thr,
    cal_fun   = cal_fun,
    oof       = data.frame(obs = oof$obs, p = p_cal),
    data_ext  = dext,
    p_ext     = p_ex,
    features  = features,
    yvar      = yvar
  )
}

# ======================================================================
# 5. permutation ΔMCC, single panel
# ======================================================================

perm_drop_mcc_simple <- function(hold) {
  dat <- hold$data_ext
  y   <- dat[[hold$yvar]]
  pos <- cfg$pos_label
  neg <- cfg$neg_label
  
  p0 <- hold$p_ext
  p0[is.na(p0)] <- 0.5
  thr0 <- hold$threshold
  if (!is.finite(thr0)) thr0 <- 0.5
  
  pred0 <- ifelse(p0 >= thr0, pos, neg)
  m0 <- mcc_from_preds(y, pred0, pos, neg)
  
  feats <- intersect(hold$features, setdiff(names(dat), hold$yvar))
  if (!length(feats)) return(tibble::tibble())
  
  purrr::map_dfr(feats, function(f) {
    dperm <- dat
    set.seed(seed_from_name(paste(f, hold$algo, hold$set, sep = "|"), offset = 999L))
    dperm[[f]] <- sample(dperm[[f]])
    
    p_mat <- predict(hold$fit, newdata = dperm, type = "prob", na.action = na.pass)
    p_raw <- if (is.vector(p_mat)) as.numeric(p_mat) else p_mat[, cfg$pos_label]
    p_raw[is.na(p_raw)] <- 0.5
    p_cal <- hold$cal_fun(p_raw)
    p_cal[is.na(p_cal)] <- 0.5
    
    pred1 <- ifelse(p_cal >= thr0, pos, neg)
    m1 <- mcc_from_preds(y, pred1, pos, neg)
    
    tibble::tibble(Feature = f, Drop_MCC = m0 - m1)
  })
}

# ======================================================================
# 6. cascade band tuning on triage OOF, NA safe
# ======================================================================

tune_band_from_triage_oof <- function(oof_df, thr_tri,
                                      grid_tl = cascade_cfg$grid_tl,
                                      grid_th = cascade_cfg$grid_th,
                                      max_defer = cascade_cfg$max_defer) {
  stopifnot(all(c("obs", "p") %in% names(oof_df)))
  
  pos <- cfg$pos_label
  neg <- cfg$neg_label
  
  p <- oof_df$p
  p[is.na(p)] <- 0.5
  
  thr_tri <- if (is.finite(thr_tri)) thr_tri else 0.5
  
  best <- NULL
  
  for (tl in grid_tl) for (th in grid_th) if (tl < th) {
    defer <- (p > tl & p < th)
    defer[is.na(defer)] <- FALSE
    
    dr <- mean(defer)
    if (!is.finite(dr)) next
    if (dr > max_defer) next
    
    keep <- !defer
    if (!any(keep)) next
    
    obs_keep <- factor(oof_df$obs[keep], levels = c(neg, pos))
    p_keep   <- p[keep]
    
    pred_keep <- ifelse(p_keep >= thr_tri, pos, neg)
    score <- mcc_from_preds(obs_keep, pred_keep, pos, neg)
    if (!is.finite(score)) next
    
    cand <- list(
      t_low  = tl,
      t_high = th,
      mcc    = score,
      defer  = dr,
      n_keep = sum(keep)
    )
    
    if (is.null(best) ||
        cand$mcc > best$mcc + 1e-12 ||
        (abs(cand$mcc - best$mcc) < 1e-12 &&
         (cand$defer < best$defer - 1e-12 ||
          (abs(cand$defer - best$defer) < 1e-12 && cand$n_keep > best$n_keep)))) {
      best <- cand
    }
  }
  
  if (is.null(best)) {
    list(t_low = 0.33, t_high = 0.66, mcc = NA_real_, defer = NA_real_, n_keep = NA_integer_)
  } else {
    best
  }
}

# ======================================================================
# 7. permutation ΔMCC, cascade
# ======================================================================

mcc_cascade_from_probs <- function(y, p_tri, thr_tri, p_full, thr_full, band) {
  pos <- cfg$pos_label
  neg <- cfg$neg_label
  
  p_tri[is.na(p_tri)] <- 0.5
  p_full[is.na(p_full)] <- 0.5
  
  thr_tri  <- if (is.finite(thr_tri))  thr_tri  else 0.5
  thr_full <- if (is.finite(thr_full)) thr_full else 0.5
  
  defer <- (p_tri > band$t_low & p_tri < band$t_high)
  defer[is.na(defer)] <- FALSE
  
  pred <- character(length(p_tri))
  pred[!defer] <- ifelse(p_tri[!defer] >= thr_tri, pos, neg)
  pred[defer]  <- ifelse(p_full[defer] >= thr_full, pos, neg)
  
  mcc_from_preds(y, pred, pos, neg)
}

perm_drop_mcc_cascade <- function(hT, hC, features_union, band) {
  dat_tri  <- hT$data_ext
  dat_full <- hC$data_ext
  stopifnot(nrow(dat_tri) == nrow(dat_full))
  
  y <- dat_full[[hC$yvar]]
  
  m0 <- mcc_cascade_from_probs(
    y        = y,
    p_tri    = hT$p_ext,
    thr_tri  = hT$threshold,
    p_full   = hC$p_ext,
    thr_full = hC$threshold,
    band     = band
  )
  
  feats <- intersect(features_union, union(hT$features, hC$features))
  feats <- intersect(feats, union(names(dat_tri), names(dat_full)))
  feats <- setdiff(feats, c(hT$yvar, hC$yvar))
  if (!length(feats)) return(tibble::tibble())
  
  purrr::map_dfr(feats, function(f) {
    dtri  <- dat_tri
    dfull <- dat_full
    set.seed(seed_from_name(f, offset = 20251010L))
    idx <- sample.int(nrow(dtri), nrow(dtri))
    
    if (f %in% names(dtri))  dtri[[f]]  <- dtri[[f]][idx]
    if (f %in% names(dfull)) dfull[[f]] <- dfull[[f]][idx]
    
    p_tri_mat <- predict(hT$fit, newdata = dtri, type = "prob", na.action = na.pass)
    p_tri_raw <- if (is.vector(p_tri_mat)) as.numeric(p_tri_mat) else p_tri_mat[, cfg$pos_label]
    p_tri_raw[is.na(p_tri_raw)] <- 0.5
    p_tri <- hT$cal_fun(p_tri_raw)
    p_tri[is.na(p_tri)] <- 0.5
    
    p_full_mat <- predict(hC$fit, newdata = dfull, type = "prob", na.action = na.pass)
    p_full_raw <- if (is.vector(p_full_mat)) as.numeric(p_full_mat) else p_full_mat[, cfg$pos_label]
    p_full_raw[is.na(p_full_raw)] <- 0.5
    p_full <- hC$cal_fun(p_full_raw)
    p_full[is.na(p_full)] <- 0.5
    
    m1 <- mcc_cascade_from_probs(
      y        = y,
      p_tri    = p_tri,
      thr_tri  = hT$threshold,
      p_full   = p_full,
      thr_full = hC$threshold,
      band     = band
    )
    
    tibble::tibble(Feature = f, Drop_MCC = m0 - m1)
  })
}

# ======================================================================
# 8. fit all models, compute importances
# ======================================================================

holds_triage <- setNames(vector("list", length(algo_order)), algo_order)
holds_full   <- setNames(vector("list", length(algo_order)), algo_order)

for (algo in algo_order) {
  message("Fitting ", algo, " Triage for ΔMCC")
  holds_triage[[algo]] <- fit_holdout_for_importance(df_train, df_test, yvar, feat_triage, algo, set_label = "Triage")
  
  message("Fitting ", algo, " Complete for ΔMCC")
  holds_full[[algo]] <- fit_holdout_for_importance(df_train, df_test, yvar, feat_full, algo, set_label = "Complete")
}

imp_tri <- purrr::map_dfr(algo_order, function(algo) {
  h <- holds_triage[[algo]]
  if (is.null(h)) return(tibble::tibble())
  perm_drop_mcc_simple(h) %>%
    dplyr::mutate(Set = "Triage", Algorithm = factor(algo, levels = algo_order))
})

imp_full <- purrr::map_dfr(algo_order, function(algo) {
  h <- holds_full[[algo]]
  if (is.null(h)) return(tibble::tibble())
  perm_drop_mcc_simple(h) %>%
    dplyr::mutate(Set = "Complete", Algorithm = factor(algo, levels = algo_order))
})

feat_union <- union(feat_full, feat_triage)

imp_cas <- purrr::map_dfr(algo_order, function(algo) {
  hT <- holds_triage[[algo]]
  hC <- holds_full[[algo]]
  if (is.null(hT) || is.null(hC) || is.null(hT$oof) || !nrow(hT$oof)) return(tibble::tibble())
  
  tuned <- tune_band_from_triage_oof(oof_df = hT$oof, thr_tri = hT$threshold)
  band <- list(t_low = tuned$t_low, t_high = tuned$t_high)
  
  perm_drop_mcc_cascade(hT, hC, features_union = feat_union, band = band) %>%
    dplyr::mutate(Set = "Cascade", Algorithm = factor(algo, levels = algo_order))
})

imp_all_raw <- dplyr::bind_rows(imp_full, imp_tri, imp_cas)
if (!nrow(imp_all_raw)) stop("No importance results were produced.")

# ======================================================================
# 9. pretty feature labels and ordering
# ======================================================================

pretty_feature_label <- function(x) {
  base <- x
  base <- gsub("_raw_rank$", "", base)
  base <- gsub("_rank$", "", base)
  base <- gsub("_raw$", "", base)
  
  map_special <- c(
    "age"        = "Age",
    "sex"        = "Sex (M=1)",
    "spo2"       = "SpO\u2082",
    "heart_rate" = "Heart rate",
    "resp_rate"  = "Respiratory rate",
    "temp"       = "Temperature",
    "sbp"        = "Systolic BP",
    "dbp"        = "Diastolic BP",
    "map"        = "MAP",
    "wbc"        = "White blood cells",
    "crp"        = "CRP",
    "urea"       = "Urea",
    "inr"        = "INR",
    "alk_phos"   = "Alkaline phosphatase",
    "ph_venous"  = "Venous pH"
  )
  
  lab <- unname(map_special[match(base, names(map_special))])
  idx <- is.na(lab)
  if (any(idx)) {
    lab[idx] <- base[idx]
    lab[idx] <- gsub("_", " ", lab[idx], fixed = TRUE)
    lab[idx] <- tools::toTitleCase(lab[idx])
  }
  ifelse(is.na(lab), base, lab)
}

imp_all <- imp_all_raw %>%
  dplyr::mutate(
    Feature_label = pretty_feature_label(.data$Feature),
    Set = factor(.data$Set, levels = c("Complete", "Triage", "Cascade"))
  )

tri_labels_vec   <- unique(pretty_feature_label(feat_triage))
other_labels_vec <- setdiff(unique(pretty_feature_label(feat_full)), tri_labels_vec)

tri_rev   <- sort(tri_labels_vec, decreasing = TRUE)
other_rev <- sort(other_labels_vec, decreasing = TRUE)

levels_full_like <- c(other_rev, tri_rev)
levels_tri       <- rev(tri_rev)
tri_labels       <- tri_labels_vec

# ======================================================================
# 10. grid, heatmap
# ======================================================================

grid_full <- tidyr::expand_grid(
  Set = factor("Complete", levels = c("Complete", "Triage", "Cascade")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_full_like, levels = levels_full_like)
)

grid_tri <- tidyr::expand_grid(
  Set = factor("Triage", levels = c("Complete", "Triage", "Cascade")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_tri, levels = levels_tri)
)

grid_cas <- tidyr::expand_grid(
  Set = factor("Cascade", levels = c("Complete", "Triage", "Cascade")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_full_like, levels = levels_full_like)
)

df_full <- dplyr::left_join(
  grid_full,
  imp_all %>% dplyr::filter(.data$Set == "Complete") %>% dplyr::select(.data$Set, .data$Algorithm, .data$Feature_label, .data$Drop_MCC),
  by = c("Set", "Algorithm", "Feature_label")
)

df_tri <- dplyr::left_join(
  grid_tri,
  imp_all %>% dplyr::filter(.data$Set == "Triage") %>% dplyr::select(.data$Set, .data$Algorithm, .data$Feature_label, .data$Drop_MCC),
  by = c("Set", "Algorithm", "Feature_label")
)

df_cas <- dplyr::left_join(
  grid_cas,
  imp_all %>% dplyr::filter(.data$Set == "Cascade") %>% dplyr::select(.data$Set, .data$Algorithm, .data$Feature_label, .data$Drop_MCC),
  by = c("Set", "Algorithm", "Feature_label")
)

df_plot <- dplyr::bind_rows(df_full, df_tri, df_cas)

lim <- stats::quantile(abs(df_plot$Drop_MCC), 0.95, na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- max(abs(df_plot$Drop_MCC), na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- 1e-6

super_y_labels <- function(labs, tri = tri_labels, mark = "\u2020") {
  lapply(labs, function(lab) if (lab %in% tri) bquote(.(lab)^.(mark)) else bquote(.(lab)))
}

p_heat <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$Algorithm, y = .data$Feature_label, fill = .data$Drop_MCC)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::scale_fill_gradient2(
    low = "#D55E00",
    mid = "white",
    high = "#0072B2",
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish,
    name = expression(Delta * "MCC"),
    na.value = "white"
  ) +
  ggplot2::facet_grid(rows = vars(Set), scales = "free_y", space = "free_y") +
  ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
  ggplot2::scale_y_discrete(labels = super_y_labels) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 8) +
  ggplot2::theme(
    panel.grid         = ggplot2::element_blank(),
    panel.border       = ggplot2::element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.text.x.top    = ggplot2::element_text(size = 7, face = "plain"),
    axis.text.y        = ggplot2::element_text(size = 7),
    strip.text.y.right = ggplot2::element_text(size = 8, face = "plain"),
    legend.title       = ggplot2::element_text(size = 7, face = "plain"),
    legend.text        = ggplot2::element_text(size = 6)
  ) +
  ggplot2::guides(fill = ggplot2::guide_colorbar(
    barheight = grid::unit(18, "mm"),
    barwidth  = grid::unit(2.5, "mm"),
    ticks = TRUE
  ))

print(p_heat)
