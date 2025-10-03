# ============================================================
# FULL PIPELINE + TABLE FACTORY, group-filtered + printed tables
# - Filters to Group %in% {CTRL_noCOVID, COVID}
# - Trains Full vs Triage models with outer CV + holdouts
# - Builds Tables 1–10 + Supplement, prints them
# - Stable calibration (ridge Platt), robust subgroup CIs, clean Table 1/2
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(yardstick)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))

set.seed(123)

# ---------------- config ----------------
cfg <- list(
  # I/O
  file_path = "biomarkers_acuteCOVID_meta.xlsx",
  sheet     = "meta",
  outcome   = "Hospital_ID",  # original outcome; derived to factor {Yes, No}
  
  # modeling
  pos_label = "Yes",
  neg_label = "No",
  fast_mode = TRUE,
  cv_k      = 10,
  cv_R      = 1,
  inner_k   = 5,
  tune_len  = 10,
  seed_cv   = 444,
  
  # bootstrap
  n_boot_holdout = 1000,
  n_boot_subgrp  = 1000,
  
  # misc
  print_fold_diag = TRUE
)

if (cfg$fast_mode) {
  cfg$cv_k           <- 3
  cfg$inner_k        <- 3
  cfg$tune_len       <- 3
  cfg$n_boot_holdout <- 250
  cfg$n_boot_subgrp  <- 300
}

# exact predictors
feat_full <- c(
  "Diagnosis","severity_admission","Age","Gender",
  "SpO2_admission",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc"
)
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- data (FILTER Group ∈ {CTRL_noCOVID, COVID}) ----------------
load_data <- function(path, sheet, outcome, keep) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s\nPlease set cfg$file_path to your .xlsx", path))
  }
  df <- readxl::read_excel(path, sheet = sheet)
  df <- as.data.frame(df)
  
  # harmonize column name then filter to desired groups
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found, filtering cannot be applied")
  names(df)[names(df) == grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  # simple row id (for diagnostics)
  df$.rid <- seq_len(nrow(df))
  
  # cast likely categoricals
  guess_factor <- intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))
  for (v in guess_factor) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  # normalize outcome to factor {Yes, No} with Yes as FIRST level (caret's twoClassSummary expects eventLevel="first")
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' not in data", outcome))
  y <- df[[outcome]]
  to_YN <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) x <- ifelse(x, "Yes", "No")
    if (is.numeric(x)) {
      if (!all(x %in% c(0, 1), na.rm = TRUE)) stop("Outcome numeric but not 0/1")
      out <- ifelse(x == 1, cfg$pos_label, cfg$neg_label)
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    if (is.character(x)) {
      s <- trimws(tolower(x))
      map_yes <- c("1","yes","y","true","pos","positive")
      map_no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% map_yes, cfg$pos_label,
                    ifelse(s %in% map_no,  cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome must be Yes/No, Y/N, TRUE/FALSE, or 0/1 after trimming")
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    stop("Outcome type unsupported for coercion")
  }
  df[[outcome]] <- to_YN(y)
  
  # keep only desired columns
  keepx <- unique(c(keep, outcome, "data_split", "group", ".rid"))
  keepx <- intersect(keepx, names(df))
  df <- df[, keepx, drop = FALSE]
  df
}

# ---------------- recipe ----------------
make_recipe <- function(dat, yvar) {
  rec <- recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c("data_split", ".rid"), names(dat))
  if (length(ign)) rec <- rec |> update_role(all_of(ign), new_role = "ignore")
  
  rec <- rec |>
    step_impute_median(all_numeric_predictors()) |>
    step_impute_mode(all_nominal_predictors()) |>
    step_novel(all_nominal_predictors()) |>
    step_other(all_nominal_predictors(), threshold = 0.01) |>
    step_dummy(all_nominal_predictors()) |>
    step_zv(all_predictors())
  
  # Yeo–Johnson if available
  if (exists("step_YeoJohnson", where = asNamespace("recipes"), inherits = FALSE)) {
    rec <- rec |> step_YeoJohnson(all_numeric_predictors())
  } else if (exists("step_yeojohnson", where = asNamespace("recipes"), inherits = FALSE)) {
    rec <- rec |> step_yeojohnson(all_numeric_predictors())
  }
  
  # explicit z-score to match Table 2
  rec |> step_normalize(all_numeric_predictors())
}

# ---------------- threshold selection ----------------
mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(NA_real_)
  num/den
}

best_threshold_mcc <- function(obs, p_pos, pos = "Yes", neg = "No",
                               grid = seq(0.01, 0.99, by = 0.001)) {
  y <- factor(obs, levels = c(neg, pos))
  best_t <- 0.5
  best_m <- -Inf
  for (t in grid) {
    pred <- factor(ifelse(p_pos >= t, pos, neg), levels = c(neg, pos))
    tab <- table(y, pred)
    TP <- as_num(tab[pos, pos] %||% 0); TN <- as_num(tab[neg, neg] %||% 0)
    FP <- as_num(tab[neg, pos] %||% 0); FN <- as_num(tab[pos, neg] %||% 0)
    m  <- mcc_from_counts(TP, FP, FN, TN)
    if (is.finite(m) && m > best_m) { best_m <- m; best_t <- t }
  }
  list(t = best_t, mcc = best_m)
}

# ---------------- binary metrics ----------------
binary_auc <- function(obs, p_pos, pos_label = "Yes") {
  y <- as.integer(obs == pos_label)
  if (length(unique(y)) < 2) return(NA_real_)
  suppressMessages(suppressWarnings(
    as.numeric(pROC::auc(y, p_pos, quiet = TRUE)) # let pROC auto-detect direction
  ))
}

compute_metrics_binary <- function(obs, pred, p_pos, pos_label = "Yes", neg_label = "No") {
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
  c(MCC = mcc, F1 = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec, AUC = aucv)
}

# ---------------- model registry ----------------
model_specs <- function(tune_len, algos = c("LR","RF","SVM","k-NN","C4.5")) {
  avail <- names(caret::getModelInfo())
  all <- list(
    "LR"   = list(method = "glmnet", tuneLength = tune_len,
                  grid = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 0, length.out = tune_len))),
    "RF"   = list(method = "rf", tuneLength = tune_len),
    "SVM"  = list(method = "svmRadial", tuneLength = tune_len),
    "k-NN" = list(method = "kknn", tuneLength = tune_len),
    "C4.5" = list(method = "J48", tuneLength = tune_len) # requires RWeka installed
  )
  # keep only requested and available
  all[names(all) %in% algos & vapply(all, function(s) s$method %in% avail, logical(1))]
}

# ---------------- CV utilities ----------------
choose_k_for_task <- function(y, k_desired) {
  if (is.null(y) || !length(y) || anyNA(y)) return(k_desired)
  tab <- table(y); max_k <- max(2, min(tab))
  kk  <- min(as.integer(k_desired), as.integer(max_k))
  kk
}

build_cv_splits <- function(y, R, k_desired, seed = 999) {
  set.seed(seed)
  splits <- list()
  for (r in seq_len(R)) {
    kk <- choose_k_for_task(y, k_desired)
    idx <- caret::createFolds(y, k = kk, list = TRUE, returnTrain = FALSE)
    for (i in seq_along(idx)) {
      splits[[paste0("r", r, "_f", i)]] <- list(va = idx[[i]])
    }
  }
  splits
}

# ---------------- training wrappers ----------------
train_inner <- function(dat_train, yvar, spec, inner_k) {
  rec <- make_recipe(dat_train, yvar)
  tr_ctrl <- caret::trainControl(
    method = "cv", number = inner_k,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", allowParallel = TRUE
  )
  if (!is.null(spec$grid)) {
    fit <- caret::train(
      rec, data = dat_train, method = spec$method,
      trControl = tr_ctrl, metric = "ROC", tuneGrid = spec$grid
    )
  } else {
    fit <- caret::train(
      rec, data = dat_train, method = spec$method,
      trControl = tr_ctrl, metric = "ROC", tuneLength = spec$tuneLength
    )
  }
  fit
}

# ---------------- single fold eval ----------------
eval_fold_binary <- function(dat_tr, dat_va, spec, yvar, pos_label, algo_name, fold_tag, inner_k) {
  fit <- train_inner(dat_tr, yvar, spec, inner_k)
  p_va <- as.numeric(predict(fit, newdata = dat_va, type = "prob")[, pos_label])
  thr  <- best_threshold_mcc(dat_va[[yvar]], p_va, pos = pos_label, neg = setdiff(levels(dat_va[[yvar]]), pos_label)[1])
  pred <- factor(ifelse(p_va >= thr$t, pos_label, setdiff(levels(dat_va[[yvar]]), pos_label)[1]),
                 levels = levels(dat_va[[yvar]]))
  mets <- compute_metrics_binary(dat_va[[yvar]], pred, p_va, pos_label = pos_label)
  cat(sprintf("[Fold %s] %s | %s | tuned.t=%.3f | pos=%s, neg=%s\n",
              fold_tag, yvar, algo_name, thr$t, pos_label, setdiff(levels(dat_va[[yvar]]), pos_label)[1]))
  tibble::tibble(
    Fold = fold_tag,
    Algorithm = algo_name,
    Metric = names(mets),
    Value  = as.numeric(mets)
  )
}

# ---------------- unified CV runner ----------------
run_cv <- function(df, yvar, features, algos, type = "binary",
                   pos_label = "Yes", inner_k = cfg$inner_k,
                   k_desired = cfg$cv_k, R = cfg$cv_R) {
  dat <- df[, unique(c(features, yvar, "data_split")), drop = FALSE]
  dat <- dat[complete.cases(dat[[yvar]]), , drop = FALSE]
  y <- dat[[yvar]]
  splits <- build_cv_splits(y, R, k_desired, seed = cfg$seed_cv)
  specs  <- model_specs(cfg$tune_len, algos = algos)
  out <- list()
  for (algo in names(specs)) {
    spec <- specs[[algo]]
    for (nm in names(splits)) {
      idx_va <- splits[[nm]]$va
      dat_va <- dat[idx_va, , drop = FALSE]
      dat_tr <- dat[-idx_va, , drop = FALSE]
      out[[paste0(algo, "_", nm)]] <- eval_fold_binary(dat_tr, dat_va, spec, yvar, pos_label, algo, nm, inner_k)
    }
  }
  dplyr::bind_rows(out)
}

# ---------------- CIs & helpers ----------------
binom_ci_wilson <- function(x, n, conf.level = 0.95) {
  if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf.level)/2)
  p <- x/n
  denom <- 1 + z^2/n
  center <- (p + z^2/(2*n)) / denom
  half   <- z * sqrt((p*(1-p) + z^2/(4*n)) / n) / denom
  c(center - half, center + half)
}

auc_delong_ci <- function(obs, p_pos, conf.level = 0.95) {
  y <- as.integer(obs == cfg$pos_label)
  if (length(unique(y)) < 2) return(c(NA_real_, NA_real_))
  suppressMessages(suppressWarnings({
    r <- pROC::roc(y, p_pos, quiet = TRUE)
    as.numeric(pROC::ci.auc(r, conf.level = conf.level))
  }))
}

mk_cm <- function(obs, pred, pos = cfg$pos_label) {
  lev <- levels(obs); neg <- setdiff(lev, pos)[1]
  tab <- table(obs, pred)
  TP <- as_num(tab[pos, pos] %||% 0); TN <- as_num(tab[neg, neg] %||% 0)
  FP <- as_num(tab[neg, pos] %||% 0); FN <- as_num(tab[pos, neg] %||% 0)
  n  <- sum(tab)
  sens <- if ((TP+FN)==0) NA_real_ else TP/(TP+FN)
  spec <- if ((TN+FP)==0) NA_real_ else TN/(TN+FP)
  prec <- if ((TP+FP)==0) NA_real_ else TP/(TP+FP)
  acc  <- (TP+TN)/n
  list(TP=TP, FP=FP, FN=FN, TN=TN, Sensitivity=sens, Specificity=spec, Precision=prec, Accuracy=acc, N=n)
}

# ---- fit on train, predict on holdout (robust to missing 'group') ----
run_holdout <- function(df, yvar, features, algo_name,
                        split_train = "train", split_test = c("test","external")) {
  present_splits <- intersect(split_test, unique(as.character(df$data_split)))
  if (!length(present_splits)) return(NULL)
  specs <- model_specs(cfg$tune_len, algos = algo_name)
  if (!length(specs)) return(NULL)
  spec  <- specs[[algo_name]]
  base_cols <- unique(c(features, yvar, intersect("group", names(df))))
  base_cols <- intersect(base_cols, names(df))
  dtrain <- df[df$data_split == split_train, base_cols, drop = FALSE]
  if (!nrow(dtrain)) return(NULL)
  fit <- train_inner(dtrain, yvar, spec, cfg$inner_k)
  out <- list()
  for (sp in present_splits) {
    dtest <- df[df$data_split == sp, base_cols, drop = FALSE]
    if (!nrow(dtest)) next
    p_tr  <- as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label])
    thr   <- best_threshold_mcc(dtrain[[yvar]], p_tr, pos = cfg$pos_label, neg = cfg$neg_label)
    p_te  <- as.numeric(predict(fit, newdata = dtest, type = "prob")[, cfg$pos_label])
    pred  <- factor(ifelse(p_te >= thr$t, cfg$pos_label, cfg$neg_label),
                    levels = c(cfg$neg_label, cfg$pos_label))
    out[[sp]] <- list(
      obs = factor(dtest[[yvar]], levels = c(cfg$neg_label, cfg$pos_label)),
      p   = p_te,
      pred = pred,
      threshold = thr$t,
      model = fit,
      holdout_df = dtest
    )
  }
  out
}

# ---- calibration metrics (ridge Platt) ----
safe_logit <- function(p, eps = 1e-4) qlogis(pmin(pmax(p, eps), 1 - eps))

calibration_metrics <- function(obs, p, bins = 10, conf.level = 0.95,
                                n_boot = cfg$n_boot_holdout, lambda = 0.1) {
  y <- as.integer(obs == cfg$pos_label)
  if (length(unique(y)) < 2 || all(!is.finite(p))) {
    return(tibble::tibble(
      Metric = character(), Point = double(), CI_low = double(), CI_high = double(),
      Phase = character()
    ))
  }
  p <- pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)
  logitp <- safe_logit(p)
  
  # PRE: raw metrics
  brier <- mean((y - p)^2)
  brks <- stats::quantile(p, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE, type = 1)
  brks[1] <- 0; brks[length(brks)] <- 1
  bin_id <- cut(p, breaks = unique(brks), include.lowest = TRUE, labels = FALSE)
  ece <- {
    dfb <- tibble::tibble(bin = bin_id, y = y, p = p) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(n = dplyr::n(), obs = mean(y), pred = mean(p), .groups = "drop") |>
      dplyr::mutate(w = n / sum(n))
    sum(dfb$w * abs(dfb$obs - dfb$pred))
  }
  
  # ridge slope-only Platt with robust coef extraction
  fit_slope <- try(glmnet::glmnet(x = matrix(logitp, ncol = 1), y = y,
                                  family = "binomial", alpha = 0,
                                  lambda = lambda, standardize = TRUE), silent = TRUE)
  if (inherits(fit_slope, "try-error")) {
    a_slope <- NA_real_
    b1 <- NA_real_
  } else {
    cf <- as.numeric(glmnet::coef(fit_slope))  # (Intercept), x
    a_slope <- cf[1]
    b1 <- cf[2]
  }
  
  # intercept at fixed slope 1 via offset (bias-reduced if brglm2 is available)
  fit_int <- suppressWarnings(try({
    if (requireNamespace("brglm2", quietly = TRUE)) {
      brglm2::brglm(y ~ 1 + offset(logitp), family = binomial("logit"), type = "AS_median")
    } else {
      glm(y ~ 1 + offset(logitp), family = binomial())
    }
  }, silent = TRUE))
  a0 <- if (inherits(fit_int, "try-error")) NA_real_ else as.numeric(coef(fit_int)[1])
  
  # bootstrap PRE metrics
  boot_once <- function(idx) {
    yy <- y[idx]; pp <- p[idx]; lp <- safe_logit(pp)
    b <- mean((yy - pp)^2)
    br <- stats::quantile(pp, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE, type = 1)
    br[1] <- 0; br[length(br)] <- 1
    bi <- cut(pp, breaks = unique(br), include.lowest = TRUE, labels = FALSE)
    e <- {
      dfb <- tibble::tibble(bin = bi, y = yy, p = pp) |>
        dplyr::group_by(bin) |>
        dplyr::summarise(n = dplyr::n(), obs = mean(y), pred = mean(p), .groups = "drop") |>
        dplyr::mutate(w = n / sum(n))
      sum(dfb$w * abs(dfb$obs - dfb$pred))
    }
    fs <- try(glmnet::glmnet(x = matrix(lp, ncol = 1), y = yy, family = "binomial",
                             alpha = 0, lambda = lambda, standardize = TRUE), silent = TRUE)
    if (inherits(fs, "try-error")) {
      bb1 <- NA_real_
    } else {
      bb1 <- as.numeric(glmnet::coef(fs))[2]
    }
    fi <- suppressWarnings(try(glm(yy ~ 1 + offset(lp), family = binomial()), silent = TRUE))
    aa0 <- if (inherits(fi, "try-error")) NA_real_ else as.numeric(coef(fi))
    c(b, e, aa0, bb1)
  }
  nB <- min(n_boot, max(200, ceiling(length(y) * 2)))
  set.seed(1L)
  idx_list <- replicate(nB, sample.int(length(y), replace = TRUE), simplify = FALSE)
  boots <- purrr::map_dfr(idx_list, ~{
    res <- boot_once(.x)
    tibble::tibble(Brier = res[1], ECE = res[2], Intercept = res[3], Slope = res[4])
  })
  ci <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA, NA))
    as.numeric(stats::quantile(v, probs = c((1 - conf.level)/2, 1 - (1 - conf.level)/2),
                               na.rm = TRUE, names = FALSE))
  }
  out_pre <- tibble::tibble(
    Metric = c("Brier","Intercept","Slope","ECE"),
    Point  = c(brier, a0, b1, ece),
    CI_low = c(ci(boots$Brier)[1], ci(boots$Intercept)[1], ci(boots$Slope)[1], ci(boots$ECE)[1]),
    CI_high= c(ci(boots$Brier)[2], ci(boots$Intercept)[2], ci(boots$Slope)[2], ci(boots$ECE)[2]),
    Phase  = "Pre"
  )
  
  # POST: ridge Platt calibrated
  p_cal <- if (is.finite(a_slope) && is.finite(b1)) plogis(a_slope + b1 * logitp) else p
  p_cal <- pmin(pmax(p_cal, 1e-6), 1 - 1e-6)
  brier2 <- mean((y - p_cal)^2)
  brks2 <- stats::quantile(p_cal, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE, type = 1)
  brks2[1] <- 0; brks2[length(brks2)] <- 1
  bin_id2 <- cut(p_cal, breaks = unique(brks2), include.lowest = TRUE, labels = FALSE)
  ece2 <- {
    dfb2 <- tibble::tibble(bin = bin_id2, y = y, p = p_cal) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(n = dplyr::n(), obs = mean(y), pred = mean(p), .groups = "drop") |>
      dplyr::mutate(w = n / sum(n))
    sum(dfb2$w * abs(dfb2$obs - dfb2$pred))
  }
  
  boots2 <- purrr::map_dfr(idx_list, ~{
    ii <- .x
    yy <- y[ii]; pp <- p[ii]; lp <- safe_logit(pp)
    
    # Try ridge (glmnet) first
    fs <- try(glmnet::glmnet(
      x = matrix(lp, ncol = 1), y = yy,
      family = "binomial", alpha = 0, lambda = lambda, standardize = TRUE
    ), silent = TRUE)
    
    # Fallback: vanilla glm if glmnet errors or gives degenerate coefs
    f2 <- suppressWarnings(try(
      glm(yy ~ lp, family = binomial(), control = list(maxit = 50)),
      silent = TRUE
    ))
    
    # Coefficients: prefer glmnet; fallback to glm; final fallback = identity (aa=0, bb=1)
    aa <- if (inherits(fs, "try-error")) {
      if (inherits(f2, "try-error")) 0 else as.numeric(coef(f2)[1])
    } else {
      as.numeric(glmnet::coef(fs))[1]
    }
    
    bb <- if (inherits(fs, "try-error")) {
      if (inherits(f2, "try-error")) 1 else as.numeric(coef(f2)[2])
    } else {
      as.numeric(glmnet::coef(fs))[2]
    }
    
    # Calibrated probabilities for this bootstrap sample
    pc <- plogis(aa + bb * lp)
    pc <- pmin(pmax(pc, 1e-6), 1 - 1e-6)
    
    # Metrics on calibrated probs
    b <- mean((yy - pc)^2)
    br <- stats::quantile(pc, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE, type = 1)
    br[1] <- 0; br[length(br)] <- 1
    bi <- cut(pc, breaks = unique(br), include.lowest = TRUE, labels = FALSE)
    e <- {
      dfb <- tibble::tibble(bin = bi, y = yy, p = pc) |>
        dplyr::group_by(bin) |>
        dplyr::summarise(n = dplyr::n(), obs = mean(y), pred = mean(p), .groups = "drop") |>
        dplyr::mutate(w = n / sum(n))
      sum(dfb$w * abs(dfb$obs - dfb$pred))
    }
    
    # Intercept-at-fixed-slope for reporting (same as your original)
    fi <- suppressWarnings(try(glm(yy ~ 1 + offset(safe_logit(pc)), family = binomial()), silent = TRUE))
    aa0 <- if (inherits(fi, "try-error")) NA_real_ else as.numeric(coef(fi))
    
    # Slope-only ridge on calibrated logit for reporting
    fs2 <- try(glmnet::glmnet(
      x = matrix(safe_logit(pc), ncol = 1), y = yy,
      family = "binomial", alpha = 0, lambda = lambda, standardize = TRUE
    ), silent = TRUE)
    bb1 <- if (inherits(fs2, "try-error")) NA_real_ else as.numeric(glmnet::coef(fs2))[2]
    
    tibble::tibble(Brier = b, ECE = e, Intercept = aa0, Slope = bb1)
  })
  
  
  out_post <- tibble::tibble(
    Metric = c("Brier","Intercept","Slope","ECE"),
    Point  = c(brier2, NA_real_, NA_real_, ece2),
    CI_low = c(ci(boots2$Brier)[1], ci(boots2$Intercept)[1], ci(boots2$Slope)[1], ci(boots2$ECE)[1]),
    CI_high= c(ci(boots2$Brier)[2], ci(boots2$Intercept)[2], ci(boots2$Slope)[2], ci(boots2$ECE)[2]),
    Phase  = "Post"
  )
  
  dplyr::bind_rows(out_pre, out_post)
}

# ---- decision curve numerics ----
decision_curve_table <- function(obs, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  y <- as.integer(obs == cfg$pos_label)
  N <- length(y); prev <- mean(y)
  purrr::map_dfr(thresholds, function(pt) {
    pred <- as.integer(p >= pt)
    TP <- sum(pred == 1 & y == 1)
    FP <- sum(pred == 1 & y == 0)
    NB <- TP/N - FP/N * (pt/(1-pt))
    NB_all <- prev - (1 - prev) * (pt/(1-pt))
    Precision <- if ((TP + sum(pred == 1 & y == 0)) == 0) NA_real_ else TP / sum(pred == 1)
    tibble::tibble(
      threshold = pt,
      Net_Benefit = NB,
      NB_TreatAll = NB_all,
      NB_TreatNone = 0,
      Precision = Precision,
      NNE = ifelse(is.finite(1/(Precision - pt/(1-pt))), 1/(Precision - pt/(1-pt)), NA_real_)
    )
  })
}

# ---- subgroup performance ----
subgroup_metrics <- function(df_eval, obs, p, pred,
                             vars = c("group","Gender","Diagnosis","severity_admission"),
                             pos_label = cfg$pos_label,
                             min_n = 8, conf.level = 0.95, n_boot = cfg$n_boot_subgrp) {
  if (!length(vars)) return(tibble::tibble())
  dat <- tibble::tibble(obs = obs, pred = pred, p = p) %>%
    dplyr::bind_cols(df_eval[intersect(vars, names(df_eval))])
  each_var <- function(v) {
    if (!(v %in% names(dat))) return(tibble::tibble())
    levs <- levels(factor(dat[[v]]))
    purrr::map_dfr(levs, function(lv) {
      dd <- dat[dat[[v]] == lv & !is.na(dat[[v]]), , drop = FALSE]
      n  <- nrow(dd)
      if (n < min_n || length(unique(dd$obs)) < 2) {
        return(tibble::tibble(Subgroup = v, Level = as.character(lv),
                              Metric = c("MCC","F1","Accuracy","Precision","Sensitivity","Specificity","AUC"),
                              Point = NA_real_, CI_low = NA_real_, CI_high = NA_real_, N = n))
      }
      m <- compute_metrics_binary(dd$obs, dd$pred, dd$p, pos_label = pos_label)
      boot_once <- function(idx) {
        mm <- compute_metrics_binary(dd$obs[idx], dd$pred[idx], dd$p[idx], pos_label = pos_label)
        unname(mm[c("MCC","F1","Accuracy","Precision","Sensitivity","Specificity","AUC")])
      }
      nB <- min(n_boot, max(300, n*5))
      set.seed(2L)
      idxs <- replicate(nB, sample.int(n, replace=TRUE), simplify = FALSE)
      B <- purrr::map_dfr(idxs, ~{
        vals <- boot_once(.x)
        tibble::tibble(MCC=vals[1], F1=vals[2], Accuracy=vals[3], Precision=vals[4], Sensitivity=vals[5], Specificity=vals[6], AUC=vals[7])
      })
      mkci <- function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(c(NA, NA))
        stats::quantile(x, probs = c((1-conf.level)/2, 1-(1-conf.level)/2), na.rm=TRUE, names=FALSE)
      }
      out <- tibble::tibble(
        Metric = names(m),
        Point  = as.numeric(m),
        CI_low = c(mkci(B$MCC)[1], mkci(B$F1)[1], mkci(B$Accuracy)[1], mkci(B$Precision)[1], mkci(B$Sensitivity)[1], mkci(B$Specificity)[1], mkci(B$AUC)[1]),
        CI_high= c(mkci(B$MCC)[2], mkci(B$F1)[2], mkci(B$Accuracy)[2], mkci(B$Precision)[2], mkci(B$Sensitivity)[2], mkci(B$Specificity)[2], mkci(B$AUC)[2])
      )
      out$Subgroup <- v; out$Level <- as.character(lv); out$N <- n
      out[, c("Subgroup","Level","Metric","Point","CI_low","CI_high","N")]
    })
  }
  purrr::map_dfr(vars, each_var)
}

# ---- missingness and availability ----
missingness_table <- function(df, predictors, split_col = "data_split") {
  present <- predictors[predictors %in% names(df)]
  base <- tibble::tibble(
    Predictor = present,
    Overall_N = vapply(present, function(v) length(df[[v]]), integer(1)),
    Overall_Missing_pct = vapply(present, function(v) 100*mean(is.na(df[[v]])), numeric(1))
  )
  if (split_col %in% names(df)) {
    spl <- levels(factor(df[[split_col]]))
    for (sp in spl) {
      idx <- df[[split_col]] == sp
      base[[paste0("Missing_", sp)]] <- vapply(present, function(v) 100*mean(is.na(df[[v]][idx])), numeric(1))
    }
  }
  base
}

# ---- baseline characteristics ----
baseline_table <- function(df, yvar, split_col = "data_split") {
  df2 <- df
  drop_cols <- c(".rid", yvar, "Hosp_Bin", cfg$outcome, split_col)
  drop_cols <- intersect(drop_cols, names(df2))
  df2 <- df2[, setdiff(names(df2), drop_cols), drop = FALSE]
  has_split <- split_col %in% names(df)
  
  num_cols <- names(df2)[sapply(df2, is.numeric)]
  fac_cols <- names(df2)[sapply(df2, is.factor)]
  
  num_tbl <- if (length(num_cols)) {
    agg <- function(z) {
      c(N = sum(is.finite(z)),
        Mean = mean(z, na.rm = TRUE),
        SD = stats::sd(z, na.rm = TRUE),
        Median = stats::median(z, na.rm = TRUE),
        Q1 = stats::quantile(z, 0.25, na.rm = TRUE),
        Q3 = stats::quantile(z, 0.75, na.rm = TRUE))
    }
    if (has_split) {
      sp <- split(df[[split_col]], df[[split_col]])
      out <- lapply(names(sp), function(s) {
        idx <- df[[split_col]] == s
        as.data.frame(t(vapply(df2[idx, num_cols, drop=FALSE], agg, numeric(6))))
      })
      names(out) <- names(sp)
      out <- Map(function(x, nm) { colnames(x) <- paste0(colnames(x), "_", nm); x },
                 out, names(out))
      tbl <- Reduce(function(a,b) cbind(a,b), out)
      tibble::tibble(Variable = rownames(tbl)) %>% cbind(tbl, row.names = NULL)
    } else {
      tbl <- as.data.frame(t(vapply(df2[, num_cols, drop=FALSE], agg, numeric(6))))
      tibble::tibble(Variable = rownames(tbl)) %>% cbind(tbl, row.names = NULL)
    }
  } else tibble::tibble()
  
  fac_tbl <- if (length(fac_cols)) {
    count_one <- function(var) {
      if (has_split) {
        by(df[[var]], df[[split_col]], function(x) round(100*prop.table(table(x)),1)) |>
          lapply(function(pct) tibble::tibble(Level = names(pct), Pct = as.numeric(pct))) |>
          purrr::imap(~dplyr::rename(.x, !!paste0("pct_", .y) := Pct)) |>
          purrr::reduce(dplyr::left_join, by = "Level") |>
          dplyr::mutate(Variable = var, .before = 1)
      } else {
        pct <- round(100*prop.table(table(df[[var]])),1)
        tibble::tibble(Variable = var, Level = names(pct), pct_overall = as.numeric(pct))
      }
    }
    purrr::map_dfr(fac_cols, count_one)
  } else tibble::tibble()
  
  dplyr::bind_rows(num_tbl, fac_tbl)
}

# ---- TRIPOD-AI skeleton (Table 10) ----
tripod_ai_skeleton <- function() {
  tibble::tibble(
    Item = c("Title","Abstract","Background","Objectives","Source of data","Study dates",
             "Participants","Outcome definition","Predictors","Sample size and missing data",
             "Statistical analysis methods","Model development","Internal validation",
             "External validation","Model performance","Calibration",
             "Model updating/recalibration","Clinical utility (decision curve)",
             "Subgroup analyses","Sensitivity analyses","Interpretation","Limitations",
             "Implications for practice","Funding","Ethical approval","Data/Code availability"),
    Manuscript_location = "",
    Figure_or_Table = "",
    Notes = ""
  )
}

# ---- sensitivity analyses ledger (Table 9) ----
sensitivity_ledger <- function(cv_tbl_full, cv_tbl_triage, primary_algo = "LR") {
  best_full <- cv_tbl_full %>%
    dplyr::filter(Metric %in% c("MCC","AUC")) %>%
    dplyr::group_by(Algorithm, Metric) %>%
    dplyr::summarise(Mean = mean(Value, na.rm = TRUE), .groups="drop") %>%
    dplyr::group_by(Metric) %>%
    dplyr::slice_max(order_by = Mean, n = 1, with_ties = FALSE)
  deltas <- full_join(
    cv_tbl_full %>% dplyr::filter(Metric %in% c("MCC","AUC")) %>%
      dplyr::group_by(Metric) %>% dplyr::summarise(Full = mean(Value, na.rm = TRUE), .groups="drop"),
    cv_tbl_triage %>% dplyr::filter(Metric %in% c("MCC","AUC")) %>%
      dplyr::group_by(Metric) %>% dplyr::summarise(Triage = mean(Value, na.rm = TRUE), .groups="drop"),
    by = "Metric"
  ) %>% dplyr::mutate(Delta = Triage - Full)
  delta_text <- paste0("ΔMCC=", sprintf("%.3f", deltas$Delta[deltas$Metric=="MCC"]), ", ",
                       "ΔAUC=", sprintf("%.3f", deltas$Delta[deltas$Metric=="AUC"]))
  tibble::tibble(
    Analysis = c("Triage feature set vs Full feature set",
                 sprintf("Best algorithm (%s) vs LR", best_full$Algorithm[best_full$Metric=="AUC"] %||% "LR"),
                 "Threshold at 0.5 vs MCC-tuned"),
    Change   = c("Feature set reduced to triage from full",
                 sprintf("Algorithm switched to %s from LR", best_full$Algorithm[best_full$Metric=="AUC"] %||% "LR"),
                 "Classification threshold fixed at 0.5"),
    Headline = c(delta_text,
                 "Compare CV means across algorithms for AUC and MCC",
                 "Re-compute holdout confusion and rates at t=0.5")
  )
}

# ============================================================
# RUN PIPELINE
# ============================================================

df <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
yvar <- cfg$outcome

# ensure train, test, external splits if missing
if (!("data_split" %in% names(df)) || all(is.na(df$data_split))) {
  set.seed(cfg$seed_cv)
  n <- nrow(df)
  idx <- sample(n)
  n_tr <- floor(0.7*n); n_te <- floor(0.15*n)
  df$data_split <- factor(c(rep("train", n_tr), rep("test", n_te), rep("external", n - n_tr - n_te))[order(idx)])
}

# sanity
if ("group" %in% names(df)) {
  cat("\n[Sanity] group counts after filter:\n"); print(table(df$group, useNA = "ifany"))
}
cat("\n[Sanity] data_split counts:\n"); print(table(df$data_split, useNA = "ifany"))

algos <- c("C4.5","k-NN","SVM","RF","LR")

full_folds <- run_cv(df, yvar, features = feat_full, algos = algos,
                     type = "binary", pos_label = cfg$pos_label,
                     inner_k = cfg$inner_k, k_desired = cfg$cv_k, R = cfg$cv_R)

full_agg <- full_folds

triage_folds <- run_cv(df, yvar, features = feat_triage, algos = algos,
                       type = "binary", pos_label = cfg$pos_label,
                       inner_k = cfg$inner_k, k_desired = cfg$cv_k, R = cfg$cv_R)

triage_agg <- triage_folds

# ============================================================
# TABLE FACTORY
# ============================================================

# 1) Table 1
tbl1_baseline <- baseline_table(df, yvar, split_col = "data_split")

# 2) Table 2
mk_predictor_row <- function(predictors, set_name) {
  tibble::tibble(
    Set = set_name,
    Predictor = predictors
  ) |>
    dplyr::mutate(
      Type = vapply(Predictor, function(v) {
        if (v %in% names(df)) {
          if (is.numeric(df[[v]])) "numeric"
          else if (is.factor(df[[v]])) "factor"
          else "other"
        } else "unknown"
      }, character(1)),
      Preprocessing = dplyr::case_when(
        Type == "numeric" ~ "impute_median; yeojohnson; zscore",
        Type == "factor"  ~ "impute_mode; novel; other; one-hot",
        TRUE              ~ "—"
      )
    )
}
tbl2_predictors <- tibble::tibble(
  Predictor = sort(unique(c(feat_full, feat_triage)))
) %>%
  dplyr::mutate(
    In_Triage = Predictor %in% feat_triage,
    In_Full   = Predictor %in% feat_full
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Type = {
      if (Predictor %in% names(df)) {
        if (is.numeric(df[[Predictor]])) "numeric"
        else if (is.factor(df[[Predictor]])) "factor"
        else "other"
      } else "unknown"
    },
    Preprocessing = dplyr::case_when(
      Type == "numeric" ~ "impute_median; yeojohnson; zscore",
      Type == "factor"  ~ "impute_mode; novel; other; one-hot",
      TRUE              ~ "—"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    missingness_table(
      df,
      unique(c(feat_full, feat_triage)),
      split_col = "data_split"
    ) %>%
      dplyr::select(Predictor, Overall_N, Overall_Missing_pct, dplyr::starts_with("Missing_")),
    by = "Predictor"
  ) %>%
  dplyr::arrange(Predictor)

# 3) Table 3 (CV + holdouts) with CI clamping
agg_to_ci_table <- function(folds_tbl) {
  if (!nrow(folds_tbl)) return(tibble::tibble())
  folds_tbl %>%
    dplyr::group_by(Algorithm, Metric) %>%
    dplyr::summarise(
      Mean = mean(Value, na.rm = TRUE),
      CI_low = quantile(Value, probs = 0.025, na.rm = TRUE),
      CI_high= quantile(Value, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CI_low = pmax(ifelse(Metric %in% c("AUC","F1","Accuracy","Precision","Sensitivity","Specificity"), 0, -Inf), CI_low),
      CI_high= pmin(ifelse(Metric %in% c("AUC","F1","Accuracy","Precision","Sensitivity","Specificity"), 1, Inf), CI_high)
    )
}
cv_full_ci <- if (exists("full_folds")) agg_to_ci_table(full_folds) else tibble::tibble()
cv_tri_ci  <- if (exists("triage_folds")) agg_to_ci_table(triage_folds) else tibble::tibble()
cv_full_ci <- cv_full_ci %>% dplyr::mutate(Split = "Train-CV", Feature_Set = "Full")
cv_tri_ci  <- cv_tri_ci  %>% dplyr::mutate(Split = "Train-CV", Feature_Set = "Triage")

eval_all_holdouts_for_set <- function(feature_set, set_name, algos) {
  holder <- lapply(algos, function(a) run_holdout(df, yvar, features = feature_set, algo_name = a))
  names(holder) <- algos
  list(set_name = set_name, objects = holder)
}
hold_full  <- eval_all_holdouts_for_set(feat_full, "Full", algos)
hold_tri   <- eval_all_holdouts_for_set(feat_triage, "Triage", algos)

to_perf_table <- function(holder, split_name) {
  objs <- holder$objects
  out <- list()
  for (algo in names(objs)) {
    o <- objs[[algo]][[split_name]]
    if (is.null(o)) next
    auc  <- binary_auc(o$obs, o$p, pos_label = cfg$pos_label)
    aucC <- auc_delong_ci(o$obs, o$p, conf.level = 0.95)
    mets <- compute_metrics_binary(o$obs, o$pred, o$p, pos_label = cfg$pos_label)
    tbl <- tibble::tibble(
      Algorithm = algo,
      Feature_Set = holder$set_name,
      Split = tools::toTitleCase(split_name),
      Metric = c("AUC","MCC","F1","Accuracy","Precision","Sensitivity","Specificity","Threshold"),
      Point = c(auc, mets["MCC"], mets["F1"], mets["Accuracy"], mets["Precision"], mets["Sensitivity"], mets["Specificity"], o$threshold),
      CI_low = c(aucC[1], NA, NA, NA, NA, NA, NA, NA),
      CI_high= c(aucC[3], NA, NA, NA, NA, NA, NA, NA)
    )
    out[[algo]] <- tbl
  }
  dplyr::bind_rows(out)
}

test_full_tbl <- to_perf_table(hold_full, "test")
ext_full_tbl  <- to_perf_table(hold_full, "external")
test_tri_tbl  <- to_perf_table(hold_tri,  "test")
ext_tri_tbl   <- to_perf_table(hold_tri,  "external")

tbl3_perf_by_split <- dplyr::bind_rows(
  cv_full_ci %>% dplyr::rename(Point = Mean) %>% dplyr::select(Algorithm, Feature_Set, Split, Metric, Point, CI_low, CI_high),
  cv_tri_ci  %>% dplyr::rename(Point = Mean) %>% dplyr::select(Algorithm, Feature_Set, Split, Metric, Point, CI_low, CI_high),
  test_full_tbl, ext_full_tbl, test_tri_tbl, ext_tri_tbl
)

# 4) Table 4 Calibration pre/post
calib_one <- function(o) {
  if (is.null(o)) return(tibble::tibble())
  prepost <- calibration_metrics(o$obs, o$p, bins = 10, conf.level = 0.95, n_boot = cfg$n_boot_holdout)
  prepost
}
collect_calibration <- function(holder, split_name) {
  objs <- holder$objects
  out <- list()
  for (algo in names(objs)) {
    o <- objs[[algo]][[split_name]]
    if (is.null(o)) next
    tbl <- calib_one(o)
    if (!nrow(tbl)) next
    tbl$Algorithm <- algo; tbl$Feature_Set <- holder$set_name; tbl$Split <- tools::toTitleCase(split_name)
    out[[algo]] <- tbl
  }
  dplyr::bind_rows(out)
}
tbl4_calibration <- dplyr::bind_rows(
  collect_calibration(hold_full, "test"),
  collect_calibration(hold_full, "external"),
  collect_calibration(hold_tri,  "test"),
  collect_calibration(hold_tri,  "external")
)

# 5) Table 5 Confusion matrices (with Wilson CIs incl. precision)
confusion_table_one <- function(o) {
  if (is.null(o)) return(tibble::tibble())
  cm <- mk_cm(o$obs, o$pred, pos = cfg$pos_label)
  sens_ci <- binom_ci_wilson(cm$TP, cm$TP + cm$FN)
  spec_ci <- binom_ci_wilson(cm$TN, cm$TN + cm$FP)
  acc_ci  <- binom_ci_wilson(cm$TP + cm$TN, cm$N)
  prec_ci <- binom_ci_wilson(cm$TP, cm$TP + cm$FP)
  tibble::tibble(
    TP = cm$TP, FP = cm$FP, FN = cm$FN, TN = cm$TN,
    Sensitivity = cm$Sensitivity, Sens_CI_low = sens_ci[1], Sens_CI_high = sens_ci[2],
    Specificity = cm$Specificity, Spec_CI_low = spec_ci[1], Spec_CI_high = spec_ci[2],
    Precision = cm$Precision, Prec_CI_low = prec_ci[1], Prec_CI_high = prec_ci[2],
    Accuracy = cm$Accuracy, Acc_CI_low = acc_ci[1], Acc_CI_high = acc_ci[2],
    Threshold = o$threshold
  )
}
confusion_collect <- function(holder, split_name) {
  objs <- holder$objects
  out <- list()
  for (algo in names(objs)) {
    o <- objs[[algo]][[split_name]]
    if (is.null(o)) next
    tbl <- confusion_table_one(o)
    tbl$Algorithm <- algo; tbl$Feature_Set <- holder$set_name; tbl$Split <- tools::toTitleCase(split_name)
    out[[algo]] <- tbl
  }
  dplyr::bind_rows(out)
}
tbl5_confusion <- dplyr::bind_rows(
  confusion_collect(hold_full, "test"),
  confusion_collect(hold_full, "external"),
  confusion_collect(hold_tri,  "test"),
  confusion_collect(hold_tri,  "external")
)

# 6) Table 6 DCA numerics
dca_collect <- function(holder, split_name) {
  objs <- holder$objects
  out <- list()
  for (algo in names(objs)) {
    o <- objs[[algo]][[split_name]]
    if (is.null(o)) next
    dca <- decision_curve_table(o$obs, o$p)
    dca$Algorithm <- algo; dca$Feature_Set <- holder$set_name; dca$Split <- tools::toTitleCase(split_name)
    out[[algo]] <- dca
  }
  dplyr::bind_rows(out)
}
tbl6_decision_curve <- dplyr::bind_rows(
  dca_collect(hold_full, "test"),
  dca_collect(hold_full, "external"),
  dca_collect(hold_tri,  "test"),
  dca_collect(hold_tri,  "external")
)

# 7) Table 7 Subgroups
sub_collect <- function(holder, split_name) {
  objs <- holder$objects
  out <- list()
  for (algo in names(objs)) {
    o <- objs[[algo]][[split_name]]
    if (is.null(o)) next
    sg <- subgroup_metrics(o$holdout_df, o$obs, o$p, o$pred,
                           vars = c("group","Gender","Diagnosis","severity_admission"),
                           pos_label = cfg$pos_label)
    if (!nrow(sg)) next
    sg$Algorithm <- algo; sg$Feature_Set <- holder$set_name; sg$Split <- tools::toTitleCase(split_name)
    out[[algo]] <- sg
  }
  dplyr::bind_rows(out)
}
tbl7_subgroups <- dplyr::bind_rows(
  sub_collect(hold_full, "test"),
  sub_collect(hold_full, "external"),
  sub_collect(hold_tri,  "test"),
  sub_collect(hold_tri,  "external")
)

# 8) Table 8 Missingness
tbl8_missingness <- missingness_table(df, unique(c(feat_full, feat_triage)), split_col = "data_split") %>%
  dplyr::arrange(Predictor)

# 9) Table 9 Sensitivity
tbl9_sensitivity <- sensitivity_ledger(
  cv_tbl_full   = if (exists("full_agg")) full_agg else tibble::tibble(),
  cv_tbl_triage = if (exists("triage_agg")) triage_agg else tibble::tibble(),
  primary_algo  = "LR"
)

# 10) Table 10 TRIPOD-AI
tbl10_tripod_ai <- tripod_ai_skeleton()

# Supplement: hyperparameters and session info
collect_hparams <- function(holder, split_name) {
  objs <- holder$objects
  out <- list()
  for (algo in names(objs)) {
    o <- objs[[algo]][[split_name]]
    if (is.null(o)) next
    bt <- try(o$model$bestTune, silent = TRUE)
    if (inherits(bt, "try-error") || is.null(bt)) next
    kv <- tibble::as_tibble(bt) %>% mutate(across(everything(), as.character)) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "Param", values_to = "Value")
    kv$Algorithm <- algo; kv$Feature_Set <- holder$set_name; kv$Split <- tools::toTitleCase(split_name)
    kv$Threshold <- o$threshold; kv$Seed_CV <- cfg$seed_cv
    out[[algo]] <- kv[, c("Algorithm","Feature_Set","Split","Param","Value","Threshold","Seed_CV")]
  }
  dplyr::bind_rows(out)
}
supp_hparams <- dplyr::bind_rows(
  collect_hparams(hold_full, "test"),
  collect_hparams(hold_full, "external"),
  collect_hparams(hold_tri,  "test"),
  collect_hparams(hold_tri,  "external")
)

supp_session <- {
  pkgs <- sessionInfo()$otherPkgs
  if (is.null(pkgs)) tibble::tibble(Package = character(), Version = character())
  else tibble::tibble(
    Package = names(pkgs),
    Version = vapply(pkgs, function(x) x$Version, character(1))
  ) %>% dplyr::arrange(Package)
}

# ============================================================
# PRINT TABLES + SANITY CHECKS
# ============================================================

print_table <- function(title, x, n = 12) {
  cat(sprintf("\n[%s] n=%d x %d\n", title, nrow(x), ncol(x)))
  print(head(x, n))
}

check_tbl3 <- function(tbl) {
  if (!nrow(tbl)) { cat("\n[Sanity] Table 3 empty\n"); return(invisible()) }
  b <- tbl %>% dplyr::mutate(
    low_ok  = is.na(CI_low)  | (Metric %in% c("Threshold")) | (CI_low  >= 0),
    high_ok = is.na(CI_high) | (Metric %in% c("Threshold")) | (CI_high <= 1),
    order_ok= is.na(CI_low) | is.na(CI_high) | (CI_low <= CI_high)
  )
  if (all(b$low_ok) && all(b$high_ok) && all(b$order_ok)) cat("\n[Sanity] Table 3 CI bounds and ordering OK\n")
  else cat("\n[Sanity] Table 3 has CI issues\n")
}

cat(sprintf("\n[Sanity] bootstrap reps used for holdouts: %d\n", cfg$n_boot_holdout))

print_table("Table 1 Baseline characteristics", tbl1_baseline)
print_table("Table 2 Predictors and preprocessing", tbl2_predictors)
print_table("Table 3 Performance by split with CIs", tbl3_perf_by_split); check_tbl3(tbl3_perf_by_split)
print_table("Table 4 Calibration pre/post", tbl4_calibration)
print_table("Table 5 Confusion matrices", tbl5_confusion)
print_table("Table 6 Decision-curve numerics", tbl6_decision_curve)
print_table("Table 7 Subgroup performance", tbl7_subgroups)
print_table("Table 8 Missingness & availability", tbl8_missingness)
print_table("Table 9 Sensitivity analyses", tbl9_sensitivity)
print_table("Table 10 TRIPOD-AI checklist", tbl10_tripod_ai)
print_table("Supp hyperparameters", supp_hparams)
print_table("Supp session info", supp_session)

# Objects created:
# tbl1_baseline, tbl2_predictors, tbl3_perf_by_split, tbl4_calibration,
# tbl5_confusion, tbl6_decision_curve, tbl7_subgroups, tbl8_missingness,
# tbl9_sensitivity, tbl10_tripod_ai, supp_hparams, supp_session
