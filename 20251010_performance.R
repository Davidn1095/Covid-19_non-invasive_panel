# ===============================================================
# External-only heatmap, one panel
# Top = Complete, middle = Triage, bottom = Cascade (Triage → Complete if uncertain)
# Metrics shown with 95% percentile CIs
# TP TN FP FN added as rows with 95% CIs
# Split boundary = 2020-04-15
# Tuned on MCC in CV, τ chosen by MCC with BACC then median-p tie-break
# Calibration = Platt ridge fitted on train CV OOF, applied external
# Inner CV uses time blocked 5 folds, repeats disabled under time blocks
# 'group' kept for cohort filtering then ignored for modeling
# Engineered features: NLR MLR SIRI NMR CAR CLR PNI DLR
# ENN plus BLSMOTE applied within CV for all algorithms, optional oversample refit on full train
# ENN removals are logged per algorithm and feature set
# Leakage fix: WHO score dropped for Hospital_ID
# Recipes: NA indicators and Yeo-Johnson for non tree models
# Wider tuning for LR with ridge to lasso path
# C4.5 uses Laplace smoothing with robust probability columns
# Adds: Cascade model, triage then complete when in doubt, with band tuning on TRIAGE OOF
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(kknn)
  library(UBL)
  library(NoiseFiltersR)
  library(RWeka)
  library(randomForest)
  library(kernlab)
  library(ggplot2)
})

# ---------------- config ----------------
cfg <- list(
  file_path       = "biomarkers_acuteCOVID_meta.xlsx",
  sheet           = "meta",
  outcome         = "Hospital_ID",
  pos_label       = "Yes",
  neg_label       = "No",
  boundary_date   = as.Date("2020-04-15"),
  inner_k         = 5,
  inner_repeats   = 1,
  seed_cv         = 123,
  tune_len        = 10,
  B_boot          = 2000,
  calibration          = "platt_ridge",
  platt_ridge_lambda   = 1e-3,
  use_defaults         = FALSE,
  lr_solver            = "ridge",
  cv_blocked_by_time   = TRUE,
  adasyn_refit_on_train= TRUE,
  sampler              = "blsmote",  # "blsmote" or "adasyn"
  # ---- Cascade config ----
  cascade = list(
    enabled         = TRUE,
    tune_band       = TRUE,
    t_low           = 0.20,
    t_high          = 0.80,
    grid_tl         = seq(0.05, 0.45, by = 0.05),
    grid_th         = seq(0.55, 0.95, by = 0.05),
    max_defer       = 1.00
  )
)

# ---------------- feature sets ----------------
feat_triage <- c("Diagnosis","WHO_score_admission_mod","Age","Gender","SpO2_admission")
feat_full <- c(
  feat_triage,
  "CRP","D_Dimer","albumin",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc",
  "NLR","MLR","SIRI","NMR","CAR","CLR","PNI","DLR"
)

# ---- leak free admission modeling: drop WHO score when predicting Hospital_ID ----
if (identical(cfg$outcome, "Hospital_ID")) {
  feat_triage <- setdiff(feat_triage, "WHO_score_admission_mod")
  feat_full   <- setdiff(feat_full,   "WHO_score_admission_mod")
  stopifnot(!"WHO_score_admission_mod" %in% feat_triage,
            !"WHO_score_admission_mod" %in% feat_full)
  message("[leakage] WHO_score_admission_mod removed from feature sets for Hospital_ID.")
}

# ---------------- order ----------------
metric_order <- c("MCC","AUC-ROC","F1-score","Accuracy","Precision","Sensitivity","Specificity")
algo_order   <- c("C4.5","k-NN","SVM","RF","LR")

# ---------------- reproducibility ----------------
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(cfg$seed_cv %||% 123)

# ---------------- helpers ----------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))
sanitize_name <- function(x) gsub("[^A-Za-z0-9]", "", tolower(x))
clip01 <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)
effective_repeats <- function() if (isTRUE(cfg$cv_blocked_by_time)) 1 else (cfg$inner_repeats %||% 1)
seed_for <- function(algo, set_label, base = cfg$seed_cv) base + sum(utf8ToInt(paste(algo, set_label, sep = "_")))

make_caret_seeds <- function(num_resamples, num_tunes, base_seed) {
  num_tunes <- max(1L, as.integer(num_tunes %||% 1L))
  set.seed(base_seed)
  out <- vector("list", num_resamples + 1L)
  for (i in seq_len(num_resamples)) out[[i]] <- sample.int(999999, num_tunes)
  out[[num_resamples + 1L]] <- sample.int(999999, 1L)
  out
}

to_YN <- function(x, pos = cfg$pos_label, neg = cfg$neg_label) {
  make_fac <- function(z) factor(z, levels = c(neg, pos))
  if (is.factor(x)) return(make_fac(as.character(x)))
  if (is.logical(x)) return(make_fac(ifelse(x, pos, neg)))
  if (is.numeric(x)) {
    stopifnot(all(x %in% c(0,1), na.rm = TRUE))
    return(make_fac(ifelse(x == 1, pos, neg)))
  }
  if (is.character(x)) {
    s <- trimws(tolower(x))
    map_yes <- c("1","yes","y","true","pos","positive")
    map_no  <- c("0","no","n","false","neg","negative")
    out <- ifelse(s %in% map_yes, pos, ifelse(s %in% map_no, neg, NA_character_))
    if (any(is.na(out))) stop("Outcome not coercible to Yes/No")
    return(make_fac(out))
  }
  stop("Unsupported outcome type")
}
enforce_y_levels <- function(df, yvar) {
  df[[yvar]] <- factor(df[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  df
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

# ---------------- ENN logging ----------------
._enn_log <- new.env(parent = emptyenv())
._enn_log$events <- tibble::tibble(
  key = character(),
  before = integer(),
  after = integer(),
  removed = integer()
)
enn_log_summary <- function() {
  ev <- ._enn_log$events
  if (!nrow(ev)) return(tibble::tibble())
  ev %>%
    tidyr::separate(key, into = c("Algorithm","Set"), sep = "\\|", fill = "right") %>%
    dplyr::group_by(Algorithm, Set) %>%
    dplyr::summarise(
      calls = dplyr::n(),
      total_removed = sum(removed, na.rm = TRUE),
      mean_removed  = mean(removed, na.rm = TRUE),
      median_removed= stats::median(removed, na.rm = TRUE),
      min_removed   = min(removed, na.rm = TRUE),
      max_removed   = max(removed, na.rm = TRUE),
      frac_removed  = sum(removed, na.rm = TRUE) / sum(before, na.rm = TRUE),
      .groups = "drop"
    ) %>% dplyr::arrange(Algorithm, Set)
}

# ---------------- data loader plus cohort filter ----------------
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
  df <- enforce_y_levels(df, outcome)
  if ("WHO_score_admission_mod" %in% names(df)) {
    if (!is.numeric(df$WHO_score_admission_mod)) {
      df$WHO_score_admission_mod <- as_num(df$WHO_score_admission_mod)
    }
    if (!identical(outcome, "Hospital_ID")) {
      idx_ctrl_na <- with(df, is.na(WHO_score_admission_mod) & group == "CTRL_noCOVID")
      df$WHO_score_admission_mod[idx_ctrl_na] <- 0
    }
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

# ---------------- engineered features ----------------
engineer_features <- function(d) {
  eps <- 1e-6
  has <- function(...) all(c(...) %in% names(d))
  div_safe <- function(a, b) a / pmax(b, eps)
  albumin_gdl <- if (has("albumin")) d$albumin else NA_real_
  albumin_gl  <- if (has("albumin")) d$albumin * 10 else NA_real_
  lymph_mm3   <- if (has("lymphocyte_abs_number")) d$lymphocyte_abs_number * 1000 else NA_real_
  d$NLR  <- if (has("neutrophil_abs_number","lymphocyte_abs_number")) div_safe(d$neutrophil_abs_number, d$lymphocyte_abs_number) else NA_real_
  d$MLR  <- if (has("monocyte_abs_number","lymphocyte_abs_number"))    div_safe(d$monocyte_abs_number,    d$lymphocyte_abs_number) else NA_real_
  d$SIRI <- if (has("neutrophil_abs_number","monocyte_abs_number","lymphocyte_abs_number"))
    div_safe(d$neutrophil_abs_number * d$monocyte_abs_number, d$lymphocyte_abs_number) else NA_real_
  d$NMR  <- if (has("neutrophil_abs_number","monocyte_abs_number"))    div_safe(d$neutrophil_abs_number,  d$monocyte_abs_number)  else NA_real_
  d$CAR  <- if (has("CRP","albumin"))                                  div_safe(d$CRP, albumin_gl)       else NA_real_
  d$CLR  <- if (has("CRP","lymphocyte_abs_number"))                    div_safe(d$CRP, d$lymphocyte_abs_number) else NA_real_
  d$PNI  <- if (has("albumin","lymphocyte_abs_number"))                10 * albumin_gdl + 0.005 * lymph_mm3 else NA_real_
  d$DLR  <- if (has("D_Dimer","lymphocyte_abs_number"))                div_safe(d$D_Dimer, d$lymphocyte_abs_number) else NA_real_
  d
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

# ---------------- recipes ----------------
make_recipe <- function(dat, yvar, method) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
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
      step_indicate_na(all_numeric_predictors())  %>%
      step_novel(all_nominal_predictors())        %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_dummy(all_nominal_predictors())        %>%
      step_zv(all_predictors())                   %>%
      step_YeoJohnson(all_numeric_predictors())   %>%
      step_normalize(all_numeric_predictors())
  }
}

# ---------------- metrics ----------------
mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(0)
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
  TN <- as_num(tab[neg_label,  neg_label] %||% 0)
  FP <- as_num(tab[neg_label,  pos_label] %||% 0)
  FN <- as_num(tab[pos_label,  neg_label] %||% 0)
  preci <- if ((TP+FP)==0) NA_real_ else TP/(TP+FP)
  sens  <- if ((TP+FN)==0) NA_real_ else TP/(TP+FN)
  spec  <- if ((TN+FP)==0) NA_real_ else TN/(TN+FP)
  acc   <- (TP+TN)/sum(tab)
  f1    <- if (is.na(preci) || is.na(sens) || (preci+sens)==0) NA_real_ else 2*preci*sens/(preci+sens)
  mcc   <- mcc_from_counts(TP, FP, FN, TN)
  aucv  <- binary_auc(y, as_num(p_pos), pos_label = cfg$pos_label)
  c(MCC = mcc, AUC = aucv, `F1` = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec)
}

# ---------------- τ rule ----------------
best_threshold_mcc_bacc_med <- function(obs, p_pos,
                                        grid = seq(0.01, 0.99, by = 0.001),
                                        tol = 1e-12) {
  pos <- cfg$pos_label; neg <- cfg$neg_label
  y <- factor(obs, levels = c(neg, pos))
  nG <- length(grid)
  TP <- TN <- FP <- FN <- numeric(nG)
  for (i in seq_len(nG)) {
    t  <- grid[i]
    pr <- factor(ifelse(p_pos >= t, pos, neg), levels = c(neg, pos))
    tab <- table(y, pr)
    TP[i] <- as.numeric(tab[pos, pos] %||% 0)
    TN[i] <- as.numeric(tab[neg, neg] %||% 0)
    FP[i] <- as.numeric(tab[neg, pos] %||% 0)
    FN[i] <- as.numeric(tab[pos, neg] %||% 0)
  }
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- ifelse(den == 0, -Inf, num/den)
  sens <- ifelse((TP+FN)==0, NA_real_, TP/(TP+FN))
  spec <- ifelse((TN+FP)==0, NA_real_, TN/(TN+FP))
  bacc <- (sens + spec)/2; bacc[is.na(bacc)] <- -Inf
  idx <- which(mcc >= max(mcc, na.rm = TRUE) - tol)
  idx <- idx[bacc[idx] == max(bacc[idx], na.rm = TRUE)]
  medp <- stats::median(p_pos, na.rm = TRUE)
  dthr <- abs(grid[idx] - medp)
  idx  <- idx[dthr == min(dthr, na.rm = TRUE)]
  t_star <- min(grid[idx], na.rm = TRUE)
  list(t = t_star, mcc = max(mcc, na.rm = TRUE))
}

# ---------------- Weka safety helper ----------------
make_weka_safe <- function(d) {
  d2 <- as.data.frame(d)
  is_num <- vapply(d2, is.numeric, logical(1))
  if (any(is_num)) {
    d2[is_num] <- lapply(d2[is_num], function(v) { v[!is.finite(v)] <- NA_real_; v })
  }
  is_fac <- vapply(d2, is.factor, logical(1))
  if (any(is_fac)) {
    d2[is_fac] <- lapply(d2[is_fac], function(f) {
      f <- droplevels(f)
      if (nlevels(f) > 50) {
        tab <- sort(table(f), decreasing = TRUE)
        keep <- names(tab)[seq_len(50)]
        f <- forcats::fct_other(f, keep = keep, other_level = "other")
        f <- droplevels(f)
      }
      f
    })
  }
  all_na <- vapply(d2, function(col) all(is.na(col)), logical(1))
  if (any(all_na)) d2 <- d2[, !all_na, drop = FALSE]
  names(d2) <- make.names(names(d2), unique = TRUE)
  d2
}

# ---------------- custom caret model: J48 with Laplace and robust prob ----------------
get_J48_laplace_model <- function() {
  list(
    label = "C4.5 with Laplace",
    library = "RWeka",
    type = c("Classification"),
    parameters = data.frame(
      parameter = c("C","M"),
      class     = c("numeric","integer"),
      label     = c("Confidence Factor","Min Instances")
    ),
    grid = function(x, y, len = NULL, search = "grid") {
      expand.grid(C = c(0.01, 0.05, 0.10, 0.25, 0.50),
                  M = c(2, 5, 10))
    },
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      dat <- if (inherits(x, "recipe")) recipes::juice(recipes::prep(x)) else as.data.frame(x)
      dat$.y <- y
      dat <- make_weka_safe(dat)
      if (!is.factor(dat$.y)) dat$.y <- factor(dat$.y, levels = lev)
      ctrl <- RWeka::Weka_control(C = param$C, M = param$M, A = TRUE)
      m <- RWeka::J48(.y ~ ., data = dat, control = ctrl, na.action = na.omit, ...)
      m$obsLevels <- lev
      m
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      nd <- make_weka_safe(newdata)
      pred <- predict(modelFit, newdata = nd)
      factor(pred, levels = modelFit$obsLevels)
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      nd <- make_weka_safe(newdata)
      pr <- as.data.frame(predict(modelFit, newdata = nd, type = "probability"))
      lev <- modelFit$obsLevels
      missing_cols <- setdiff(lev, colnames(pr))
      if (length(missing_cols)) for (cc in missing_cols) pr[[cc]] <- 0
      pr <- pr[, lev, drop = FALSE]
      colnames(pr) <- lev
      pr
    },
    levels = function(x) x$obsLevels,
    tags = c("Tree-Based Model"),
    sort = function(x) x[order(x$C, x$M), ]
  )
}

# ---------------- model registry ----------------
algo_requires <- function(method) {
  switch(method,
         "glmnet"    = "glmnet",
         "rf"        = "randomForest",
         "svmRadial" = "kernlab",
         "kknn"      = "kknn",
         NULL)
}
method_available <- function(method) {
  pkg <- algo_requires(method)
  if (is.null(pkg)) return(TRUE)
  requireNamespace(pkg, quietly = TRUE)
}
model_specs <- function(tune_len, algos = c("C4.5","k-NN","SVM","RF","LR")) {
  all <- list(
    "LR"   = if (isTRUE(cfg$use_defaults) && identical(cfg$lr_solver, "glm")) list(
      name   = "LR",
      method = "glm",
      args   = list(family = binomial(), control = glm.control(maxit = 200))
    ) else list(
      name   = "LR",
      method = "glmnet",
      grid   = expand.grid(
        alpha  = c(0, 0.05, 0.25, 0.5, 0.75, 1),
        lambda = 10^seq(-5, 0.5, length.out = 30)
      ),
      args   = list(family = "binomial", standardize = TRUE)
    ),
    "RF"   = list(name = "RF",   method = "rf",        tuneLength = tune_len),
    "SVM"  = list(name = "SVM",  method = "svmRadial", tuneLength = tune_len),
    "k-NN" = list(name = "k-NN", method = "kknn",      tuneLength = tune_len),
    "C4.5" = list(name = "C4.5", method = get_J48_laplace_model())
  )
  out <- all[names(all)[order(match(names(all), algo_order))]]
  out
}

# ---------------- CV folds ----------------
make_cv_indices <- function(y, k, seed, repeats = 1) {
  set.seed(seed)
  if (repeats > 1) {
    idx_tr <- caret::createMultiFolds(y, k = k, times = repeats)
  } else {
    idx_tr <- caret::createFolds(y, k = k, returnTrain = TRUE)
  }
  idx_te <- lapply(idx_tr, function(tr) setdiff(seq_along(y), tr))
  list(index = idx_tr, indexOut = idx_te)
}
make_timeblock_indices <- function(dat_train_with_date, k) {
  dcol <- resolve_date_col(dat_train_with_date)
  ord  <- order(dat_train_with_date[[dcol]])
  n    <- length(ord)
  fold_id <- cut(seq_len(n), breaks = k, labels = FALSE)
  idx_tr <- vector("list", k)
  idx_te <- vector("list", k)
  for (i in seq_len(k)) {
    te <- ord[fold_id == i]
    tr <- setdiff(seq_len(nrow(dat_train_with_date)), te)
    idx_tr[[i]] <- tr
    idx_te[[i]] <- te
  }
  list(index = idx_tr, indexOut = idx_te)
}

# ---------------- MCC summary for tuning ----------------
mccSummaryAdaptive <- function(data, lev = NULL, model = NULL) {
  pos <- cfg$pos_label
  if (!(pos %in% colnames(data)) && !is.null(lev) && length(lev)) {
    cand <- lev[length(lev)]
    if (cand %in% colnames(data)) pos <- cand
  }
  if (!(pos %in% colnames(data))) return(c(MCC = NA_real_))
  obs <- factor(data$obs, levels = c(cfg$neg_label, cfg$pos_label))
  p   <- as.numeric(data[[pos]])
  if (length(stats::na.omit(unique(obs))) < 2) {
    pred <- factor(ifelse(p >= 0.5, cfg$pos_label, cfg$neg_label), levels = levels(obs))
    tab <- table(obs, pred)
    TP <- as.numeric(tab[cfg$pos_label, cfg$pos_label] %||% 0)
    TN <- as.numeric(tab[cfg$neg_label, cfg$neg_label] %||% 0)
    FP <- as.numeric(tab[cfg$neg_label, cfg$pos_label] %||% 0)
    FN <- as.numeric(tab[cfg$pos_label, cfg$neg_label] %||% 0)
    mcc <- mcc_from_counts(TP, FP, FN, TN)
    return(c(MCC = mcc))
  }
  br <- best_threshold_mcc_bacc_med(obs, p)
  c(MCC = as.numeric(br$mcc))
}

# ---------------- ENN plus BLSMOTE sampler ----------------
apply_ENN <- function(df, yname = ".y", k = 3) {
  df[[yname]] <- factor(df[[yname]], levels = c(cfg$neg_label, cfg$pos_label))
  before_n <- nrow(df)
  out <- try(NoiseFiltersR::ENN(stats::as.formula(paste(yname, "~ .")), df, k = k), silent = TRUE)
  out <- if (!inherits(out, "try-error") && !is.null(out$cleanData)) out$cleanData else df
  ok <- nrow(out) > 0 && length(unique(out[[yname]])) == length(levels(df[[yname]]))
  if (!ok) out <- df
  out[[yname]] <- factor(out[[yname]], levels = c(cfg$neg_label, cfg$pos_label))
  after_n <- nrow(out)
  key <- getOption("enn_log_key", default = "unknown|unknown")
  ._enn_log$events <- dplyr::bind_rows(
    ._enn_log$events,
    tibble::tibble(key = key, before = before_n, after = after_n, removed = max(0L, before_n - after_n))
  )
  out
}
noise_adasyn_sampler <- function(x, y) {
  y <- factor(y, levels = c(cfg$neg_label, cfg$pos_label))
  if (length(unique(y)) < 2) return(list(x = x, y = y))
  df <- as.data.frame(x)
  df$.y <- y
  df2 <- apply_ENN(df, ".y", k = 3)
  ad <- try(UBL::ADASYN(.y ~ ., df2, C.perc = "balance", k = 5, dist = "HEOM"), silent = TRUE)
  res <- if (!inherits(ad, "try-error") && length(unique(df2$.y)) > 1) ad else {
    up <- try(UBL::RandOverClassif(.y ~ ., df2, C.perc = "balance"), silent = TRUE)
    if (!inherits(up, "try-error")) up else df2
  }
  res$.y <- factor(res$.y, levels = c(cfg$neg_label, cfg$pos_label))
  list(x = res[, setdiff(names(res), ".y"), drop = FALSE], y = res$.y)
}
blsmote_on_df <- function(dat, yvar) {
  df <- enforce_y_levels(dat, yvar)
  f  <- stats::as.formula(paste(yvar, "~ ."))
  out <- try(UBL::BLSMOTE(f, df, C.perc = "balance", k = 5), silent = TRUE)
  if (inherits(out, "try-error")) df else as.data.frame(out)
}
noise_blsmote_sampler <- function(x, y) {
  y <- factor(y, levels = c(cfg$neg_label, cfg$pos_label))
  if (length(unique(y)) < 2) return(list(x = x, y = y))
  df <- as.data.frame(x); df.y <- y
  names(df)[names(df) == "。.y"] <- ".y"
  if (!".y" %in% names(df)) df$.y <- y else df$.y <- df[[".y"]]
  df2 <- apply_ENN(df, ".y", k = 3)
  sm  <- try(UBL::BLSMOTE(.y ~ ., df2, C.perc = "balance", k = 5), silent = TRUE)
  res <- if (!inherits(sm, "try-error")) sm else df2
  res$.y <- factor(res$.y, levels = c(cfg$neg_label, cfg$pos_label))
  list(x = res[, setdiff(names(res), ".y"), drop = FALSE], y = res$.y)
}
oversample_on_df <- function(dat, yvar) {
  if (identical(cfg$sampler, "blsmote")) blsmote_on_df(dat, yvar) else adasyn_on_df(dat, yvar)
}
adasyn_on_df <- function(dat, yvar) {
  df <- enforce_y_levels(dat, yvar)
  f <- stats::as.formula(paste(yvar, "~ ."))
  out <- try(UBL::ADASYN(f, df, C.perc = "balance", k = 5, dist = "HEOM"), silent = TRUE)
  if (inherits(out, "try-error")) df else as.data.frame(out)
}

# ---------------- calibrators ----------------
get_oof_bestTune <- function(fit, pos_label = cfg$pos_label) {
  pd <- try(fit$pred, silent = TRUE)
  if (inherits(pd, "try-error") || is.null(pd) || !nrow(pd)) return(NULL)
  bt <- fit$bestTune; if (is.null(bt) || !nrow(bt)) return(NULL)
  for (nm in names(bt)) pd <- pd[pd[[nm]] == bt[[nm]][1], , drop = FALSE]
  lvl_in_pred <- intersect(levels(pd$obs), names(pd))
  pcol <- if (pos_label %in% names(pd)) pos_label else tail(lvl_in_pred, 1)
  if (is.null(pcol) || !(pcol %in% names(pd))) return(NULL)
  tibble::tibble(obs = pd$obs, p = as_num(pd[[pcol]])) %>% dplyr::filter(!is.na(p))
}
fit_platt_from_oof <- function(y, p, pos = cfg$pos_label) {
  y01 <- as.integer(factor(y, levels = c(cfg$neg_label, pos))) - 1L
  x   <- qlogis(clip01(p))
  m   <- try(stats::glm(y01 ~ x, family = binomial()), silent = TRUE)
  if (inherits(m, "try-error")) return(function(z) clip01(z))
  function(p_new) stats::predict(m, newdata = data.frame(x = qlogis(clip01(p_new))), type = "response") %>% as_num()
}
fit_platt_ridge_from_oof <- function(y, p, pos = cfg$pos_label,
                                     lambda = cfg$platt_ridge_lambda %||% 1e-4) {
  y01 <- as.integer(factor(y, levels = c(cfg$neg_label, cfg$pos_label))) - 1L
  X   <- matrix(qlogis(clip01(p)), ncol = 1)
  g   <- try(glmnet::glmnet(X, y01, family = "binomial", alpha = 0,
                            lambda = lambda, standardize = FALSE), silent = TRUE)
  if (inherits(g, "try-error")) return(function(z) clip01(z))
  function(p_new) {
    pn <- matrix(qlogis(clip01(p_new)), ncol = 1)
    as.numeric(glmnet::predict(g, pn, type = "response")[,1])
  }
}
get_calibrator_from_oof <- function(oof) {
  if (is.null(oof) || !nrow(oof)) return(function(z) z)
  switch(tolower(cfg$calibration),
         "platt"        = fit_platt_from_oof(oof$obs, oof$p),
         "platt_ridge"  = fit_platt_ridge_from_oof(oof$obs, oof$p),
         "isotonic"     = {
           y01 <- as.integer(factor(oof$obs, levels = c(cfg$neg_label, cfg$pos_label))) - 1L
           o   <- order(oof$p); xs <- as_num(oof$p[o]); ys <- y01[o]
           iso <- stats::isoreg(xs, ys)
           xs2 <- as_num(iso$x); ys2 <- as_num(iso$yf)
           function(p_new) clip01(approx(xs2, ys2, xout = p_new, rule = 2, ties = mean)$y)
         },
         function(z) z)
}

# ---------------- training and holdout ----------------
ensure_two_classes <- function(d, y) {
  if (length(unique(d[[y]])) < 2) oversample_on_df(d, y) else d
}

train_inner <- function(dat_train, yvar, spec, inner_k, cv_idx, seed_base) {
  dat_train <- enforce_y_levels(dat_train, yvar)
  is_j48_custom <- is.list(spec$method) && identical(spec$name, "C4.5")
  rec <- if (is_j48_custom) {
    recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat_train) %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_impute_mode(all_nominal_predictors())  %>%
      step_novel(all_nominal_predictors())        %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_zv(all_predictors())
  } else {
    make_recipe(dat_train, yvar, method = spec$method)
  }
  num_tunes <- if (!is.null(spec$grid)) {
    nrow(spec$grid)
  } else if (!is.null(spec$tuneLength)) {
    spec$tuneLength
  } else if (is_j48_custom && is.list(spec$method) && !is.null(spec$method$grid)) {
    tmpg <- try(spec$method$grid(x = NULL, y = NULL, len = cfg$tune_len %||% 10, search = "grid"), silent = TRUE)
    if (inherits(tmpg, "try-error")) 15L else max(1L, nrow(tmpg))
  } else {
    1L
  }
  sampler_fun <- if (identical(cfg$sampler, "blsmote")) noise_blsmote_sampler else noise_adasyn_sampler
  tr_ctrl <- caret::trainControl(
    method = if (isTRUE(cfg$cv_blocked_by_time) || effective_repeats() == 1) "cv" else "repeatedcv",
    number = inner_k,
    classProbs = TRUE,
    summaryFunction = mccSummaryAdaptive,
    savePredictions = "final",
    allowParallel = TRUE,
    index = cv_idx$index, indexOut = cv_idx$indexOut,
    seeds = make_caret_seeds(length(cv_idx$index), num_tunes, seed_base),
    sampling = sampler_fun
  )
  base_args <- list(
    x = rec, data = dat_train, method = spec$method,
    trControl = tr_ctrl, metric = "MCC"
  )
  if (!is.null(spec$grid)) base_args$tuneGrid <- spec$grid else base_args$tuneLength <- spec$tuneLength
  if (!is.null(spec$args)) base_args <- c(base_args, spec$args)
  fit0 <- tryCatch(
    do.call(caret::train, base_args),
    error = function(e) {
      message(sprintf("[C4.5] training error caught: %s", conditionMessage(e)))
      structure(e, class = c("try-error","error","condition"))
    }
  )
  if (inherits(fit0, "try-error") || is.null(fit0$bestTune)) return(fit0)
  tr_ctrl2 <- caret::trainControl(method = "none", classProbs = TRUE, summaryFunction = mccSummaryAdaptive)
  base_dat <- if (isTRUE(cfg$adasyn_refit_on_train)) oversample_on_df(dat_train, yvar) else dat_train
  final_dat <- ensure_two_classes(base_dat, yvar)
  rec2 <- if (is_j48_custom) {
    recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = final_dat) %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_impute_mode(all_nominal_predictors())  %>%
      step_novel(all_nominal_predictors())        %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_zv(all_predictors())
  } else {
    make_recipe(final_dat, yvar, method = spec$method)
  }
  base2 <- list(
    x = rec2, data = final_dat, method = spec$method,
    trControl = tr_ctrl2, metric = "MCC",
    tuneGrid = fit0$bestTune
  )
  if (!is.null(spec$args)) base2 <- c(base2, spec$args)
  tryCatch(
    do.call(caret::train, base2),
    error = function(e) {
      message(sprintf("[C4.5] refit error caught: %s", conditionMessage(e)))
      structure(e, class = c("try-error","error","condition"))
    }
  )
}

run_holdout <- function(df, yvar, features, algo_name, set_label) {
  specs <- model_specs(cfg$tune_len, algos = algo_name)
  if (!length(specs)) return(NULL)
  spec  <- specs[[algo_name]]
  row_train <- which(df$data_split == "train")
  row_ext   <- which(df$data_split == "external")
  model_cols <- unique(c(features, yvar))
  model_cols <- intersect(model_cols, names(df))
  dtrain <- df[row_train, model_cols, drop = FALSE]
  dext   <- df[row_ext,   model_cols, drop = FALSE]
  if (!nrow(dtrain) || !nrow(dext)) return(NULL)
  dtrain_with_date <- df[row_train, , drop = FALSE]
  seed_base <- seed_for(algo_name, set_label)
  cv_idx <- if (isTRUE(cfg$cv_blocked_by_time)) {
    make_timeblock_indices(dtrain_with_date, k = cfg$inner_k)
  } else {
    make_cv_indices(dtrain[[yvar]], k = cfg$inner_k, seed = seed_base, repeats = effective_repeats())
  }
  old_key <- getOption("enn_log_key", NULL)
  options(enn_log_key = paste(algo_name, set_label, sep = "|"))
  on.exit(options(enn_log_key = old_key), add = TRUE)
  fit <- try(train_inner(dtrain, yvar, spec, cfg$inner_k, cv_idx, seed_base), silent = TRUE)
  if (inherits(fit, "try-error")) {
    message(sprintf("[skip] %s returned error for %s: %s",
                    algo_name, set_label, attr(fit, "condition")$message %||% as.character(fit)))
    return(NULL)
  }
  if (is.null(fit) || is.null(fit$finalModel)) {
    message(sprintf("[skip] %s produced no finalModel for %s", algo_name, set_label))
    return(NULL)
  }
  oof <- get_oof_bestTune(fit, pos_label = cfg$pos_label)
  cal_fun <- get_calibrator_from_oof(oof)
  thr_info <- if (!is.null(oof) && nrow(oof)) {
    best_threshold_mcc_bacc_med(oof$obs, cal_fun(oof$p))
  } else {
    p_tr_raw <- try(as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label]), silent = TRUE)
    if (inherits(p_tr_raw, "try-error")) return(NULL)
    p_tr <- cal_fun(p_tr_raw)
    best_threshold_mcc_bacc_med(dtrain[[yvar]], p_tr)
  }
  thr <- thr_info$t
  p_ex_raw <- try(as.numeric(predict(fit, newdata = dext, type = "prob")[, cfg$pos_label]), silent = TRUE)
  if (inherits(p_ex_raw, "try-error") || all(is.na(p_ex_raw))) return(NULL)
  p_ex <- cal_fun(p_ex_raw)
  if (sd(p_ex, na.rm = TRUE) < .Machine$double.eps) {
    set.seed(seed_base + 4242)
    p_ex <- pmin(pmax(p_ex + rnorm(length(p_ex), sd = 1e-9), 0), 1)
  }
  neg  <- cfg$neg_label
  pred <- factor(ifelse(p_ex >= thr, cfg$pos_label, neg), levels = c(neg, cfg$pos_label))
  list(
    obs = factor(dext[[yvar]], levels = c(neg, cfg$pos_label)),
    p = p_ex, pred = pred, threshold = thr,
    oof = oof
  )
}

# ---------------- bootstrap CIs ----------------
bootstrap_ci_all_metrics <- function(obs, pred, p, B = cfg$B_boot, level = 0.95, seed = cfg$seed_cv + 99) {
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
  alpha <- (1 - level) / 2
  lwr <- apply(M, 2, quantile, probs = alpha,  na.rm = TRUE)
  upr <- apply(M, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  sdv <- apply(M, 2, sd, na.rm = TRUE)
  list(lwr = lwr, upr = upr, sd = sdv)
}
bootstrap_ci_counts <- function(obs, pred, B = cfg$B_boot, level = 0.95, seed = cfg$seed_cv + 777) {
  set.seed(seed)
  n <- length(obs)
  labs <- c("TP","TN","FP","FN")
  M <- matrix(NA_real_, nrow = B, ncol = 4, dimnames = list(NULL, labs))
  pos <- cfg$pos_label; neg <- cfg$neg_label
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    yb  <- factor(obs[idx],  levels = c(neg, pos))
    pb  <- factor(as.character(pred[idx]), levels = c(neg, pos))
    tab <- table(yb, pb, useNA = "no")
    TP <- as_num(tab[pos, pos] %||% 0)
    TN <- as_num(tab[neg, neg] %||% 0)
    FP <- as_num(tab[neg, pos] %||% 0)
    FN <- as_num(tab[pos, neg] %||% 0)
    M[b,] <- c(TP, TN, FP, FN)
  }
  alpha <- (1 - level) / 2
  lwr <- apply(M, 2, quantile, probs = alpha,  na.rm = TRUE)
  upr <- apply(M, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  list(lwr = lwr, upr = upr)
}

# -------------------- Cascade helpers --------------------
# Tune uncertainty band [t_low, t_high] on TRIAGE OOF by maximizing MCC
tune_band_from_triage_oof <- function(oof_df,
                                      thr_tri,
                                      grid_tl   = cfg$cascade$grid_tl,
                                      grid_th   = cfg$cascade$grid_th,
                                      max_defer = cfg$cascade$max_defer) {
  stopifnot(all(c("obs","p") %in% names(oof_df)))
  pos <- cfg$pos_label; neg <- cfg$neg_label
  best <- NULL
  for (tl in grid_tl) for (th in grid_th) if (tl < th) {
    defer <- (oof_df$p > tl & oof_df$p < th)
    defer_rate <- mean(defer)
    if (defer_rate > max_defer + 1e-12) next
    keep <- !defer
    if (!any(keep)) next
    pred_keep <- factor(ifelse(oof_df$p[keep] >= thr_tri, pos, neg), levels = c(neg, pos))
    obs_keep  <- factor(oof_df$obs[keep], levels = c(neg, pos))
    m <- compute_metrics_binary(obs_keep, pred_keep, p_pos = oof_df$p[keep], pos_label = pos, neg_label = neg)
    score <- as_num(m["MCC"])
    cand <- list(t_low = tl, t_high = th, mcc = score, defer = defer_rate, n_keep = sum(keep))
    if (is.null(best) ||
        cand$mcc > best$mcc + 1e-12 ||
        (abs(cand$mcc - best$mcc) < 1e-12 && (cand$defer < best$defer - 1e-12 ||
                                              (abs(cand$defer - best$defer) < 1e-12 && cand$n_keep > best$n_keep)))) best <- cand
  }
  if (is.null(best)) {
    list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high, mcc = NA_real_, defer = NA_real_, n_keep = NA_integer_)
  } else best
}

# Apply cascade on EXTERNAL: outside band use TRIAGE, inside band use COMPLETE
apply_cascade_external <- function(tri_res, comp_res, t_low, t_high) {
  stopifnot(length(tri_res$obs) == length(comp_res$obs))
  obs <- tri_res$obs
  pos <- cfg$pos_label; neg <- cfg$neg_label
  defer <- (tri_res$p > t_low & tri_res$p < t_high)
  use_tri <- !defer
  pred_tri <- factor(ifelse(tri_res$p >= tri_res$threshold, pos, neg), levels = c(neg, pos))
  pred_com <- factor(ifelse(comp_res$p >= comp_res$threshold, pos, neg), levels = c(neg, pos))
  pred_final <- ifelse(use_tri, as.character(pred_tri), as.character(pred_com)) %>%
    factor(levels = c(neg, pos))
  p_final <- ifelse(use_tri, tri_res$p, comp_res$p)
  list(
    obs = obs,
    pred = pred_final,
    p = p_final,
    defer_rate = mean(defer),
    defer_n = sum(defer),
    n_external = length(defer)
  )
}

# ---------------- External summary for heatmap, Complete Triage Cascade ----------------
external_summary_for_heatmap_with_cascade <- function(df, yvar, feat_full, feat_triage) {
  specs <- model_specs(cfg$tune_len)
  get_counts <- function(obs, pred) {
    pos <- cfg$pos_label; neg <- cfg$neg_label
    y  <- factor(obs, levels = c(neg, pos))
    pr <- factor(as.character(pred), levels = c(neg, pos))
    tab <- table(y, pr, useNA = "no")
    c(
      TP = as_num(tab[pos, pos] %||% 0),
      TN = as_num(tab[neg, neg] %||% 0),
      FP = as_num(tab[neg, pos] %||% 0),
      FN = as_num(tab[pos, neg] %||% 0)
    )
  }
  build_all_three <- function(algo) {
    oC <- run_holdout(df, yvar, feat_full,  algo, "Complete")
    oT <- run_holdout(df, yvar, feat_triage, algo, "Triage")
    if (is.null(oC) || is.null(oT)) {
      message(sprintf("[skip] %s returned NULL for %s", algo, ifelse(is.null(oC), "Complete", "Triage")))
      return(NULL)
    }
    # 1) Complete plus Triage
    vC <- compute_metrics_binary(oC$obs, oC$pred, oC$p, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
    vT <- compute_metrics_binary(oT$obs, oT$pred, oT$p, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
    ciC <- bootstrap_ci_all_metrics(oC$obs, oC$pred, oC$p)
    ciT <- bootstrap_ci_all_metrics(oT$obs, oT$pred, oT$p)
    cntC <- get_counts(oC$obs, oC$pred)
    cntT <- get_counts(oT$obs, oT$pred)
    ciCntC <- bootstrap_ci_counts(oC$obs, oC$pred)
    ciCntT <- bootstrap_ci_counts(oT$obs, oT$pred)
    # 2) Cascade band tuning on TRIAGE OOF
    if (isTRUE(cfg$cascade$enabled)) {
      if (!is.null(oT$oof) && nrow(oT$oof) > 0 && is.finite(oT$threshold)) {
        tuned <- tune_band_from_triage_oof(oT$oof, thr_tri = oT$threshold)
        t_low  <- tuned$t_low; t_high <- tuned$t_high
      } else {
        t_low  <- cfg$cascade$t_low; t_high <- cfg$cascade$t_high
      }
      cas <- apply_cascade_external(oT, oC, t_low, t_high)
      vS  <- compute_metrics_binary(cas$obs, cas$pred, cas$p, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
      ciS <- bootstrap_ci_all_metrics(cas$obs, cas$pred, cas$p)
      cntS <- get_counts(cas$obs, cas$pred)
      ciCntS <- bootstrap_ci_counts(cas$obs, cas$pred)
      Defer_N     <- cas$defer_n
      Defer_Rate  <- cas$defer_rate
      N_external  <- cas$n_external
    } else {
      vS <- c(MCC=NA, AUC=NA, F1=NA, Accuracy=NA, Precision=NA, Sensitivity=NA, Specificity=NA)
      ciS <- list(lwr = setNames(rep(NA,7), c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")),
                  upr = setNames(rep(NA,7), c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")))
      cntS <- c(TP=NA, TN=NA, FP=NA, FN=NA)
      ciCntS <- list(lwr = setNames(rep(NA,4), c("TP","TN","FP","FN")),
                     upr = setNames(rep(NA,4), c("TP","TN","FP","FN")))
      Defer_N    <- NA_integer_
      Defer_Rate <- NA_real_
      N_external <- NA_integer_
    }
    tibble::tibble(
      Algorithm   = algo,
      Metric_raw  = c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity"),
      # Complete
      Mean_C      = as.numeric(c(vC["MCC"], vC["AUC"], vC["F1"], vC["Accuracy"], vC["Precision"], vC["Sensitivity"], vC["Specificity"])),
      L_C         = as.numeric(ciC$lwr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      U_C         = as.numeric(ciC$upr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      # Triage
      Mean_T      = as.numeric(c(vT["MCC"], vT["AUC"], vT["F1"], vT["Accuracy"], vT["Precision"], vT["Sensitivity"], vT["Specificity"])),
      L_T         = as.numeric(ciT$lwr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      U_T         = as.numeric(ciT$upr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      # Cascade
      Mean_S      = as.numeric(c(vS["MCC"], vS["AUC"], vS["F1"], vS["Accuracy"], vS["Precision"], vS["Sensitivity"], vS["Specificity"])),
      L_S         = as.numeric(ciS$lwr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      U_S         = as.numeric(ciS$upr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      # Counts plus CIs
      TP_C = cntC["TP"], TN_C = cntC["TN"], FP_C = cntC["FP"], FN_C = cntC["FN"],
      TP_T = cntT["TP"], TN_T = cntT["TN"], FP_T = cntT["FP"], FN_T = cntT["FN"],
      TP_S = cntS["TP"], TN_S = cntS["TN"], FP_S = cntS["FP"], FN_S = cntS["FN"],
      TP_C_L = ciCntC$lwr["TP"], TP_C_U = ciCntC$upr["TP"],
      TN_C_L = ciCntC$lwr["TN"], TN_C_U = ciCntC$upr["TN"],
      FP_C_L = ciCntC$lwr["FP"], FP_C_U = ciCntC$upr["FP"],
      FN_C_L = ciCntC$lwr["FN"], FN_C_U = ciCntC$upr["FN"],
      TP_T_L = ciCntT$lwr["TP"], TP_T_U = ciCntT$upr["TP"],
      TN_T_L = ciCntT$lwr["TN"], TN_T_U = ciCntT$upr["TN"],
      FP_T_L = ciCntT$lwr["FP"], FP_T_U = ciCntT$upr["FP"],
      FN_T_L = ciCntT$lwr["FN"], FN_T_U = ciCntT$upr["FN"],
      TP_S_L = ciCntS$lwr["TP"], TP_S_U = ciCntS$upr["TP"],
      TN_S_L = ciCntS$lwr["TN"], TN_S_U = ciCntS$upr["TN"],
      FP_S_L = ciCntS$lwr["FP"], FP_S_U = ciCntS$upr["FP"],
      FN_S_L = ciCntS$lwr["FN"], FN_S_U = ciCntS$upr["FN"],
      # New: cascade defer summary
      Defer_N    = Defer_N,
      Defer_Rate = Defer_Rate,
      N_external = N_external
    )
  }
  out <- lapply(names(specs), build_all_three)
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) stop("No algorithms produced results, see [skip] messages above")
  dfw <- dplyr::bind_rows(out)
  name_map <- c("MCC" = "MCC", "AUC" = "AUC-ROC", "F1" = "F1-score",
                "Accuracy" = "Accuracy", "Precision" = "Precision",
                "Sensitivity" = "Sensitivity", "Specificity" = "Specificity")
  dfw %>%
    dplyr::mutate(
      Metric    = factor(name_map[Metric_raw], levels = rev(metric_order), ordered = TRUE),
      Algorithm = factor(Algorithm, levels = algo_order, ordered = TRUE),
      MaxMean   = pmax(Mean_C, Mean_T, Mean_S)
    ) %>%
    dplyr::arrange(Algorithm, Metric)
}

# ---------------- figure, rows: Complete Triage Cascade ----------------
plot_external_heatmap <- function(dfw, base_size = 12, family = "sans", digits_label = 4) {
  fmt_val <- function(m, l, u) ifelse(is.na(m), "NA", sprintf("%.4f [%.4f, %.4f]", m, l, u))
  fmt_cnt <- function(n, l, u) sprintf("%d [%d, %d]", as.integer(n), as.integer(round(l)), as.integer(round(u)))
  dfw <- dfw %>% dplyr::mutate(Algorithm = factor(Algorithm, levels = algo_order, ordered = TRUE))
  base_levels <- levels(dfw$Algorithm)
  df_round <- dfw %>% dplyr::mutate(
    MeanC_r = round(Mean_C, digits_label),
    MeanT_r = round(Mean_T, digits_label),
    MeanS_r = round(Mean_S, digits_label)
  )
  row_max <- df_round %>%
    dplyr::group_by(Metric) %>%
    dplyr::summarise(
      MaxC_r = max(MeanC_r, na.rm = TRUE),
      MaxT_r = max(MeanT_r, na.rm = TRUE),
      MaxS_r = max(MeanS_r, na.rm = TRUE),
      .groups = "drop"
    )
  df_metrics <- df_round %>%
    dplyr::left_join(row_max, by = "Metric") %>%
    dplyr::mutate(
      Label_C = fmt_val(Mean_C, L_C, U_C),
      Label_T = fmt_val(Mean_T, L_T, U_T),
      Label_S = fmt_val(Mean_S, L_S, U_S),
      Bold_C  = !is.na(MeanC_r) & (MeanC_r == MaxC_r),
      Bold_T  = !is.na(MeanT_r) & (MeanT_r == MaxT_r),
      Bold_S  = !is.na(MeanS_r) & (MeanS_r == MaxS_r),
      Block   = "metrics"
    )
  desired_top_order <- c("TP","TN","FP","FN")
  counts_levels <- rev(desired_top_order)
  counts_core <- dfw %>%
    dplyr::distinct(Algorithm,
                    TP_C, TN_C, FP_C, FN_C,
                    TP_T, TN_T, FP_T, FN_T,
                    TP_S, TN_S, FP_S, FN_S,
                    TP_C_L, TP_C_U, TN_C_L, TN_C_U, FP_C_L, FP_C_U, FN_C_L, FN_C_U,
                    TP_T_L, TP_T_U, TN_T_L, TN_T_U, FP_T_L, FP_T_U, FN_T_L, FN_T_U,
                    TP_S_L, TP_S_U, TN_S_L, TN_S_U, FP_S_L, FP_S_U, FN_S_L, FN_S_U
    ) %>% dplyr::mutate(Algorithm_chr = as.character(Algorithm)) %>% dplyr::select(-Algorithm)
  mk_count_row <- function(which_label, nC, lC, uC, nT, lT, uT, nS, lS, uS) {
    tidyr::crossing(Algorithm_chr = base_levels) %>%
      dplyr::left_join(counts_core, by = "Algorithm_chr") %>%
      dplyr::transmute(
        Algorithm = factor(Algorithm_chr, levels = base_levels, ordered = TRUE),
        Metric    = factor(which_label, levels = counts_levels, ordered = TRUE),
        Label_C   = fmt_cnt(.data[[nC]], .data[[lC]], .data[[uC]]),
        Label_T   = fmt_cnt(.data[[nT]], .data[[lT]], .data[[uT]]),
        Label_S   = fmt_cnt(.data[[nS]], .data[[lS]], .data[[uS]]),
        Bold_C    = FALSE, Bold_T = FALSE, Bold_S = FALSE,
        Block     = "counts"
      )
  }
  df_counts <- dplyr::bind_rows(
    mk_count_row("TP","TP_C","TP_C_L","TP_C_U","TP_T","TP_T_L","TP_T_U","TP_S","TP_S_L","TP_S_U"),
    mk_count_row("TN","TN_C","TN_C_L","TN_C_U","TN_T","TN_T_L","TN_T_U","TN_S","TN_S_L","TN_S_U"),
    mk_count_row("FP","FP_C","FP_C_L","FP_C_U","FP_T","FP_T_L","FP_T_U","FP_S","FP_S_L","FP_S_U"),
    mk_count_row("FN","FN_C","FN_C_L","FN_C_U","FN_T","FN_T_L","FN_T_U","FN_S","FN_S_L","FN_S_U")
  )
  metric_levels <- levels(dfw$Metric) %||% unique(as.character(dfw$Metric))
  all_levels    <- c(counts_levels, metric_levels)
  df_metrics <- df_metrics %>% dplyr::mutate(Metric = factor(as.character(Metric), levels = all_levels, ordered = TRUE))
  df_counts  <- df_counts  %>% dplyr::mutate(Metric = factor(as.character(Metric), levels = all_levels, ordered = TRUE))
  stub_df <- dplyr::bind_rows(df_counts %>% dplyr::select(Metric),
                              df_metrics %>% dplyr::select(Metric)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Metric    = factor(as.character(Metric), levels = all_levels, ordered = TRUE),
      Algorithm = factor("Feature set", levels = c("Feature set", base_levels), ordered = TRUE),
      StubTop   = "Complete",
      StubMid   = "Triage",
      StubBottom= "Cascade"
    )
  x_levels <- c("Feature set", base_levels)
  y_idx <- sort(unique(as.numeric(stub_df$Metric)))
  y_top <- max(y_idx) + 0.5
  y_bot <- min(y_idx) - 0.5
  seg_dark  <- tibble::tibble(y = c(head(y_idx, -1) + 0.5, y_top, y_bot), x = 0.5, xend = length(x_levels) + 0.5)
  seg_white <- tibble::tibble(y = y_idx, x = 1.5, xend = length(x_levels) + 0.5)
  ggplot2::ggplot(NULL) +
    ggplot2::geom_tile(data = df_metrics, ggplot2::aes(x = Algorithm, y = Metric),
                       fill = "white", color = "white", linewidth = 0.25) +
    ggplot2::geom_tile(data = df_counts,  ggplot2::aes(x = Algorithm, y = Metric),
                       fill = "white", color = "white", linewidth = 0.25) +
    ggplot2::geom_tile(data = stub_df, inherit.aes = FALSE,
                       ggplot2::aes(x = Algorithm, y = Metric),
                       fill = "white", color = "white", linewidth = 0.25) +
    ggplot2::geom_segment(data = seg_dark,  ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          linewidth = 0.6, color = "grey20") +
    ggplot2::geom_segment(data = seg_white, ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          linewidth = 0.6, color = "white") +
    ggplot2::geom_text(data = df_metrics,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_C,
                                    fontface = ifelse(Bold_C, "bold", "plain")),
                       nudge_y = +0.28, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_metrics,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_T,
                                    fontface = ifelse(Bold_T, "bold", "plain")),
                       nudge_y = 0.00, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_metrics,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_S,
                                    fontface = ifelse(Bold_S, "bold", "plain")),
                       nudge_y = -0.28, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_C),
                       nudge_y = +0.28, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_T),
                       nudge_y = 0.00, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_S),
                       nudge_y = -0.28, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubTop),
                       nudge_y = +0.28, size = 3.1, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubMid),
                       nudge_y = 0.00, size = 3.1, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubBottom),
                       nudge_y = -0.28, size = 3.1, lineheight = 0.95, family = family) +
    ggplot2::scale_x_discrete(limits = c("Feature set", base_levels), position = "top") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = base_size, base_family = family) +
    ggplot2::theme(
      plot.title  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      axis.text.x.bottom = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(margin = ggplot2::margin(b = 6), face = "bold"),
      axis.text.y  = ggplot2::element_text(margin = ggplot2::margin(r = 6), face = "bold"),
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 14, 8, 10)
    )
}

# ========================= DRIVER =========================
set.seed(123)

ld  <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
df0 <- engineer_features(ld$df)
df  <- enforce_temporal_split(df0, admit_col = ld$admit_col, boundary_date = cfg$boundary_date)
yvar <- cfg$outcome

ext_df <- external_summary_for_heatmap_with_cascade(df, yvar, feat_full, feat_triage)
p_ext  <- plot_external_heatmap(ext_df)
print(p_ext)

# ENN summary
enn_sum <- enn_log_summary()
print(enn_sum)

# ---- cascade counts, how many are sent to Complete per algorithm ----
cascade_counts <- ext_df %>%
  dplyr::distinct(Algorithm, N_external, Defer_N, Defer_Rate) %>%
  dplyr::arrange(Algorithm) %>%
  dplyr::mutate(Percent = sprintf("%.1f%%", 100 * Defer_Rate))

print(cascade_counts)

# session info
si <- utils::sessionInfo()
print(si)
