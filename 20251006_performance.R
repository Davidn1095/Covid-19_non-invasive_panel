# ===============================================================
# External-only heatmap, one panel
# Top = Complete, bottom = Triage
# Metrics shown with 95% percentile CIs
# TP TN FP FN added as rows with 95% CIs
# Split boundary = 2020-04-15
# Tuned on MCC in CV, τ chosen by MCC with BACC then median-p tie-break
# WHO_score_admission_mod: controls' NA -> 0
# Calibration = Platt-ridge fitted on train CV-OOF, applied external
# Inner CV uses time-blocked 5 folds, repeats disabled under time blocks
# 'group' kept for cohort filtering then ignored for modeling
# Engineered features added: NLR MLR SIRI NMR CAR CLR PNI DLR
# ENN + ADASYN applied within CV for all algorithms, optional ADASYN refit on full train
# ENN removals are logged per algorithm and feature set
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
  tune_len        = 5,
  B_boot          = 2000,
  calibration          = "platt_ridge",
  platt_ridge_lambda   = 1e-3,
  use_defaults         = FALSE,
  lr_solver            = "ridge",
  cv_blocked_by_time   = TRUE,
  adasyn_refit_on_train= TRUE
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
  set.seed(base_seed)
  out <- vector("list", num_resamples + 1)
  for (i in seq_len(num_resamples)) out[[i]] <- sample.int(999999, num_tunes)
  out[[num_resamples + 1]] <- sample.int(999999, 1)
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
  df <- enforce_y_levels(df, outcome)
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
  TN <- as_num(tab[neg_label, neg_label] %||% 0)
  FP <- as_num(tab[neg_label, pos_label] %||% 0)
  FN <- as_num(tab[pos_label, neg_label] %||% 0)
  preci <- if ((TP+FP)==0) NA_real_ else TP/(TP+FP)
  sens  <- if ((TP+FN)==0) NA_real_ else TP/(TP+FN)
  spec  <- if ((TN+FP)==0) NA_real_ else TN/(TN+FP)
  acc   <- (TP+TN)/sum(tab)
  f1    <- if (is.na(preci) || is.na(sens) || (preci+sens)==0) NA_real_ else 2*preci*sens/(preci+sens)
  mcc   <- mcc_from_counts(TP, FP, FN, TN)
  aucv  <- binary_auc(y, as_num(p_pos), pos_label = cfg$pos_label)
  c(MCC = mcc, AUC = aucv, `F1` = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec)
}

# ---------------- τ rule (robust) ----------------
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
  bacc <- (sens + spec)/2
  bacc[is.na(bacc)] <- -Inf
  
  idx <- which(mcc >= max(mcc, na.rm = TRUE) - tol)
  idx <- idx[bacc[idx] == max(bacc[idx], na.rm = TRUE)]
  
  medp <- stats::median(p_pos, na.rm = TRUE)
  dthr <- abs(grid[idx] - medp)
  idx  <- idx[dthr == min(dthr, na.rm = TRUE)]
  
  t_star <- min(grid[idx], na.rm = TRUE)
  list(t = t_star, mcc = max(mcc, na.rm = TRUE))
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
model_specs <- function(tune_len, algos = c("C4.5","k-NN","SVM","RF","LR")) {
  all <- list(
    "LR"   = if (isTRUE(cfg$use_defaults) && identical(cfg$lr_solver, "glm")) list(
      name   = "LR",
      method = "glm",
      args   = list(family = binomial(), control = glm.control(maxit = 200))
    ) else if (isTRUE(cfg$use_defaults) && identical(cfg$lr_solver, "ridge")) list(
      name   = "LR",
      method = "glmnet",
      args   = list(family = "binomial", standardize = TRUE),
      grid   = data.frame(alpha = 0, lambda = 1e-4)
    ) else list(
      name   = "LR",
      method = "glmnet",
      grid   = expand.grid(
        alpha  = c(0, 0.05, 0.25, 0.5, 0.75, 1),
        lambda = 10^seq(-4, 1, length.out = 15)
      ),
      args   = list(family = "binomial", standardize = TRUE)
    ),
    "RF"   = list(name = "RF",   method = "rf",        tuneLength = tune_len),
    "SVM"  = list(name = "SVM",  method = "svmRadial", tuneLength = tune_len),
    "k-NN" = list(name = "k-NN", method = "kknn",      tuneLength = tune_len),
    "C4.5" = list(
      name   = "C4.5",
      method = "J48",
      grid   = expand.grid(
        C = c(0.01, 0.05, 0.10, 0.25, 0.50),
        M = c(2, 5, 10)
      )
    )
  )
  avail <- names(caret::getModelInfo())
  keep  <- names(all)[vapply(all, \(s) is.null(s$method) || s$method %in% avail, logical(1))]
  out   <- all[keep]
  out[names(out)[order(match(names(out), algo_order))]]
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

# ---------------- MCC summary for tuning (robust) ----------------
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
    pred <- factor(ifelse(p >= 0.5, cfg$pos_label, cfg$neg_label),
                   levels = levels(obs))
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

# ---------------- ENN + ADASYN sampler (with logging) ----------------
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
  if (!(pos_label %in% names(pd))) return(NULL)
  tibble::tibble(obs = pd$obs, p = as_num(pd[[pos_label]])) %>% dplyr::filter(!is.na(p))
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
  if (length(unique(d[[y]])) < 2) adasyn_on_df(d, y) else d
}

train_inner <- function(dat_train, yvar, spec, inner_k, cv_idx, seed_base) {
  dat_train <- enforce_y_levels(dat_train, yvar)
  rec <- make_recipe(dat_train, yvar, method = spec$method)
  
  if (isTRUE(cfg$use_defaults)) {
    tg <- NULL
    if (identical(spec$method, "rf")) {
      pp <- recipes::prep(rec, training = dat_train, retain = TRUE)
      X  <- recipes::juice(pp)
      p  <- ncol(X) - 1
      tg <- data.frame(mtry = max(1, floor(sqrt(max(1, p)))))
    } else if (identical(spec$method, "svmRadial")) {
      pp <- recipes::prep(rec, training = dat_train, retain = TRUE)
      Xj <- recipes::juice(pp)
      X  <- Xj[, setdiff(names(Xj), yvar), drop = FALSE]
      sg <- tryCatch(as.numeric(kernlab::sigest(as.matrix(X), frac = 1))[1],
                     error = function(e) 1 / max(1, ncol(X)))
      tg <- data.frame(C = 1, sigma = sg)
    } else if (identical(spec$method, "kknn")) {
      tg <- data.frame(kmax = 7, distance = 2, kernel = "optimal")
    } else if (identical(spec$method, "J48")) {
      tg <- data.frame(C = 0.25, M = 2)
    } else if (identical(spec$method, "glmnet")) {
      tg <- data.frame(alpha = 0, lambda = 1e-4)
    }
    tr_ctrl2 <- caret::trainControl(method = "none", classProbs = TRUE, summaryFunction = mccSummaryAdaptive)
    base_dat <- if (isTRUE(cfg$adasyn_refit_on_train)) adasyn_on_df(dat_train, yvar) else dat_train
    final_dat <- ensure_two_classes(base_dat, yvar)
    base2 <- list(
      x = make_recipe(final_dat, yvar, method = spec$method),
      data = final_dat, method = spec$method,
      trControl = tr_ctrl2, metric = "MCC"
    )
    if (!is.null(tg)) base2$tuneGrid <- tg
    if (!is.null(spec$args)) base2 <- c(base2, spec$args)
    return(try(do.call(caret::train, base2), silent = TRUE))
  }
  
  num_tunes   <- if (!is.null(spec$grid)) nrow(spec$grid) else spec$tuneLength
  tr_ctrl <- caret::trainControl(
    method = if (isTRUE(cfg$cv_blocked_by_time) || effective_repeats() == 1) "cv" else "repeatedcv",
    number = inner_k,
    classProbs = TRUE,
    summaryFunction = mccSummaryAdaptive,
    savePredictions = "final",
    allowParallel = TRUE,
    index = cv_idx$index, indexOut = cv_idx$indexOut,
    seeds = make_caret_seeds(length(cv_idx$index), num_tunes, seed_base),
    sampling = noise_adasyn_sampler
  )
  base_args <- list(
    x = rec, data = dat_train, method = spec$method,
    trControl = tr_ctrl, metric = "MCC"
  )
  if (!is.null(spec$grid)) base_args$tuneGrid <- spec$grid else base_args$tuneLength <- spec$tuneLength
  if (!is.null(spec$args)) base_args <- c(base_args, spec$args)
  fit0 <- try(do.call(caret::train, base_args), silent = TRUE)
  if (inherits(fit0, "try-error") || is.null(fit0$bestTune)) return(fit0)
  
  tr_ctrl2 <- caret::trainControl(method = "none", classProbs = TRUE, summaryFunction = mccSummaryAdaptive)
  base_dat <- if (isTRUE(cfg$adasyn_refit_on_train)) adasyn_on_df(dat_train, yvar) else dat_train
  final_dat <- ensure_two_classes(base_dat, yvar)
  base2 <- list(
    x = make_recipe(final_dat, yvar, method = spec$method),
    data = final_dat, method = spec$method,
    trControl = tr_ctrl2, metric = "MCC",
    tuneGrid = fit0$bestTune
  )
  if (!is.null(spec$args)) base2 <- c(base2, spec$args)
  try(do.call(caret::train, base2), silent = TRUE)
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
    make_cv_indices(dtrain[[yvar]], k = cfg$inner_k, seed = seed_base,
                    repeats = effective_repeats())
  }
  
  old_key <- getOption("enn_log_key", NULL)
  options(enn_log_key = paste(algo_name, set_label, sep = "|"))
  on.exit(options(enn_log_key = old_key), add = TRUE)
  
  fit <- try(train_inner(dtrain, yvar, spec, cfg$inner_k, cv_idx, seed_base), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
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
    p = p_ex, pred = pred, threshold = thr
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
    pb  <- factor(pred[idx], levels = c(neg, pos))
    tab <- table(yb, pb)
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

# ---------------- External summary for heatmap ----------------
external_summary_for_heatmap <- function(df, yvar, feat_full, feat_triage) {
  specs <- model_specs(cfg$tune_len)
  get_counts <- function(obs, pred) {
    pos <- cfg$pos_label; neg <- cfg$neg_label
    y <- factor(obs, levels = c(neg, pos))
    p <- factor(pred, levels = c(neg, pos))
    tab <- table(y, p)
    c(
      TP = as_num(tab[pos, pos] %||% 0),
      TN = as_num(tab[neg, neg] %||% 0),
      FP = as_num(tab[neg, pos] %||% 0),
      FN = as_num(tab[pos, neg] %||% 0)
    )
  }
  build_both <- function(algo) {
    oC <- run_holdout(df, yvar, feat_full,  algo, "Complete")
    oT <- run_holdout(df, yvar, feat_triage, algo, "Triage")
    if (is.null(oC) || is.null(oT)) {
      message(sprintf("[skip] %s returned NULL for %s", algo, ifelse(is.null(oC), "Complete", "Triage")))
      return(NULL)
    }
    vC <- compute_metrics_binary(oC$obs, oC$pred, oC$p, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
    vT <- compute_metrics_binary(oT$obs, oT$pred, oT$p, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
    ciC <- bootstrap_ci_all_metrics(oC$obs, oC$pred, oC$p)
    ciT <- bootstrap_ci_all_metrics(oT$obs, oT$pred, oT$p)
    cntC <- get_counts(oC$obs, oC$pred)
    cntT <- get_counts(oT$obs, oT$pred)
    ciCntC <- bootstrap_ci_counts(oC$obs, oC$pred)
    ciCntT <- bootstrap_ci_counts(oT$obs, oT$pred)
    tibble::tibble(
      Algorithm   = algo,
      Metric_raw  = c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity"),
      Mean_C      = as.numeric(c(vC["MCC"], vC["AUC"], vC["F1"], vC["Accuracy"], vC["Precision"], vC["Sensitivity"], vC["Specificity"])),
      L_C         = as.numeric(ciC$lwr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      U_C         = as.numeric(ciC$upr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      Mean_T      = as.numeric(c(vT["MCC"], vT["AUC"], vT["F1"], vT["Accuracy"], vT["Precision"], vT["Sensitivity"], vT["Specificity"])),
      L_T         = as.numeric(ciT$lwr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      U_T         = as.numeric(ciT$upr[c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")]),
      TP_C = cntC["TP"], TN_C = cntC["TN"], FP_C = cntC["FP"], FN_C = cntC["FN"],
      TP_T = cntT["TP"], TN_T = cntT["TN"], FP_T = cntT["FP"], FN_T = cntT["FN"],
      TP_C_L = ciCntC$lwr["TP"], TP_C_U = ciCntC$upr["TP"],
      TN_C_L = ciCntC$lwr["TN"], TN_C_U = ciCntC$upr["TN"],
      FP_C_L = ciCntC$lwr["FP"], FP_C_U = ciCntC$upr["FP"],
      FN_C_L = ciCntC$lwr["FN"], FN_C_U = ciCntC$upr["FN"],
      TP_T_L = ciCntT$lwr["TP"], TP_T_U = ciCntT$upr["TP"],
      TN_T_L = ciCntT$lwr["TN"], TN_T_U = ciCntT$upr["TN"],
      FP_T_L = ciCntT$lwr["FP"], FP_T_U = ciCntT$upr["FP"],
      FN_T_L = ciCntT$lwr["FN"], FN_T_U = ciCntT$upr["FN"]
    )
  }
  out <- lapply(names(specs), build_both)
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
      MaxMean   = pmax(Mean_C, Mean_T)
    ) %>%
    dplyr::arrange(Algorithm, Metric)
}

# ---------------- figure ----------------
# Bold highest value per row across algorithms, separately for Complete and Triage
plot_external_heatmap <- function(dfw, base_size = 12, family = "sans") {
  fmt_val  <- function(m, l, u) ifelse(is.na(m), "NA", sprintf("%.4f [%.4f, %.4f]", m, l, u))
  fmt_cnt  <- function(n, l, u) sprintf("%d [%d, %d]", as.integer(n), as.integer(round(l)), as.integer(round(u)))
  
  row_max_all <- dfw %>%
    dplyr::group_by(Metric) %>%
    dplyr::summarise(RowMax_C = max(Mean_C, na.rm = TRUE),
                     RowMax_T = max(Mean_T, na.rm = TRUE), .groups = "drop")
  
  df_metrics <- dfw %>%
    dplyr::left_join(row_max_all, by = "Metric") %>%
    dplyr::mutate(
      Label_C = fmt_val(Mean_C, L_C, U_C),
      Label_T = fmt_val(Mean_T, L_T, U_T),
      Bold_C  = !is.na(Mean_C) & (RowMax_C - Mean_C) <= 1e-6,
      Bold_T  = !is.na(Mean_T) & (RowMax_T - Mean_T) <= 1e-6,
      Algorithm = factor(as.character(Algorithm), levels = levels(Algorithm), ordered = FALSE),
      Block = "metrics"
    )
  
  counts_levels <- c("TP","TN","FP","FN")
  counts_core <- dfw %>%
    dplyr::distinct(Algorithm,
                    TP_C, TN_C, FP_C, FN_C, TP_T, TN_T, FP_T, FN_T,
                    TP_C_L, TP_C_U, TN_C_L, TN_C_U, FP_C_L, FP_C_U, FN_C_L, FN_C_U,
                    TP_T_L, TP_T_U, TN_T_L, TN_T_U, FP_T_L, FP_T_U, FN_T_L, FN_T_U
    ) %>% dplyr::mutate(Algorithm_chr = as.character(Algorithm)) %>% dplyr::select(-Algorithm)
  
  base_levels <- levels(df_metrics$Algorithm)
  mk_count_row <- function(which_label, nC, lC, uC, nT, lT, uT) {
    tidyr::crossing(Algorithm_chr = base_levels) %>%
      dplyr::left_join(counts_core, by = "Algorithm_chr") %>%
      dplyr::transmute(
        Algorithm = factor(Algorithm_chr, levels = base_levels, ordered = FALSE),
        Metric    = factor(which_label, levels = counts_levels, ordered = TRUE),
        Label_C   = fmt_cnt(.data[[nC]], .data[[lC]], .data[[uC]]),
        Label_T   = fmt_cnt(.data[[nT]], .data[[lT]], .data[[uT]]),
        Bold_C    = FALSE, Bold_T = FALSE,
        Block     = "counts"
      )
  }
  df_counts <- dplyr::bind_rows(
    mk_count_row("TP","TP_C","TP_C_L","TP_C_U","TP_T","TP_T_L","TP_T_U"),
    mk_count_row("TN","TN_C","TN_C_L","TN_C_U","TN_T","TN_T_L","TN_T_U"),
    mk_count_row("FP","FP_C","FP_C_L","FP_C_U","FP_T","FP_T_L","FP_T_U"),
    mk_count_row("FN","FN_C","FN_C_L","FN_C_U","FN_T","FN_T_L","FN_T_U")
  )
  
  metric_levels <- levels(df_metrics$Metric)
  all_levels    <- c(counts_levels, metric_levels)
  df_metrics <- df_metrics %>% dplyr::mutate(Metric = factor(as.character(Metric), levels = all_levels, ordered = TRUE))
  df_counts  <- df_counts  %>% dplyr::mutate(Metric = factor(as.character(Metric), levels = all_levels, ordered = TRUE))
  
  stub_df <- dplyr::bind_rows(df_counts %>% dplyr::select(Metric),
                              df_metrics %>% dplyr::select(Metric)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Metric    = factor(as.character(Metric), levels = all_levels, ordered = TRUE),
      Algorithm = factor("Feature set", levels = c("Feature set", base_levels), ordered = FALSE),
      StubTop   = "Complete",
      StubBottom= "Triage"
    )
  
  x_levels <- c("Feature set", base_levels)
  y_idx <- sort(unique(as.numeric(stub_df$Metric)))
  y_top <- max(y_idx) + 0.5
  y_bot <- min(y_idx) - 0.5
  seg_dark <- tibble::tibble(y = c(head(y_idx, -1) + 0.5, y_top, y_bot), x = 0.5, xend = length(x_levels) + 0.5)
  seg_white<- tibble::tibble(y = y_idx, x = 1.5, xend = length(x_levels) + 0.5)
  
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
                       nudge_y = +0.18, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_metrics,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_T,
                                    fontface = ifelse(Bold_T, "bold", "plain")),
                       nudge_y = -0.18, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_C),
                       nudge_y = +0.18, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_T),
                       nudge_y = -0.18, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubTop),
                       nudge_y = +0.18, size = 3.1, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubBottom),
                       nudge_y = -0.18, size = 3.1, lineheight = 0.95, family = family) +
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

# ---------------- paired external comparison ----------------
plot_external_heatmap <- function(dfw, base_size = 12, family = "sans", digits_label = 4) {
  fmt_val <- function(m, l, u) ifelse(is.na(m), "NA", sprintf("%.4f [%.4f, %.4f]", m, l, u))
  fmt_cnt <- function(n, l, u) sprintf("%d [%d, %d]", as.integer(n), as.integer(round(l)), as.integer(round(u)))
  
  # enforce algorithm order for stable layout
  dfw <- dfw %>%
    dplyr::mutate(Algorithm = factor(Algorithm, levels = algo_order, ordered = TRUE))
  base_levels <- levels(dfw$Algorithm)
  
  # bold all ties per metric-feature set using label rounding
  df_round <- dfw %>%
    dplyr::mutate(MeanC_r = round(Mean_C, digits_label),
                  MeanT_r = round(Mean_T, digits_label))
  
  row_max <- df_round %>%
    dplyr::group_by(Metric) %>%
    dplyr::summarise(MaxC_r = max(MeanC_r, na.rm = TRUE),
                     MaxT_r = max(MeanT_r, na.rm = TRUE),
                     .groups = "drop")
  
  df_metrics <- df_round %>%
    dplyr::left_join(row_max, by = "Metric") %>%
    dplyr::mutate(
      Label_C = fmt_val(Mean_C, L_C, U_C),
      Label_T = fmt_val(Mean_T, L_T, U_T),
      Bold_C  = !is.na(MeanC_r) & MeanC_r == MaxC_r,
      Bold_T  = !is.na(MeanT_r) & MeanT_r == MaxT_r,
      Block   = "metrics"
    )
  
  # counts block, order wanted top-to-bottom: TP, TN, FP, FN
  desired_top_order <- c("TP","TN","FP","FN")
  counts_levels <- rev(desired_top_order)  # factor levels bottom-to-top become TP, TN, FP, FN
  
  counts_core <- dfw %>%
    dplyr::distinct(Algorithm,
                    TP_C, TN_C, FP_C, FN_C, TP_T, TN_T, FP_T, FN_T,
                    TP_C_L, TP_C_U, TN_C_L, TN_C_U, FP_C_L, FP_C_U, FN_C_L, FN_C_U,
                    TP_T_L, TP_T_U, TN_T_L, TN_T_U, FP_T_L, FP_T_U, FN_T_L, FN_T_U
    ) %>%
    dplyr::mutate(Algorithm_chr = as.character(Algorithm)) %>%
    dplyr::select(-Algorithm)
  
  mk_count_row <- function(which_label, nC, lC, uC, nT, lT, uT) {
    tidyr::crossing(Algorithm_chr = base_levels) %>%
      dplyr::left_join(counts_core, by = "Algorithm_chr") %>%
      dplyr::transmute(
        Algorithm = factor(Algorithm_chr, levels = base_levels, ordered = TRUE),
        Metric    = factor(which_label, levels = counts_levels, ordered = TRUE),
        Label_C   = fmt_cnt(.data[[nC]], .data[[lC]], .data[[uC]]),
        Label_T   = fmt_cnt(.data[[nT]], .data[[lT]], .data[[uT]]),
        Bold_C    = FALSE, Bold_T = FALSE,
        Block     = "counts"
      )
  }
  
  df_counts <- dplyr::bind_rows(
    mk_count_row("TP","TP_C","TP_C_L","TP_C_U","TP_T","TP_T_L","TP_T_U"),
    mk_count_row("TN","TN_C","TN_C_L","TN_C_U","TN_T","TN_T_L","TN_T_U"),
    mk_count_row("FP","FP_C","FP_C_L","FP_C_U","FP_T","FP_T_L","FP_T_U"),
    mk_count_row("FN","FN_C","FN_C_L","FN_C_U","FN_T","FN_T_L","FN_T_U")
  )
  
  # vertical order: counts then metrics
  metric_levels <- levels(dfw$Metric) %||% unique(as.character(dfw$Metric))
  all_levels    <- c(counts_levels, metric_levels)
  df_metrics <- df_metrics %>% dplyr::mutate(Metric = factor(as.character(Metric), levels = all_levels, ordered = TRUE))
  df_counts  <- df_counts  %>% dplyr::mutate(Metric = factor(as.character(Metric), levels = all_levels, ordered = TRUE))
  
  # stub column with “Complete” top, “Triage” bottom
  stub_df <- dplyr::bind_rows(df_counts %>% dplyr::select(Metric),
                              df_metrics %>% dplyr::select(Metric)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Metric    = factor(as.character(Metric), levels = all_levels, ordered = TRUE),
      Algorithm = factor("Feature set", levels = c("Feature set", base_levels), ordered = TRUE),
      StubTop   = "Complete",
      StubBottom= "Triage"
    )
  
  # grid separators
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
                       nudge_y = +0.18, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_metrics,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_T,
                                    fontface = ifelse(Bold_T, "bold", "plain")),
                       nudge_y = -0.18, size = 2.9, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_C),
                       nudge_y = +0.18, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = df_counts,
                       ggplot2::aes(x = Algorithm, y = Metric, label = Label_T),
                       nudge_y = -0.18, size = 3.0, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubTop),
                       nudge_y = +0.18, size = 3.1, lineheight = 0.95, family = family) +
    ggplot2::geom_text(data = stub_df, ggplot2::aes(x = Algorithm, y = Metric, label = StubBottom),
                       nudge_y = -0.18, size = 3.1, lineheight = 0.95, family = family) +
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

paired_external_comparison <- function(df, yvar, feat_full, feat_triage,
                                       B = 2000, seed = cfg$seed_cv + 202) {
  set.seed(seed)
  algos <- names(model_specs(cfg$tune_len))
  
  res <- lapply(algos, function(algo) {
    oC <- run_holdout(df, yvar, feat_full,  algo, "Complete")
    oT <- run_holdout(df, yvar, feat_triage, algo, "Triage")
    if (is.null(oC) || is.null(oT)) return(NULL)
    
    y  <- oC$obs
    stopifnot(identical(y, oT$obs))
    pC <- as_num(oC$p); pT <- as_num(oT$p)
    prC <- oC$pred;    prT <- oT$pred
    neg <- cfg$neg_label; pos <- cfg$pos_label
    
    mC <- compute_metrics_binary(y, prC, pC)
    mT <- compute_metrics_binary(y, prT, pT)
    dAUC <- as_num(mT["AUC"]) - as_num(mC["AUC"])
    dMCC <- as_num(mT["MCC"]) - as_num(mC["MCC"])
    dF1  <- as_num(mT["F1"])  - as_num(mC["F1"])
    
    rC <- try(pROC::roc(response = factor(y, levels = c(neg, pos)), predictor = pC, quiet = TRUE), silent = TRUE)
    rT <- try(pROC::roc(response = factor(y, levels = c(neg, pos)), predictor = pT, quiet = TRUE), silent = TRUE)
    p_delong <- if (!inherits(rC, "try-error") && !inherits(rT, "try-error"))
      as.numeric(pROC::roc.test(rC, rT, method = "delong", paired = TRUE)$p.value) else NA_real_
    
    accC <- prC == y; accT <- prT == y
    tb   <- table(accC, accT)
    p_mcnemar <- try(as.numeric(stats::mcnemar.test(tb, correct = TRUE)$p.value), silent = TRUE)
    if (inherits(p_mcnemar, "try-error")) p_mcnemar <- NA_real_
    
    n <- length(y)
    boot_m <- matrix(NA_real_, nrow = B, ncol = 3,
                     dimnames = list(NULL, c("dAUC","dMCC","dF1")))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      mbC <- compute_metrics_binary(y[idx], prC[idx], pC[idx])
      mbT <- compute_metrics_binary(y[idx], prT[idx], pT[idx])
      boot_m[b,"dAUC"] <- as_num(mbT["AUC"]) - as_num(mbC["AUC"])
      boot_m[b,"dMCC"] <- as_num(mbT["MCC"]) - as_num(mbC["MCC"])
      boot_m[b,"dF1"]  <- as_num(mbT["F1"])  - as_num(mbC["F1"])
    }
    q <- function(v) stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
    ci_auc <- q(boot_m[,"dAUC"])
    ci_mcc <- q(boot_m[,"dMCC"])
    ci_f1  <- q(boot_m[,"dF1"])
    
    tibble::tibble(
      Algorithm = algo,
      delta_AUC = dAUC,  AUC_CI_L = ci_auc[1], AUC_CI_U = ci_auc[2],  p_AUC_DeLong = p_delong,
      delta_MCC = dMCC,  MCC_CI_L = ci_mcc[1], MCC_CI_U = ci_mcc[2],
      delta_F1  = dF1,   F1_CI_L  = ci_f1[1],  F1_CI_U  = ci_f1[2],
      p_McNemar = p_mcnemar
    )
  })
  
  dplyr::bind_rows(Filter(Negate(is.null), res))
}


# ========================= DRIVER =========================
set.seed(123)

ld  <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
df0 <- engineer_features(ld$df)
df  <- enforce_temporal_split(df0, admit_col = ld$admit_col, boundary_date = cfg$boundary_date)
yvar <- cfg$outcome

ext_df <- external_summary_for_heatmap(df, yvar, feat_full, feat_triage)
p_ext  <- plot_external_heatmap(ext_df)
print(p_ext)

paired_df <- paired_external_comparison(df, yvar, feat_full, feat_triage)
print(paired_df)

enn_sum <- enn_log_summary()
print(enn_sum)

si <- utils::sessionInfo()
print(si)
