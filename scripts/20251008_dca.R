# ===============================================================
# Read data -> engineer features -> temporal split -> train models
# -> predict external probabilities -> DCA with bootstrap ribbons
# Facet order: C4.5, k-NN, SVM, RF, LR
# Palettes: Complete = orange (#D55E00), Triage = blue (#0072B2)
# No discretization, WHO kept. C4.5 encouraged to use several features.
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(kernlab)
  library(kknn)
  library(UBL)
  library(NoiseFiltersR)
  library(RWeka)
  library(randomForest)
  library(ggplot2)
  library(patchwork)
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
  seed_cv         = 123,
  tune_len        = 5,
  calibration          = "platt_ridge",
  platt_ridge_lambda   = 1e-3,
  cv_blocked_by_time   = TRUE
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

# desired algorithm order in facets
algo_order <- c("C4.5","k-NN","SVM","RF","LR")

# ---------------- helpers ----------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))
sanitize_name <- function(x) gsub("[^A-Za-z0-9]", "", tolower(x))
clip01 <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

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

# ---------------- data loader ----------------
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

# ---------------- recipes and model specs ----------------
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

model_specs <- function(tune_len, algos = c("C4.5","k-NN","SVM","RF","LR")) {
  list(
    "C4.5" = list(
      name   = "C4.5",
      method = "J48",
      grid   = expand.grid(
        C = c(0.25, 0.35, 0.45, 0.50),
        M = c(1, 2, 5)
      ),
      args   = list(control = RWeka::Weka_control(A = TRUE, B = TRUE, S = TRUE))
    ),
    "k-NN" = list(name = "k-NN", method = "kknn",      tuneLength = tune_len),
    "SVM"  = list(name = "SVM",  method = "svmRadial", tuneLength = tune_len),
    "RF"   = list(name = "RF",   method = "rf",        tuneLength = tune_len),
    "LR"   = list(
      name   = "LR",
      method = "glmnet",
      grid   = expand.grid(
        alpha  = 0,
        lambda = 10^seq(-7, 0, length.out = 60)
      ),
      args   = list(family = "binomial", standardize = FALSE)
    )
  )[algos]
}

# ---------------- CV folds and tuning metric ----------------
mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(0)
  num/den
}
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
    return(c(MCC = mcc_from_counts(TP, FP, FN, TN)))
  }
  br <- best_threshold_mcc_bacc_med(obs, p)
  c(MCC = as.numeric(br$mcc))
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

# ---------------- ENN + ADASYN sampler ----------------
apply_ENN <- function(df, yname = ".y", k = 3) {
  df[[yname]] <- factor(df[[yname]], levels = c(cfg$neg_label, cfg$pos_label))
  out <- try(NoiseFiltersR::ENN(stats::as.formula(paste(yname, "~ .")), df, k = k), silent = TRUE)
  if (!inherits(out, "try-error") && !is.null(out$cleanData)) out$cleanData else df
}
noise_adasyn_sampler <- function(x, y) {
  y <- factor(y, levels = c(cfg$neg_label, cfg$pos_label))
  if (length(unique(y)) < 2) return(list(x = x, y = y))
  df <- as.data.frame(x); df$.y <- y
  df2 <- apply_ENN(df, ".y", k = 3)
  ad <- try(UBL::ADASYN(.y ~ ., df2, C.perc = "balance", k = 5, dist = "HEOM"), silent = TRUE)
  res <- if (!inherits(ad, "try-error") && length(unique(df2$.y)) > 1) ad else df2
  res$.y <- factor(res$.y, levels = c(cfg$neg_label, cfg$pos_label))
  list(x = res[, setdiff(names(res), ".y"), drop = FALSE], y = res$.y)
}

# ---------------- calibration from CV-OOF ----------------
get_oof_bestTune <- function(fit, pos_label = cfg$pos_label) {
  pd <- try(fit$pred, silent = TRUE)
  if (inherits(pd, "try-error") || is.null(pd) || !nrow(pd)) return(NULL)
  bt <- fit$bestTune; if (is.null(bt) || !nrow(bt)) return(NULL)
  for (nm in names(bt)) pd <- pd[pd[[nm]] == bt[[nm]][1], , drop = FALSE]
  if (!(pos_label %in% names(pd))) return(NULL)
  tibble::tibble(obs = pd$obs, p = as_num(pd[[pos_label]])) %>% dplyr::filter(!is.na(p))
}
fit_platt_ridge_from_oof <- function(y, p, lambda = cfg$platt_ridge_lambda %||% 1e-3) {
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
  fit_platt_ridge_from_oof(oof$obs, oof$p)
}

# ---------------- training and external predictions ----------------
train_inner <- function(dat_train, yvar, spec, inner_k, cv_idx, seed_base) {
  dat_train <- enforce_y_levels(dat_train, yvar)
  rec <- make_recipe(dat_train, yvar, method = spec$method)
  tr_ctrl <- caret::trainControl(
    method = "cv",
    number = inner_k,
    classProbs = TRUE,
    summaryFunction = mccSummaryAdaptive,
    savePredictions = "final",
    allowParallel = TRUE,
    index = cv_idx$index,
    indexOut = cv_idx$indexOut,
    sampling = noise_adasyn_sampler
  )
  base_args <- list(
    x = rec, data = dat_train, method = spec$method,
    trControl = tr_ctrl, metric = "MCC"
  )
  if (!is.null(spec$grid)) base_args$tuneGrid <- spec$grid else base_args$tuneLength <- spec$tuneLength
  if (!is.null(spec$args)) base_args <- c(base_args, spec$args)
  set.seed(seed_base)
  do.call(caret::train, base_args)
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
  seed_base <- cfg$seed_cv + sum(utf8ToInt(paste(algo_name, set_label, sep = "_")))
  cv_idx <- make_timeblock_indices(dtrain_with_date, k = cfg$inner_k)
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
  neg  <- cfg$neg_label
  pred <- factor(ifelse(p_ex >= thr, cfg$pos_label, neg), levels = c(neg, cfg$pos_label))
  list(
    obs = factor(dext[[yvar]], levels = c(neg, cfg$pos_label)),
    p = p_ex, pred = pred, algo = algo_name, set = set_label
  )
}

# ---------------- DCA helpers ----------------
decision_curve_table <- function(obs, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  y_raw <- as.integer(obs == cfg$pos_label)
  keep  <- is.finite(as.numeric(p)) & !is.na(y_raw)
  y <- y_raw[keep]; p <- as.numeric(p)[keep]
  if (!length(y)) return(tibble(threshold=thresholds, Net_Benefit=NA_real_, NB_TreatAll=NA_real_, NB_TreatNone=0, prevalence=NA_real_))
  N <- length(y); prev <- mean(y)
  purrr::map_dfr(thresholds, function(pt){
    pred <- as.integer(p >= pt)
    TP <- sum(pred==1 & y==1)
    FP <- sum(pred==1 & y==0)
    NB <- TP/N - FP/N * (pt/(1-pt))
    NB_all <- prev - (1 - prev) * (pt/(1-pt))
    tibble(threshold=pt, Net_Benefit=NB, NB_TreatAll=NB_all, NB_TreatNone=0, prevalence=prev)
  })
}

dca_boot_band <- function(obj, B = 500){
  if (is.null(obj)) return(NULL)
  y <- as.integer(obj$obs == cfg$pos_label)
  p <- as.numeric(obj$p)
  th <- seq(0.01, 0.99, by = 0.01)
  N  <- length(y)
  nb_mat <- matrix(NA_real_, nrow = B, ncol = length(th))
  for (b in seq_len(B)){
    idx <- sample.int(N, N, replace = TRUE)
    yb <- y[idx]; pb <- p[idx]
    nb_mat[b,] <- vapply(th, function(pt){
      pred <- as.integer(pb >= pt)
      TP <- sum(pred == 1 & yb == 1)
      FP <- sum(pred == 1 & yb == 0)
      TP/N - FP/N * (pt/(1-pt))
    }, numeric(1))
  }
  qs <- apply(nb_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  tibble(
    threshold = th,
    NB_L = qs[1,], NB_U = qs[2,],
    set = obj$set, algo = obj$algo
  )
}

# ---------------- DRIVER: read, split, train, DCA ----------------
set.seed(cfg$seed_cv)

ld  <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
df0 <- engineer_features(ld$df)
df  <- enforce_temporal_split(df0, admit_col = ld$admit_col, boundary_date = cfg$boundary_date)
yvar <- cfg$outcome

algos_avail <- names(model_specs(cfg$tune_len))
get_hold <- function(a, set_feats, set_label) run_holdout(df, yvar, set_feats, a, set_label)

holds_external <- purrr::compact(unlist(lapply(algos_avail, function(a) {
  list(
    full  = get_hold(a, feat_full,   "Complete feature set"),
    tri   = get_hold(a, feat_triage, "Triage feature set")
  )
}), recursive = FALSE))

if (!length(holds_external)) stop("No external predictions available for DCA")
by_algo <- split(holds_external, vapply(holds_external, function(x) x$algo, character(1)))

combo <- lapply(names(by_algo), function(a) {
  hx <- by_algo[[a]]
  full <- purrr::detect(hx, ~ .x$set == "Complete feature set")
  tri  <- purrr::detect(hx, ~ .x$set == "Triage feature set")
  d_full <- if (!is.null(full)) decision_curve_table(full$obs, full$p) %>% dplyr::mutate(set = full$set, algo = a) else NULL
  d_tri  <- if (!is.null(tri))  decision_curve_table(tri$obs,  tri$p)  %>% dplyr::mutate(set = tri$set,  algo = a) else NULL
  ref <- if (!is.null(full)) full else tri
  t_df <- decision_curve_table(ref$obs, ref$p) %>% dplyr::transmute(algo = a, threshold, NB_TreatAll, NB_TreatNone)
  list(dca = dplyr::bind_rows(d_full, d_tri), treat = t_df)
})

dca_df   <- dplyr::bind_rows(purrr::map(combo, "dca"))
treat_df <- dplyr::bind_rows(purrr::map(combo, "treat"))

bands_df <- holds_external %>%
  purrr::map(dca_boot_band) %>%
  purrr::compact() %>%
  { if (!length(.)) NULL else dplyr::bind_rows(.) }

valid_algos <- algo_order
clean_algo <- function(x) x %>% dplyr::filter(!is.na(algo), algo %in% valid_algos)
dca_df   <- clean_algo(dca_df)
treat_df <- clean_algo(treat_df)
if (!is.null(bands_df)) bands_df <- clean_algo(bands_df)

lvl <- valid_algos[valid_algos %in% unique(as.character(dca_df$algo))]
dca_df$algo   <- factor(as.character(dca_df$algo),   levels = lvl, ordered = TRUE)
treat_df$algo <- factor(as.character(treat_df$algo), levels = lvl, ordered = TRUE)
if (!is.null(bands_df)) bands_df$algo <- factor(as.character(bands_df$algo), levels = lvl, ordered = TRUE)
dca_df$set    <- factor(dca_df$set,    levels = c("Complete feature set","Triage feature set"))
if (!is.null(bands_df)) bands_df$set  <- factor(bands_df$set,  levels = levels(dca_df$set))

y_cap <- max(c(dca_df$Net_Benefit, treat_df$NB_TreatAll %||% 0, 0), na.rm = TRUE)
if (!is.finite(y_cap) || y_cap <= 0) y_cap <- 0.05

# palettes
cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"

# --- axis breaks and theme ---
axis_breaks_y <- function(y_cap) {
  if (y_cap <= 0.06) seq(0, y_cap, by = 0.01) else scales::breaks_extended()(c(0, y_cap))
}
theme_ref <- function(base_size = 15) {
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

# --- CI ribbons drawn safely with fixed fills ---
y0 <- -0.02
safe_ribbon <- function(df_b) {
  if (is.null(df_b) || !nrow(df_b)) return(list())
  df_b <- df_b %>% dplyr::filter(is.finite(NB_L), is.finite(NB_U), !is.na(set))
  layers <- list()
  for (s in c("Complete feature set","Triage feature set")) {
    d <- dplyr::filter(df_b, set == s)
    if (nrow(d)) {
      layers[[length(layers) + 1]] <- ggplot2::geom_ribbon(
        data = d,
        ggplot2::aes(threshold, ymin = NB_L, ymax = NB_U),
        fill = if (s == "Complete feature set") cb_orange else cb_blue,
        alpha = 0.18,
        inherit.aes = FALSE
      )
    }
  }
  layers
}

# --- single-panel builder (no facets) ---
dca_panel <- function(a, show_y = FALSE) {
  df_d <- dplyr::filter(dca_df, algo == a, is.finite(Net_Benefit))
  df_t <- dplyr::filter(treat_df, algo == a, is.finite(NB_TreatAll), is.finite(NB_TreatNone))
  df_b <- if (!is.null(bands_df)) dplyr::filter(bands_df, algo == a) else NULL
  
  ggplot2::ggplot() +
    safe_ribbon(df_b) +
    ggplot2::geom_line(
      data = df_t,
      ggplot2::aes(threshold, NB_TreatAll, linetype = "Treat-all"),
      linewidth = 0.6, color = "black", na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = df_t,
      ggplot2::aes(threshold, NB_TreatNone, linetype = "Treat-none"),
      linewidth = 0.6, color = "grey40", na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = df_d,
      ggplot2::aes(threshold, Net_Benefit, color = set),
      linewidth = 0.9, na.rm = TRUE
    ) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0, y = y0, yend = y_cap),
                          inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = y0, yend = y0),
                          inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    ggplot2::scale_color_manual(values = c("Complete feature set" = cb_orange,
                                           "Triage feature set"   = cb_blue), name = "") +
    ggplot2::scale_linetype_manual(values = c("Treat-all" = "solid", "Treat-none" = "dashed"), name = "") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    ggplot2::scale_y_continuous(breaks = axis_breaks_y(y_cap),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::coord_cartesian(ylim = c(y0, y_cap)) +
    ggplot2::labs(x = "Threshold probability", y = if (show_y) "Net benefit" else NULL, title = a) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA)) +
    theme_ref()
}

# --- assemble: 3 top, 2 bottom, legend in bottom-right tile ---
algos_in_data <- levels(dca_df$algo)

top <- list(
  dca_panel(algos_in_data[1], show_y = TRUE)  + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
  dca_panel(algos_in_data[2], show_y = FALSE) + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
  dca_panel(algos_in_data[3], show_y = FALSE) + ggplot2::theme(axis.title.x = ggplot2::element_blank())
)
bottom <- list(
  dca_panel(algos_in_data[4], show_y = TRUE),
  dca_panel(algos_in_data[5], show_y = FALSE)
)

p_dca_grid <-
  wrap_plots(
    plotlist = c(top, bottom, list(guide_area())),
    ncol = 3
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "center", legend.direction = "horizontal")

print(p_dca_grid)

# save PNG
ggsave(
  filename = "dca_grid.png",
  plot     = p_dca_grid,
  width    = 10,
  height   = 7.2,
  units    = "in",
  dpi      = 320,
  bg       = "white"
)
