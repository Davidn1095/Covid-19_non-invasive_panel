# ============================================================
# CONCHI BENCHMARK (EXTERNAL RESULTS ONLY) — CASCADE POLICY — 20260107
# - No file writing, prints tables to stdout
# - Development used for OOF threshold selection only
# - External used once at the end for evaluation
# - Algorithms: LR, RF, SVMRBF, kNN, C4.5
# - Cascade stages: Non-invasive, Laboratory augmented
# - kNN fix: store training data and call kknn::kknn at predict time
# ============================================================

suppressWarnings(options(warn = -1))
options(na.action = "na.pass")

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(dplyr)
  library(pROC)
  library(randomForest)
  library(kernlab)
  library(kknn)
  library(rpart)
})

source(file.path("R", "common_utils.R"))
source(file.path("R", "utils_metrics.R"))
source(file.path("R", "utils_thresholds.R"))
source(file.path("R", "cascade_pipeline.R"))

# =========================
# USER SETTINGS
# =========================
dev_csv <- Sys.getenv("DEV_CSV", unset = "20260101_development_mimic_raw_rowcomplete60.csv")
ext_csv <- Sys.getenv("EXT_CSV", unset = "20260101_external_mcmed_raw_rowcomplete60.csv")

seed_cv   <- 123
k_folds   <- 5L
repeats   <- 2L

dev_max_n <- 10000L
ext_max_n <- 10000L
B_boot    <- 500L

# Global row completeness filter, applied after coercion, before sampling
row_thr_global <- 0.80

# Fixed preprocessing + calibration
clip_q      <- 0.01
standardise <- TRUE
use_weights <- TRUE
filter_rate <- 0.00
calibration <- "isotonic"  # "none" or "isotonic"

# Threshold search grid
thr_grid <- seq(0.01, 0.99, by = 0.01)

# Threshold selection constraints on DEV OOF
sens_min_target <- -Inf
spec_min_target <- -Inf
max_def_rate_target <- 0.30

# If TRUE: enforce Sensitivity >= Specificity when selecting cascade params on DEV OOF
require_sens_ge_spec <- TRUE

set.seed(seed_cv)
stopifnot(file.exists(dev_csv))
stopifnot(file.exists(ext_csv))

df_dev_full <- readr::read_csv(dev_csv, show_col_types = FALSE) |> prep_df_any("Development")
df_ext_full <- readr::read_csv(ext_csv, show_col_types = FALSE) |> prep_df_any("External")

# =========================
# 3) GLOBAL ROW COMPLETENESS FILTER, THEN STRATIFIED SAMPLE
# =========================
df_dev <- apply_global_filter_then_sample(df_dev_full, row_thr = row_thr_global, n_max = dev_max_n, seed = seed_cv + 10001L)
df_ext <- apply_global_filter_then_sample(df_ext_full, row_thr = row_thr_global, n_max = ext_max_n, seed = seed_cv + 20001L)

if (nrow(df_dev) < 100L) stop("Too few development rows after global filter and sampling")
if (nrow(df_ext) < 100L) stop("Too few external rows after global filter and sampling")

# =========================
# 4) CV SPLITS ON DEVELOPMENT ONLY
# =========================
folds <- make_repeated_stratified_folds(df_dev$outcome, k = k_folds, repeats = repeats, seed = seed_cv)

# =========================
# 12) FIT ON FULL DEV AND PREDICT EXTERNAL, SINGLE PANEL
# =========================
fit_full_and_predict_panel <- function(df_dev, df_ext, algo, predictors,
                                       clip_q, standardise, use_weights, filter_rate, calibration,
                                       flip_from_dev = FALSE) {
  xy_dev <- make_Xy(df_dev, predictors)
  xy_ext <- make_Xy(df_ext, predictors)

  X_dev_raw <- xy_dev$X
  y_dev01 <- xy_dev$y01
  X_ext_raw <- xy_ext$X

  pp <- prep_fit(X_dev_raw, clip_q = clip_q, standardise = standardise)
  X_dev <- prep_apply(X_dev_raw, pp)
  X_ext <- prep_apply(X_ext_raw, pp)

  w_case <- if (isTRUE(use_weights)) make_case_weights(y_dev01) else NULL

  nf <- apply_noise_filter(X_dev, y_dev01, filter_rate = filter_rate, seed = seed_cv + 777)
  X_dev2 <- nf$X
  y_dev2 <- nf$y01
  w_case2 <- if (!is.null(w_case)) w_case[nf$keep] else NULL

  mod <- fit_model(algo, X_dev2, y_dev2, w_case = w_case2)
  p_ext_raw <- predict_prob(mod, X_ext)

  if (calibration == "isotonic") {
    p_dev_raw <- predict_prob(mod, X_dev2)
    cal <- fit_calibrator(p_dev_raw, y_dev2, method = calibration)
    p_ext <- apply_calibrator(cal, p_ext_raw)
  } else {
    p_ext <- p_ext_raw
  }

  if (isTRUE(flip_from_dev)) p_ext <- 1 - p_ext
  pmin(pmax(as.numeric(p_ext), 1e-6), 1 - 1e-6)
}

# =========================
# 13) DEVELOPMENT: OOF PROBS FOR STAGE 1 + STAGE 2, THEN 3 POLICIES
# =========================
algos <- c("LR", "RF", "SVMRBF", "kNN", "C4.5")
stage1_name <- "Non-invasive"
stage2_name <- "Laboratory augmented"

dev_params_rows <- list()
dev_policy_rows <- list()
dev_cache <- list()

for (algo in algos) {
  message("DEV OOF: ", algo, " (Stage 1: ", stage1_name, ")")
  preds1 <- get_panel_predictors(df_dev, stage1_name)
  res1 <- run_oof_single_panel(
    df_dev = df_dev, algo = algo, predictors = preds1,
    clip_q = clip_q, standardise = standardise, use_weights = use_weights,
    filter_rate = filter_rate, calibration = calibration,
    folds = folds
  )

  message("DEV OOF: ", algo, " (Stage 2: ", stage2_name, ")")
  preds2 <- get_panel_predictors(df_dev, stage2_name)
  res2 <- run_oof_single_panel(
    df_dev = df_dev, algo = algo, predictors = preds2,
    clip_q = clip_q, standardise = standardise, use_weights = use_weights,
    filter_rate = filter_rate, calibration = calibration,
    folds = folds
  )

  stopifnot(length(res1$y01) == length(res2$y01))

  # Policy 1: Stage 1 only
  thr1_sel <- best_thr_mcc_at_sens(
    y01 = res1$y01, p = res1$p_oof, thr_grid = thr_grid,
    sens_min = sens_min_target, spec_min = spec_min_target
  )
  met1 <- compute_metrics_at_thr(res1$y01, res1$p_oof, thr1_sel$thr)

  dev_policy_rows[[length(dev_policy_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    Policy = "Non-invasive",
    tau_low = NA_real_, tau_high = NA_real_,
    thr = thr1_sel$thr,
    def_rate = 0.0,
    hit_constraint = thr1_sel$hit_constraint,
    TP = met1$TP, TN = met1$TN, FP = met1$FP, FN = met1$FN,
    MCC = met1$MCC, AUC = met1$AUC, F1 = met1$F1,
    Accuracy = met1$Accuracy, Precision = met1$Precision,
    Sensitivity = met1$Sensitivity, Specificity = met1$Specificity
  )

  # Policy 2: Stage 2 only
  thr2_sel <- best_thr_mcc_at_sens(
    y01 = res2$y01, p = res2$p_oof, thr_grid = thr_grid,
    sens_min = sens_min_target, spec_min = spec_min_target
  )
  met2 <- compute_metrics_at_thr(res2$y01, res2$p_oof, thr2_sel$thr)

  dev_policy_rows[[length(dev_policy_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    Policy = "Laboratory augmented",
    tau_low = NA_real_, tau_high = NA_real_,
    thr = thr2_sel$thr,
    def_rate = 0.0,
    hit_constraint = thr2_sel$hit_constraint,
    TP = met2$TP, TN = met2$TN, FP = met2$FP, FN = met2$FN,
    MCC = met2$MCC, AUC = met2$AUC, F1 = met2$F1,
    Accuracy = met2$Accuracy, Precision = met2$Precision,
    Sensitivity = met2$Sensitivity, Specificity = met2$Specificity
  )

  # Policy 3: Cascade
  sel <- best_cascade_params(
    y01 = res1$y01, p1 = res1$p_oof, p2 = res2$p_oof,
    tau_grid = thr_grid, thr2_grid = thr_grid,
    sens_min = sens_min_target, spec_min = spec_min_target,
    max_def_rate = max_def_rate_target,
    require_sens_ge_spec = require_sens_ge_spec
  )

  cas <- cascade_apply(res1$p_oof, res2$p_oof, sel$tau_low, sel$tau_high, sel$thr2)
  metc <- compute_metrics_from_pred(res1$y01, cas$pred01, score = cas$score)
  def_rate <- mean(cas$defer, na.rm = TRUE)

  dev_policy_rows[[length(dev_policy_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    Policy = "Cascade",
    tau_low = sel$tau_low, tau_high = sel$tau_high,
    thr = sel$thr2,
    def_rate = def_rate,
    hit_constraint = sel$hit_constraint,
    TP = metc$TP, TN = metc$TN, FP = metc$FP, FN = metc$FN,
    MCC = metc$MCC, AUC = metc$AUC, F1 = metc$F1,
    Accuracy = metc$Accuracy, Precision = metc$Precision,
    Sensitivity = metc$Sensitivity, Specificity = metc$Specificity
  )

  dev_params_rows[[length(dev_params_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    thr_stage1 = thr1_sel$thr,
    thr_stage2 = thr2_sel$thr,
    tau_low = sel$tau_low,
    tau_high = sel$tau_high,
    thr_cascade = sel$thr2
  )

  dev_cache[[algo]] <- list(
    preds1 = preds1, flip1 = res1$flip,
    preds2 = preds2, flip2 = res2$flip
  )
}

dev_params_tbl <- dplyr::bind_rows(dev_params_rows) |>
  dplyr::arrange(.data$Algorithm)

dev_policy_tbl <- dplyr::bind_rows(dev_policy_rows) |>
  dplyr::mutate(Policy = factor(Policy, levels = c("Non-invasive","Laboratory augmented","Cascade"))) |>
  dplyr::arrange(.data$Algorithm, .data$Policy)

print(dev_policy_tbl, n = Inf, width = Inf)

# =========================
# 14) EXTERNAL EVALUATION + BOOTSTRAP CIs, 3 POLICIES
# =========================
fmt_ci_num <- function(est, lo, hi, digits = 4) {
  if (!is.finite(est) || is.na(est)) return("NA")
  if (!is.finite(lo) || is.na(lo) || !is.finite(hi) || is.na(hi)) return(sprintf(paste0("%.", digits, "f"), est))
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), est, lo, hi)
}

fmt_ci_int <- function(est, lo, hi) {
  if (!is.finite(est) || is.na(est)) return("NA")
  if (!is.finite(lo) || is.na(lo) || !is.finite(hi) || is.na(hi)) return(sprintf("%.0f", est))
  sprintf("%.0f [%.0f, %.0f]", est, lo, hi)
}

fmt_ci_rate <- function(est, lo, hi, digits = 4) {
  if (!is.finite(est) || is.na(est)) return("NA")
  if (!is.finite(lo) || is.na(lo) || !is.finite(hi) || is.na(hi)) return(sprintf(paste0("%.", digits, "f"), est))
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), est, lo, hi)
}

boot_metrics_external_single <- function(y01, p, thr, B, seed) {
  set.seed(seed)

  met0 <- compute_metrics_at_thr(y01, p, thr)
  cc0  <- confusion_counts(as.integer(y01), as.integer(as.numeric(p) >= thr))

  stat <- matrix(NA_real_, nrow = B, ncol = 7)
  colnames(stat) <- c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")

  cnt <- matrix(NA_real_, nrow = B, ncol = 4)
  colnames(cnt) <- c("TP","TN","FP","FN")

  n <- length(y01)

  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    yb <- y01[ii]
    pb <- p[ii]

    metb <- compute_metrics_at_thr(yb, pb, thr)
    ccb  <- confusion_counts(as.integer(yb), as.integer(as.numeric(pb) >= thr))

    stat[b, ] <- c(metb$MCC, metb$AUC, metb$F1, metb$Accuracy, metb$Precision, metb$Sensitivity, metb$Specificity)
    cnt[b, ]  <- c(ccb["TP"], ccb["TN"], ccb["FP"], ccb["FN"])
  }

  lo_m <- apply(stat, 2, function(x) suppressWarnings(stats::quantile(x, 0.025, na.rm = TRUE)))
  hi_m <- apply(stat, 2, function(x) suppressWarnings(stats::quantile(x, 0.975, na.rm = TRUE)))
  lo_c <- apply(cnt, 2, function(x) suppressWarnings(stats::quantile(x, 0.025, na.rm = TRUE)))
  hi_c <- apply(cnt, 2, function(x) suppressWarnings(stats::quantile(x, 0.975, na.rm = TRUE)))

  list(met0 = met0, cc0 = cc0, lo_m = lo_m, hi_m = hi_m, lo_c = lo_c, hi_c = hi_c)
}

boot_metrics_external_cascade <- function(y01, p1, p2, tau_low, tau_high, thr2, B, seed) {
  set.seed(seed)

  cas0 <- cascade_apply(p1, p2, tau_low, tau_high, thr2)
  met0 <- compute_metrics_from_pred(y01, cas0$pred01, score = cas0$score)
  def0 <- mean(cas0$defer, na.rm = TRUE)

  stat <- matrix(NA_real_, nrow = B, ncol = 8)
  colnames(stat) <- c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity","Deferral")

  cnt <- matrix(NA_real_, nrow = B, ncol = 4)
  colnames(cnt) <- c("TP","TN","FP","FN")

  n <- length(y01)

  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    yb <- y01[ii]
    p1b <- p1[ii]
    p2b <- p2[ii]

    casb <- cascade_apply(p1b, p2b, tau_low, tau_high, thr2)
    metb <- compute_metrics_from_pred(yb, casb$pred01, score = casb$score)
    defb <- mean(casb$defer, na.rm = TRUE)

    stat[b, ] <- c(metb$MCC, metb$AUC, metb$F1, metb$Accuracy, metb$Precision, metb$Sensitivity, metb$Specificity, defb)
    cnt[b, ]  <- c(metb$TP,  metb$TN,  metb$FP,  metb$FN)
  }

  lo_m <- apply(stat, 2, function(x) suppressWarnings(stats::quantile(x, 0.025, na.rm = TRUE)))
  hi_m <- apply(stat, 2, function(x) suppressWarnings(stats::quantile(x, 0.975, na.rm = TRUE)))
  lo_c <- apply(cnt, 2, function(x) suppressWarnings(stats::quantile(x, 0.025, na.rm = TRUE)))
  hi_c <- apply(cnt, 2, function(x) suppressWarnings(stats::quantile(x, 0.975, na.rm = TRUE)))

  list(met0 = met0, def0 = def0, lo_m = lo_m, hi_m = hi_m, lo_c = lo_c, hi_c = hi_c)
}

# =========================
# 15) EXTERNAL RESULTS TABLE, 3 POLICIES
# =========================
ext_rows <- list()

for (algo in algos) {
  cache <- dev_cache[[algo]]
  rowp <- dev_params_tbl[dev_params_tbl$Algorithm == algo, , drop = FALSE]
  stopifnot(nrow(rowp) == 1)

  p1_ext <- fit_full_and_predict_panel(
    df_dev = df_dev, df_ext = df_ext,
    algo = algo, predictors = cache$preds1,
    clip_q = clip_q, standardise = standardise,
    use_weights = use_weights, filter_rate = filter_rate,
    calibration = calibration,
    flip_from_dev = cache$flip1
  )

  p2_ext <- fit_full_and_predict_panel(
    df_dev = df_dev, df_ext = df_ext,
    algo = algo, predictors = cache$preds2,
    clip_q = clip_q, standardise = standardise,
    use_weights = use_weights, filter_rate = filter_rate,
    calibration = calibration,
    flip_from_dev = cache$flip2
  )

  # Stage 1 only
  bm1 <- boot_metrics_external_single(
    y01 = df_ext$outcome, p = p1_ext, thr = rowp$thr_stage1,
    B = B_boot, seed = seed_cv + 3100 + match(algo, algos)
  )

  ext_rows[[length(ext_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    Policy = "Non-invasive",
    tau_low = NA_real_, tau_high = NA_real_, thr = rowp$thr_stage1,
    Deferral = "0.0000",
    TP = fmt_ci_int(bm1$cc0["TP"], bm1$lo_c["TP"], bm1$hi_c["TP"]),
    TN = fmt_ci_int(bm1$cc0["TN"], bm1$lo_c["TN"], bm1$hi_c["TN"]),
    FP = fmt_ci_int(bm1$cc0["FP"], bm1$lo_c["FP"], bm1$hi_c["FP"]),
    FN = fmt_ci_int(bm1$cc0["FN"], bm1$lo_c["FN"], bm1$hi_c["FN"]),
    MCC = fmt_ci_num(bm1$met0$MCC, bm1$lo_m["MCC"], bm1$hi_m["MCC"], 4),
    AUC = fmt_ci_num(bm1$met0$AUC, bm1$lo_m["AUC"], bm1$hi_m["AUC"], 4),
    F1  = fmt_ci_num(bm1$met0$F1,  bm1$lo_m["F1"],  bm1$hi_m["F1"],  4),
    Accuracy    = fmt_ci_num(bm1$met0$Accuracy,    bm1$lo_m["Accuracy"],    bm1$hi_m["Accuracy"],    4),
    Precision   = fmt_ci_num(bm1$met0$Precision,   bm1$lo_m["Precision"],   bm1$hi_m["Precision"],   4),
    Sensitivity = fmt_ci_num(bm1$met0$Sensitivity, bm1$lo_m["Sensitivity"], bm1$hi_m["Sensitivity"], 4),
    Specificity = fmt_ci_num(bm1$met0$Specificity, bm1$lo_m["Specificity"], bm1$hi_m["Specificity"], 4)
  )

  # Stage 2 only
  bm2 <- boot_metrics_external_single(
    y01 = df_ext$outcome, p = p2_ext, thr = rowp$thr_stage2,
    B = B_boot, seed = seed_cv + 3200 + match(algo, algos)
  )

  ext_rows[[length(ext_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    Policy = "Laboratory augmented",
    tau_low = NA_real_, tau_high = NA_real_, thr = rowp$thr_stage2,
    Deferral = "0.0000",
    TP = fmt_ci_int(bm2$cc0["TP"], bm2$lo_c["TP"], bm2$hi_c["TP"]),
    TN = fmt_ci_int(bm2$cc0["TN"], bm2$lo_c["TN"], bm2$hi_c["TN"]),
    FP = fmt_ci_int(bm2$cc0["FP"], bm2$lo_c["FP"], bm2$hi_c["FP"]),
    FN = fmt_ci_int(bm2$cc0["FN"], bm2$lo_c["FN"], bm2$hi_c["FN"]),
    MCC = fmt_ci_num(bm2$met0$MCC, bm2$lo_m["MCC"], bm2$hi_m["MCC"], 4),
    AUC = fmt_ci_num(bm2$met0$AUC, bm2$lo_m["AUC"], bm2$hi_m["AUC"], 4),
    F1  = fmt_ci_num(bm2$met0$F1,  bm2$lo_m["F1"],  bm2$hi_m["F1"],  4),
    Accuracy    = fmt_ci_num(bm2$met0$Accuracy,    bm2$lo_m["Accuracy"],    bm2$hi_m["Accuracy"],    4),
    Precision   = fmt_ci_num(bm2$met0$Precision,   bm2$lo_m["Precision"],   bm2$hi_m["Precision"],   4),
    Sensitivity = fmt_ci_num(bm2$met0$Sensitivity, bm2$lo_m["Sensitivity"], bm2$hi_m["Sensitivity"], 4),
    Specificity = fmt_ci_num(bm2$met0$Specificity, bm2$lo_m["Specificity"], bm2$hi_m["Specificity"], 4)
  )

  # Cascade
  bm3 <- boot_metrics_external_cascade(
    y01 = df_ext$outcome, p1 = p1_ext, p2 = p2_ext,
    tau_low = rowp$tau_low, tau_high = rowp$tau_high, thr2 = rowp$thr_cascade,
    B = B_boot, seed = seed_cv + 3300 + match(algo, algos)
  )

  ext_rows[[length(ext_rows) + 1L]] <- tibble::tibble(
    Algorithm = algo,
    Policy = "Cascade",
    tau_low = rowp$tau_low, tau_high = rowp$tau_high, thr = rowp$thr_cascade,
    Deferral = fmt_ci_rate(bm3$def0, bm3$lo_m["Deferral"], bm3$hi_m["Deferral"], 4),
    TP = fmt_ci_int(bm3$met0$TP, bm3$lo_c["TP"], bm3$hi_c["TP"]),
    TN = fmt_ci_int(bm3$met0$TN, bm3$lo_c["TN"], bm3$hi_c["TN"]),
    FP = fmt_ci_int(bm3$met0$FP, bm3$lo_c["FP"], bm3$hi_c["FP"]),
    FN = fmt_ci_int(bm3$met0$FN, bm3$lo_c["FN"], bm3$hi_c["FN"]),
    MCC = fmt_ci_num(bm3$met0$MCC, bm3$lo_m["MCC"], bm3$hi_m["MCC"], 4),
    AUC = fmt_ci_num(bm3$met0$AUC, bm3$lo_m["AUC"], bm3$hi_m["AUC"], 4),
    F1  = fmt_ci_num(bm3$met0$F1,  bm3$lo_m["F1"],  bm3$hi_m["F1"],  4),
    Accuracy    = fmt_ci_num(bm3$met0$Accuracy,    bm3$lo_m["Accuracy"],    bm3$hi_m["Accuracy"],    4),
    Precision   = fmt_ci_num(bm3$met0$Precision,   bm3$lo_m["Precision"],   bm3$hi_m["Precision"],   4),
    Sensitivity = fmt_ci_num(bm3$met0$Sensitivity, bm3$lo_m["Sensitivity"], bm3$hi_m["Sensitivity"], 4),
    Specificity = fmt_ci_num(bm3$met0$Specificity, bm3$lo_m["Specificity"], bm3$hi_m["Specificity"], 4)
  )
}

external_policy_tbl <- dplyr::bind_rows(ext_rows) |>
  dplyr::mutate(Policy = factor(Policy, levels = c("Non-invasive","Laboratory augmented","Cascade"))) |>
  dplyr::arrange(.data$Algorithm, .data$Policy)

print(external_policy_tbl, n = Inf, width = Inf)
