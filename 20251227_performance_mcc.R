# ============================================================
# FINAL BENCHMARK (RAW CSVs) + MCC TUNING + GRID SEARCH (t_low, t_high)
# - Development: MIMIC  (development_*.csv)
# - External:    MC-MED (external_*.csv)
# Panels:
#   - Non-invasive
#   - Laboratory augmented
#   - Cascade: GRID-OPT (t_low, t_high) on TRAIN OOF + MCC-opt cutoffs (thrT, thrC)
# Outputs:
#   - External test metrics with bootstrap 95% percentile confidence intervals (NOW incl. Specificity)
#   - Best (t_low, t_high) per algorithm (tuned on TRAIN OOF, no leakage)
#   - Optional: External paired bootstrap deltas among the 3 policies for all metrics
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  
  library(caret)
  library(glmnet)
  library(pROC)
  library(kknn)
  library(randomForest)
  library(kernlab)
  library(RWeka)
})

options(warn = -1)
options(na.action = "na.pass")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# =========================
# USER SETTINGS
# =========================
max_n_per_dataset <- 5000L   # set 500L for quick checks
B_boot            <- 500L    # bootstrap replicates (external metrics CIs)
use_class_weights <- TRUE
use_calibration   <- FALSE  # placeholder; not used by default

cascade_auc_convention <- "mixture"  # "mixture" or "complete"

# ---- GRID SEARCH settings for (t_low, t_high)
grid_low_min  <- 0.05
grid_low_max  <- 0.45
grid_high_min <- 0.55
grid_high_max <- 0.95
grid_step     <- 0.01

min_gap   <- 0.05     # enforce t_high - t_low >= min_gap
defer_max <- 0.30     # cap defer rate, set NULL for no constraint

cfg <- list(
  pos_label = "Yes",
  neg_label = "No",
  seed_cv   = 123
)

set.seed(cfg$seed_cv)

message(
  ">>> FINAL BENCHMARK (RAW): max_n_per_dataset=", max_n_per_dataset,
  "; bootstrap B=", B_boot,
  "; class_weights=", use_class_weights,
  "; calibration=", use_calibration,
  "; cascade AUC=", cascade_auc_convention,
  "; GRID (t_low,t_high) step=", grid_step,
  "; min_gap=", min_gap,
  "; defer_max=", (defer_max %||% "None"),
  " <<<"
)

# =========================
# 0) FIND INPUT FILES (NO GUESSING)
# =========================
csvs <- list.files(getwd(), pattern = "\\.csv$", ignore.case = TRUE)
if (length(csvs) == 0L) stop("No .csv files were found in current working directory: ", getwd())

pick_single <- function(hits, label) {
  if (length(hits) == 0L) {
    stop(
      "Could not find ", label, " CSV in: ", getwd(),
      "\nAvailable CSVs:\n- ", paste(csvs, collapse = "\n- ")
    )
  }
  if (length(hits) > 1L) {
    stop(
      "Ambiguous matches for ", label, ":\n- ", paste(hits, collapse = "\n- "),
      "\nRename files so exactly one matches."
    )
  }
  hits[[1]]
}

dev_path <- pick_single(
  grep("development.*mimic|mimic.*development", csvs, value = TRUE, ignore.case = TRUE),
  "development (MIMIC)"
)
ext_path <- pick_single(
  grep("external.*mcmed|mcmed.*external", csvs, value = TRUE, ignore.case = TRUE),
  "external (MC-MED)"
)

message("Using dev_path = ", dev_path)
message("Using ext_path = ", ext_path)

# =========================
# 1) LOAD + BASIC PREP
# =========================
dev_raw <- readr::read_csv(dev_path, show_col_types = FALSE)
ext_raw <- readr::read_csv(ext_path, show_col_types = FALSE)

prep_outcome <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  
  if (is.numeric(x) || is.integer(x)) {
    x <- ifelse(is.na(x), NA_character_, ifelse(x == 1, cfg$pos_label, cfg$neg_label))
  } else {
    x0 <- trimws(tolower(as.character(x)))
    x  <- ifelse(
      x0 %in% c("1","yes","y","true","t","pos","positive","case","covid","covid19"),
      cfg$pos_label,
      ifelse(
        x0 %in% c("0","no","n","false","f","neg","negative","control","ctrl","nocovid"),
        cfg$neg_label,
        NA_character_
      )
    )
  }
  
  # IMPORTANT: POS first (caret summary uses lev[1] as event)
  factor(x, levels = c(cfg$pos_label, cfg$neg_label))
}

coerce_numeric_robust <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x1 <- as.character(x)
  x1 <- trimws(x1)
  x1 <- gsub(",", ".", x1, fixed = TRUE)
  x1 <- gsub("[^0-9eE\\+\\-\\.]", "", x1)
  suppressWarnings(as.numeric(x1))
}

coerce_sex_to_01 <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x0 <- trimws(tolower(as.character(x)))
  out <- rep(NA_real_, length(x0))
  out[x0 %in% c("m","male","man","1")] <- 1
  out[x0 %in% c("f","female","woman","0")] <- 0
  out
}

harmonize_units <- function(df, label) {
  out <- df
  
  # Temperature: Fahrenheit-like values -> Celsius if most > 60
  if ("temp_raw" %in% names(out)) {
    x <- coerce_numeric_robust(out$temp_raw)
    frac_gt60 <- mean(is.finite(x) & !is.na(x) & x > 60, na.rm = TRUE)
    if (is.finite(frac_gt60) && frac_gt60 > 0.70) {
      out$temp_raw <- (x - 32) * 5/9
      message("[", label, "] temp_raw: Fahrenheit-like values detected (", round(100*frac_gt60, 1),
              "% > 60), converted to Celsius.")
    } else {
      out$temp_raw <- x
    }
  }
  
  # Creatinine: µmol/L -> mg/dL if values look like 40–200 range
  if ("creatinine_raw" %in% names(out)) {
    x <- coerce_numeric_robust(out$creatinine_raw)
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    p95 <- suppressWarnings(stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE))
    if (is.finite(med) && is.finite(p95) && med > 20 && p95 > 40) {
      out$creatinine_raw <- x / 88.4
      message("[", label, "] creatinine_raw: median=", round(med, 2), ", p95=", round(p95, 2),
              ", micromole/L suspected, converted to mg/dL by /88.4.")
    } else {
      out$creatinine_raw <- x
    }
  }
  
  # Bilirubin: µmol/L -> mg/dL if values look like 2–30
  if ("bilirubin_raw" %in% names(out)) {
    x <- coerce_numeric_robust(out$bilirubin_raw)
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    p95 <- suppressWarnings(stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE))
    if (is.finite(med) && is.finite(p95) && med > 2 && p95 > 4) {
      out$bilirubin_raw <- x / 17.1
      message("[", label, "] bilirubin_raw: median=", round(med, 2), ", p95=", round(p95, 2),
              ", micromole/L suspected, converted to mg/dL by /17.1.")
    } else {
      out$bilirubin_raw <- x
    }
  }
  
  out
}

prep_df <- function(df, label) {
  if (!("outcome" %in% names(df))) stop("Missing column 'outcome' in ", label)
  
  names(df) <- trimws(names(df))
  
  df$outcome <- prep_outcome(df$outcome)
  if (any(is.na(df$outcome))) stop("Parsed NA in outcome for ", label, ". Inspect raw 'outcome' values.")
  
  raw_cols <- c(
    "age_raw","sex_raw","spo2_raw","heart_rate_raw","resp_rate_raw","temp_raw",
    "sbp_raw","dbp_raw","map_raw",
    "wbc_raw","neutrophils_raw","lymphocytes_raw","monocytes_raw","platelets_raw",
    "hemoglobin_raw","albumin_raw","creatinine_raw","urea_raw",
    "sodium_raw","potassium_raw","chloride_raw","anion_gap_raw","glucose_raw",
    "bilirubin_raw","alk_phos_raw","crp_raw","d_dimer_raw","ast_raw","alt_raw","lactate_raw","troponin_raw"
  )
  
  for (nm in intersect(raw_cols, names(df))) {
    if (nm == "sex_raw") df[[nm]] <- coerce_sex_to_01(df[[nm]]) else df[[nm]] <- coerce_numeric_robust(df[[nm]])
  }
  
  harmonize_units(df, label)
}

df_dev_full <- prep_df(dev_raw, "Development")
df_ext_full <- prep_df(ext_raw, "External")

# =========================
# 1b) RANDOM SUBSAMPLING
# =========================
sample_rows <- function(df, n, seed) {
  set.seed(seed)
  if (n >= nrow(df)) return(df)
  idx <- sample.int(nrow(df), size = n, replace = FALSE)
  df[idx, , drop = FALSE]
}

n_dev <- min(max_n_per_dataset, nrow(df_dev_full))
n_ext <- min(max_n_per_dataset, nrow(df_ext_full))

df_train <- sample_rows(df_dev_full, n_dev, cfg$seed_cv + 1001)
df_test  <- sample_rows(df_ext_full, n_ext, cfg$seed_cv + 2001)

df_train$outcome <- factor(as.character(df_train$outcome), levels = c(cfg$pos_label, cfg$neg_label))
df_test$outcome  <- factor(as.character(df_test$outcome),  levels = c(cfg$pos_label, cfg$neg_label))

if (length(unique(df_test$outcome)) < 2L) {
  stop("External test had fewer than 2 outcome levels after sampling. Increase max_n_per_dataset.")
}

# =========================
# 2) ALIGN + IMPUTE
# =========================
align_test_to_train <- function(train, test) {
  miss <- setdiff(names(train), names(test))
  for (m in miss) test[[m]] <- NA
  test <- test[, names(train), drop = FALSE]
  test
}

median_impute_numeric <- function(train, test) {
  num_cols <- names(train)[vapply(train, is.numeric, logical(1))]
  for (nm in num_cols) {
    med <- suppressWarnings(stats::median(train[[nm]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0.0
    if (anyNA(train[[nm]])) train[[nm]][is.na(train[[nm]])] <- med
    if (anyNA(test[[nm]]))  test[[nm]][is.na(test[[nm]])]   <- med
  }
  list(train = train, test = test)
}

fill_remaining_numeric_na <- function(train, test, fill_value = 0.0) {
  num_cols <- names(train)[vapply(train, is.numeric, logical(1))]
  for (nm in num_cols) {
    train[[nm]][is.na(train[[nm]])] <- fill_value
    test[[nm]][is.na(test[[nm]])]   <- fill_value
  }
  list(train = train, test = test)
}

df_test <- align_test_to_train(df_train, df_test)
imp1 <- median_impute_numeric(df_train, df_test); df_train <- imp1$train; df_test <- imp1$test
fix1 <- fill_remaining_numeric_na(df_train, df_test, 0.0); df_train <- fix1$train; df_test <- fix1$test

# =========================
# 3) ENGINEERED FEATURES
# =========================
safe_ratio <- function(num, den, min_den = 1e-6) {
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & !is.na(num) & !is.na(den) & (abs(den) >= min_den)
  out[ok] <- num[ok] / den[ok]
  out
}

safe_diff <- function(a, b) {
  out <- rep(NA_real_, length(a))
  ok  <- is.finite(a) & is.finite(b) & !is.na(a) & !is.na(b)
  out[ok] <- a[ok] - b[ok]
  out
}

safe_log1p_pos <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & !is.na(x) & (x >= 0)
  out[ok] <- log1p(x[ok])
  out
}

add_features <- function(df) {
  created <- character(0)
  skipped <- character(0)
  
  names(df) <- trimws(names(df))
  
  pick_ci <- function(...) {
    cands <- as.character(c(...))
    idx <- which(tolower(names(df)) %in% tolower(cands))
    if (length(idx) == 0L) return(NULL)
    names(df)[idx[1]]
  }
  
  add_one <- function(new_name, expr) {
    if (new_name %in% names(df)) return(invisible(NULL))
    v <- tryCatch(expr, error = function(e) rep(NA_real_, nrow(df)))
    if (all(is.na(v))) {
      skipped <<- c(skipped, new_name)
    } else {
      df[[new_name]] <<- as.numeric(v)
      created <<- c(created, new_name)
    }
    invisible(NULL)
  }
  
  neut <- pick_ci("neutrophils_raw","neutrophils")
  lymph <- pick_ci("lymphocytes_raw","lymphocytes")
  hr   <- pick_ci("heart_rate_raw","heartrate_raw","hr_raw","heart_rate","hr")
  rr   <- pick_ci("resp_rate_raw","resprate_raw","rr_raw","resp_rate","rr")
  spo2 <- pick_ci("spo2_raw","spo2")
  sbp  <- pick_ci("sbp_raw","sbp")
  dbp  <- pick_ci("dbp_raw","dbp")
  mapc <- pick_ci("map_raw","map")
  
  urea  <- pick_ci("urea_raw","urea","bun_raw","bun")
  creat <- pick_ci("creatinine_raw","creatinine","creat_raw","creat")
  na_   <- pick_ci("sodium_raw","sodium","na_raw","na")
  k_    <- pick_ci("potassium_raw","potassium","k_raw","k")
  cl_   <- pick_ci("chloride_raw","chloride","cl_raw","cl")
  ag_   <- pick_ci("anion_gap_raw","anion_gap","ag_raw","ag")
  glu   <- pick_ci("glucose_raw","glucose","glu_raw","glu")
  bili  <- pick_ci("bilirubin_raw","bilirubin","bili_raw","bili")
  alk   <- pick_ci("alk_phos_raw","alk_phos","alkphos_raw","alp_raw","alp")
  alb   <- pick_ci("albumin_raw","albumin","alb_raw","alb")
  
  add_one("NLR",           if (!is.null(neut) && !is.null(lymph)) safe_ratio(df[[neut]], df[[lymph]]) else rep(NA_real_, nrow(df)))
  add_one("ShockIndex",    if (!is.null(hr) && !is.null(sbp))     safe_ratio(df[[hr]], df[[sbp]])     else rep(NA_real_, nrow(df)))
  add_one("PulsePressure", if (!is.null(sbp) && !is.null(dbp))    safe_diff(df[[sbp]], df[[dbp]])     else rep(NA_real_, nrow(df)))
  add_one("MAP_est",       if (!is.null(sbp) && !is.null(dbp))    (df[[sbp]] + 2*df[[dbp]])/3         else rep(NA_real_, nrow(df)))
  add_one("MAP_delta",     if (!is.null(mapc) && ("MAP_est" %in% names(df))) safe_diff(df[[mapc]], df[["MAP_est"]]) else rep(NA_real_, nrow(df)))
  add_one("RR_SpO2",       if (!is.null(rr) && !is.null(spo2))    safe_ratio(df[[rr]], df[[spo2]])     else rep(NA_real_, nrow(df)))
  
  add_one("UreaCreatRatio", if (!is.null(urea) && !is.null(creat)) safe_ratio(df[[urea]], df[[creat]]) else rep(NA_real_, nrow(df)))
  add_one("NaKRatio",       if (!is.null(na_)  && !is.null(k_))    safe_ratio(df[[na_]],  df[[k_]])    else rep(NA_real_, nrow(df)))
  add_one("ClNaRatio",      if (!is.null(cl_)  && !is.null(na_))   safe_ratio(df[[cl_]],  df[[na_]])   else rep(NA_real_, nrow(df)))
  add_one("AGClRatio",      if (!is.null(ag_)  && !is.null(cl_))   safe_ratio(df[[ag_]],  df[[cl_]])   else rep(NA_real_, nrow(df)))
  add_one("GluUreaRatio",   if (!is.null(glu)  && !is.null(urea))  safe_ratio(df[[glu]],  df[[urea]])  else rep(NA_real_, nrow(df)))
  add_one("BiliAlkRatio",   if (!is.null(bili) && !is.null(alk))   safe_ratio(df[[bili]], df[[alk]])   else rep(NA_real_, nrow(df)))
  
  add_one("PNI", if (!is.null(alb) && !is.null(lymph)) df[[alb]] + 5 * df[[lymph]] else rep(NA_real_, nrow(df)))
  
  if (!is.null(creat)) add_one("log1p_creatinine", safe_log1p_pos(df[[creat]])) else skipped <- c(skipped, "log1p_creatinine")
  if (!is.null(urea))  add_one("log1p_urea",       safe_log1p_pos(df[[urea]]))  else skipped <- c(skipped, "log1p_urea")
  if (!is.null(bili))  add_one("log1p_bilirubin",  safe_log1p_pos(df[[bili]]))  else skipped <- c(skipped, "log1p_bilirubin")
  if (!is.null(glu))   add_one("log1p_glucose",    safe_log1p_pos(df[[glu]]))   else skipped <- c(skipped, "log1p_glucose")
  if (!is.null(ag_))   add_one("log1p_anion_gap",  safe_log1p_pos(df[[ag_]]))   else skipped <- c(skipped, "log1p_anion_gap")
  
  attr(df, "engineered_created") <- unique(created)
  attr(df, "engineered_skipped") <- unique(skipped)
  df
}

before_cols <- names(df_train)
df_train <- add_features(df_train)
df_test  <- add_features(df_test)

eng_created <- setdiff(names(df_train), before_cols)
eng_created <- setdiff(eng_created, "outcome")
eng_skipped <- attr(df_train, "engineered_skipped") %||% character(0)

message("Engineered created (present): ", if (length(eng_created)) paste(eng_created, collapse = ", ") else "None")
message("Engineered skipped (not computable or all-NA): ", if (length(eng_skipped)) paste(eng_skipped, collapse = ", ") else "None")

df_test <- align_test_to_train(df_train, df_test)
imp2 <- median_impute_numeric(df_train, df_test); df_train <- imp2$train; df_test <- imp2$test
fix2 <- fill_remaining_numeric_na(df_train, df_test, 0.0); df_train <- fix2$train; df_test <- fix2$test

# =========================
# 4) FEATURE SETS
# =========================
all_cols <- setdiff(names(df_train), "outcome")

noninv_bases <- c(
  "age_raw","sex_raw","spo2_raw","heart_rate_raw","resp_rate_raw","temp_raw","sbp_raw","dbp_raw","map_raw"
)
feat_noninv <- noninv_bases[noninv_bases %in% all_cols]
if (length(feat_noninv) != length(noninv_bases)) {
  miss <- setdiff(noninv_bases, feat_noninv)
  stop("Missing non-invasive columns: ", paste(miss, collapse = ", "))
}

lab_map <- list(
  wbc         = c("wbc_raw","wbc"),
  neutrophils = c("neutrophils_raw","neutrophils"),
  lymphocytes = c("lymphocytes_raw","lymphocytes"),
  monocytes   = c("monocytes_raw","monocytes"),
  platelets   = c("platelets_raw","platelets"),
  hemoglobin  = c("hemoglobin_raw","hemoglobin"),
  albumin     = c("albumin_raw","albumin"),
  creatinine  = c("creatinine_raw","creatinine"),
  urea        = c("urea_raw","urea","bun_raw","bun"),
  sodium      = c("sodium_raw","sodium"),
  potassium   = c("potassium_raw","potassium"),
  chloride    = c("chloride_raw","chloride"),
  anion_gap   = c("anion_gap_raw","anion_gap"),
  glucose     = c("glucose_raw","glucose"),
  bilirubin   = c("bilirubin_raw","bilirubin"),
  alk_phos    = c("alk_phos_raw","alk_phos","alp_raw","alp"),
  crp         = c("crp_raw","crp"),
  d_dimer     = c("d_dimer_raw","d_dimer"),
  ast         = c("ast_raw","ast"),
  alt         = c("alt_raw","alt"),
  lactate     = c("lactate_raw","lactate"),
  troponin    = c("troponin_raw","troponin")
)

resolve_present <- function(map, all_cols) {
  cols <- character(0)
  found <- character(0)
  missing <- character(0)
  
  for (nm in names(map)) {
    cands <- map[[nm]]
    hit <- cands[cands %in% all_cols]
    if (length(hit) > 0L) {
      cols <- c(cols, hit[1])
      found <- c(found, nm)
    } else {
      missing <- c(missing, nm)
    }
  }
  list(cols = unique(cols), found = found, missing = missing)
}

lab_resolved <- resolve_present(lab_map, all_cols)
feat_labs <- lab_resolved$cols

message("Labs found: ", if (length(lab_resolved$found)) paste(lab_resolved$found, collapse = ", ") else "None")
message("Labs missing: ", if (length(lab_resolved$missing)) paste(lab_resolved$missing, collapse = ", ") else "None")

feat_full <- unique(c(feat_noninv, feat_labs, eng_created))
feat_full <- feat_full[feat_full %in% names(df_train)]

only_numeric <- function(df, feats) {
  feats <- feats[feats %in% names(df)]
  feats[vapply(df[, feats, drop = FALSE], is.numeric, logical(1))]
}

feat_noninv <- only_numeric(df_train, feat_noninv)
feat_full   <- only_numeric(df_train, feat_full)

if (length(feat_noninv) < 2L) stop("Non-invasive numeric predictors are fewer than 2 after coercion.")
if (length(feat_full)   < 2L) stop("Laboratory augmented numeric predictors are fewer than 2 after coercion.")

drop_nzv <- function(train, feat_names) {
  x <- train[, feat_names, drop = FALSE]
  nz <- caret::nearZeroVar(x)
  dropped <- character(0)
  keep <- feat_names
  if (length(nz) > 0L) {
    dropped <- feat_names[nz]
    keep <- setdiff(feat_names, dropped)
  }
  list(keep = keep, dropped = dropped)
}

nz1 <- drop_nzv(df_train, feat_noninv)
feat_noninv <- nz1$keep
dropped_noninv <- nz1$dropped

nz2 <- drop_nzv(df_train, feat_full)
feat_full <- nz2$keep
dropped_full <- nz2$dropped

eng_used <- intersect(eng_created, feat_full)
eng_dropped_nzv <- intersect(eng_created, dropped_full)

# =========================
# 5) METRICS + THRESHOLDS
# =========================
confusion_counts <- function(obs, pred, pos = "Yes", neg = "No") {
  obs <- as.character(obs)
  pred <- as.character(pred)
  TP <- sum(obs == pos & pred == pos, na.rm = TRUE)
  TN <- sum(obs == neg & pred == neg, na.rm = TRUE)
  FP <- sum(obs == neg & pred == pos, na.rm = TRUE)
  FN <- sum(obs == pos & pred == neg, na.rm = TRUE)
  list(TP = TP, TN = TN, FP = FP, FN = FN)
}

mcc_from_counts <- function(TP, TN, FP, FN) {
  TP <- as.numeric(TP); TN <- as.numeric(TN); FP <- as.numeric(FP); FN <- as.numeric(FN)
  num <- (TP * TN) - (FP * FN)
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (!is.finite(den) || den == 0) return(NA_real_)
  num / den
}

auc_of_probs <- function(obs, prob, pos = "Yes", neg = "No") {
  p <- suppressWarnings(as.numeric(prob))
  y <- as.character(obs)
  ok <- is.finite(p) & !is.na(p) & y %in% c(pos, neg)
  if (sum(ok) < 10) return(NA_real_)
  rocobj <- pROC::roc(response = y[ok], predictor = p[ok], levels = c(neg, pos), direction = "<", quiet = TRUE)
  as.numeric(pROC::auc(rocobj))
}

compute_metrics <- function(obs, pred, prob, pos = "Yes", neg = "No") {
  cc <- confusion_counts(obs, pred, pos, neg)
  TP <- cc$TP; TN <- cc$TN; FP <- cc$FP; FN <- cc$FN
  mcc <- mcc_from_counts(TP, TN, FP, FN)
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1        <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA_real_
  
  acc <- (TP + TN) / (TP + TN + FP + FN)
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  auc <- auc_of_probs(obs, prob, pos, neg)
  
  list(
    TP = TP, TN = TN, FP = FP, FN = FN,
    MCC = mcc, AUC = auc, F1 = f1,
    Accuracy = acc, Precision = precision,
    Sensitivity = recall, Specificity = spec
  )
}

best_mcc_and_thresh <- function(obs, prob, pos = "Yes", neg = "No", fallback = 0.5) {
  p <- suppressWarnings(as.numeric(prob))
  y <- as.character(obs)
  ok <- is.finite(p) & !is.na(p) & y %in% c(pos, neg)
  if (sum(ok) < 10) return(list(thr = fallback, mcc = NA_real_))
  
  p <- p[ok]
  y <- y[ok]
  y_pos <- (y == pos)
  
  o <- order(p, decreasing = TRUE)
  p <- p[o]
  y_pos <- y_pos[o]
  
  total_pos <- sum(y_pos)
  total_neg <- length(y_pos) - total_pos
  
  TP <- 0; FP <- 0
  TN <- total_neg; FN <- total_pos
  best_mcc <- mcc_from_counts(TP, TN, FP, FN)
  best_thr <- max(p) + 1e-12
  
  i <- 1L
  n <- length(p)
  while (i <= n) {
    j <- i
    while (j <= n && p[j] == p[i]) j <- j + 1L
    
    grp_pos <- sum(y_pos[i:(j-1)])
    grp_neg <- (j - i) - grp_pos
    
    TP <- TP + grp_pos
    FP <- FP + grp_neg
    FN <- total_pos - TP
    TN <- total_neg - FP
    
    mcc <- mcc_from_counts(TP, TN, FP, FN)
    if (is.finite(mcc) && (!is.finite(best_mcc) || mcc > best_mcc)) {
      best_mcc <- mcc
      best_thr <- p[i]
    }
    i <- j
  }
  
  list(thr = best_thr, mcc = best_mcc)
}

best_thresh_mcc <- function(obs, prob, pos = "Yes", neg = "No", fallback = 0.5) {
  best_mcc_and_thresh(obs, prob, pos, neg, fallback)$thr
}

cascade_prob_mixture <- function(pT, pC, t_low, t_high) {
  ifelse(pT <= t_low | pT >= t_high, pT, pC)
}

# =========================
# 5b) SHARED CV SPLITS
# =========================
make_cv_indices <- function(y, k = 5, repeats = 2, seed = 123) {
  set.seed(seed)
  idx <- vector("list", k * repeats)
  idx_out <- vector("list", k * repeats)
  ctr <- 1L
  for (r in seq_len(repeats)) {
    folds <- caret::createFolds(y, k = k, returnTrain = TRUE)
    for (i in seq_len(k)) {
      idx[[ctr]] <- folds[[i]]
      idx_out[[ctr]] <- setdiff(seq_along(y), folds[[i]])
      ctr <- ctr + 1L
    }
  }
  list(index = idx, indexOut = idx_out)
}

cv_idx <- make_cv_indices(df_train$outcome, k = 5, repeats = 2, seed = cfg$seed_cv)

mcc_bestthr_summary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev) || length(lev) < 2) return(c(MCC = NA_real_))
  pos <- lev[1]; neg <- lev[2]
  
  if (!(pos %in% names(data))) return(c(MCC = NA_real_))
  prob <- suppressWarnings(as.numeric(data[[pos]]))
  obs  <- as.character(data$obs)
  ok <- is.finite(prob) & !is.na(prob) & obs %in% c(pos, neg)
  if (sum(ok) < 10) return(c(MCC = NA_real_))
  
  bm <- best_mcc_and_thresh(obs[ok], prob[ok], pos = pos, neg = neg, fallback = 0.5)
  c(MCC = bm$mcc)
}

ctrl <- caret::trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 2,
  index           = cv_idx$index,
  indexOut        = cv_idx$indexOut,
  classProbs      = TRUE,
  summaryFunction = mcc_bestthr_summary,
  savePredictions = "final",
  allowParallel   = TRUE
)

make_case_weights <- function(y) {
  y <- factor(as.character(y), levels = c(cfg$pos_label, cfg$neg_label))
  tab <- table(y)
  if (any(tab == 0L)) {
    stop("Training outcome is missing a class: ", paste(names(tab)[tab == 0L], collapse = ", "))
  }
  w_class <- sum(tab) / (2 * tab)
  ifelse(y == cfg$pos_label, w_class[[cfg$pos_label]], w_class[[cfg$neg_label]]) |> as.numeric()
}

extract_oof <- function(fit, pos_label = "Yes", neg_label = "No") {
  pred <- fit$pred
  if (is.null(pred) || nrow(pred) == 0L) stop("OOF predictions were not available from caret::train (fit$pred empty).")
  
  bt <- fit$bestTune
  if (!is.null(bt) && ncol(bt) > 0L) {
    for (nm in names(bt)) {
      if (nm %in% names(pred)) pred <- pred[pred[[nm]] == bt[[nm]], , drop = FALSE]
    }
  }
  
  prob <- NULL
  if (pos_label %in% names(pred)) {
    prob <- suppressWarnings(as.numeric(pred[[pos_label]]))
  } else if (neg_label %in% names(pred)) {
    prob <- 1 - suppressWarnings(as.numeric(pred[[neg_label]]))
  } else {
    skip <- c("pred","obs","rowIndex","Resample")
    tune_nms <- names(bt %||% data.frame())
    cand <- setdiff(names(pred), c(skip, tune_nms))
    num_cand <- cand[vapply(pred[, cand, drop = FALSE], is.numeric, logical(1))]
    looks_prob <- num_cand[sapply(num_cand, function(nm) {
      x <- pred[[nm]]
      ok <- is.finite(x) & !is.na(x)
      if (sum(ok) < 50) return(FALSE)
      mean(x[ok] >= 0 & x[ok] <= 1) > 0.98
    })]
    if (length(looks_prob) == 1L) {
      prob <- as.numeric(pred[[looks_prob]])
    } else {
      stop("Could not locate probability column in fit$pred. Columns are: ", paste(names(pred), collapse = ", "))
    }
  }
  
  tibble::tibble(
    rowIndex = pred$rowIndex,
    Resample = pred$Resample %||% NA_character_,
    obs      = as.character(pred$obs),
    pred     = as.character(pred$pred),
    prob     = prob
  )
}

# =========================
# 6) TRAIN ONE ALGORITHM
# =========================
run_algo <- function(algo, feat_names, df_train, df_test) {
  x_train <- df_train[, feat_names, drop = FALSE]
  y_train <- factor(as.character(df_train$outcome), levels = c(cfg$pos_label, cfg$neg_label))
  x_test  <- df_test[, feat_names, drop = FALSE]
  
  w <- NULL
  if (isTRUE(use_class_weights)) w <- make_case_weights(y_train)
  
  if (algo == "LR") {
    tune <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-4, 1, length.out = 30))
    fit <- caret::train(
      x = x_train, y = y_train,
      method = "glmnet",
      metric = "MCC",
      trControl = ctrl,
      tuneGrid = tune,
      weights = w
    )
  } else if (algo == "RF") {
    fit <- caret::train(
      x = x_train, y = y_train,
      method = "rf",
      metric = "MCC",
      trControl = ctrl,
      tuneLength = 5,
      weights = w
    )
  } else if (algo == "SVM") {
    fit <- caret::train(
      x = x_train, y = y_train,
      method = "svmRadial",
      metric = "MCC",
      trControl = ctrl,
      tuneLength = 10,
      weights = w
    )
  } else if (algo == "k-NN") {
    fit <- caret::train(
      x = x_train, y = y_train,
      method = "kknn",
      metric = "MCC",
      trControl = ctrl,
      tuneLength = 15,
      weights = w
    )
  } else if (algo == "C4.5") {
    fit <- caret::train(
      x = x_train, y = y_train,
      method = "J48",
      metric = "MCC",
      trControl = ctrl,
      tuneLength = 10,
      weights = w
    )
  } else {
    stop("Unknown algo: ", algo)
  }
  
  prob_test <- as.numeric(predict(fit, newdata = x_test, type = "prob")[[cfg$pos_label]])
  oof <- extract_oof(fit, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
  
  list(fit = fit, prob_test = prob_test, oof = oof)
}

models <- c("LR","RF","SVM","k-NN","C4.5")

noninv_res <- setNames(vector("list", length(models)), models)
full_res   <- setNames(vector("list", length(models)), models)

for (algo in models) {
  message("Fitting ", algo, " (non-invasive)...")
  set.seed(cfg$seed_cv + match(algo, models) * 10 + 1)
  noninv_res[[algo]] <- run_algo(algo, feat_noninv, df_train, df_test)
  
  message("Fitting ", algo, " (laboratory augmented)...")
  set.seed(cfg$seed_cv + match(algo, models) * 10 + 2)
  full_res[[algo]]   <- run_algo(algo, feat_full, df_train, df_test)
}

# =========================
# 6b) OPTIMIZE thrT, thrC BY MCC (TRAIN OOF)
# =========================
opt_tbl <- setNames(vector("list", length(models)), models)

for (algo in models) {
  oT <- noninv_res[[algo]]$oof %>% dplyr::rename(probT = .data$prob)
  oC <- full_res[[algo]]$oof   %>% dplyr::rename(probC = .data$prob)
  
  join_keys <- intersect(c("rowIndex","Resample"), intersect(names(oT), names(oC)))
  o <- dplyr::inner_join(oT, oC, by = join_keys, suffix = c("_T","_C"))
  
  bT <- best_mcc_and_thresh(o$obs_T, o$probT, pos = cfg$pos_label, neg = cfg$neg_label, fallback = 0.5)
  bC <- best_mcc_and_thresh(o$obs_T, o$probC, pos = cfg$pos_label, neg = cfg$neg_label, fallback = 0.5)
  
  opt_tbl[[algo]] <- list(thrT = bT$thr, thrC = bC$thr, oof_mcc_T = bT$mcc, oof_mcc_C = bC$mcc)
  
  message(sprintf(
    "Cutoffs (OOF MCC-opt) %s: thrT=%.3f (MCC=%.4f), thrC=%.3f (MCC=%.4f)",
    algo, bT$thr, bT$mcc, bC$thr, bC$mcc
  ))
}

# =========================
# 6c) GRID SEARCH (t_low, t_high) BY MCC (TRAIN OOF)
# =========================
grid_low  <- seq(grid_low_min,  grid_low_max,  by = grid_step)
grid_high <- seq(grid_high_min, grid_high_max, by = grid_step)

band_grid <- expand.grid(t_low = grid_low, t_high = grid_high, KEEP.OUT.ATTRS = FALSE) %>%
  dplyr::filter(.data$t_high - .data$t_low >= min_gap) %>%
  dplyr::arrange(.data$t_low, .data$t_high)

oof_band_score <- function(obs, probT, probC, thrT, thrC, t_low, t_high,
                           pos = "Yes", neg = "No") {
  pred <- ifelse(
    probT <= t_low | probT >= t_high,
    ifelse(probT >= thrT, pos, neg),
    ifelse(probC >= thrC, pos, neg)
  )
  defer <- mean(probT > t_low & probT < t_high)
  cc  <- confusion_counts(obs, pred, pos = pos, neg = neg)
  mcc <- mcc_from_counts(cc$TP, cc$TN, cc$FP, cc$FN)
  list(mcc = mcc, defer = defer)
}

tune_band_for_algo <- function(algo) {
  oT <- noninv_res[[algo]]$oof %>% dplyr::rename(probT = .data$prob)
  oC <- full_res[[algo]]$oof   %>% dplyr::rename(probC = .data$prob)
  
  join_keys <- intersect(c("rowIndex","Resample"), intersect(names(oT), names(oC)))
  o <- dplyr::inner_join(oT, oC, by = join_keys, suffix = c("_T","_C"))
  
  obs <- as.character(o$obs_T)
  pT  <- as.numeric(o$probT)
  pC  <- as.numeric(o$probC)
  
  thrT <- opt_tbl[[algo]]$thrT
  thrC <- opt_tbl[[algo]]$thrC
  
  scores <- band_grid %>%
    dplyr::mutate(mcc = NA_real_, defer = NA_real_)
  
  for (i in seq_len(nrow(scores))) {
    s <- oof_band_score(
      obs = obs, probT = pT, probC = pC,
      thrT = thrT, thrC = thrC,
      t_low = scores$t_low[i], t_high = scores$t_high[i],
      pos = cfg$pos_label, neg = cfg$neg_label
    )
    scores$mcc[i]   <- s$mcc
    scores$defer[i] <- s$defer
  }
  
  if (!is.null(defer_max)) {
    scores <- dplyr::filter(scores, .data$defer <= defer_max)
    if (nrow(scores) == 0L) stop("No (t_low,t_high) satisfied defer_max for algo=", algo,
                                 ". Relax defer_max or grid.")
  }
  
  best <- scores %>%
    dplyr::arrange(dplyr::desc(.data$mcc), .data$defer, .data$t_low, .data$t_high) %>%
    dplyr::slice(1)
  
  tibble::tibble(
    Algorithm = algo,
    thrT = thrT,
    thrC = thrC,
    t_low = best$t_low,
    t_high = best$t_high,
    OOF_Cascade_MCC = best$mcc,
    OOF_DeferRate = best$defer
  )
}

band_tbl <- dplyr::bind_rows(lapply(models, tune_band_for_algo))

cat("\n=== BEST (t_low, t_high) PER ALGORITHM (TRAIN OOF, MCC objective) ===\n")
print(band_tbl %>% dplyr::arrange(dplyr::desc(.data$OOF_Cascade_MCC)), n = Inf, width = Inf)

# =========================
# 7) BOOTSTRAP CONFIDENCE INTERVALS (EXTERNAL)
# =========================
boot_ci_counts <- function(obs, pred, B = 500, seed0 = 1) {
  set.seed(seed0)
  n <- length(obs)
  out <- matrix(NA_real_, nrow = B, ncol = 4)
  colnames(out) <- c("TP","TN","FP","FN")
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    cc <- confusion_counts(obs[idx], pred[idx], pos = cfg$pos_label, neg = cfg$neg_label)
    out[b, ] <- c(cc$TP, cc$TN, cc$FP, cc$FN)
  }
  
  qs <- apply(out, 2, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  qs <- t(qs)
  colnames(qs) <- c("lo","hi")
  qs
}

# FIXED: now includes Specificity
boot_ci_metrics <- function(obs, pred, prob, B = 500, seed0 = 1) {
  set.seed(seed0)
  n <- length(obs)
  
  stat <- matrix(NA_real_, nrow = B, ncol = 7)
  colnames(stat) <- c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    m <- compute_metrics(obs[idx], pred[idx], prob[idx], pos = cfg$pos_label, neg = cfg$neg_label)
    stat[b, ] <- c(m$MCC, m$AUC, m$F1, m$Accuracy, m$Precision, m$Sensitivity, m$Specificity)
  }
  
  qs <- apply(stat, 2, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  qs <- t(qs)
  colnames(qs) <- c("lo","hi")
  qs
}

fmt_ci <- function(est, lo, hi, digits = 4) {
  if (!is.finite(est)) return("NA")
  if (!is.finite(lo) || !is.finite(hi)) return(sprintf(paste0("%.", digits, "f"), est))
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), est, lo, hi)
}

fmt_ci_int <- function(est, lo, hi) {
  if (!is.finite(est)) return("NA")
  if (!is.finite(lo) || !is.finite(hi)) return(sprintf("%.0f", est))
  sprintf("%.0f [%.0f, %.0f]", est, lo, hi)
}

# =========================
# 8) EVALUATE ON EXTERNAL
# =========================
metrics_row <- function(algo, panel, t_low, t_high, defer, thrT, thrC, m, ci_counts, ci_metrics) {
  tibble::tibble(
    Algorithm   = algo,
    Panel       = panel,
    t_low       = t_low,
    t_high      = t_high,
    DeferRate   = defer,
    thrT        = thrT,
    thrC        = thrC,
    TP          = fmt_ci_int(m$TP, ci_counts["TP","lo"], ci_counts["TP","hi"]),
    TN          = fmt_ci_int(m$TN, ci_counts["TN","lo"], ci_counts["TN","hi"]),
    FP          = fmt_ci_int(m$FP, ci_counts["FP","lo"], ci_counts["FP","hi"]),
    FN          = fmt_ci_int(m$FN, ci_counts["FN","lo"], ci_counts["FN","hi"]),
    MCC         = fmt_ci(m$MCC, ci_metrics["MCC","lo"], ci_metrics["MCC","hi"], digits = 4),
    AUC         = fmt_ci(m$AUC, ci_metrics["AUC","lo"], ci_metrics["AUC","hi"], digits = 4),
    F1          = fmt_ci(m$F1,  ci_metrics["F1","lo"],  ci_metrics["F1","hi"],  digits = 4),
    Accuracy    = fmt_ci(m$Accuracy, ci_metrics["Accuracy","lo"], ci_metrics["Accuracy","hi"], digits = 4),
    Precision   = fmt_ci(m$Precision, ci_metrics["Precision","lo"], ci_metrics["Precision","hi"], digits = 4),
    Sensitivity = fmt_ci(m$Sensitivity, ci_metrics["Sensitivity","lo"], ci_metrics["Sensitivity","hi"], digits = 4),
    # FIXED: now prints CI for specificity
    Specificity = fmt_ci(m$Specificity, ci_metrics["Specificity","lo"], ci_metrics["Specificity","hi"], digits = 4)
  )
}

eval_algo <- function(algo) {
  y  <- as.character(df_test$outcome)
  pT <- as.numeric(noninv_res[[algo]]$prob_test)
  pC <- as.numeric(full_res[[algo]]$prob_test)
  
  thrT <- opt_tbl[[algo]]$thrT
  thrC <- opt_tbl[[algo]]$thrC
  
  t_low  <- band_tbl$t_low[band_tbl$Algorithm == algo][1]
  t_high <- band_tbl$t_high[band_tbl$Algorithm == algo][1]
  if (!is.finite(t_low) || !is.finite(t_high)) stop("Missing tuned band for algo=", algo)
  
  # Non-invasive
  predT <- ifelse(pT >= thrT, cfg$pos_label, cfg$neg_label)
  mT <- compute_metrics(y, predT, pT, pos = cfg$pos_label, neg = cfg$neg_label)
  ciCT <- boot_ci_counts(y, predT, B = B_boot, seed0 = cfg$seed_cv + 10000 + match(algo, models))
  ciMT <- boot_ci_metrics(y, predT, pT, B = B_boot, seed0 = cfg$seed_cv + 20000 + match(algo, models))
  rowT <- metrics_row(algo, "Non-invasive", NA_real_, NA_real_, NA_real_, thrT, NA_real_, mT, ciCT, ciMT)
  
  # Laboratory augmented
  predC <- ifelse(pC >= thrC, cfg$pos_label, cfg$neg_label)
  mC <- compute_metrics(y, predC, pC, pos = cfg$pos_label, neg = cfg$neg_label)
  ciCC <- boot_ci_counts(y, predC, B = B_boot, seed0 = cfg$seed_cv + 30000 + match(algo, models))
  ciMC <- boot_ci_metrics(y, predC, pC, B = B_boot, seed0 = cfg$seed_cv + 40000 + match(algo, models))
  rowC <- metrics_row(algo, "Laboratory augmented", NA_real_, NA_real_, NA_real_, NA_real_, thrC, mC, ciCC, ciMC)
  
  # Cascade
  defer <- mean(pT > t_low & pT < t_high)
  
  predX <- ifelse(
    pT <= t_low | pT >= t_high,
    ifelse(pT >= thrT, cfg$pos_label, cfg$neg_label),
    ifelse(pC >= thrC, cfg$pos_label, cfg$neg_label)
  )
  
  if (cascade_auc_convention == "mixture") {
    p_mix <- cascade_prob_mixture(pT, pC, t_low, t_high)
    mX <- compute_metrics(y, predX, p_mix, pos = cfg$pos_label, neg = cfg$neg_label)
    prob_for_boot <- p_mix
  } else {
    deferred <- (pT > t_low & pT < t_high)
    prob_for_auc <- ifelse(deferred, pC, NA_real_)
    mX <- compute_metrics(y, predX, prob_for_auc, pos = cfg$pos_label, neg = cfg$neg_label)
    prob_for_boot <- prob_for_auc
  }
  
  ciCX <- boot_ci_counts(y, predX, B = B_boot, seed0 = cfg$seed_cv + 50000 + match(algo, models))
  ciMX <- boot_ci_metrics(y, predX, prob_for_boot, B = B_boot, seed0 = cfg$seed_cv + 60000 + match(algo, models))
  rowX <- metrics_row(
    algo,
    "Cascade (GRID-OPT t_low,t_high + MCC-opt cutoffs)",
    t_low, t_high, defer, thrT, thrC, mX, ciCX, ciMX
  )
  
  dplyr::bind_rows(rowT, rowC, rowX)
}

res_tbl <- dplyr::bind_rows(lapply(models, eval_algo))

cat("\n=== RESULTS (bootstrap 95% percentile confidence intervals on EXTERNAL test) ===\n")
print(res_tbl, n = Inf, width = Inf)

cat("\n=== OUTCOME COUNTS ===\n")
cat("Development:\n"); print(table(df_train$outcome))
cat("External:\n");     print(table(df_test$outcome))

cat("\n=== FEATURE SETS (USED) ===\n")
cat("Non-invasive: ", paste(feat_noninv, collapse = ", "), "\n", sep = "")
cat("Laboratory augmented: ", paste(feat_full, collapse = ", "), "\n", sep = "")

cat("\n=== ENGINEERED FEATURES ===\n")
cat("Created (present): ", if (length(eng_created)) paste(eng_created, collapse = ", ") else "None", "\n", sep = "")
cat("Used (in lab-augmented model): ", if (length(eng_used)) paste(eng_used, collapse = ", ") else "None", "\n", sep = "")
cat("Dropped by NZV: ", if (length(eng_dropped_nzv)) paste(eng_dropped_nzv, collapse = ", ") else "None", "\n", sep = "")

if (length(dropped_noninv)) cat("\nDropped NZV (non-invasive): ", paste(dropped_noninv, collapse = ", "), "\n", sep = "")
if (length(dropped_full))   cat("Dropped NZV (lab-augmented): ", paste(dropped_full, collapse = ", "), "\n", sep = "")

# ============================================================
# 9) OPTIONAL: EXTERNAL paired bootstrap comparisons (A minus B)
# among the 3 policies for all metrics
# ============================================================

safe_metrics_vec <- function(y, pred, prob, pos, neg) {
  out <- tryCatch(
    compute_metrics(y, pred, prob, pos = pos, neg = neg),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(c(MCC=NA_real_, AUC=NA_real_, F1=NA_real_, Accuracy=NA_real_,
             Precision=NA_real_, Sensitivity=NA_real_, Specificity=NA_real_))
  }
  c(MCC=out$MCC, AUC=out$AUC, F1=out$F1, Accuracy=out$Accuracy,
    Precision=out$Precision, Sensitivity=out$Sensitivity, Specificity=out$Specificity)
}

paired_boot_delta_all_metrics <- function(y, predA, probA, predB, probB,
                                          B = 2000, seed = 1, pos, neg) {
  set.seed(seed)
  n <- length(y)
  
  estA <- safe_metrics_vec(y, predA, probA, pos, neg)
  estB <- safe_metrics_vec(y, predB, probB, pos, neg)
  delta_hat <- estA - estB
  
  deltas <- matrix(NA_real_, nrow = B, ncol = length(delta_hat))
  colnames(deltas) <- names(delta_hat)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    mA <- safe_metrics_vec(y[idx], predA[idx], probA[idx], pos, neg)
    mB <- safe_metrics_vec(y[idx], predB[idx], probB[idx], pos, neg)
    deltas[b, ] <- mA - mB
  }
  
  ci <- apply(deltas, 2, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  ci <- t(ci); colnames(ci) <- c("lo","hi")
  
  p_raw <- vapply(seq_along(delta_hat), function(j) {
    x <- deltas[, j]
    x <- x[is.finite(x) & !is.na(x)]
    if (length(x) < 50) return(NA_real_)
    p2 <- 2 * min(mean(x <= 0), mean(x >= 0))
    min(1, p2)
  }, numeric(1))
  
  tibble::tibble(
    metric = names(delta_hat),
    est_A  = as.numeric(estA),
    est_B  = as.numeric(estB),
    delta  = as.numeric(delta_hat),
    lo     = ci[, "lo"],
    hi     = ci[, "hi"],
    p_raw  = as.numeric(p_raw)
  )
}

build_three_external <- function(algo) {
  y <- as.character(df_test$outcome)
  pos <- cfg$pos_label
  neg <- cfg$neg_label
  
  pT <- as.numeric(noninv_res[[algo]]$prob_test)
  pC <- as.numeric(full_res[[algo]]$prob_test)
  
  thrT <- opt_tbl[[algo]]$thrT
  thrC <- opt_tbl[[algo]]$thrC
  
  t_low  <- band_tbl$t_low[band_tbl$Algorithm == algo][1]
  t_high <- band_tbl$t_high[band_tbl$Algorithm == algo][1]
  
  pred_noninv <- ifelse(pT >= thrT, pos, neg)
  prob_noninv <- pT
  
  pred_aug <- ifelse(pC >= thrC, pos, neg)
  prob_aug <- pC
  
  pred_cas <- ifelse(
    pT <= t_low | pT >= t_high,
    ifelse(pT >= thrT, pos, neg),
    ifelse(pC >= thrC, pos, neg)
  )
  
  if (cascade_auc_convention == "mixture") {
    prob_cas <- cascade_prob_mixture(pT, pC, t_low, t_high)
  } else {
    deferred <- (pT > t_low & pT < t_high)
    prob_cas <- ifelse(deferred, pC, NA_real_)
  }
  
  policy_info <- tibble::tibble(
    Algorithm = algo,
    thrT = thrT,
    thrC = thrC,
    t_low = t_low,
    t_high = t_high,
    DeferRate_External = mean(pT > t_low & pT < t_high)
  )
  
  list(
    y = y, pos = pos, neg = neg,
    policy_info = policy_info,
    Non_invasive = list(pred = pred_noninv, prob = prob_noninv),
    Augmented    = list(pred = pred_aug,    prob = prob_aug),
    Cascade      = list(pred = pred_cas,    prob = prob_cas)
  )
}

pairwise_compare_three <- function(three, B = 2000, seed = 1) {
  y <- three$y
  pos <- three$pos
  neg <- three$neg
  
  pairs <- list(
    c("Augmented", "Non_invasive"),
    c("Cascade",   "Augmented"),
    c("Cascade",   "Non_invasive")
  )
  
  out <- list()
  for (pp in pairs) {
    A <- pp[1]; Bname <- pp[2]
    tab <- paired_boot_delta_all_metrics(
      y = y,
      predA = three[[A]]$pred, probA = three[[A]]$prob,
      predB = three[[Bname]]$pred, probB = three[[Bname]]$prob,
      B = B, seed = seed + sum(utf8ToInt(A)) + sum(utf8ToInt(Bname)),
      pos = pos, neg = neg
    ) %>%
      dplyr::mutate(Comparison = paste(A, "minus", Bname))
    out[[paste(A, Bname, sep = "_")]] <- tab
  }
  
  dplyr::bind_rows(out)
}

# Run pairwise comparisons for each algorithm
B_pairwise <- 2000L
all_policy <- vector("list", length(models))
all_pairs  <- vector("list", length(models))

for (i in seq_along(models)) {
  algo <- models[i]
  three <- build_three_external(algo)
  
  all_policy[[i]] <- three$policy_info
  all_pairs[[i]] <- pairwise_compare_three(three, B = B_pairwise, seed = cfg$seed_cv + 90000 + i) %>%
    dplyr::mutate(Algorithm = algo, .before = 1)
}

external_cascade_policy <- dplyr::bind_rows(all_policy)
external_pairwise_all_metrics <- dplyr::bind_rows(all_pairs) %>%
  dplyr::select(Algorithm, Comparison, metric, est_A, est_B, delta, lo, hi, p_raw)

# Holm correction within each metric across all (Algorithm x Comparison)
external_pairwise_all_metrics <- external_pairwise_all_metrics %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(p_holm_within_metric = p.adjust(p_raw, method = "holm")) %>%
  dplyr::ungroup()

cat("\n=== External: cascade policy and defer rate (frozen from Development tuning) ===\n")
print(external_cascade_policy, n = Inf, width = Inf)

cat("\n=== External: paired bootstrap deltas for all metrics (A minus B) ===\n")
print(external_pairwise_all_metrics, n = Inf, width = Inf)



