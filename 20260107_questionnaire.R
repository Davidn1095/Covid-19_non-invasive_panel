options(warn = -1)
options(na.action = "na.pass")

# ============================================================
# QUESTIONNAIRE (NON-INVASIVE) — PERFORMANCE VERSION
# - Repeated stratified k-fold OOF on development
# - Threshold chosen on DEV OOF (max MCC)
# - Final model trained on full development, evaluated on external
# - No file writing
# ============================================================

pos_label <- "Yes"
neg_label <- "No"
use_class_weights <- TRUE

source(file.path("R", "common_utils.R"))


risk_low  <- 0.33
risk_high <- 0.66
# binary_mcc | binary_low | binary_high | three_zone
decision_mode <- "binary_mcc"

# -------------------------
# USER SETTINGS
# -------------------------
dev_csv <- Sys.getenv("DEV_CSV", unset = "")
ext_csv <- Sys.getenv("EXT_CSV", unset = "")

seed_cv <- 123
k_folds <- 5L
repeats <- 2L

dev_max_n <- as.integer(Sys.getenv("DEV_MAX_N", unset = "10000"))  # 0 = no cap  # default 10k
ext_max_n <- as.integer(Sys.getenv("EXT_MAX_N", unset = "10000"))  # 0 = no cap  # default 10k

# -------------------------
# 0) find CSV paths
# -------------------------
find_one_csv <- function(pattern) {
  csvs <- list.files(".", pattern = "\\.csv$", ignore.case = TRUE)
  if (!length(csvs)) stop("No .csv files found in: ", getwd())
  hits <- grep(pattern, csvs, value = TRUE, ignore.case = TRUE)
  if (!length(hits)) return(NA_character_)
  if (length(hits) > 1L) stop("More than one CSV match for pattern: ", pattern, "\n- ", paste(hits, collapse = "\n- "))
  hits[[1]]
}

if (nzchar(dev_csv)) {
  if (!file.exists(dev_csv)) stop("DEV_CSV does not exist: ", dev_csv)
  dev_path <- dev_csv
} else {
  dev_path <- find_one_csv("development.*mimic|mimic.*development")
  if (!nzchar(dev_path)) stop("No development MIMIC CSV match in current folder.")
}

if (nzchar(ext_csv)) {
  if (!file.exists(ext_csv)) stop("EXT_CSV does not exist: ", ext_csv)
  ext_path <- ext_csv
} else {
  ext_path <- find_one_csv("external.*mcmed|mcmed.*external|external")
  if (!nzchar(ext_path)) ext_path <- NA_character_
}

cat("Using dev_path = ", dev_path, "\n", sep = "")
if (nzchar(ext_path)) cat("Using ext_path = ", ext_path, "\n", sep = "") else cat("No external CSV found, external eval skipped.\n")

# -------------------------
# 1) helpers
# -------------------------
prep_outcome <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) || is.integer(x)) {
    x <- ifelse(is.na(x), NA_character_, ifelse(x == 1, pos_label, neg_label))
  } else {
    x0 <- trimws(tolower(as.character(x)))
    x  <- ifelse(
      x0 %in% c("1","yes","y","true","t","pos","positive","case","covid","covid19"),
      pos_label,
      ifelse(
        x0 %in% c("0","no","n","false","f","neg","negative","control","ctrl","nocovid"),
        neg_label,
        NA_character_
      )
    )
  }
  factor(x, levels = c(neg_label, pos_label))
}

clamp_na <- function(x, lo, hi, zero_is_na = FALSE) {
  x <- as.numeric(x)
  if (zero_is_na) x[x == 0] <- NA_real_
  x[!is.finite(x)] <- NA_real_
  x[x < lo | x > hi] <- NA_real_
  x
}

cut3_fixed <- function(x, breaks, labs) {
  z <- cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE, labels = labs)
  factor(z, levels = labs)
}

mode01 <- function(x) {
  x1 <- x[is.finite(x) & !is.na(x)]
  if (!length(x1)) return(0)
  if (mean(x1 == 1) >= 0.5) 1 else 0
}

auc_fast <- function(score, y01) {
  ok <- is.finite(score) & !is.na(score) & is.finite(y01) & !is.na(y01)
  score <- score[ok]
  y01 <- y01[ok]
  n1 <- sum(y01 == 1)
  n0 <- sum(y01 == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

best_cut_mcc <- function(score, y01) {
  ok <- is.finite(score) & !is.na(score) & is.finite(y01) & !is.na(y01)
  score <- score[ok]
  y01 <- y01[ok]
  if (!length(score)) return(list(cut = NA_real_, mcc = NA_real_))
  
  cuts <- sort(unique(score))
  best_mcc <- -Inf
  best_cut <- cuts[1]
  
  for (ct in cuts) {
    pred1 <- as.integer(score >= ct)
    TP <- sum(pred1 == 1 & y01 == 1)
    TN <- sum(pred1 == 0 & y01 == 0)
    FP <- sum(pred1 == 1 & y01 == 0)
    FN <- sum(pred1 == 0 & y01 == 1)
    m <- mcc_from_counts(c(TP = TP, TN = TN, FP = FP, FN = FN))
    if (is.finite(m) && m > best_mcc) {
      best_mcc <- m
      best_cut <- ct
    }
  }
  list(cut = best_cut, mcc = best_mcc)
}

metrics_from_prob <- function(prob, y01, thr) {
  pred1 <- as.integer(prob >= thr)
  TP <- sum(pred1 == 1 & y01 == 1)
  TN <- sum(pred1 == 0 & y01 == 0)
  FP <- sum(pred1 == 1 & y01 == 0)
  FN <- sum(pred1 == 0 & y01 == 1)
  
  mcc <- mcc_from_counts(c(TP = TP, TN = TN, FP = FP, FN = FN))
  prec <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  f1   <- if (!is.finite(prec) || !is.finite(sens) || (prec + sens) == 0) NA_real_ else 2 * prec * sens / (prec + sens)
  auc  <- auc_fast(prob, y01)
  
  data.frame(
    thr = thr,
    TP = TP, TN = TN, FP = FP, FN = FN,
    MCC = mcc, AUC = auc, F1 = f1,
    Accuracy = acc, Precision = prec, Sensitivity = sens, Specificity = spec,
    stringsAsFactors = FALSE
  )
}

make_stratified_folds <- function(y01, k, seed) {
  set.seed(seed)
  idx1 <- which(y01 == 1)
  idx0 <- which(y01 == 0)
  idx1 <- sample(idx1)
  idx0 <- sample(idx0)
  fold <- integer(length(y01))
  fold[idx1] <- rep(seq_len(k), length.out = length(idx1))
  fold[idx0] <- rep(seq_len(k), length.out = length(idx0))
  fold
}

# -------------------------
# 2) fixed questionnaire bins
# -------------------------
age_breaks  <- c(-Inf, 50, 68, Inf)
spo2_breaks <- c(-Inf, 92, 96, Inf)
hr_breaks   <- c(-Inf, 90, 110, Inf)
rr_breaks   <- c(-Inf, 20, 30, Inf)
t_breaks    <- c(-Inf, 36, 38, Inf)
sbp_breaks  <- c(-Inf, 100, 140, Inf)
dbp_breaks  <- c(-Inf, 60, 90, Inf)
map_breaks  <- c(-Inf, 65, 100, Inf)

# -------------------------
# 3) load and basic checks
# -------------------------
need <- c(
  "outcome","age_raw","sex_raw","spo2_raw","heart_rate_raw","resp_rate_raw","temp_raw","sbp_raw","dbp_raw","map_raw"
)

read_and_standardize <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- trimws(names(df))
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))
  df
}

dev <- read_and_standardize(dev_path)

# cap dev rows if requested
if (is.finite(dev_max_n) && dev_max_n > 0L && nrow(dev) > dev_max_n) {
  set.seed(seed_cv)
  dev <- dev[sample.int(nrow(dev), dev_max_n), , drop = FALSE]
  cat("DEV capped to n = ", nrow(dev), "\n", sep = "")
}

prep_numeric_and_impute <- function(df, impute = NULL) {
  out <- df
  
  out$outcome <- prep_outcome(out$outcome)
  if (any(is.na(out$outcome))) stop("Parsed NA in outcome, inspect raw outcome values.")
  
  out$age_raw        <- coerce_numeric_robust(out$age_raw)
  out$sex_raw        <- coerce_sex_to_01(out$sex_raw)
  out$spo2_raw       <- coerce_numeric_robust(out$spo2_raw)
  out$heart_rate_raw <- coerce_numeric_robust(out$heart_rate_raw)
  out$resp_rate_raw  <- coerce_numeric_robust(out$resp_rate_raw)
  out$temp_raw       <- coerce_numeric_robust(out$temp_raw)
  out$sbp_raw        <- coerce_numeric_robust(out$sbp_raw)
  out$dbp_raw        <- coerce_numeric_robust(out$dbp_raw)
  out$map_raw        <- coerce_numeric_robust(out$map_raw)
  
  tm <- suppressWarnings(stats::median(out$temp_raw, na.rm = TRUE))
  if (is.finite(tm) && tm > 60) out$temp_raw <- (out$temp_raw - 32) * 5/9
  
  out$age_raw        <- clamp_na(out$age_raw,          0, 120)
  out$spo2_raw       <- clamp_na(out$spo2_raw,        50, 100, zero_is_na = TRUE)
  out$heart_rate_raw <- clamp_na(out$heart_rate_raw,  20, 250)
  out$resp_rate_raw  <- clamp_na(out$resp_rate_raw,    5,  80)
  out$temp_raw       <- clamp_na(out$temp_raw,        30,  43)
  out$sbp_raw        <- clamp_na(out$sbp_raw,         50, 250)
  out$dbp_raw        <- clamp_na(out$dbp_raw,         30, 150)
  out$map_raw        <- clamp_na(out$map_raw,         40, 200)
  
  if (is.null(impute)) {
    sex_mode <- mode01(out$sex_raw)
    meds <- list()
    for (nm in c("age_raw","spo2_raw","heart_rate_raw","resp_rate_raw","temp_raw","sbp_raw","dbp_raw","map_raw")) {
      med <- suppressWarnings(stats::median(out[[nm]], na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      meds[[nm]] <- med
    }
    impute <- list(sex_mode = sex_mode, medians = meds)
  }
  
  out$sex_raw[is.na(out$sex_raw)] <- impute$sex_mode
  for (nm in names(impute$medians)) {
    out[[nm]][is.na(out[[nm]])] <- impute$medians[[nm]]
  }
  
  list(df = out, impute = impute)
}

make_bins <- function(df) {
  df$age_b  <- cut3_fixed(df$age_raw,        age_breaks,  c("<50", "50 to <68", ">=68"))
  df$sex_b  <- factor(ifelse(df$sex_raw == 1, "Male", "Female"), levels = c("Female","Male"))
  df$spo2_b <- cut3_fixed(df$spo2_raw,       spo2_breaks, c("<92", "92 to <96", ">=96"))
  df$hr_b   <- cut3_fixed(df$heart_rate_raw, hr_breaks,   c("<90", "90 to <110", ">=110"))
  df$rr_b   <- cut3_fixed(df$resp_rate_raw,  rr_breaks,   c("<20", "20 to <30", ">=30"))
  df$temp_b <- cut3_fixed(df$temp_raw,       t_breaks,    c("<36", "36 to <38", ">=38"))
  df$sbp_b  <- cut3_fixed(df$sbp_raw,        sbp_breaks,  c("<100", "100 to <140", ">=140"))
  df$dbp_b  <- cut3_fixed(df$dbp_raw,        dbp_breaks,  c("<60", "60 to <90", ">=90"))
  df$map_b  <- cut3_fixed(df$map_raw,        map_breaks,  c("<65", "65 to <100", ">=100"))
  df
}

# -------------------------
# 2b) Monotonic bins: build ordinal codes per variable from training data
#     - For each binned variable, order its 3 levels by observed risk in training
#     - Encode as 0, 1, 2 and fit a single slope so risk is monotonic in that order
# -------------------------
make_ord_map <- function(df_bins, y01, vars) {
  maps <- list()
  overall <- mean(y01, na.rm = TRUE)
  for (v in vars) {
    levs <- levels(df_bins[[v]])
    m <- tapply(y01, df_bins[[v]], mean, na.rm = TRUE)
    m <- m[levs]
    m[is.na(m)] <- overall
    # low risk -> high risk
    ord_levels <- levs[order(m, decreasing = FALSE)]
    maps[[v]] <- setNames(as.integer(seq_along(ord_levels) - 1L), ord_levels)
  }
  maps
}

apply_ord_map <- function(x_factor, map_named_int) {
  z <- as.character(x_factor)
  out <- unname(map_named_int[z])
  # if a level is missing, fall back to the middle code
  if (anyNA(out)) {
    mid <- as.integer(round(stats::median(unname(map_named_int), na.rm = TRUE)))
    out[is.na(out)] <- mid
  }
  as.numeric(out)
}

fit_questionnaire_model <- function(train_raw) {
  p  <- prep_numeric_and_impute(train_raw, impute = NULL)
  tr <- make_bins(p$df)

  y01 <- as.integer(tr$outcome == pos_label)

  w <- rep(1, nrow(tr))
  if (isTRUE(use_class_weights)) {
    tab <- table(tr$outcome)
    w_class <- sum(tab) / (2 * tab)
    w <- ifelse(tr$outcome == pos_label,
                as.numeric(w_class[[pos_label]]),
                as.numeric(w_class[[neg_label]]))
  }

  bin_vars <- c("age_b","spo2_b","hr_b","rr_b","temp_b","sbp_b","dbp_b","map_b")
  ord_map  <- make_ord_map(tr, y01, bin_vars)

  # ordinal encodings, 0..2 where 2 means higher observed risk in training
  tr$sex_o  <- ifelse(tr$sex_b == "Male", 1, 0)
  tr$age_o  <- apply_ord_map(tr$age_b,  ord_map$age_b)
  tr$spo2_o <- apply_ord_map(tr$spo2_b, ord_map$spo2_b)
  tr$hr_o   <- apply_ord_map(tr$hr_b,   ord_map$hr_b)
  tr$rr_o   <- apply_ord_map(tr$rr_b,   ord_map$rr_b)
  tr$temp_o <- apply_ord_map(tr$temp_b, ord_map$temp_b)
  tr$sbp_o  <- apply_ord_map(tr$sbp_b,  ord_map$sbp_b)
  tr$dbp_o  <- apply_ord_map(tr$dbp_b,  ord_map$dbp_b)
  tr$map_o  <- apply_ord_map(tr$map_b,  ord_map$map_b)

  # flip any variable whose slope goes negative, then refit
  flip_var <- function(tr_df, ord_map, bin_name, ord_name) {
    m <- max(unname(ord_map[[bin_name]]), na.rm = TRUE)
    ord_map[[bin_name]] <- m - ord_map[[bin_name]]
    tr_df[[ord_name]]   <- m - tr_df[[ord_name]]
    list(tr = tr_df, ord_map = ord_map)
  }

  fit <- NULL
  for (it in 1:5) {
    fit <- glm(
      y01 ~ sex_o + age_o + spo2_o + hr_o + rr_o + temp_o + sbp_o + dbp_o + map_o,
      data = tr,
      family = binomial(),
      weights = w
    )
    b <- stats::coef(fit)

    chk <- c(
      age_o  = b[["age_o"]],
      spo2_o = b[["spo2_o"]],
      hr_o   = b[["hr_o"]],
      rr_o   = b[["rr_o"]],
      temp_o = b[["temp_o"]],
      sbp_o  = b[["sbp_o"]],
      dbp_o  = b[["dbp_o"]],
      map_o  = b[["map_o"]]
    )
    to_flip <- names(chk)[is.finite(chk) & chk < 0]

    if (!length(to_flip)) break

    for (v in to_flip) {
      bin <- sub("_o$", "_b", v)
      tmp <- flip_var(tr, ord_map, bin_name = bin, ord_name = v)
      tr <- tmp$tr
      ord_map <- tmp$ord_map
    }
  }

  list(
    fit = fit,
    impute = p$impute,
    ord_map = ord_map
  )
}


predict_prob_questionnaire <- function(model, new_raw) {
  p  <- prep_numeric_and_impute(new_raw, impute = model$impute)
  nx <- make_bins(p$df)

  nx$sex_o  <- ifelse(nx$sex_b == "Male", 1, 0)
  nx$age_o  <- apply_ord_map(nx$age_b,  model$ord_map$age_b)
  nx$spo2_o <- apply_ord_map(nx$spo2_b, model$ord_map$spo2_b)
  nx$hr_o   <- apply_ord_map(nx$hr_b,   model$ord_map$hr_b)
  nx$rr_o   <- apply_ord_map(nx$rr_b,   model$ord_map$rr_b)
  nx$temp_o <- apply_ord_map(nx$temp_b, model$ord_map$temp_b)
  nx$sbp_o  <- apply_ord_map(nx$sbp_b,  model$ord_map$sbp_b)
  nx$dbp_o  <- apply_ord_map(nx$dbp_b,  model$ord_map$dbp_b)
  nx$map_o  <- apply_ord_map(nx$map_b,  model$ord_map$map_b)

  as.numeric(stats::predict(model$fit, newdata = nx, type = "response"))
}


# -------------------------
# 4) DEV: repeated stratified k-fold OOF probabilities
# -------------------------
dev0 <- prep_numeric_and_impute(dev, impute = NULL)$df
dev0$outcome <- prep_outcome(dev0$outcome)
y01_dev <- as.integer(dev0$outcome == pos_label)

oof_mat <- matrix(NA_real_, nrow(dev0), repeats)

for (r in seq_len(repeats)) {
  fold_id <- make_stratified_folds(y01_dev, k = k_folds, seed = seed_cv + 1000L * r)
  
  for (k in seq_len(k_folds)) {
    te <- which(fold_id == k)
    tr <- which(fold_id != k)
    
    mod <- fit_questionnaire_model(dev0[tr, , drop = FALSE])
    oof_mat[te, r] <- predict_prob_questionnaire(mod, dev0[te, , drop = FALSE])
  }
}

dev_oof_prob <- rowMeans(oof_mat, na.rm = TRUE)

# choose threshold on DEV OOF for max MCC
cut_oof <- best_cut_mcc(dev_oof_prob, y01_dev)
thr_prob <- cut_oof$cut

cat("\n=== DEV OOF THRESHOLD (max MCC) ===\n")
cat("thr_prob = ", format(thr_prob, digits = 6), "\n", sep = "")

cat("\n=== DEV OOF METRICS (threshold above) ===\n")
print(metrics_from_prob(dev_oof_prob, y01_dev, thr_prob), row.names = FALSE)

# -------------------------
# 5) Train final model on full DEV, evaluate EXT (if present)
# -------------------------
final_mod <- fit_questionnaire_model(dev0)

if (nzchar(ext_path)) {
  ext <- read_and_standardize(ext_path)
  
  if (is.finite(ext_max_n) && ext_max_n > 0L && nrow(ext) > ext_max_n) {
    set.seed(seed_cv + 999)
    ext <- ext[sample.int(nrow(ext), ext_max_n), , drop = FALSE]
    cat("EXT capped to n = ", nrow(ext), "\n", sep = "")
  }
  
  ext0 <- prep_numeric_and_impute(ext, impute = final_mod$impute)$df
  ext0$outcome <- prep_outcome(ext0$outcome)
  if (any(is.na(ext0$outcome))) stop("Parsed NA in external outcome, inspect raw values in: ", ext_path)
  y01_ext <- as.integer(ext0$outcome == pos_label)
  
  ext_prob <- predict_prob_questionnaire(final_mod, ext0)
  
  cat("\n=== EXTERNAL METRICS (fixed thr_prob from DEV OOF) ===\n")
  print(metrics_from_prob(ext_prob, y01_ext, thr_prob), row.names = FALSE)
}

# -------------------------
# 6) Optional: print bin spec used by this model
# -------------------------
cat("\n=== BIN THRESHOLDS (fixed) ===\n")
cat("Age breaks: ", paste(age_breaks, collapse = " , "), "\n", sep = "")
cat("SpO2 breaks: ", paste(spo2_breaks, collapse = " , "), "\n", sep = "")
cat("Heart rate breaks: ", paste(hr_breaks, collapse = " , "), "\n", sep = "")
cat("Respiratory rate breaks: ", paste(rr_breaks, collapse = " , "), "\n", sep = "")
cat("Temperature breaks, Celsius: ", paste(t_breaks, collapse = " , "), "\n", sep = "")
cat("SBP breaks: ", paste(sbp_breaks, collapse = " , "), "\n", sep = "")
cat("DBP breaks: ", paste(dbp_breaks, collapse = " , "), "\n", sep = "")
cat("MAP breaks: ", paste(map_breaks, collapse = " , "), "\n", sep = "")
cat("Sex mapping: Female=0, Male=1\n")

# -------------------------
# 6b) POINTS TABLE (derived from final DEV model coefficients)
# -------------------------
build_points_table <- function(model, sf_target = 10) {
  fit <- model$fit
  b <- stats::coef(fit)
  b <- b[names(b) != "(Intercept)"]
  b <- b[is.finite(b)]
  if (!length(b)) stop("No finite coefficients found for points table")

  coef_or0 <- function(nm) {
    if (!nm %in% names(b)) return(0)
    as.numeric(b[[nm]])
  }

  c_sex  <- coef_or0("sex_o")
  c_age  <- coef_or0("age_o")
  c_spo2 <- coef_or0("spo2_o")
  c_hr   <- coef_or0("hr_o")
  c_rr   <- coef_or0("rr_o")
  c_temp <- coef_or0("temp_o")
  c_sbp  <- coef_or0("sbp_o")
  c_dbp  <- coef_or0("dbp_o")
  c_map  <- coef_or0("map_o")

  used <- c(c_sex, c_age, c_spo2, c_hr, c_rr, c_temp, c_sbp, c_dbp, c_map)
  used <- used[is.finite(used) & used != 0]
  if (!length(used)) {
    stop("No non-zero coefficients found for points table. Coefs are: ",
         paste(names(b), collapse = ", "))
  }

  sf <- sf_target / stats::median(abs(used))

  make_rows_ord <- function(attr, coef1, ord_map_named_int) {
    levs <- names(ord_map_named_int)
    codes <- as.integer(ord_map_named_int[levs])
    ord <- order(codes, decreasing = FALSE)  # low risk -> high risk
    levs <- levs[ord]
    codes <- codes[ord]
    pts <- as.integer(round(coef1 * sf * codes))
    pts <- pts - min(pts, na.rm = TRUE)
    data.frame(Attribute = rep(attr, length(levs)),
               Level = levs,
               Points = pts,
               stringsAsFactors = FALSE)
  }

  sex_tbl <- data.frame(
    Attribute = c("Sex, male=1", "Sex, male=1"),
    Level = c("Female", "Male"),
    Points = as.integer(round(c_sex * sf * c(0, 1))),
    stringsAsFactors = FALSE
  )

  out <- rbind(
    sex_tbl,
    make_rows_ord("Age", c_age, model$ord_map$age_b),
    make_rows_ord("SpO2", c_spo2, model$ord_map$spo2_b),
    make_rows_ord("Heart rate", c_hr, model$ord_map$hr_b),
    make_rows_ord("Respiratory rate", c_rr, model$ord_map$rr_b),
    make_rows_ord("Temperature, Celsius", c_temp, model$ord_map$temp_b),
    make_rows_ord("SBP", c_sbp, model$ord_map$sbp_b),
    make_rows_ord("DBP", c_dbp, model$ord_map$dbp_b),
    make_rows_ord("MAP", c_map, model$ord_map$map_b)
  )

  # rebase each attribute to start at 0
  for (a in unique(out$Attribute)) {
    idx <- out$Attribute == a
    out$Points[idx] <- out$Points[idx] - min(out$Points[idx], na.rm = TRUE)
  }

  attr(out, "sf") <- sf
  out
}


pt <- build_points_table(final_mod, sf_target = 10)
cat("\n=== POINTS TABLE (derived from final DEV model, monotonic) ===\n")
cat("Scale factor sf = ", format(attr(pt, "sf"), digits = 6), "\n", sep = "")
print(pt, row.names = FALSE)

max_score <- sum(tapply(pt$Points, pt$Attribute, max))
cat("Max total score = ", max_score, "\n", sep = "")



# -------------------------
# 6c) SCORE THRESHOLDS IN POINTS + QUESTIONNAIRE SPEC
# -------------------------
score_points_from_bins <- function(df_bins, pt_tbl) {
  pt_map <- setNames(as.integer(pt_tbl$Points), paste(pt_tbl$Attribute, pt_tbl$Level, sep = "|"))
  take <- function(attr, level) {
    k <- paste(attr, as.character(level), sep = "|")
    v <- pt_map[[k]]
    if (is.null(v) || is.na(v)) 0L else as.integer(v)
  }

  vapply(df_bins$sex_b,  function(z) take("Sex, male=1", z), integer(1)) +
    vapply(df_bins$age_b,  function(z) take("Age", z), integer(1)) +
    vapply(df_bins$spo2_b, function(z) take("SpO2", z), integer(1)) +
    vapply(df_bins$hr_b,   function(z) take("Heart rate", z), integer(1)) +
    vapply(df_bins$rr_b,   function(z) take("Respiratory rate", z), integer(1)) +
    vapply(df_bins$temp_b, function(z) take("Temperature, Celsius", z), integer(1)) +
    vapply(df_bins$sbp_b,  function(z) take("SBP", z), integer(1)) +
    vapply(df_bins$dbp_b,  function(z) take("DBP", z), integer(1)) +
    vapply(df_bins$map_b,  function(z) take("MAP", z), integer(1))
}

metrics_from_score <- function(score, y01, thr) {
  pred1 <- as.integer(score >= thr)
  TP <- sum(pred1 == 1 & y01 == 1)
  TN <- sum(pred1 == 0 & y01 == 0)
  FP <- sum(pred1 == 1 & y01 == 0)
  FN <- sum(pred1 == 0 & y01 == 1)

  mcc <- mcc_from_counts(c(TP = TP, TN = TN, FP = FP, FN = FN))
  prec <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  f1   <- if (!is.finite(prec) || !is.finite(sens) || (prec + sens) == 0) NA_real_ else 2 * prec * sens / (prec + sens)
  auc  <- auc_fast(score, y01)

  data.frame(
    thr = thr, TP = TP, TN = TN, FP = FP, FN = FN,
    MCC = mcc, AUC = auc, F1 = f1,
    Accuracy = acc, Precision = prec, Sensitivity = sens, Specificity = spec,
    stringsAsFactors = FALSE
  )
}

# DEV: compute points from the same fixed bins used everywhere
dev_bins_full <- make_bins(dev0)
dev_score_points <- score_points_from_bins(dev_bins_full, pt)

# For the 3-zone thresholds we map points -> mean predicted risk (from final monotonic model)
dev_prob_final <- predict_prob_questionnaire(final_mod, dev0)
agg_pts <- stats::aggregate(dev_prob_final ~ dev_score_points, FUN = mean)
agg_pts <- agg_pts[order(agg_pts$dev_score_points), , drop = FALSE]

s_no  <- suppressWarnings(max(agg_pts$dev_score_points[agg_pts$dev_prob_final < risk_low]))
s_yes <- suppressWarnings(min(agg_pts$dev_score_points[agg_pts$dev_prob_final >= risk_high]))
if (!is.finite(s_no))  s_no  <- min(agg_pts$dev_score_points)
if (!is.finite(s_yes)) s_yes <- max(agg_pts$dev_score_points)

# Binary MCC threshold directly in points space
cut_mcc_pts <- best_cut_mcc(dev_score_points, y01_dev)

cat("\n=== SCORE THRESHOLDS (points) ===\n")
cat("Three zone, No if score <= ", s_no, ", Yes if score >= ", s_yes, "\n", sep = "")
cat("Binary MCC, Yes if score >= ", cut_mcc_pts$cut, "\n", sep = "")
cat("Decision mode used: ", decision_mode, "\n", sep = "")

questionnaire_spec <- list(
  breaks = list(
    age  = age_breaks,
    spo2 = spo2_breaks,
    hr   = hr_breaks,
    rr   = rr_breaks,
    temp = t_breaks,
    sbp  = sbp_breaks,
    dbp  = dbp_breaks,
    map  = map_breaks
  ),
  points = pt,
  thresholds = list(
    three_zone_no  = s_no,
    three_zone_yes = s_yes,
    binary_mcc_yes = cut_mcc_pts$cut
  ),
  impute = list(
    sex_mode = final_mod$impute$sex_mode,
    medians  = final_mod$impute$medians
  ),
  decision_mode = decision_mode,
  labels = list(pos = pos_label, neg = neg_label)
)

apply_questionnaire <- function(new_df, spec = questionnaire_spec) {
  x <- new_df

  # coerce
  x$age_raw        <- coerce_numeric_robust(x$age_raw)
  x$sex_raw        <- coerce_sex_to_01(x$sex_raw)
  x$spo2_raw       <- coerce_numeric_robust(x$spo2_raw)
  x$heart_rate_raw <- coerce_numeric_robust(x$heart_rate_raw)
  x$resp_rate_raw  <- coerce_numeric_robust(x$resp_rate_raw)
  x$temp_raw       <- coerce_numeric_robust(x$temp_raw)
  x$sbp_raw        <- coerce_numeric_robust(x$sbp_raw)
  x$dbp_raw        <- coerce_numeric_robust(x$dbp_raw)
  x$map_raw        <- coerce_numeric_robust(x$map_raw)

  tm <- suppressWarnings(stats::median(x$temp_raw, na.rm = TRUE))
  if (is.finite(tm) && tm > 60) x$temp_raw <- (x$temp_raw - 32) * 5/9

  # clamp
  x$age_raw        <- clamp_na(x$age_raw,          0, 120)
  x$spo2_raw       <- clamp_na(x$spo2_raw,        50, 100, zero_is_na = TRUE)
  x$heart_rate_raw <- clamp_na(x$heart_rate_raw,  20, 250)
  x$resp_rate_raw  <- clamp_na(x$resp_rate_raw,    5,  80)
  x$temp_raw       <- clamp_na(x$temp_raw,        30,  43)
  x$sbp_raw        <- clamp_na(x$sbp_raw,         50, 250)
  x$dbp_raw        <- clamp_na(x$dbp_raw,         30, 150)
  x$map_raw        <- clamp_na(x$map_raw,         40, 200)

  # impute using DEV
  x$sex_raw[is.na(x$sex_raw)] <- spec$impute$sex_mode
  for (nm in names(spec$impute$medians)) {
    x[[nm]][is.na(x[[nm]])] <- spec$impute$medians[[nm]]
  }

  # bins
  x$age_b  <- cut3_fixed(x$age_raw,        spec$breaks$age,  c("<50", "50 to <68", ">=68"))
  x$sex_b  <- factor(ifelse(x$sex_raw == 1, "Male", "Female"), levels = c("Female","Male"))
  x$spo2_b <- cut3_fixed(x$spo2_raw,       spec$breaks$spo2, c("<92", "92 to <96", ">=96"))
  x$hr_b   <- cut3_fixed(x$heart_rate_raw, spec$breaks$hr,   c("<90", "90 to <110", ">=110"))
  x$rr_b   <- cut3_fixed(x$resp_rate_raw,  spec$breaks$rr,   c("<20", "20 to <30", ">=30"))
  x$temp_b <- cut3_fixed(x$temp_raw,       spec$breaks$temp, c("<36", "36 to <38", ">=38"))
  x$sbp_b  <- cut3_fixed(x$sbp_raw,        spec$breaks$sbp,  c("<100", "100 to <140", ">=140"))
  x$dbp_b  <- cut3_fixed(x$dbp_raw,        spec$breaks$dbp,  c("<60", "60 to <90", ">=90"))
  x$map_b  <- cut3_fixed(x$map_raw,        spec$breaks$map,  c("<65", "65 to <100", ">=100"))

  x$score_points <- score_points_from_bins(x, spec$points)

  neg <- spec$labels$neg
  pos <- spec$labels$pos
  mode <- spec$decision_mode

  if (mode == "three_zone") {
    x$decision <- ifelse(
      x$score_points <= spec$thresholds$three_zone_no, neg,
      ifelse(x$score_points >= spec$thresholds$three_zone_yes, pos, "Indeterminate")
    )
  } else if (mode == "binary_low") {
    x$decision <- ifelse(x$score_points <= spec$thresholds$three_zone_no, neg, pos)
  } else if (mode == "binary_high") {
    x$decision <- ifelse(x$score_points >= spec$thresholds$three_zone_yes, pos, neg)
  } else {
    x$decision <- ifelse(x$score_points >= spec$thresholds$binary_mcc_yes, pos, neg)
  }

  x
}

# DEV: decisions + binary metrics using the points MCC cut
dev_decision <- ifelse(dev_score_points >= cut_mcc_pts$cut, pos_label, neg_label)

cat("\n=== DEV DECISION COUNTS (points) ===\n")
print(table(dev_decision, useNA = "ifany"))

cat("\n=== DEV METRICS (binary by points, MCC cut) ===\n")
print(metrics_from_score(dev_score_points, y01_dev, cut_mcc_pts$cut), row.names = FALSE)

if (nzchar(ext_path)) {
  ext_bins_full <- make_bins(ext0)
  ext_score_points <- score_points_from_bins(ext_bins_full, pt)
  ext_decision <- ifelse(ext_score_points >= cut_mcc_pts$cut, pos_label, neg_label)

  cat("\n=== EXTERNAL DECISION COUNTS (points) ===\n")
  print(table(ext_decision, useNA = "ifany"))

  cat("\n=== EXTERNAL METRICS (binary by points, MCC cut) ===\n")
  print(metrics_from_score(ext_score_points, y01_ext, cut_mcc_pts$cut), row.names = FALSE)
}
