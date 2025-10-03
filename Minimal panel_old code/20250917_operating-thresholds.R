# ================================
# Confusion matrices @ MCC-tuned threshold (time split)
# KEEP ONLY: Group %in% {"CTRL_noCOVID","COVID"}
# BOUNDARY column: "Sampling _date" (no other column used)
# PRINT-ONLY
# ================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))

# ---------- CONFIG ----------
FILE_PATH <- "biomarkers_acuteCOVID_meta.xlsx"
SHEET     <- "meta"           # set to 1 if needed
OUTCOME   <- "Hospital_ID"
POS       <- "Yes"
NEG       <- "No"
BOUNDARY  <- as.Date("2020-04-21")   # <-- boundary date for time split
N_BOOT    <- 1000                    # bootstrap reps for F1/MCC CIs
set.seed(444)

# Feature sets (albumin/CRP/D_Dimer removed from Full)
feat_full <- c(
  "Diagnosis","severity_admission","Age","Gender",
  "SpO2_admission",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc"
)
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------- HELPERS ----------
stop_if_missing <- function(df, nm) if (!(nm %in% names(df))) stop(sprintf("Column '%s' not found.", nm))

parse_excel_date <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x))         return(as.Date(x, origin = "1899-12-30"))  # Excel epoch
  if (is.character(x))       return(as.Date(x))
  as.Date(x)
}

to_YN <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(factor(ifelse(x, POS, NEG), levels = c(POS, NEG)))
  if (is.numeric(x)) {
    if (!all(x %in% c(0,1), na.rm = TRUE)) stop("Outcome numeric but not 0/1.")
    return(factor(ifelse(x==1, POS, NEG), levels = c(POS, NEG)))
  }
  if (is.character(x)) {
    s <- trimws(tolower(x))
    map_yes <- c("1","yes","y","true","pos","positive")
    map_no  <- c("0","no","n","false","neg","negative")
    out <- ifelse(s %in% map_yes, POS, ifelse(s %in% map_no, NEG, NA_character_))
    if (any(is.na(out))) stop("Outcome cannot be coerced to Yes/No.")
    return(factor(out, levels = c(POS, NEG)))
  }
  stop("Unsupported outcome type.")
}

mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(NA_real_)
  num/den
}

best_threshold_mcc <- function(y, p_pos, pos = POS, neg = NEG, grid = seq(0.01,0.99,by=0.001)) {
  yy <- factor(y, levels = c(neg, pos))
  best_t <- 0.5; best_m <- -Inf
  for (t in grid) {
    pred <- factor(ifelse(p_pos >= t, pos, neg), levels = c(neg, pos))
    tab <- table(yy, pred)
    TP <- as_num(tab[pos, pos] %||% 0); TN <- as_num(tab[neg, neg] %||% 0)
    FP <- as_num(tab[neg, pos] %||% 0); FN <- as_num(tab[pos, neg] %||% 0)
    m  <- mcc_from_counts(TP, FP, FN, TN)
    if (is.finite(m) && m > best_m) { best_m <- m; best_t <- t }
  }
  list(t = best_t, mcc = best_m)
}

binom_ci_wilson <- function(x, n, conf = 0.95) {
  if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1-conf)/2)
  p <- x/n
  denom <- 1 + z^2/n
  center <- (p + z^2/(2*n)) / denom
  half   <- z * sqrt((p*(1-p) + z^2/(4*n)) / n) / denom
  c(center - half, center + half)
}

compute_counts_and_rates <- function(obs, pred, p_pos, pos = POS, neg = NEG) {
  y    <- factor(obs,  levels = c(neg, pos))
  pr   <- factor(pred, levels = c(neg, pos))
  tab  <- table(y, pr)
  TP <- as_num(tab[pos, pos] %||% 0)
  TN <- as_num(tab[neg, neg] %||% 0)
  FP <- as_num(tab[neg, pos] %||% 0)
  FN <- as_num(tab[pos, neg] %||% 0)
  n  <- TP + TN + FP + FN
  prec <- if ((TP+FP)==0) NA_real_ else TP/(TP+FP)
  sens <- if ((TP+FN)==0) NA_real_ else TP/(TP+FN)
  spec <- if ((TN+FP)==0) NA_real_ else TN/(TN+FP)
  acc  <- (TP+TN)/n
  f1   <- if (is.na(prec) || is.na(sens) || (prec+sens)==0) NA_real_ else 2*prec*sens/(prec+sens)
  mcc  <- mcc_from_counts(TP, FP, FN, TN)
  # AUC
  y01  <- as.integer(y == pos)
  aucv <- if (length(unique(y01)) < 2 || all(!is.finite(p_pos))) NA_real_ else suppressMessages(as.numeric(auc(y01, p_pos, quiet = TRUE)))
  list(TP=TP, FP=FP, FN=FN, TN=TN, Precision=prec, Sensitivity=sens, Specificity=spec, Accuracy=acc, F1=f1, MCC=mcc, AUC=aucv, N=n)
}

auc_delong_ci <- function(obs, p_pos, conf.level = 0.95) {
  y <- as.integer(obs == POS)
  if (length(unique(y)) < 2 || all(!is.finite(p_pos))) return(c(NA_real_, NA_real_, NA_real_))
  suppressMessages(as.numeric(ci.auc(y, p_pos, conf.level = conf.level)))
}

boot_metric_ci <- function(obs, pred, p_pos, metric = c("F1","MCC"), n_boot = N_BOOT, conf.level = 0.95) {
  metric <- match.arg(metric)
  n <- length(obs)
  if (n <= 1) return(c(NA_real_, NA_real_))
  getm <- function(o, pr, pp) {
    cr <- compute_counts_and_rates(o, pr, pp)
    if (metric == "F1") cr$F1 else cr$MCC
  }
  vals <- replicate(n_boot, {
    idx <- sample.int(n, replace = TRUE)
    getm(obs[idx], pred[idx], p_pos[idx])
  })
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(NA_real_, NA_real_))
  qs <- quantile(vals, probs = c((1-conf.level)/2, 1-(1-conf.level)/2), names = FALSE, na.rm = TRUE)
  as.numeric(qs)
}

fmt_ci <- function(lo, hi, digits = 3) {
  paste0("[", ifelse(is.finite(lo), sprintf(paste0("%.",digits,"f"), lo), "NA"),
         ", ", ifelse(is.finite(hi), sprintf(paste0("%.",digits,"f"), hi), "NA"), "]")
}

# --- recipe (avoid pipes with RHS blocks) ---
make_recipe <- function(dat, yvar) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c(".rid","data_split"), names(dat))
  if (length(ign)) rec <- recipes::update_role(rec, tidyselect::all_of(ign), new_role = "ignore")
  rec <- recipes::step_impute_median(rec, recipes::all_numeric_predictors())
  rec <- recipes::step_impute_mode(rec,   recipes::all_nominal_predictors())
  rec <- recipes::step_novel(rec,         recipes::all_nominal_predictors())
  rec <- recipes::step_other(rec,         recipes::all_nominal_predictors(), threshold = 0.01)
  rec <- recipes::step_dummy(rec,         recipes::all_nominal_predictors())
  rec <- recipes::step_zv(rec,            recipes::all_predictors())
  # Yeo–Johnson if available
  if ("step_YeoJohnson" %in% ls(getNamespace("recipes"))) {
    rec <- recipes::step_YeoJohnson(rec, recipes::all_numeric_predictors())
  } else if ("step_yeojohnson" %in% ls(getNamespace("recipes"))) {
    rec <- recipes::step_yeojohnson(rec, recipes::all_numeric_predictors())
  }
  rec <- recipes::step_normalize(rec, recipes::all_numeric_predictors())
  rec
}

# --- model availability guard ---
method_ready <- function(method) {
  req <- switch(method,
                "J48"       = "RWeka",
                "kknn"      = "kknn",
                "svmRadial" = "kernlab",
                "rf"        = "randomForest",
                "glmnet"    = "glmnet",
                NULL)
  if (is.null(req)) TRUE else requireNamespace(req, quietly = TRUE)
}

# --- train + eval on test, threshold tuned on train by MCC ---
train_and_eval <- function(train, test, yvar, features, algo_name) {
  method <- switch(algo_name,
                   "C4.5" = "J48",
                   "k-NN" = "kknn",
                   "SVM"  = "svmRadial",
                   "RF"   = "rf",
                   "LR"   = "glmnet",
                   stop("Unknown algo: ", algo_name))
  if (!method_ready(method)) {
    message(sprintf("[skip] %s (%s) not available; missing package.", algo_name, method))
    return(NULL)
  }
  
  rec <- make_recipe(train[, unique(c(features, yvar, "data_split")), drop = FALSE], yvar)
  
  tr_ctrl <- caret::trainControl(
    method = "cv", number = 5,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", allowParallel = TRUE
  )
  
  spec <- list(method = method, tuneLength = 5)
  if (identical(method, "glmnet")) {
    spec$grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 0, length.out = 5))
  }
  
  fit <- if (!is.null(spec$grid)) {
    caret::train(rec, data = train, method = spec$method,
                 trControl = tr_ctrl, metric = "ROC", tuneGrid = spec$grid)
  } else {
    caret::train(rec, data = train, method = spec$method,
                 trControl = tr_ctrl, metric = "ROC", tuneLength = spec$tuneLength)
  }
  
  p_tr <- as.numeric(predict(fit, newdata = train, type = "prob")[, POS])
  thr  <- best_threshold_mcc(train[[yvar]], p_tr, pos = POS, neg = NEG)$t
  
  p_te <- as.numeric(predict(fit, newdata = test, type = "prob")[, POS])
  pred_te <- factor(ifelse(p_te >= thr, POS, NEG), levels = c(POS, NEG))
  list(prob = p_te, pred = pred_te, threshold = thr)
}

# --- build one result row with metrics + CIs ---
row_from_eval <- function(y_true, eval, algo, set_name) {
  cr <- compute_counts_and_rates(y_true, eval$pred, eval$prob)
  
  # Wilson CIs
  sens_ci <- binom_ci_wilson(cr$TP, cr$TP + cr$FN)
  spec_ci <- binom_ci_wilson(cr$TN, cr$TN + cr$FP)
  acc_ci  <- binom_ci_wilson(cr$TP + cr$TN, cr$N)
  prec_ci <- binom_ci_wilson(cr$TP, cr$TP + cr$FP)
  
  # Bootstrap CIs
  f1_ci   <- boot_metric_ci(y_true, eval$pred, eval$prob, metric = "F1")
  mcc_ci  <- boot_metric_ci(y_true, eval$pred, eval$prob, metric = "MCC")
  
  # AUC + DeLong CI
  auc_ci  <- auc_delong_ci(y_true, eval$prob)
  
  tibble::tibble(
    TP = cr$TP, FP = cr$FP, FN = cr$FN, TN = cr$TN,
    
    # metrics in required order with CI strings
    MCC        = cr$MCC,
    MCC_CI     = fmt_ci(mcc_ci[1], mcc_ci[2]),
    
    AUC_ROC    = cr$AUC,
    AUC_ROC_CI = fmt_ci(auc_ci[1], auc_ci[3]),
    
    F1         = cr$F1,
    F1_CI      = fmt_ci(f1_ci[1], f1_ci[2]),
    
    Accuracy   = cr$Accuracy,
    Accuracy_CI= fmt_ci(acc_ci[1], acc_ci[2]),
    
    Precision  = cr$Precision,
    Precision_CI = fmt_ci(prec_ci[1], prec_ci[2]),
    
    Sensitivity  = cr$Sensitivity,
    Sensitivity_CI = fmt_ci(sens_ci[1], sens_ci[2]),
    
    Specificity  = cr$Specificity,
    Specificity_CI = fmt_ci(spec_ci[1], spec_ci[2]),
    
    Threshold = eval$threshold,
    Algorithm = algo,
    Feature_Set = set_name,
    Split = "Test"
  )
}

# ---------- LOAD & FILTER ----------
df <- readxl::read_excel(FILE_PATH, sheet = SHEET) |> as.data.frame()

stop_if_missing(df, "Group")
stop_if_missing(df, "Sampling _date")
stop_if_missing(df, OUTCOME)

# keep EXACT groups (no remapping)
df <- df[df$Group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
df$Group <- droplevels(factor(df$Group))

# parse ONLY boundary column from 'Sampling _date'
df$Sampling_date_boundary <- parse_excel_date(df[["Sampling _date"]])
drop_n <- sum(is.na(df$Sampling_date_boundary))
if (drop_n > 0) {
  message("[warn] Dropping ", drop_n, " rows with missing Sampling _date.")
  df <- df[!is.na(df$Sampling_date_boundary), , drop = FALSE]
}

# split by boundary
df$data_split <- ifelse(df$Sampling_date_boundary < BOUNDARY, "train", "test")
df$data_split <- factor(df$data_split, levels = c("train","test"))

# normalize outcome to {Yes, No}, with Yes FIRST level (caret expects eventLevel="first")
df[[OUTCOME]] <- to_YN(df[[OUTCOME]])

cat("\n[Sanity] Using ONLY boundary column: 'Sampling _date'\n")
cat("[Sanity] Kept groups exactly: CTRL_noCOVID, COVID\n")
cat("[Sanity] Split counts (by Sampling _date < ", format(BOUNDARY), "):\n", sep = "")
print(table(df$data_split))

# ---------- TRAIN (pre) / TEST (post) ----------
dtr <- df[df$data_split == "train", , drop = FALSE]
dte <- df[df$data_split == "test",  , drop = FALSE]
if (!nrow(dtr) || !nrow(dte)) stop("Empty train/test after boundary split.")

algos <- c("C4.5","k-NN","SVM","RF","LR")
sets  <- list(Full = feat_full, Triage = feat_triage)

res_rows <- list()
set.seed(444)

for (set_name in names(sets)) {
  feats <- sets[[set_name]]
  feats <- feats[feats %in% names(df)]  # use only available predictors
  if (!length(feats)) { message("[warn] No predictors found for set: ", set_name); next }
  for (algo in algos) {
    ok <- try({
      ev <- train_and_eval(dtr[, unique(c(feats, OUTCOME, "data_split")), drop = FALSE],
                           dte[, unique(c(feats, OUTCOME, "data_split")), drop = FALSE],
                           yvar = OUTCOME, features = feats, algo_name = algo)
      if (!is.null(ev)) {
        res_rows[[paste(set_name, algo, sep = "::")]] <- row_from_eval(dte[[OUTCOME]], ev, algo, set_name)
      }
    }, silent = TRUE)
    if (inherits(ok, "try-error")) message("[skip] ", algo, " failed on set ", set_name, ": ", as.character(ok))
  }
}

out_tbl <- dplyr::bind_rows(res_rows)

# Ensure required column order:
desired_cols <- c(
  "TP","FP","FN","TN",
  "MCC","MCC_CI",
  "AUC_ROC","AUC_ROC_CI",
  "F1","F1_CI",
  "Accuracy","Accuracy_CI",
  "Precision","Precision_CI",
  "Sensitivity","Sensitivity_CI",
  "Specificity","Specificity_CI",
  "Threshold","Algorithm","Feature_Set","Split"
)
out_tbl <- out_tbl[, intersect(desired_cols, names(out_tbl)), drop = FALSE]

# ---------- PRINT ----------
cat("\n[Results] Confusion matrices at MCC-tuned thresholds (time split by Sampling _date)\n\n")
if (nrow(out_tbl)) {
  out_tbl <- out_tbl %>%
    mutate(Algorithm = factor(Algorithm, levels = algos[algos %in% Algorithm])) %>%
    arrange(Feature_Set, Algorithm)
  print(out_tbl, n = nrow(out_tbl))
} else {
  cat("(no rows)\n")
}
