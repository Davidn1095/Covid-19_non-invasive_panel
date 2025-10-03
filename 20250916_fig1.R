# ===============================================================
# COMPLETE: Load data -> temporal split -> CV -> Heatmap + Tables
# Prints heatmap plot, sensitivity (Δ all metrics), and full perf table
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(kknn)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))

# ---------------- config ----------------
cfg <- list(
  file_path     = "biomarkers_acuteCOVID_meta.xlsx",
  sheet         = "meta",
  outcome       = "Hospital_ID",
  pos_label     = "Yes",
  neg_label     = "No",
  boundary_date = as.Date("2020-04-21"),
  cv_k          = 3,
  inner_k       = 3,
  tune_len      = 3,
  seed_cv       = 444
)

# ---------------- features ----------------
feat_full   <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission",
                 "monocyte_abs_number","monocytes_perc",
                 "lymphocyte_abs_number","lymphocytes_perc",
                 "neutrophil_abs_number","neutrophils_perc")
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- labels, orders, palette ----------------
metric_label_map <- c(
  "MCC"         = "MCC",
  "AUC"         = "AUC-ROC",
  "F1"          = "F1-score",
  "Accuracy"    = "Accuracy",
  "Precision"   = "Precision",
  "Sensitivity" = "Sensitivity",
  "Specificity" = "Specificity"
)
metric_order <- c("MCC","AUC-ROC","F1-score","Accuracy","Precision","Sensitivity","Specificity")
algo_order   <- c("C4.5","k-NN","SVM","RF","LR")
pal_red_green <- c("#F4CCCC", "#D9EAD3")

# ---------------- helpers: coercions & dates ----------------
to_YN <- function(x, pos = cfg$pos_label, neg = cfg$neg_label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(factor(ifelse(x, pos, neg), levels = c(pos, neg)))
  if (is.numeric(x)) {
    stopifnot(all(x %in% c(0,1), na.rm = TRUE))
    return(factor(ifelse(x == 1, pos, neg), levels = c(pos, neg)))
  }
  if (is.character(x)) {
    s <- trimws(tolower(x))
    map_yes <- c("1","yes","y","true","pos","positive")
    map_no  <- c("0","no","n","false","neg","negative")
    out <- ifelse(s %in% map_yes, pos, ifelse(s %in% map_no, neg, NA_character_))
    if (any(is.na(out))) stop("Outcome not coercible to Yes/No")
    return(factor(out, levels = c(pos, neg)))
  }
  stop("Unsupported outcome type")
}

sanitize_name <- function(x) gsub("[^A-Za-z0-9]", "", tolower(x))

parse_excel_date <- function(v) {
  if (inherits(v, "Date")) return(v)
  if (inherits(v, "POSIXct") || inherits(v, "POSIXt")) return(as.Date(v))
  if (is.numeric(v)) return(as.Date(v, origin = "1899-12-30"))
  if (is.character(v)) {
    tryers <- c("%d/%m/%Y","%Y-%m-%d","%m/%d/%Y","%d-%m-%Y","%d.%m.%Y")
    for (fmt in tryers) {
      d <- suppressWarnings(as.Date(v, format = fmt))
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
  stop("Could not find a date column (sampling/admission/date)")
}

# ---------------- data loader + filter + outcome ----------------
load_data <- function(path, sheet, outcome, keep_predictors) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet = sheet) |> as.data.frame()
  
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df) == grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  for (v in intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' missing", outcome))
  df[[outcome]] <- to_YN(df[[outcome]])
  
  dcol <- resolve_date_col(df)
  message(sprintf("[date] Using '%s' as the admission/sampling date.", dcol))
  df[[dcol]] <- parse_excel_date(df[[dcol]])
  
  keepx <- unique(c(keep_predictors, outcome, "group", dcol, "data_split"))
  keepx <- intersect(keepx, names(df))
  list(df = df[, keepx, drop = FALSE], admit_col = dcol)
}

# ---------------- temporal split ----------------
enforce_temporal_split <- function(df, admit_col, boundary_date = cfg$boundary_date) {
  stopifnot(admit_col %in% names(df))
  ad <- parse_excel_date(df[[admit_col]])
  df$data_split <- factor(ifelse(ad <= boundary_date, "train", "external"), levels = c("train","external"))
  message("\n[Temporal split]\n  Boundary date: ", format(boundary_date), "\n")
  print(table(df$data_split, useNA = "ifany"))
  df
}

# ---------------- recipe ----------------
make_recipe <- function(dat, yvar) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c("data_split"), names(dat))
  if (length(ign)) rec <- rec |> update_role(all_of(ign), new_role = "ignore")
  rec |>
    step_impute_median(all_numeric_predictors()) |>
    step_impute_mode(all_nominal_predictors())  |>
    step_novel(all_nominal_predictors())        |>
    step_other(all_nominal_predictors(), threshold = 0.01) |>
    step_dummy(all_nominal_predictors())        |>
    step_zv(all_predictors())                   |>
    step_normalize(all_numeric_predictors())
}

# ---------------- metrics ----------------
mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(NA_real_)
  num/den
}
get_neg_level <- function(y) { levs <- levels(y); nl <- setdiff(levs, cfg$pos_label); if (length(nl)) nl[1] else cfg$neg_label }
best_threshold_mcc <- function(obs, p_pos, pos = cfg$pos_label, neg = cfg$neg_label, grid = seq(0.01,0.99,by=0.001)) {
  y <- factor(obs, levels = c(neg, pos))
  best_t <- 0.5; best_m <- -Inf
  for (t in grid) {
    pred <- factor(ifelse(p_pos >= t, pos, neg), levels = c(neg,pos))
    tab <- table(y,pred)
    TP <- as_num(tab[pos,pos] %||% 0); TN <- as_num(tab[neg,neg] %||% 0)
    FP <- as_num(tab[neg,pos] %||% 0); FN <- as_num(tab[pos,neg] %||% 0)
    m  <- mcc_from_counts(TP, FP, FN, TN)
    if (is.finite(m) && m > best_m) { best_m <- m; best_t <- t }
  }
  list(t = best_t, mcc = best_m)
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
  aucv  <- binary_auc(y, as_num(p_pos), pos_label)
  c(MCC = mcc, AUC = aucv, `F1` = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec)
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
model_specs <- function(tune_len, algos = c("LR","RF","SVM","k-NN","C4.5")) {
  all <- list(
    "LR"   = list(name="LR",   method = "glmnet",    tuneLength = tune_len),
    "RF"   = list(name="RF",   method = "rf",        tuneLength = tune_len),
    "SVM"  = list(name="SVM",  method = "svmRadial", tuneLength = tune_len),
    "k-NN" = list(name="k-NN", method = "kknn",      tuneLength = tune_len),
    "C4.5" = list(name="C4.5", method = "J48",       tuneLength = tune_len)
  )
  avail <- names(caret::getModelInfo())
  keep <- names(all)[vapply(all, \(s) s$method %in% avail && method_available(s$method), logical(1))]
  out  <- all[keep]
  out[names(out)[order(match(names(out), algo_order))]]
}

# ---------------- CV utilities ----------------
choose_k_for_task <- function(y, k_desired) {
  if (is.null(y) || !length(y) || anyNA(y)) return(k_desired)
  tab <- table(y); max_k <- max(2, min(tab))
  min(as.integer(k_desired), as.integer(max_k))
}
build_cv_splits <- function(y, R = 1, k_desired = cfg$cv_k, seed = cfg$seed_cv) {
  set.seed(seed)
  kk <- choose_k_for_task(y, k_desired)
  caret::createFolds(y, k = kk, list = TRUE, returnTrain = FALSE)
}
train_inner <- function(dat_train, yvar, spec, inner_k) {
  rec <- make_recipe(dat_train, yvar)
  tr_ctrl <- caret::trainControl(
    method = "cv", number = inner_k,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", allowParallel = TRUE
  )
  if (!is.null(spec$grid)) {
    caret::train(rec, data = dat_train, method = spec$method,
                 trControl = tr_ctrl, metric = "ROC", tuneGrid = spec$grid)
  } else {
    caret::train(rec, data = dat_train, method = spec$method,
                 trControl = tr_ctrl, metric = "ROC", tuneLength = spec$tuneLength)
  }
}
eval_fold_binary <- function(dat_tr, dat_va, spec, yvar, pos_label, algo_name, fold_tag, inner_k) {
  fit <- train_inner(dat_tr, yvar, spec, inner_k)
  p_va <- as.numeric(predict(fit, newdata = dat_va, type = "prob")[, pos_label])
  neg_label <- get_neg_level(dat_va[[yvar]])
  thr  <- best_threshold_mcc(dat_va[[yvar]], p_va, pos = pos_label, neg = neg_label)
  pred <- factor(ifelse(p_va >= thr$t, pos_label, neg_label), levels = levels(dat_va[[yvar]]))
  mets <- compute_metrics_binary(dat_va[[yvar]], pred, p_va, pos_label = pos_label, neg_label = neg_label)
  tibble(
    Fold = fold_tag, Algorithm = algo_name,
    Metric = names(mets), Value = as.numeric(mets)
  )
}
run_cv <- function(df, yvar, features, algos,
                   pos_label = cfg$pos_label, inner_k = cfg$inner_k,
                   k_desired = cfg$cv_k) {
  use_cols <- unique(c(features, yvar, "data_split"))
  dat <- df[df$data_split == "train", intersect(use_cols, names(df)), drop = FALSE]
  dat <- dat[complete.cases(dat[[yvar]]), , drop = FALSE]
  y <- dat[[yvar]]
  splits <- build_cv_splits(y, k_desired = k_desired)
  specs  <- model_specs(cfg$tune_len, algos = algos)
  out <- list()
  for (algo in names(specs)) {
    spec <- specs[[algo]]
    i <- 1
    for (idx_va in splits) {
      dat_va <- dat[idx_va, , drop = FALSE]
      dat_tr <- dat[-idx_va, , drop = FALSE]
      out[[paste0(algo, "_", i)]] <- eval_fold_binary(dat_tr, dat_va, spec, yvar, pos_label, algo, paste0("Fold", i), inner_k)
      i <- i + 1
    }
  }
  bind_rows(out)
}

# ---------------- holdout (external) ----------------
run_holdout <- function(df, yvar, features, algo_name) {
  specs <- model_specs(cfg$tune_len, algos = algo_name)
  if (!length(specs)) return(NULL)
  spec  <- specs[[algo_name]]
  use_cols <- unique(c(features, yvar, "group", "data_split"))
  use_cols <- intersect(use_cols, names(df))
  dtrain <- df[df$data_split == "train", use_cols, drop = FALSE]
  dext   <- df[df$data_split == "external", use_cols, drop = FALSE]
  if (!nrow(dtrain) || !nrow(dext)) return(NULL)
  fit <- train_inner(dtrain, yvar, spec, cfg$inner_k)
  p_tr <- as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label])
  thr  <- best_threshold_mcc(dtrain[[yvar]], p_tr, pos = cfg$pos_label, neg = cfg$neg_label)$t
  p_ex <- as.numeric(predict(fit, newdata = dext, type = "prob")[, cfg$pos_label])
  neg  <- cfg$neg_label
  pred <- factor(ifelse(p_ex >= thr, cfg$pos_label, neg), levels = c(neg, cfg$pos_label))
  list(obs = factor(dext[[yvar]], levels = c(neg, cfg$pos_label)),
       p = p_ex, pred = pred, threshold = thr)
}

# ---------------- summaries ----------------
summarise_mean_sd <- function(folds_tbl, set_name){
  folds_tbl %>%
    filter(Metric %in% names(metric_label_map)) %>%
    mutate(Metric = metric_label_map[Metric],
           Algorithm = factor(Algorithm, levels = algo_order, ordered = TRUE)) %>%
    group_by(Algorithm, Metric, .drop = FALSE) %>%
    summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Set = ifelse(set_name == "Full", "Complete feature set", "Triage feature set"))
}
agg_to_ci_table <- function(folds_tbl) {
  folds_tbl %>%
    group_by(Algorithm, Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      CI_low = quantile(Value, probs = 0.025, na.rm = TRUE),
      CI_high= quantile(Value, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# ===== helpers for CIs on the external split =====
binom_ci <- function(k, n, conf = 0.95){
  if (is.na(k) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  suppressWarnings(stats::binom.test(k, n, conf.level = conf)$conf.int[1:2])
}
counts_from_cm <- function(obs, pred, pos = cfg$pos_label, neg = cfg$neg_label){
  tab <- table(obs, pred)
  TP <- as_num(tab[pos, pos] %||% 0); TN <- as_num(tab[neg, neg] %||% 0)
  FP <- as_num(tab[neg, pos] %||% 0); FN <- as_num(tab[pos, neg] %||% 0)
  list(TP = TP, TN = TN, FP = FP, FN = FN)
}
bootstrap_ci_metric <- function(obs, pred, p, metric_name, B = 2000, conf = 0.95, seed = cfg$seed_cv + 37){
  set.seed(seed)
  n <- length(obs)
  if (n <= 1) return(c(NA_real_, NA_real_))
  vals <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    m <- compute_metrics_binary(obs[idx], pred[idx], p[idx],
                                pos_label = cfg$pos_label, neg_label = cfg$neg_label)
    as_num(m[metric_name])
  })
  stats::quantile(vals, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE, names = FALSE)
}

# ---------------- Heatmap (vertical: Complete on top, Triage below) ----------------
plot_cv_heatmap <- function(cv_full_ci, cv_tri_ci) {
  stopifnot(nrow(cv_full_ci) + nrow(cv_tri_ci) > 0)
  
  df <- dplyr::bind_rows(
    cv_full_ci %>% dplyr::mutate(Set = "Complete feature set"),
    cv_tri_ci  %>% dplyr::mutate(Set = "Triage feature set")
  ) %>%
    dplyr::mutate(
      Metric    = factor(Metric, levels = rev(metric_order), ordered = TRUE),
      Set       = factor(Set, levels = c("Complete feature set","Triage feature set"), ordered = TRUE),
      Algorithm = factor(Algorithm, levels = algo_order, ordered = TRUE),
      LabelText = sprintf("%.4f \u00B1 %.4f", pmax(pmin(Mean, 1), 0), SD)
    ) %>%
    dplyr::group_by(Set, Metric) %>%
    dplyr::mutate(is_max = Mean == max(Mean, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(df, ggplot2::aes(Algorithm, Metric, fill = Mean)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = LabelText,
                                    fontface = ifelse(is_max, "bold", "plain")),
                       size = 3) +
    ggplot2::scale_fill_gradient(limits = c(0,1),
                                 low = pal_red_green[1],
                                 high = pal_red_green[2],
                                 name = "Score") +
    ggplot2::facet_wrap(~Set, ncol = 1) +   # <- single column, vertical stacking
    ggplot2::labs(x = "Algorithm", y = "Metric") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title   = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      strip.text   = ggplot2::element_text(face = "bold"),
      panel.spacing.y = grid::unit(1.0, "lines")
    )
}


# ---------------- Sensitivity Δ all metrics ----------------
sensitivity_deltas_all <- function(cv_full_ci, cv_tri_ci){
  full_mean <- cv_full_ci %>% group_by(Metric) %>% summarise(MeanFull = mean(Mean, na.rm=TRUE), .groups="drop")
  tri_mean  <- cv_tri_ci  %>% group_by(Metric) %>% summarise(MeanTri  = mean(Mean, na.rm=TRUE), .groups="drop")
  full_join(full_mean, tri_mean, by = "Metric") %>%
    mutate(Delta = MeanTri - MeanFull) %>%
    mutate(Metric = factor(Metric, levels = metric_order, ordered = TRUE)) %>%
    arrange(Metric)
}

# ---------------- Performance-by-split with CIs (all rows) ----------------
performance_by_split_with_cis <- function(df, yvar, full_folds, tri_folds) {
  mk_cv <- function(folds) {
    agg_to_ci_table(folds) %>%
      mutate(
        Metric = metric_label_map[Metric],
        Feature_Set = "Complete feature set",
        Split = "Train-CV"
      )
  }
  cv_full <- mk_cv(full_folds) %>% mutate(Feature_Set = "Complete feature set")
  cv_tri  <- mk_cv(tri_folds)  %>% mutate(Feature_Set = "Triage feature set")
  
  make_hold <- function(set_feats, set_name){
    specs <- model_specs(cfg$tune_len)
    out <- list()
    for (algo in names(specs)) {
      o <- run_holdout(df, yvar, set_feats, algo)
      if (is.null(o)) next
      
      mets <- compute_metrics_binary(o$obs, o$pred, o$p, pos_label = cfg$pos_label, neg_label = cfg$neg_label)
      lev  <- levels(o$obs); neg <- lev[1]; pos <- lev[2]
      cnt  <- counts_from_cm(o$obs, o$pred, pos = pos, neg = neg)
      
      r    <- suppressMessages(pROC::roc(o$obs, o$p, levels = c(neg, pos), quiet = TRUE, direction = "<"))
      aucv <- as.numeric(pROC::auc(r))
      aucC <- suppressMessages(pROC::ci.auc(r, conf.level = 0.95, method = "delong"))
      auc_low  <- as.numeric(aucC[1]); auc_high <- as.numeric(aucC[3])
      
      acc_ci <- binom_ci(cnt$TP + cnt$TN, cnt$TP + cnt$TN + cnt$FP + cnt$FN)
      sen_ci <- binom_ci(cnt$TP, cnt$TP + cnt$FN)
      spe_ci <- binom_ci(cnt$TN, cnt$TN + cnt$FP)
      pre_ci <- if ((cnt$TP + cnt$FP) > 0) binom_ci(cnt$TP, cnt$TP + cnt$FP) else c(NA_real_, NA_real_)
      
      f1_ci  <- bootstrap_ci_metric(o$obs, o$pred, o$p, "F1",  B = 2000)
      mcc_ci <- bootstrap_ci_metric(o$obs, o$pred, o$p, "MCC", B = 2000)
      
      rows <- tibble::tibble(
        Algorithm   = algo,
        Feature_Set = set_name,
        Split       = "External",
        Metric      = c("MCC","AUC-ROC","F1-score","Accuracy","Precision","Sensitivity","Specificity"),
        Point       = c(as_num(mets["MCC"]), aucv, as_num(mets["F1"]),
                        as_num(mets["Accuracy"]), as_num(mets["Precision"]),
                        as_num(mets["Sensitivity"]), as_num(mets["Specificity"])),
        CI_low      = c(mcc_ci[1], auc_low, f1_ci[1], acc_ci[1], pre_ci[1], sen_ci[1], spe_ci[1]),
        CI_high     = c(mcc_ci[2], auc_high, f1_ci[2], acc_ci[2], pre_ci[2], sen_ci[2], spe_ci[2])
      )
      out[[algo]] <- rows
    }
    dplyr::bind_rows(out)
  }
  
  ext_full <- make_hold(feat_full,  "Complete feature set")
  ext_tri  <- make_hold(feat_triage,"Triage feature set")
  
  train_cv <- dplyr::bind_rows(cv_full, cv_tri) %>%
    dplyr::transmute(Algorithm, Feature_Set, Split, Metric, Point = Mean, CI_low, CI_high)
  
  all <- dplyr::bind_rows(train_cv, ext_full, ext_tri) %>%
    dplyr::mutate(
      Algorithm   = factor(Algorithm, levels = algo_order, ordered = TRUE),
      Metric      = factor(Metric, levels = metric_order, ordered = TRUE),
      Feature_Set = factor(Feature_Set, levels = c("Complete feature set","Triage feature set"), ordered = TRUE),
      Split       = factor(Split, levels = c("Train-CV","External"), ordered = TRUE)
    ) %>%
    dplyr::arrange(Feature_Set, Split, Metric, Algorithm)
  all
}

#png("20250922_performance.png", width = 2500, height = 2000, res = 300, bg = "transparent")
#print(p_hm)
#dev.off()



# ========================= DRIVER =========================
set.seed(cfg$seed_cv)

ld <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
df <- enforce_temporal_split(ld$df, admit_col = ld$admit_col, boundary_date = cfg$boundary_date)
yvar <- cfg$outcome

algos <- names(model_specs(cfg$tune_len, algos = c("C4.5","k-NN","SVM","RF","LR")))
message("\n[algorithms] ", paste(algos, collapse = ", "))

full_folds <- run_cv(df, yvar, features = feat_full,   algos = algos)
tri_folds  <- run_cv(df, yvar, features = feat_triage, algos = algos)

cv_full_ms <- summarise_mean_sd(full_folds, set_name = "Full")
cv_tri_ms  <- summarise_mean_sd(tri_folds,  set_name = "Triage")

cat("\n=== Sensitivity (Δ across ALL metrics) ===\n")
sens_tbl <- sensitivity_deltas_all(
  cv_full_ci = full_folds %>% filter(Metric %in% names(metric_label_map)) %>%
    mutate(Metric = metric_label_map[Metric]) %>%
    group_by(Algorithm, Metric) %>%
    summarise(Mean = mean(Value, na.rm=TRUE), .groups = "drop"),
  cv_tri_ci  = tri_folds %>% filter(Metric %in% names(metric_label_map)) %>%
    mutate(Metric = metric_label_map[Metric]) %>%
    group_by(Algorithm, Metric) %>%
    summarise(Mean = mean(Value, na.rm=TRUE), .groups = "drop")
)
print(sens_tbl)

cat("\n=== Heatmap (CV mean ± SD, 0–1) ===\n")
p_hm <- plot_cv_heatmap(cv_full_ms, cv_tri_ms)
print(p_hm)

cat("\n=== Performance by split with CIs — ALL ROWS ===\n")
perf_tbl <- performance_by_split_with_cis(df, yvar, full_folds, tri_folds)
print(perf_tbl, n = nrow(perf_tbl))

