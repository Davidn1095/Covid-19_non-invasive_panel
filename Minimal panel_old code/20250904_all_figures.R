# ============================================================
# COVID-19 Triage, end-to-end pipeline and figure printer
# - Full pipeline is run once, then figures are printed
# - No files are saved, all plots print to the device
# - Robust to earlier errors, PR arg, glmnet 1-col, subgroup reorder, sparse splits
# - Feature change: albumin, CRP, D_Dimer removed from feat_full
# - Metrics included everywhere: MCC, AUC-ROC, F1-score, Accuracy, Precision, Sensitivity, Specificity
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(broom)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))
set.seed(123)

# ---------------- config ----------------
cfg <- list(
  file_path = "biomarkers_acuteCOVID_meta.xlsx",
  sheet     = "meta",
  outcome   = "Hospital_ID",
  pos_label = "Yes",
  neg_label = "No",
  fast_mode = TRUE,
  cv_k      = 10,
  cv_R      = 1,
  inner_k   = 5,
  tune_len  = 10,
  seed_cv   = 444,
  n_boot_holdout = 250,
  n_boot_subgrp  = 300
)

if (cfg$fast_mode) {
  cfg$cv_k <- 3; cfg$inner_k <- 3; cfg$tune_len <- 3
  cfg$n_boot_holdout <- 250; cfg$n_boot_subgrp <- 300
}

# ---------------- predictors ----------------
# albumin, CRP, D_Dimer removed
feat_full <- c(
  "Diagnosis","severity_admission","Age","Gender",
  "SpO2_admission",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc"
)
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- data loader ----------------
load_data <- function(path, sheet, outcome, keep) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet = sheet) |> as.data.frame()
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df)==grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  df$.rid <- seq_len(nrow(df))
  for (v in intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' not found", outcome))
  to_YN <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) return(factor(ifelse(x, cfg$pos_label, cfg$neg_label), levels = c(cfg$pos_label, cfg$neg_label)))
    if (is.numeric(x)) {
      stopifnot(all(x %in% c(0,1), na.rm=TRUE))
      return(factor(ifelse(x==1,cfg$pos_label,cfg$neg_label), levels = c(cfg$pos_label, cfg$neg_label)))
    }
    if (is.character(x)) {
      s <- trimws(tolower(x))
      map_yes <- c("1","yes","y","true","pos","positive")
      map_no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% map_yes, cfg$pos_label, ifelse(s %in% map_no, cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome cannot be coerced to Yes/No")
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    stop("Unsupported outcome type")
  }
  df[[outcome]] <- to_YN(df[[outcome]])
  keepx <- unique(c(keep, outcome, "data_split","group",".rid"))
  keepx <- intersect(keepx, names(df))
  df[, keepx, drop = FALSE]
}

# ---------------- recipe ----------------
make_recipe <- function(dat, yvar) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c("data_split",".rid"), names(dat))
  if (length(ign)) rec <- rec |> update_role(all_of(ign), new_role = "ignore")
  rec <- rec |>
    step_impute_median(all_numeric_predictors()) |>
    step_impute_mode(all_nominal_predictors())  |>
    step_novel(all_nominal_predictors())        |>
    step_other(all_nominal_predictors(), threshold = 0.01) |>
    step_dummy(all_nominal_predictors())        |>
    step_zv(all_predictors())
  if ("step_YeoJohnson" %in% ls(getNamespace("recipes"))) {
    rec <- rec |> step_YeoJohnson(all_numeric_predictors())
  } else if ("step_yeojohnson" %in% ls(getNamespace("recipes"))) {
    rec <- rec |> step_yeojohnson(all_numeric_predictors())
  }
  rec |> step_normalize(all_numeric_predictors())
}

# ---------------- metrics ----------------
mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(NA_real_)
  num/den
}

get_neg_level <- function(y) {
  levs <- levels(y)
  nl <- setdiff(levs, cfg$pos_label)
  if (length(nl)) nl[1] else cfg$neg_label
}

best_threshold_mcc <- function(obs, p_pos, pos = cfg$pos_label, neg = cfg$neg_label,
                               grid = seq(0.01, 0.99, by = 0.001)) {
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

compute_metrics_binary <- function(obs, pred, p_pos,
                                   pos_label = cfg$pos_label, neg_label = cfg$neg_label) {
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
  acc   <- if (sum(tab)==0) NA_real_ else (TP+TN)/sum(tab)
  f1    <- if (is.na(preci) || is.na(sens) || (preci+sens)==0) NA_real_ else 2*preci*sens/(preci+sens)
  mcc   <- mcc_from_counts(TP, FP, FN, TN)
  aucv  <- binary_auc(y, as_num(p_pos), pos_label)
  c(MCC = mcc, AUC = aucv, `F1` = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec)
}

# ----- PR curve + AP -----
pr_curve_df <- function(obs, p, pos_label = cfg$pos_label) {
  y <- as.integer(obs == pos_label)
  ord <- order(p, decreasing = TRUE)
  tp <- cumsum(y[ord] == 1)
  fp <- cumsum(y[ord] == 0)
  fn <- sum(y == 1) - tp
  precision <- tp / pmax(tp + fp, 1)
  recall    <- tp / pmax(tp + fn, 1)
  tibble(thresh_rank = seq_along(ord), precision = precision, recall = recall) |>
    distinct(recall, .keep_all = TRUE) |>
    arrange(recall)
}

average_precision_manual <- function(obs, p, pos_label = cfg$pos_label) {
  pc <- pr_curve_df(obs, p, pos_label)
  r <- pc$recall; pr <- pc$precision
  if (!length(r)) return(NA_real_)
  r2 <- c(0, r, 1)
  pr2 <- c(1, pr, tail(pr, 1))
  sum(((r2[-1] - r2[-length(r2)]) * (pr2[-1] + pr2[-length(r2)]) / 2), na.rm = TRUE)
}

# ---------------- safe plot helpers ----------------
has_two <- function(x) length(unique(stats::na.omit(x))) >= 2

safe_smooth <- function(df, mapping, method = "loess") {
  if (nrow(df) > 1 && has_two(df[[rlang::as_name(mapping$x)]]) && has_two(df[[rlang::as_name(mapping$y)]])) {
    geom_smooth(mapping = mapping, data = df, method = method, se = FALSE)
  } else NULL
}

safe_line <- function(df, mapping) {
  if (nrow(df) > 1 && has_two(df[[rlang::as_name(mapping$x)]]) && has_two(df[[rlang::as_name(mapping$y)]])) {
    geom_line(mapping = mapping)
  } else if (nrow(df) >= 1) {
    geom_point(mapping = mapping)
  } else NULL
}

try_print <- function(p, label = "figure") {
  ok <- try(print(p), silent = TRUE)
  if (inherits(ok, "try-error")) message(sprintf("[skip] %s failed: %s", label, as.character(ok)))
  invisible(NULL)
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
  keep <- names(all) %in% algos & vapply(all, \(s) method_available(s$method), logical(1))
  all[keep]
}

# ---------------- CV utilities ----------------
choose_k_for_task <- function(y, k_desired) {
  if (is.null(y) || !length(y) || anyNA(y)) return(k_desired)
  tab <- table(y); max_k <- max(2, min(tab))
  min(as.integer(k_desired), as.integer(max_k))
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
  pred <- factor(ifelse(p_va >= thr$t, pos_label, neg_label),
                 levels = levels(dat_va[[yvar]]))
  mets <- compute_metrics_binary(dat_va[[yvar]], pred, p_va, pos_label = pos_label, neg_label = neg_label)
  tibble(
    Fold = fold_tag,
    Algorithm = algo_name,
    Metric = names(mets),
    Value  = as.numeric(mets)
  )
}

run_cv <- function(df, yvar, features, algos,
                   pos_label = cfg$pos_label, inner_k = cfg$inner_k,
                   k_desired = cfg$cv_k, R = cfg$cv_R) {
  use_cols <- unique(c(features, yvar, "data_split"))
  dat <- df[, intersect(use_cols, names(df)), drop = FALSE]
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

# ---------------- holdout evaluation ----------------
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
    neg_label <- get_neg_level(dtrain[[yvar]])
    thr   <- best_threshold_mcc(dtrain[[yvar]], p_tr, pos = cfg$pos_label, neg = neg_label)
    p_te  <- as.numeric(predict(fit, newdata = dtest, type = "prob")[, cfg$pos_label])
    pred  <- factor(ifelse(p_te >= thr$t, cfg$pos_label, neg_label),
                    levels = levels(dtest[[yvar]]))
    out[[sp]] <- list(
      obs = factor(dtest[[yvar]], levels = levels(dtest[[yvar]])),
      p   = p_te,
      pred = pred,
      threshold = thr$t,
      model = fit,
      holdout_df = dtest
    )
  }
  out
}

# ---------------- helper summaries ----------------
agg_to_ci_table <- function(folds_tbl) {
  if (is.null(folds_tbl) || !nrow(folds_tbl)) return(tibble())
  folds_tbl %>%
    group_by(Algorithm, Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      CI_low = quantile(Value, probs = 0.025, na.rm = TRUE),
      CI_high= quantile(Value, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      CI_low  = pmax(ifelse(Metric %in% c("AUC","F1","Accuracy","Precision","Sensitivity","Specificity"), 0, -Inf), CI_low),
      CI_high = pmin(ifelse(Metric %in% c("AUC","F1","Accuracy","Precision","Sensitivity","Specificity"), 1, Inf), CI_high)
    )
}

auc_delong_ci <- function(obs, p_pos, conf.level = 0.95) {
  y <- as.integer(obs == cfg$pos_label)
  if (length(unique(y)) < 2) return(c(NA_real_, NA_real_, NA_real_))
  suppressMessages({
    r <- pROC::roc(y, p_pos, quiet = TRUE)
    as.numeric(pROC::ci.auc(r, conf.level = conf.level))
  })
}

# ---- subgroup performance ----
subgroup_metrics <- function(df_eval, obs, p, pred,
                             vars = c("group","Gender","Diagnosis","severity_admission"),
                             pos_label = cfg$pos_label,
                             min_n = 8, conf.level = 0.95, n_boot = cfg$n_boot_subgrp) {
  if (!length(vars)) return(tibble())
  dat <- tibble(obs = obs, pred = pred, p = p) %>%
    bind_cols(df_eval[intersect(vars, names(df_eval))])
  each_var <- function(v) {
    if (!(v %in% names(dat))) return(tibble())
    levs <- levels(factor(dat[[v]]))
    map_dfr(levs, function(lv) {
      dd <- dat[dat[[v]] == lv & !is.na(dat[[v]]), , drop = FALSE]
      n  <- nrow(dd)
      if (n < min_n || length(unique(dd$obs)) < 2) {
        return(tibble(Subgroup = v, Level = as.character(lv),
                      Metric = c("AUC","MCC","F1","Accuracy","Precision","Sensitivity","Specificity"),
                      Point = NA_real_, CI_low = NA_real_, CI_high = NA_real_, N = n))
      }
      m <- compute_metrics_binary(dd$obs, dd$pred, dd$p, pos_label = pos_label)
      boot_once <- function(idx) {
        mm <- compute_metrics_binary(dd$obs[idx], dd$pred[idx], dd$p[idx], pos_label = pos_label)
        unname(mm[c("AUC","MCC","F1","Accuracy","Precision","Sensitivity","Specificity")])
      }
      nB <- min(n_boot, max(300, n*5))
      set.seed(2L)
      idxs <- replicate(nB, sample.int(n, replace=TRUE), simplify = FALSE)
      B <- map_dfr(idxs, ~{
        vals <- boot_once(.x)
        tibble(AUC=vals[1], MCC=vals[2], `F1`=vals[3], Accuracy=vals[4], Precision=vals[5], Sensitivity=vals[6], Specificity=vals[7])
      })
      mkci <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(c(NA,NA)); quantile(x, c(0.025,0.975), na.rm=TRUE, names=FALSE) }
      tibble(
        Subgroup = v, Level = as.character(lv), N = n,
        Metric  = names(m),
        Point   = as.numeric(m),
        CI_low  = c(mkci(B$AUC)[1], mkci(B$MCC)[1], mkci(B$`F1`)[1], mkci(B$Accuracy)[1], mkci(B$Precision)[1], mkci(B$Sensitivity)[1], mkci(B$Specificity)[1]),
        CI_high = c(mkci(B$AUC)[2], mkci(B$MCC)[2], mkci(B$`F1`)[2], mkci(B$Accuracy)[2], mkci(B$Precision)[2], mkci(B$Sensitivity)[2], mkci(B$Specificity)[2])
      )
    })
  }
  map_dfr(vars, each_var)
}

# ---------------- decision curve ----------------
decision_curve_table <- function(obs, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  y <- as.integer(obs == cfg$pos_label)
  N <- length(y); prev <- mean(y)
  map_dfr(thresholds, function(pt) {
    pred <- as.integer(p >= pt)
    TP <- sum(pred == 1 & y == 1)
    FP <- sum(pred == 1 & y == 0)
    NB <- TP/N - FP/N * (pt/(1-pt))
    NB_all <- prev - (1 - prev) * (pt/(1-pt))
    Precision <- if ((TP + sum(pred == 1 & y == 0)) == 0) NA_real_ else TP / sum(pred == 1)
    tibble(
      threshold = pt,
      Net_Benefit = NB,
      NB_TreatAll = NB_all,
      NB_TreatNone = 0,
      Precision = Precision,
      NNE = ifelse(is.finite(1/(Precision - pt/(1-pt))), 1/(Precision - pt/(1-pt)), NA_real_)
    )
  })
}

# ---------------- pipeline runner ----------------
run_pipeline <- function() {
  df <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))
  if (!("data_split" %in% names(df)) || all(is.na(df$data_split))) {
    set.seed(cfg$seed_cv)
    n <- nrow(df); idx <- sample(n); n_tr <- floor(0.7*n); n_te <- floor(0.15*n)
    df$data_split <- factor(c(rep("train", n_tr), rep("test", n_te), rep("external", n - n_tr - n_te))[order(idx)])
  }
  cat("\n[Sanity] data_split counts:\n\n"); print(table(df$data_split, useNA = "ifany"))
  yvar <- cfg$outcome
  algos <- c("C4.5","k-NN","SVM","RF","LR")
  
  full_folds <- run_cv(df, yvar, features = feat_full,  algos = algos)
  tri_folds  <- run_cv(df, yvar, features = feat_triage, algos = algos)
  cv_full_ci <- agg_to_ci_table(full_folds) |> mutate(Split = "Train-CV", Feature_Set = "Full")
  cv_tri_ci  <- agg_to_ci_table(tri_folds)  |> mutate(Split = "Train-CV", Feature_Set = "Triage")
  
  eval_all_holdouts_for_set <- function(feature_set, set_name, algos) {
    holder <- lapply(algos, function(a) run_holdout(df, yvar, features = feature_set, algo_name = a))
    names(holder) <- algos
    list(set_name = set_name, objects = holder)
  }
  hold_full <- eval_all_holdouts_for_set(feat_full,  "Full",   algos)
  hold_tri  <- eval_all_holdouts_for_set(feat_triage,"Triage", algos)
  
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
    sub_collect(hold_tri, "test"),
    sub_collect(hold_tri, "external")
  )
  
  list(df=df, yvar=yvar, algos=algos,
       cv_full_ci=cv_full_ci, cv_tri_ci=cv_tri_ci,
       hold_full=hold_full, hold_tri=hold_tri,
       tbl7_subgroups=tbl7_subgroups)
}

# ---------------- theme ----------------
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(color = "grey30"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

# ============================================================
# Fig1, workflow
# ============================================================
fig1_workflow <- function() {
  steps <- tibble(
    x = 1:10,
    y = 1,
    label = c(
      "Data ingestion\n(raw tables, QA)",
      "Preprocessing\n(normalize, QC)",
      "Cohort definition\n(inclusion, exclusion)",
      "Split\n(train, test, external)",
      "Modeling at triage\n(feature set)",
      "Nested CV\n(tuning, selection)",
      "Threshold selection\n(MCC target)",
      "Evaluation\n(ROC, PR, calibration, DCA)",
      "Recalibration\n(logistic, isotonic)",
      "Reporting\n(tables, figures, scorecard)"
    )
  )
  box_w <- 0.45; box_h <- 0.35
  ggplot(steps) +
    geom_rect(aes(xmin = x - box_w, xmax = x + box_w, ymin = y - box_h, ymax = y + box_h), fill = "white", color = "grey40") +
    geom_text(aes(x = x, y = y, label = label), size = 3) +
    geom_segment(data = steps |> filter(x < max(x)), aes(x = x + box_w, xend = x + 1 - box_w, y = y, yend = y),
                 arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
    coord_cartesian(xlim = c(0.5, 10.5), ylim = c(0.4, 1.6)) +
    labs(title = "Workflow") +
    theme_void() + theme_pub()
}

# ============================================================
# Fig S1, cohort flow
# ============================================================
figS1_cohort_flow <- function(df) {
  incl <- nrow(df)
  tr <- sum(df$data_split == "train", na.rm = TRUE)
  te <- sum(df$data_split == "test", na.rm = TRUE)
  ex <- sum(df$data_split == "external", na.rm = TRUE)
  counts <- list(screened = NA_integer_, excluded = NA_integer_, included = incl, train = tr, test = te, external = ex)
  nodes <- tibble(
    id = 1:5,
    label = c(
      sprintf("Screened (n=%s)", counts$screened %||% "?"),
      sprintf("Excluded (n=%s)", counts$excluded %||% "?"),
      sprintf("Included (n=%s)", counts$included %||% 0),
      sprintf("Train (n=%s)", counts$train %||% 0),
      sprintf("Test, External (n=%s, %s)", counts$test %||% 0, counts$external %||% 0)
    ),
    x = 1, y = c(5,4,3,2,1)
  )
  box_w <- 0.6; box_h <- 0.35
  segs <- tibble(x = 1, xend = 1, y = c(4.65,3.65,2.65,1.65), yend = c(4.35,3.35,2.35,1.35))
  ggplot(nodes) +
    geom_rect(aes(xmin = x - box_w, xmax = x + box_w, ymin = y - box_h, ymax = y + box_h), fill = "white", color = "grey40") +
    geom_text(aes(x = x, y = y, label = label), size = 3) +
    geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
    coord_cartesian(xlim = c(0.2, 1.8), ylim = c(0.4, 5.6)) +
    labs(title = "Cohort flow") +
    theme_void() + theme_pub()
}

# ============================================================
# Fig2, heatmap of CV performance, all metrics
# ============================================================
fig2_heatmap <- function(cv_full_ci, cv_tri_ci) {
  stopifnot(nrow(cv_full_ci) + nrow(cv_tri_ci) > 0)
  metrics_keep <- c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")
  recode_map <- c("MCC"="MCC","AUC-ROC"="AUC","F1-score"="F1",
                  "Accuracy"="Accuracy","Precision"="Precision",
                  "Sensitivity"="Sensitivity","Specificity"="Specificity")
  
  perf_tbl <- bind_rows(
    cv_full_ci %>% mutate(Feature_Set = "Full"),
    cv_tri_ci  %>% mutate(Feature_Set = "Triage")
  ) %>%
    filter(Metric %in% metrics_keep) %>%
    rename(Value = Mean) %>%
    mutate(
      Metric = factor(Metric, levels = metrics_keep),
      MetricLabel = forcats::fct_recode(Metric, !!!as.list(recode_map))
    )
  
  ggplot(perf_tbl, aes(Algorithm, MetricLabel, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.finite(Value), scales::percent(Value, accuracy = 0.1), "NA")), size = 3) +
    scale_fill_viridis_c(option = "C", end = 0.95, name = "Score") +
    facet_wrap(~Feature_Set, nrow = 1) +
    labs(title = "Cross-validated performance", x = "Algorithm", y = "Metric") +
    theme_pub()
}

# ============================================================
# Fig3, permutation feature importance, ΔMCC
# ============================================================
perm_drop_mcc <- function(fit, data, yvar, features, threshold, pos_label = cfg$pos_label) {
  p0 <- as.numeric(predict(fit, newdata = data, type = "prob")[, pos_label])
  neg_label <- get_neg_level(data[[yvar]])
  pred0 <- factor(ifelse(p0 >= threshold, pos_label, neg_label),
                  levels = levels(data[[yvar]]))
  m0 <- as.numeric(compute_metrics_binary(data[[yvar]], pred0, p0, pos_label = pos_label, neg_label = neg_label)["MCC"])
  map_dfr(features, function(f) {
    dperm <- data
    dperm[[f]] <- sample(dperm[[f]])
    pp <- as.numeric(predict(fit, newdata = dperm, type = "prob")[, pos_label])
    pr <- factor(ifelse(pp >= threshold, pos_label, neg_label),
                 levels = levels(data[[yvar]]))
    m1 <- as.numeric(compute_metrics_binary(data[[yvar]], pr, pp, pos_label = pos_label, neg_label = neg_label)["MCC"])
    tibble(Feature = f, Drop_MCC = m0 - m1)
  }) %>% arrange(desc(Drop_MCC))
}

fig3_feature_importance <- function(hold_full, hold_tri, cfg, feat_full, feat_triage, primary_algo = "LR", split_name = "test") {
  pick <- function(holder) {
    algos <- intersect(primary_algo, names(holder$objects))
    if (length(algos) == 0) algos <- names(holder$objects)[1]
    list(algo = algos[1], obj = holder$objects[[algos[1]]][[split_name]])
  }
  F <- pick(hold_full); T <- pick(hold_tri)
  stopifnot(!is.null(F$obj), !is.null(T$obj))
  imp_full  <- perm_drop_mcc(F$obj$model, F$obj$holdout_df, cfg$outcome, intersect(feat_full, names(F$obj$holdout_df)), F$obj$threshold)
  imp_tri   <- perm_drop_mcc(T$obj$model, T$obj$holdout_df, cfg$outcome, intersect(feat_triage, names(T$obj$holdout_df)), T$obj$threshold)
  imp_full$Set <- "Full"; imp_tri$Set <- "Triage"
  imp_all <- bind_rows(imp_full, imp_tri) %>%
    filter(is.finite(Drop_MCC)) %>%
    group_by(Set) %>% slice_max(order_by = Drop_MCC, n = 10, with_ties = FALSE) %>% ungroup()
  ggplot(imp_all, aes(Drop_MCC, forcats::fct_reorder(Feature, Drop_MCC))) +
    geom_col() + facet_wrap(~Set, scales = "free_y") +
    labs(title = sprintf("Permutation importance, algo = %s, split = %s", F$algo, tools::toTitleCase(split_name)),
         x = "ΔMCC when permuted", y = "Feature") + theme_pub()
}

# ============================================================
# Fig4, decision curves
# ============================================================
fig4_dca <- function(holder, primary_algo = "LR") {
  pick <- function(holder, split) {
    alg <- if (primary_algo %in% names(holder$objects)) primary_algo else names(holder$objects)[1]
    holder$objects[[alg]][[split]]
  }
  mk_panel <- function(o, ttl) {
    if (is.null(o)) return(ggplot() + theme_void() + ggtitle(ttl))
    dca <- decision_curve_table(o$obs, o$p)
    ggplot(dca, aes(threshold, Net_Benefit)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
      geom_line(linewidth = 0.7) +
      geom_line(aes(y = NB_TreatAll), linetype = 3) +
      geom_line(aes(y = NB_TreatNone), linetype = 3) +
      labs(title = ttl, x = "Threshold probability", y = "Net benefit") + theme_pub()
  }
  p_test <- mk_panel(pick(holder, "test"),  "Decision curve, test")
  p_ext  <- mk_panel(pick(holder, "external"), "Decision curve, external")
  p_test | p_ext
}

# ============================================================
# Fig5, ROC and PR with AP
# ============================================================
fig5_roc_pr <- function(holder, primary_algo = "LR") {
  pick <- function(holder, split) {
    alg <- if (primary_algo %in% names(holder$objects)) primary_algo else names(holder$objects)[1]
    holder$objects[[alg]][[split]]
  }
  mk_roc <- function(o, ttl) {
    if (is.null(o) || length(unique(o$obs)) < 2) return(ggplot() + theme_void() + ggtitle(ttl))
    r <- pROC::roc(o$obs, o$p, levels = rev(levels(o$obs)), quiet = TRUE)
    ci <- pROC::ci.auc(r)
    df <- tibble(x = 1 - rev(r$specificities), y = rev(r$sensitivities))
    ggplot(df, aes(x, y)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      safe_line(df, aes(x, y)) +
      coord_equal() + theme_pub() + labs(title = ttl, x = "1 − Specificity", y = "Sensitivity") +
      annotate("text", x = 0.6, y = 0.2,
               label = sprintf("AUC %.3f [%.3f, %.3f]", as.numeric(pROC::auc(r)), ci[1], ci[3]), hjust = 0)
  }
  mk_pr <- function(o, ttl) {
    if (is.null(o) || length(unique(o$obs)) < 2) return(ggplot() + theme_void() + ggtitle(ttl))
    pr <- pr_curve_df(o$obs, o$p, pos_label = cfg$pos_label)
    ap <- average_precision_manual(o$obs, o$p, pos_label = cfg$pos_label)
    ggplot(pr, aes(recall, precision)) +
      safe_line(pr, aes(recall, precision)) +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_pub() + labs(title = ttl, x = "Recall", y = "Precision") +
      annotate("text", x = 0.6, y = 0.2, label = sprintf("AP %.3f", ap %||% NA_real_), hjust = 0)
  }
  o_te <- pick(holder, "test"); o_ex <- pick(holder, "external")
  (mk_roc(o_te, "ROC, test") | mk_pr(o_te, "PR, test")) /
    (mk_roc(o_ex, "ROC, external") | mk_pr(o_ex, "PR, external"))
}

# ============================================================
# Fig6, calibration before and after, robust Platt fallback
# ============================================================
calibration_curve <- function(obs, prob, bins = 10) {
  df <- tibble(obs = as.integer(obs == levels(factor(obs))[2]), prob = prob) %>%
    mutate(bin = ntile(prob, bins)) %>%
    group_by(bin) %>% summarize(x = mean(prob, na.rm = TRUE), y = mean(obs, na.rm = TRUE), n = n(), .groups = "drop")
  df
}

platt_fit <- function(obs, prob) {
  eps <- 1e-6
  p  <- pmin(pmax(prob, eps), 1 - eps)
  y  <- as.integer(obs == levels(factor(obs))[2])
  
  fit_glm <- try(glm(y ~ qlogis(p), family = binomial(), control = list(maxit = 50)), silent = TRUE)
  bad <- inherits(fit_glm, "try-error") || any(!is.finite(coef(fit_glm)))
  if (!bad) {
    return(list(
      type = "platt",
      intercept = unname(coef(fit_glm)[1]),
      slope     = unname(coef(fit_glm)[2]),
      predict   = function(px) plogis(coef(fit_glm)[1] + coef(fit_glm)[2] * qlogis(pmin(pmax(px, eps), 1 - eps)))
    ))
  }
  iso <- stats::isoreg(p, y)
  xs <- iso$x; ys <- iso$yf
  if (length(unique(xs)) < 2) {
    mu <- mean(y)
    pred_const <- function(px) rep(mu, length(px))
    return(list(type = "isotonic", intercept = NA_real_, slope = NA_real_, predict = pred_const))
  }
  pred_iso <- function(px) {
    stats::approx(x = xs, y = ys, xout = px, method = "linear", ties = "ordered")$y %>%
      pmin(1 - eps) %>% pmax(eps)
  }
  list(type = "isotonic", intercept = NA_real_, slope = NA_real_, predict = pred_iso)
}

calib_panels <- function(o, ttl_before, ttl_after, bins = 10) {
  if (is.null(o)) return(ggplot() + theme_void() + ggtitle("No data"))
  obs <- o$obs; p0 <- o$p
  cal0 <- calibration_curve(obs, p0, bins)
  fit  <- platt_fit(obs, p0)
  p1   <- fit$predict(p0)
  cal1 <- calibration_curve(obs, p1, bins)
  y_num <- as.integer(obs == levels(factor(obs))[2])
  
  p_left <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey70") +
    geom_point(data = cal0, aes(x, y), alpha = 0.8) +
    safe_smooth(tibble(x = p0, y = y_num), aes(x, y)) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) + theme_pub() +
    labs(title = ttl_before, x = "Predicted risk", y = "Observed")
  
  p_mid <- ggplot(tibble(p = p0), aes(p)) +
    { if (sum(is.finite(p0)) >= 1) geom_histogram(bins = 30) else geom_blank() } +
    theme_pub() + labs(title = "Risk distribution", x = "Predicted risk", y = "Count")
  
  coef_lbl <- if (identical(fit$type, "platt"))
    sprintf("Intercept %.3f\nSlope %.3f", fit$intercept, fit$slope) else
      "Isotonic fit"
  
  p_right <- ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, label = coef_lbl) +
    theme_void() + theme_pub() + labs(title = if (identical(fit$type, "platt")) "Platt coefficients" else "Calibration")
  
  row1 <- p_left + p_mid + p_right + patchwork::plot_layout(widths = c(1.2, 1, 0.7))
  
  p_left2 <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey70") +
    geom_point(data = cal1, aes(x, y), alpha = 0.8, shape = 17) +
    safe_smooth(tibble(x = p1, y = y_num), aes(x, y)) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) + theme_pub() +
    labs(title = ttl_after, x = "Predicted risk", y = "Observed")
  
  p_mid2 <- ggplot(tibble(p = p1), aes(p)) +
    { if (sum(is.finite(p1)) >= 1) geom_histogram(bins = 30) else geom_blank() } +
    theme_pub() + labs(title = "Risk distribution", x = "Recalibrated risk", y = "Count")
  
  p_right2 <- ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
             label = if (identical(fit$type, "platt")) "After recalibration" else "Isotonic") +
    theme_void() + theme_pub() + labs(title = "After recalibration")
  
  row1 / (p_left2 + p_mid2 + p_right2 + patchwork::plot_layout(widths = c(1.2, 1, 0.7)))
}

fig6_calibration_triptych <- function(holder, primary_algo = "LR") {
  alg <- if (primary_algo %in% names(holder$objects)) primary_algo else names(holder$objects)[1]
  o_te <- holder$objects[[alg]][["test"]]
  o_ex <- holder$objects[[alg]][["external"]]
  calib_panels(o_te, "Calibration, test, before", "Calibration, test, after") /
    calib_panels(o_ex, "Calibration, external, before", "Calibration, external, after")
}

# ============================================================
# Fig7, parsimony curve on triage set
# ============================================================
fig7_parsimony <- function(df) {
  set.seed(123)
  dtr <- df[df$data_split == "train", , drop = FALSE]
  rec <- recipes::recipe(stats::as.formula(paste(cfg$outcome, "~", paste(feat_triage, collapse = "+"))), data = dtr) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
    recipes::step_novel(recipes::all_nominal_predictors()) |>
    recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors())
  tr_ctrl <- caret::trainControl(method = "cv", number = cfg$inner_k, classProbs = TRUE, summaryFunction = twoClassSummary)
  fit0 <- try(caret::train(rec, data = dtr, method = "glmnet", metric = "ROC", trControl = tr_ctrl, tuneLength = 10), silent = TRUE)
  order_feats <- if (!inherits(fit0, "try-error")) {
    vi <- caret::varImp(fit0)$importance %>% tibble::rownames_to_column("Feature") %>% arrange(desc(Overall))
    intersect(feat_triage, vi$Feature)
  } else feat_triage
  if (!length(order_feats)) order_feats <- feat_triage
  
  ctrl <- caret::trainControl(method = "cv", number = cfg$inner_k, classProbs = TRUE, summaryFunction = twoClassSummary)
  res <- map_dfr(seq_along(order_feats), function(k) {
    feats <- order_feats[seq_len(k)]
    rec_k <- recipes::recipe(stats::as.formula(paste(cfg$outcome, "~", paste(feats, collapse = "+"))), data = dtr) |>
      recipes::step_impute_median(recipes::all_numeric_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
      recipes::step_novel(recipes::all_nominal_predictors()) |>
      recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01) |>
      recipes::step_dummy(recipes::all_nominal_predictors()) |>
      recipes::step_zv(recipes::all_predictors()) |>
      recipes::step_normalize(recipes::all_numeric_predictors())
    
    prepped <- recipes::prep(rec_k)
    X <- recipes::bake(prepped, new_data = dtr) %>% select(-all_of(cfg$outcome))
    method_k <- if (ncol(X) < 2) "glm" else "glmnet"
    
    fit_k <- caret::train(rec_k, data = dtr, method = method_k, metric = "ROC", trControl = ctrl, tuneLength = 10)
    tibble(k = k, Score = max(fit_k$results$ROC, na.rm = TRUE))
  })
  if (!nrow(res)) return(ggplot() + theme_void() + ggtitle("Parsimony curve, no CV results"))
  ggplot(res, aes(k, Score)) +
    safe_line(res, aes(k, Score)) +
    geom_point() +
    { if (nrow(res) >= 1) geom_vline(xintercept = which.max(res$Score), linetype = 2, color = "grey50") else NULL } +
    scale_x_continuous(breaks = seq_len(max(res$k))) +
    labs(title = "Parsimony curve, triage set", x = "Number of features", y = "CV ROC") + theme_pub()
}

# ============================================================
# Fig8, subgroup forest with CIs, metric selectable
# ============================================================
fig8_subgroup_forest <- function(tbl7_subgroups, feature_set = "Triage", primary_algo = "LR", metric = "AUC") {
  d <- tbl7_subgroups %>%
    filter(Feature_Set == feature_set, Algorithm %in% c(primary_algo, unique(Algorithm)[1]), Metric == metric)
  d <- d %>% filter(is.finite(Point))
  if (!nrow(d)) return(ggplot() + theme_void() + ggtitle("No subgroup rows found for plotting"))
  d$Facet <- d$Subgroup
  d$Label <- paste0(d$Level, " (n=", d$N, ")")
  ggplot(d, aes(Point, forcats::fct_reorder(Label, Point, .na_rm = TRUE))) +
    { if (any(is.finite(d$CI_low) & is.finite(d$CI_high)))
      geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0) else NULL } +
    geom_point(size = 2) +
    geom_vline(xintercept = mean(d$Point, na.rm = TRUE), linetype = 2, color = "grey60") +
    facet_wrap(~Facet, scales = "free_y") +
    labs(title = sprintf("Subgroup performance, %s, algo=%s", feature_set, primary_algo),
         x = metric, y = NULL, subtitle = "95% bootstrap CIs") +
    theme_pub()
}

# ============================================================
# Fig9 — Scorecard, extra-robust: Firth → Firth-noCI → Ridge
# ============================================================
fig9_scorecard <- function(df) {
  dtr <- df[df$data_split == "train", , drop = FALSE]
  form <- stats::as.formula(paste(cfg$outcome, "~", paste(feat_triage, collapse = "+")))
  rec <- recipes::recipe(form, data = dtr) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
    recipes::step_novel(recipes::all_nominal_predictors()) |>
    recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors())
  
  prepped <- recipes::prep(rec)
  X <- recipes::bake(prepped, new_data = dtr)
  y <- dtr[[cfg$outcome]]
  if (!ncol(X)) return(ggplot() + theme_void() + ggtitle("No predictors available for scorecard"))
  
  coefs <- NULL
  used <- NULL
  
  # ---- Try Firth with CIs
  if (requireNamespace("logistf", quietly = TRUE)) {
    firth_fit <- try(logistf::logistf(as.formula(paste(cfg$outcome, "~ .")), data = X), silent = TRUE)
    if (!inherits(firth_fit, "try-error")) {
      # First try with CIs, if LAPACK fails, catch and retry without CIs
      tidy_ci <- try(broom::tidy(firth_fit, conf.int = TRUE, exponentiate = TRUE), silent = TRUE)
      if (!inherits(tidy_ci, "try-error")) {
        coefs <- dplyr::filter(tidy_ci, term != "(Intercept)")
        used <- "Firth logistic"
      } else {
        tidy_noci <- broom::tidy(firth_fit, conf.int = FALSE, exponentiate = TRUE)
        coefs <- dplyr::filter(tidy_noci, term != "(Intercept)") |>
          dplyr::mutate(conf.low = NA_real_, conf.high = NA_real_)
        used <- "Firth logistic, CIs unavailable"
      }
    }
  }
  
  # ---- Ridge fallback if Firth not available or failed
  if (is.null(coefs)) {
    x <- as.matrix(X[setdiff(names(X), cfg$outcome)])
    y_num <- as.integer(y == cfg$pos_label)
    cv <- glmnet::cv.glmnet(x, y_num, family = "binomial", alpha = 0, nfolds = 5)
    b  <- as.numeric(coef(cv, s = "lambda.min"))[-1]
    coefs <- tibble::tibble(term = colnames(x), estimate = exp(b),
                            conf.low = NA_real_, conf.high = NA_real_)
    used <- "Ridge logistic, CIs not shown"
  }
  
  coefs <- coefs %>%
    dplyr::mutate(ok = is.finite(estimate),
                  ok_ci = is.finite(conf.low) & is.finite(conf.high)) %>%
    dplyr::filter(ok)
  
  if (nrow(coefs) < 2) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
                 label = "Scorecard could not be drawn, insufficient stable coefficients") +
        theme_void() + theme_pub() + labs(title = "Scorecard")
    )
  }
  
  p <- ggplot(coefs, aes(estimate, forcats::fct_reorder(term, estimate, .na_rm = TRUE))) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
    { if (any(coefs$ok_ci, na.rm = TRUE))
      geom_errorbarh(aes(xmin = pmax(conf.low, .Machine$double.eps),
                         xmax = pmax(conf.high, .Machine$double.eps)), height = 0) else NULL } +
    geom_point() +
    labs(title = paste0("Scorecard, ", used), x = "Odds ratio", y = NULL) + theme_pub()
  
  if (all(is.finite(coefs$estimate)) && all(coefs$estimate > 0, na.rm = TRUE)) p <- p + scale_x_log10()
  p
}


# ============================================================
# Driver
# ============================================================
print_all_figures <- function(primary_algo = "LR") {
  message("[setup] Running pipeline to create CV summaries and holdouts...")
  pipe <- run_pipeline()
  
  message("[figs] printing Fig1–Fig9, Fig S1")
  try_print(fig1_workflow(), "Fig1 workflow")
  try_print(figS1_cohort_flow(pipe$df), "Fig S1 cohort flow")
  try_print(fig2_heatmap(pipe$cv_full_ci, pipe$cv_tri_ci), "Fig2 heatmap")
  try_print(fig3_feature_importance(pipe$hold_full, pipe$hold_tri, cfg, feat_full, feat_triage, primary_algo = primary_algo, split_name = "test"),
            "Fig3 permutation importance")
  try_print(fig4_dca(pipe$hold_tri, primary_algo = primary_algo), "Fig4 DCA")
  try_print(fig5_roc_pr(pipe$hold_tri, primary_algo = primary_algo), "Fig5 ROC and PR")
  try_print(fig6_calibration_triptych(pipe$hold_tri, primary_algo = primary_algo), "Fig6 calibration")
  try_print(fig7_parsimony(pipe$df), "Fig7 parsimony")
  try_print(fig8_subgroup_forest(pipe$tbl7_subgroups, feature_set = "Triage", primary_algo = primary_algo, metric = "AUC"),
            "Fig8 subgroup forest")
  try_print(fig9_scorecard(pipe$df), "Fig9 scorecard")
  invisible(pipe)
}

# --- RUN NOW ---
print_all_figures(primary_algo = "LR")









# ============================================================
# Save all figures — PNG only
# ============================================================

# Build all plots without printing, guard each build
make_all_plots <- function(pipe, primary_algo = "LR") {
  build <- function(expr) {
    obj <- try(expr, silent = TRUE)
    if (inherits(obj, "try-error")) NULL else obj
  }
  plots <- list(
    Fig1_Workflow        = build(fig1_workflow()),
    FigS1_CohortFlow     = build(figS1_cohort_flow(pipe$df)),
    Fig2_CVHeatmap       = build(fig2_heatmap(pipe$cv_full_ci, pipe$cv_tri_ci)),
    Fig3_PermutationMCC  = build(fig3_feature_importance(pipe$hold_full, pipe$hold_tri, cfg, feat_full, feat_triage,
                                                         primary_algo = primary_algo, split_name = "test")),
    Fig4_DecisionCurves  = build(fig4_dca(pipe$hold_tri, primary_algo = primary_algo)),
    Fig5_ROC_PR          = build(fig5_roc_pr(pipe$hold_tri, primary_algo = primary_algo)),
    Fig6_Calibration     = build(fig6_calibration_triptych(pipe$hold_tri, primary_algo = primary_algo)),
    Fig7_Parsimony       = build(fig7_parsimony(pipe$df)),
    Fig8_Subgroup_AUC    = build(fig8_subgroup_forest(pipe$tbl7_subgroups, feature_set = "Triage",
                                                      primary_algo = primary_algo, metric = "AUC")),
    Fig9_Scorecard       = build(fig9_scorecard(pipe$df))
  )
  keep <- vapply(plots, inherits, logical(1), what = "ggplot")
  plots[keep]
}

# Safe PNG saver
safe_png <- function(plot, filename, width, height, dpi = 300, bg = "white") {
  try({
    ggplot2::ggsave(filename = filename, plot = plot,
                    width = width, height = height, units = "in",
                    dpi = dpi, bg = bg, limitsize = FALSE)
  }, silent = TRUE)
}

# Per-figure sizes, fallback used if not found
.default_size <- c(w = 9, h = 6)
.figure_sizes <- list(
  Fig1_Workflow       = c(w = 11, h = 3.5),
  FigS1_CohortFlow    = c(w = 6,  h = 6.5),
  Fig2_CVHeatmap      = c(w = 10, h = 4.5),
  Fig3_PermutationMCC = c(w = 9,  h = 6),
  Fig4_DecisionCurves = c(w = 10, h = 4.5),
  Fig5_ROC_PR         = c(w = 10, h = 8),
  Fig6_Calibration    = c(w = 12, h = 8),
  Fig7_Parsimony      = c(w = 7,  h = 4.5),
  Fig8_Subgroup_AUC   = c(w = 10, h = 7),
  Fig9_Scorecard      = c(w = 7,  h = 6)
)

get_size <- function(name) {
  if (!is.null(.figure_sizes[[name]])) .figure_sizes[[name]] else .default_size
}

# Save all as PNG
save_all_figures_png <- function(primary_algo = "LR",
                                 out_dir = "figs_export",
                                 prefix = format(Sys.time(), "%Y%m%d"),
                                 dpi = 300, bg = "white") {
  pipe  <- run_pipeline()
  plots <- make_all_plots(pipe, primary_algo = primary_algo)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  mk <- function(nm, ext = "png") {
    base <- if (nzchar(prefix)) paste0(prefix, "_", nm) else nm
    file.path(out_dir, paste0(base, ".", ext))
  }
  
  for (nm in names(plots)) {
    p  <- plots[[nm]]
    sz <- get_size(nm)
    safe_png(p, mk(nm, "png"), width = sz["w"], height = sz["h"], dpi = dpi, bg = bg)
  }
  
  invisible(list(dir = out_dir, pngs = file.path(out_dir, paste0(prefix, "_", names(plots), ".png"))))
}

# =========================
# Example usage
# =========================
# Creates PNG files under figs_export
save_all_figures_png(primary_algo = "LR",
                     out_dir = "figs_export",
                     prefix = format(Sys.time(), "%Y%m%d"),
                     dpi = 300)


