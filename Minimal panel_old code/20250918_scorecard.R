# ============================================================
# Subgroup performance + Scorecard — complete corrected script
# - Reads Excel, filters to {CTRL_noCOVID, COVID}
# - Builds triage LR with CV, evaluates holdouts
# - Subgroup AUC forest with pastel palette, severity restricted to COVID, CIs clamped
# - Scorecard uses Firth logistic with explicit event coding, ridge fallback
# - Prints two figures and two tables, saves nothing
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(patchwork)
  library(forcats)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
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
  n_boot_subgrp = 300
)
if (cfg$fast_mode) {
  cfg$cv_k <- 3; cfg$inner_k <- 3; cfg$tune_len <- 3; cfg$n_boot_subgrp <- 300
}

feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- data loader ----------------
load_data <- function(path, sheet, outcome, keep) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet = sheet) |> as.data.frame()
  
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df) == grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  df$.rid <- seq_len(nrow(df))
  
  for (v in intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' not found", outcome))
  to_YN <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) return(factor(ifelse(x, cfg$pos_label, cfg$neg_label),
                                     levels = c(cfg$pos_label, cfg$neg_label)))
    if (is.numeric(x)) {
      stopifnot(all(x %in% c(0,1), na.rm = TRUE))
      return(factor(ifelse(x==1,cfg$pos_label,cfg$neg_label),
                    levels = c(cfg$pos_label,cfg$neg_label)))
    }
    if (is.character(x)) {
      s <- trimws(tolower(x))
      map_yes <- c("1","yes","y","true","pos","positive")
      map_no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% map_yes, cfg$pos_label,
                    ifelse(s %in% map_no, cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome cannot be coerced to Yes/No")
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    stop("Unsupported outcome type")
  }
  df[[outcome]] <- to_YN(df[[outcome]])
  
  # severity defined only for COVID
  if ("severity_admission" %in% names(df)) {
    keep_lv <- c("Mild","Moderate","Severe")
    df$severity_admission <- fct_drop(factor(as.character(df$severity_admission)))
    df$severity_admission <- ifelse(df$Diagnosis == "COVID" &
                                      df$severity_admission %in% keep_lv,
                                    as.character(df$severity_admission), NA_character_)
    df$severity_admission <- factor(df$severity_admission, levels = intersect(keep_lv, unique(df$severity_admission)))
  }
  
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
  nl <- setdiff(levels(y), cfg$pos_label)
  if (length(nl)) nl[1] else cfg$neg_label
}
best_threshold_mcc <- function(obs, p_pos, pos = cfg$pos_label, neg = cfg$neg_label,
                               grid = seq(0.01, 0.99, by = 0.001)) {
  y <- factor(obs, levels = c(neg, pos))
  best_t <- 0.5; best_m <- -Inf
  for (t in grid) {
    pred <- factor(ifelse(p_pos >= t, pos, neg), levels = c(neg, pos))
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
    for (i in seq_along(idx)) splits[[paste0("r", r, "_f", i)]] <- list(va = idx[[i]])
  }
  splits
}
model_specs <- function(tune_len, algos = c("LR")) {
  all <- list(
    "LR" = list(method = "glmnet", tuneLength = tune_len,
                grid = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 0, length.out = tune_len)))
  )
  avail <- names(caret::getModelInfo())
  all[names(all) %in% algos & vapply(all, function(s) s$method %in% avail, logical(1))]
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
  tibble(Fold = fold_tag, Algorithm = algo_name, Metric = names(mets), Value = as.numeric(mets))
}
run_cv <- function(df, yvar, features, algos = c("LR"),
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
run_holdout <- function(df, yvar, features, algo_name = "LR",
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
    pred  <- factor(ifelse(p_te >= thr$t, cfg$pos_label, neg_label), levels = levels(dtest[[yvar]]))
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

# ---------------- subgroup performance ----------------
subgroup_metrics <- function(df_eval, obs, p, pred,
                             vars = c("group","Gender","Diagnosis","severity_admission"),
                             pos_label = cfg$pos_label,
                             min_n = 8, conf.level = 0.95, n_boot = cfg$n_boot_subgrp) {
  if (!length(vars)) return(tibble())
  dat <- tibble(obs = obs, pred = pred, p = p) %>%
    bind_cols(df_eval[intersect(vars, names(df_eval))])
  
  each_var <- function(v) {
    if (!(v %in% names(dat))) return(tibble())
    dd_all <- if (identical(v, "severity_admission") && "Diagnosis" %in% names(dat)) {
      d <- dat[dat$Diagnosis == "COVID" & !is.na(dat$severity_admission), , drop = FALSE]
      if (!nrow(d)) return(tibble())
      d
    } else dat
    
    levs <- levels(factor(dd_all[[v]]))
    map_dfr(levs, function(lv) {
      dd <- dd_all[dd_all[[v]] == lv & !is.na(dd_all[[v]]), , drop = FALSE]
      n  <- nrow(dd)
      if (n < min_n || length(unique(dd$obs)) < 2) {
        return(tibble(Subgroup = v, Level = as.character(lv),
                      Metric = c("AUC","MCC","F1","Accuracy","Precision","Sensitivity","Specificity"),
                      Point = NA_real_, CI_low = NA_real_, CI_high = NA_real_, N = n))
      }
      m <- compute_metrics_binary(dd$obs, dd$pred, dd$p, pos_label = pos_label)
      
      nB <- min(n_boot, max(300, n*5))
      set.seed(2L)
      idxs <- replicate(nB, sample.int(n, replace = TRUE), simplify = FALSE)
      B <- map_dfr(idxs, ~{
        mm <- compute_metrics_binary(dd$obs[.x], dd$pred[.x], dd$p[.x], pos_label = pos_label)
        tibble(AUC = mm["AUC"], MCC = mm["MCC"], `F1` = mm["F1"],
               Accuracy = mm["Accuracy"], Precision = mm["Precision"],
               Sensitivity = mm["Sensitivity"], Specificity = mm["Specificity"])
      })
      
      mkci <- function(x, clamp01 = FALSE) {
        x <- x[is.finite(x)]
        if (!length(x)) return(c(NA_real_, NA_real_))
        q <- quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
        if (clamp01) q <- pmin(pmax(q, 0), 1)
        q
      }
      
      mets <- names(m)
      rows <- map_dfr(mets, function(mm) {
        cl01 <- identical(mm, "AUC")
        qs <- mkci(B[[mm]], clamp01 = cl01)
        tibble(Metric = mm, Point = as.numeric(m[[mm]]), CI_low = qs[1], CI_high = qs[2])
      })
      rows$Subgroup <- v; rows$Level <- as.character(lv); rows$N <- n
      rows[, c("Subgroup","Level","Metric","Point","CI_low","CI_high","N")]
    })
  }
  map_dfr(vars, each_var)
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

# ---------------- Fig 8: subgroup forest (AUC), pastel palette ----------------
fig8_subgroup_forest <- function(tbl7_subgroups, feature_set = "Triage",
                                 primary_algo = "LR", metric = "AUC") {
  d <- tbl7_subgroups %>%
    filter(Feature_Set == feature_set, Algorithm == primary_algo, Metric == metric) %>%
    filter(is.finite(Point))
  if (!nrow(d)) return(ggplot() + theme_void() + ggtitle("No subgroup rows found for plotting"))
  d$Facet <- d$Subgroup
  d$Label <- paste0(d$Level, " (n=", d$N, ")")
  ggplot(d, aes(Point, fct_reorder(Label, Point, .na_rm = TRUE),
                color = Split, shape = Split)) +
    { if (any(is.finite(d$CI_low) & is.finite(d$CI_high)))
      geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0,
                     position = position_dodge(width = 0.5)) else NULL } +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0.5, linetype = 3, color = "grey75") +
    facet_wrap(~Facet, scales = "free_y") +
    scale_color_brewer(palette = "Pastel1", drop = FALSE) +
    labs(title = sprintf("Subgroup performance, %s, algo=%s", feature_set, primary_algo),
         subtitle = "95% bootstrap CIs",
         x = metric, y = NULL) +
    theme_pub()
}

# ---------------- Fig 9: scorecard (Firth → Ridge), explicit event coding ----------------
fig9_scorecard <- function(df) {
  dtr <- df[df$data_split == "train", , drop = FALSE]
  if ("Diagnosis" %in% names(dtr)) dtr$Diagnosis <- factor(dtr$Diagnosis, levels = c("CTRL_noCOVID","COVID"))
  if ("Gender"    %in% names(dtr)) dtr$Gender    <- factor(dtr$Gender,    levels = c("Female","Male"))
  
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
  Xall <- recipes::juice(prepped)
  if (!(cfg$outcome %in% names(Xall))) Xall[[cfg$outcome]] <- dtr[[cfg$outcome]]
  
  coefs <- NULL
  used  <- NULL
  
  # ---- Firth with explicit event coding, no broom ----
  if (requireNamespace("logistf", quietly = TRUE)) {
    y_num <- as.integer(Xall[[cfg$outcome]] == cfg$pos_label)  # 1 = event "Yes"
    X2 <- Xall |> dplyr::select(-all_of(cfg$outcome))
    firth_fit <- try(logistf::logistf(y_num ~ ., data = X2), silent = TRUE)
    if (!inherits(firth_fit, "try-error")) {
      b  <- as.numeric(firth_fit$coefficients)
      tn <- names(firth_fit$coefficients)
      ci <- try(logistf::confint(firth_fit), silent = TRUE)
      if (!inherits(ci, "try-error")) {
        ci_df <- as.data.frame(ci)
        names(ci_df) <- tolower(names(ci_df))
        ci_df$term <- rownames(ci_df)
        coefs <- tibble(term = tn, estimate = exp(b)) |>
          left_join(tibble(term = ci_df$term,
                           conf.low = exp(ci_df$lower),
                           conf.high = exp(ci_df$upper)), by = "term")
      } else {
        coefs <- tibble(term = tn, estimate = exp(b),
                        conf.low = NA_real_, conf.high = NA_real_)
      }
      coefs <- filter(coefs, term != "(Intercept)")
      used  <- "Firth logistic"
    }
  }
  
  # ---- Ridge fallback if Firth missing or failed ----
  if (is.null(coefs)) {
    X2 <- Xall |> dplyr::select(-all_of(cfg$outcome))
    y_num <- as.integer(Xall[[cfg$outcome]] == cfg$pos_label)
    cv <- glmnet::cv.glmnet(as.matrix(X2), y_num, family = "binomial", alpha = 0, nfolds = 5)
    b  <- as.numeric(coef(cv, s = "lambda.min"))[-1]
    coefs <- tibble(term = colnames(X2), estimate = exp(b),
                    conf.low = NA_real_, conf.high = NA_real_)
    used <- "Ridge logistic, CIs not shown"
  }
  
  pretty_term <- function(t) {
    t <- gsub("\\.", "_", t)
    t <- sub("^Age$", "Age (per year)", t)
    t <- sub("^SpO2_admission$", "SpO2 at presentation (per %)", t)
    t <- sub("^Gender_Male$", "Gender = male", t)
    t <- sub("^Diagnosis_COVID$", "Diagnosis = COVID", t)
    t <- sub("^severity_admission_Mild$", "Severity_admission = Mild", t)
    t <- sub("^severity_admission_Severe$", "Severity_admission = Severe", t)
    t
  }
  
  coefs <- coefs %>%
    mutate(term = pretty_term(term),
           ok = is.finite(estimate),
           ok_ci = is.finite(conf.low) & is.finite(conf.high)) %>%
    filter(ok)
  
  p <- ggplot(coefs, aes(estimate, fct_reorder(term, estimate, .na_rm = TRUE))) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
    { if (any(coefs$ok_ci, na.rm = TRUE))
      geom_errorbarh(aes(xmin = pmax(conf.low, .Machine$double.eps),
                         xmax = pmax(conf.high, .Machine$double.eps)), height = 0) else NULL } +
    geom_point() +
    labs(title = paste0("Scorecard, ", used), x = "Odds ratio", y = NULL) + theme_pub()
  
  if (all(is.finite(coefs$estimate)) && all(coefs$estimate > 0, na.rm = TRUE)) p <- p + scale_x_log10()
  p
}

# ---------------- Ridge model specification table ----------------
final_triage_ridge_spec <- function(df) {
  stopifnot(all(c("data_split", cfg$outcome) %in% names(df)))
  feats <- intersect(feat_triage, names(df))
  if (!length(feats)) stop("No triage features found in data")
  dtr <- df[df$data_split == "train", c(cfg$outcome, feats), drop = FALSE]
  if (!nrow(dtr)) stop("No train rows found")
  
  if ("Diagnosis" %in% names(dtr)) dtr$Diagnosis <- factor(dtr$Diagnosis, levels = c("CTRL_noCOVID","COVID"))
  if ("Gender"    %in% names(dtr)) dtr$Gender    <- factor(dtr$Gender,    levels = c("Female","Male"))
  
  impute_once <- function(x) {
    if (is.numeric(x)) { m <- median(x, na.rm = TRUE); x[is.na(x)] <- m; return(x) }
    x <- as.factor(x); md <- names(sort(table(x), decreasing = TRUE))[1]
    x[is.na(x)] <- md; droplevels(x)
  }
  dtr[feats] <- lapply(dtr[feats], impute_once)
  
  X <- model.matrix(reformulate(feats), data = dtr)[, -1, drop = FALSE]
  y <- as.integer(dtr[[cfg$outcome]] == cfg$pos_label)
  
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 5, intercept = TRUE)
  co <- as.matrix(coef(cv, s = "lambda.min"))
  beta <- as.numeric(co[, 1])
  nm   <- rownames(co)
  
  prettify <- function(s) {
    if (s == "(Intercept)") return("Intercept")
    s <- gsub("\\.", "_", s)
    s <- sub("^Age$", "Age (per year)", s)
    s <- sub("^SpO2_admission$", "SpO2 at presentation (per %)", s)
    s <- sub("^GenderMale$", "Gender = male", s)
    s <- sub("^DiagnosisCOVID$", "Diagnosis = COVID", s)
    s <- sub("^severity_admissionMild$", "Severity_admission = Mild", s)
    s <- sub("^severity_admissionSevere$", "Severity_admission = Severe", s)
    s
  }
  
  out <- tibble(Term = vapply(nm, prettify, character(1)),
                Coefficient = beta) %>%
    group_by(Term) %>% slice(1) %>% ungroup()
  
  pref <- c("Intercept",
            "Age (per year)",
            "Gender = male",
            "Diagnosis = COVID",
            "Severity_admission = Mild",
            "Severity_admission = Severe",
            "SpO2 at presentation (per %)")
  out %>% arrange(match(Term, pref), Term)
}

# ---------------- pipeline ----------------
run_pipeline_section <- function() {
  df <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, feat_triage)
  
  if (!("data_split" %in% names(df)) || all(is.na(df$data_split))) {
    set.seed(cfg$seed_cv)
    n <- nrow(df); idx <- sample(n); n_tr <- floor(0.7*n); n_te <- floor(0.15*n)
    df$data_split <- factor(c(rep("train", n_tr), rep("test", n_te), rep("external", n - n_tr - n_te))[order(idx)])
  }
  
  if ("Diagnosis" %in% names(df)) df$Diagnosis <- factor(df$Diagnosis, levels = c("CTRL_noCOVID","COVID"))
  if ("Gender"    %in% names(df)) df$Gender    <- factor(df$Gender,    levels = c("Female","Male"))
  
  invisible(run_cv(df, cfg$outcome, features = feat_triage, algos = c("LR"),
                   pos_label = cfg$pos_label, inner_k = cfg$inner_k,
                   k_desired = cfg$cv_k, R = cfg$cv_R))
  
  holder <- list(objects = list(LR = run_holdout(df, cfg$outcome, features = feat_triage, algo_name = "LR")),
                 set_name = "Triage")
  
  sub_collect <- function(holder, split_name) {
    objs <- holder$objects; out <- list()
    for (algo in names(objs)) {
      o <- objs[[algo]][[split_name]]; if (is.null(o)) next
      sg <- subgroup_metrics(o$holdout_df, o$obs, o$p, o$pred,
                             vars = c("group","Gender","Diagnosis","severity_admission"),
                             pos_label = cfg$pos_label)
      if (!nrow(sg)) next
      sg$Algorithm <- algo; sg$Feature_Set <- holder$set_name; sg$Split <- tools::toTitleCase(split_name)
      out[[algo]] <- sg
    }
    dplyr::bind_rows(out)
  }
  tbl7_subgroups <- dplyr::bind_rows(sub_collect(holder, "test"),
                                     sub_collect(holder, "external"))
  
  list(df = df, tbl7_subgroups = tbl7_subgroups, holder = holder)
}

# ---------------- Driver ----------------
print_subgroup_and_scorecard <- function(primary_algo = "LR") {
  cat("[setup] running section pipeline...\n")
  pipe <- run_pipeline_section()
  
  cat("[fig] printing Fig 8, subgroup AUC forest\n")
  p8 <- fig8_subgroup_forest(pipe$tbl7_subgroups, feature_set = "Triage",
                             primary_algo = primary_algo, metric = "AUC")
  print(p8)
  
  tbl_sub_auc <- pipe$tbl7_subgroups %>%
    filter(Feature_Set == "Triage", Algorithm == primary_algo, Metric == "AUC") %>%
    mutate(CI_low = ifelse(is.finite(CI_low), pmax(0, CI_low), CI_low),
           CI_high = ifelse(is.finite(CI_high), pmin(1, CI_high), CI_high)) %>%
    arrange(Subgroup, Level, Split)
  cat("\n[Table] Subgroup performance (AUC), triage LR\n")
  print(tbl_sub_auc, n = Inf)
  
  cat("\n[fig] printing Fig 9, scorecard\n")
  p9 <- fig9_scorecard(pipe$df)
  print(p9)
  
  cat("\n[Table] Final triage model specification, coefficients on the analysis scale\n")
  coef_tbl <- final_triage_ridge_spec(pipe$df)
  print(coef_tbl, n = Inf)
  
  invisible(list(fig8 = p8, tbl_sub_auc = tbl_sub_auc, fig9 = p9, coef_tbl = coef_tbl, pipe = pipe))
}

# =========================
# RUN NOW
# =========================
print_subgroup_and_scorecard(primary_algo = "LR")
