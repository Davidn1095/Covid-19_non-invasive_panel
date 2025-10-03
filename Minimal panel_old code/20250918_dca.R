# ============================================================
# Decision analysis with temporal split, overlay Complete vs Triage
# - Train and Test: dates <= 2020-04-21
# - External: dates > 2020-04-21
# - Curves truncated where each model's net benefit reaches 0
# - Complete feature set = pastel red (#F4CCCC), Triage feature set = pastel green (#D9EAD3)
# - Treat-all = solid black, Treat-none = dotted grey
# - One shared legend on the right, no titles
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(pROC)
  library(patchwork)
  if (requireNamespace("glmnet", quietly = TRUE)) library(glmnet)
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
  boundary  = as.Date("2020-04-21"),
  pct_test  = 0.20,
  algos_req = c("C4.5","k-NN","SVM","RF","LR"),
  primary_algo = "C4.5"
)

# Predictors
feat_full <- c(
  "Diagnosis","severity_admission","Age","Gender",
  "SpO2_admission",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc"
)
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- date helpers ----------------
norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

try_parse_date <- function(x) {
  if (inherits(x, "Date"))   return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    d1 <- as.Date(x, origin = "1899-12-30")
    d2 <- as.Date(x, origin = "1970-01-01")
    inrng1 <- sum(d1 >= as.Date("2019-01-01") & d1 <= as.Date("2023-12-31"), na.rm=TRUE)
    inrng2 <- sum(d2 >= as.Date("2019-01-01") & d2 <= as.Date("2023-12-31"), na.rm=TRUE)
    return(if (inrng2 > inrng1) d2 else d1)
  }
  xs <- trimws(as.character(x))
  fmts <- c("%d/%m/%Y","%Y-%m-%d","%m/%d/%Y","%d-%m-%Y","%d.%m.%Y","%Y/%m/%d")
  mats <- lapply(fmts, function(f) suppressWarnings(as.Date(xs, format = f)))
  counts <- vapply(mats, function(d) sum(d >= as.Date("2019-01-01") & d <= as.Date("2023-12-31"), na.rm=TRUE), integer(1))
  if (all(counts == 0)) return(suppressWarnings(as.Date(xs)))
  mats[[which.max(counts)]]
}

find_sampling_date <- function(df) {
  nms_raw <- names(df)
  nms     <- norm_nm(nms_raw)
  has <- function(p) grepl(p, nms, fixed = FALSE)
  cand_idx <- which(
    (has("sampling") | has("sample") | has("collection") | has("admission") | has("draw") | has("visit")) &
      has("date")
  )
  if (!length(cand_idx)) cand_idx <- which(has("samplingdate") | has("date"))
  score <- function(col) {
    d <- try_parse_date(df[[col]])
    sum(d >= as.Date("2019-01-01") & d <= as.Date("2023-12-31"), na.rm=TRUE)
  }
  if (length(cand_idx)) {
    sc <- vapply(cand_idx, score, integer(1))
    if (max(sc, na.rm=TRUE) > 0) return(nms_raw[cand_idx[which.max(sc)]])
  }
  sc_all <- vapply(nms_raw, score, integer(1))
  if (max(sc_all, na.rm=TRUE) > 0) return(nms_raw[which.max(sc_all)])
  stop("A sampling date column was not found")
}

# ---------------- load, filter, temporal split ----------------
load_and_split <- function(path, sheet, outcome, boundary, pct_test = 0.2) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet = sheet) |> as.data.frame()
  
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df)==grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  for (v in intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' not found", outcome))
  to_YN <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) return(factor(ifelse(x, cfg$pos_label, cfg$neg_label), levels = c(cfg$pos_label, cfg$neg_label)))
    if (is.numeric(x)) {
      stopifnot(all(x %in% c(0,1), na.rm = TRUE))
      return(factor(ifelse(x==1, cfg$pos_label, cfg$neg_label), levels = c(cfg$pos_label, cfg$neg_label)))
    }
    if (is.character(x)) {
      s <- trimws(tolower(x))
      yes <- c("1","yes","y","true","pos","positive")
      no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% yes, cfg$pos_label, ifelse(s %in% no, cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome could not be coerced to Yes/No")
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    stop("Unsupported outcome type")
  }
  df[[outcome]] <- to_YN(df[[outcome]])
  
  dcol <- find_sampling_date(df)
  dts  <- try_parse_date(df[[dcol]])
  if (all(is.na(dts))) stop(sprintf("Dates in '%s' could not be parsed", dcol))
  
  pre_idx  <- which(dts <= boundary)
  post_idx <- which(dts >  boundary)
  set.seed(444)
  n_pre <- length(pre_idx)
  n_te  <- max(1L, floor(pct_test * n_pre))
  te_ids <- if (n_pre > 1) sample(pre_idx, n_te) else integer(0)
  tr_ids <- setdiff(pre_idx, te_ids)
  
  sp <- rep(NA_character_, nrow(df))
  sp[tr_ids]   <- "train"
  sp[te_ids]   <- "test"
  sp[post_idx] <- "external"
  df$data_split <- factor(sp, levels = c("train","test","external"))
  df$.rid <- seq_len(nrow(df))
  attr(df, "temporal_boundary") <- boundary
  df
}

# ---------------- modeling helpers ----------------
get_neg_level <- function(y) {
  levs <- levels(y)
  nl <- setdiff(levs, cfg$pos_label)
  if (length(nl)) nl[1] else cfg$neg_label
}

mcc_from_counts <- function(TP, FP, FN, TN) {
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (den == 0) return(NA_real_)
  num/den
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
  if (exists("step_YeoJohnson", where = asNamespace("recipes"), inherits = FALSE)) {
    rec <- rec |> step_YeoJohnson(all_numeric_predictors())
  } else if (exists("step_yeojohnson", where = asNamespace("recipes"), inherits = FALSE)) {
    rec <- rec |> step_yeojohnson(all_numeric_predictors())
  }
  rec |> step_normalize(all_numeric_predictors())
}

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

model_specs <- function(tune_len = 5, algos = cfg$algos_req) {
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

train_inner <- function(dat_train, yvar, spec, inner_k = 3) {
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

# ---- Platt scaling on train, applied to external ----
safe_logit <- function(p, eps = 1e-6) qlogis(pmin(pmax(as.numeric(p), eps), 1 - eps))
platt_fit <- function(obs, p_train) {
  y <- as.integer(obs == cfg$pos_label)
  s <- safe_logit(p_train)
  suppressWarnings(glm(y ~ s, family = binomial()))
}
platt_apply <- function(fit, p_new) {
  s <- safe_logit(p_new)
  as.numeric(plogis(drop(cbind(1, s) %*% coef(fit))))
}

run_holdout <- function(df, yvar, features, algo_name,
                        split_train = "train", split_test = c("test","external")) {
  present_splits <- intersect(split_test, unique(as.character(df$data_split)))
  if (!length(present_splits)) return(NULL)
  specs <- model_specs(algos = algo_name)
  if (!length(specs)) return(NULL)
  spec  <- specs[[algo_name]]
  base_cols <- unique(c(features, yvar, intersect("group", names(df))))
  base_cols <- intersect(base_cols, names(df))
  dtrain <- df[df$data_split == split_train, base_cols, drop = FALSE]
  if (!nrow(dtrain)) return(NULL)
  fit <- train_inner(dtrain, yvar, spec, inner_k = 3)
  
  # train probabilities for threshold and calibration
  p_tr  <- as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label])
  neg_label <- get_neg_level(dtrain[[yvar]])
  thr   <- best_threshold_mcc(dtrain[[yvar]], p_tr, pos = cfg$pos_label, neg = neg_label)
  
  # Platt calibration model from train
  platt <- tryCatch(platt_fit(dtrain[[yvar]], p_tr), error = function(e) NULL)
  
  out <- list()
  for (sp in present_splits) {
    dtest <- df[df$data_split == sp, base_cols, drop = FALSE]
    if (!nrow(dtest)) next
    p_raw <- as.numeric(predict(fit, newdata = dtest, type = "prob")[, cfg$pos_label])
    
    # Apply Platt only to external split
    p_use <- if (!is.null(platt) && identical(sp, "external")) platt_apply(platt, p_raw) else p_raw
    
    pred  <- factor(ifelse(p_use >= thr$t, cfg$pos_label, neg_label),
                    levels = levels(dtrain[[yvar]]))
    out[[sp]] <- list(
      obs = factor(dtest[[yvar]], levels = levels(dtrain[[yvar]])),
      p   = p_use, pred = pred, threshold = thr$t,
      p_raw = p_raw, p_train = p_tr, model = fit, platt = platt, holdout_df = dtest
    )
  }
  out
}

# ---------------- decision curve (robust to NAs) ----------------
decision_curve_table <- function(obs, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
  y_raw <- as.integer(obs == cfg$pos_label)
  keep  <- is.finite(as.numeric(p)) & !is.na(y_raw)
  y <- y_raw[keep]; p <- as.numeric(p)[keep]
  if (!length(y)) {
    return(tibble(
      threshold = thresholds, Net_Benefit = NA_real_,
      NB_TreatAll = NA_real_, NB_TreatNone = 0,
      Precision = NA_real_, NNE = NA_real_
    ))
  }
  N <- length(y); prev <- mean(y)
  purrr::map_dfr(thresholds, function(pt) {
    pred <- as.integer(p >= pt)
    TP <- sum(pred == 1 & y == 1, na.rm = TRUE)
    FP <- sum(pred == 1 & y == 0, na.rm = TRUE)
    NB <- TP/N - FP/N * (pt/(1-pt))
    NB_all <- prev - (1 - prev) * (pt/(1-pt))
    pos_calls <- sum(pred == 1, na.rm = TRUE)
    Precision <- if (pos_calls == 0) NA_real_ else TP / pos_calls
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

# threshold where model NB hits 0, linear interpolation
nb_zero_threshold <- function(dca_df) {
  d <- dca_df[order(dca_df$threshold), ]
  nb <- d$Net_Benefit
  th <- d$threshold
  bad <- !is.finite(nb) | !is.finite(th)
  nb <- nb[!bad]; th <- th[!bad]
  if (!length(nb)) return(NA_real_)
  i_neg <- which(nb < 0)
  if (!length(i_neg)) return(max(th, na.rm = TRUE))
  i <- i_neg[1]
  if (i == 1) return(th[1])
  t0 <- th[i - 1]; n0 <- nb[i - 1]
  t1 <- th[i];     n1 <- nb[i]
  if (!is.finite(n0) || !is.finite(n1) || t1 == t0) return(th[i])
  t_star <- t0 - n0 * (t1 - t0) / (n1 - n0)
  max(min(t_star, t1), t0)
}

# ---------------- plotting, overlay Complete vs Triage ----------------
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(color = "grey30"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

pastel_red   <- "#F4CCCC"  # Complete feature set
pastel_green <- "#D9EAD3"  # Triage feature set

overlay_dca_until_zero <- function(o_full, o_tri) {
  if (is.null(o_full) && is.null(o_tri)) return(ggplot() + theme_void())
  
  df_full <- if (!is.null(o_full)) decision_curve_table(o_full$obs, o_full$p) else NULL
  df_tri  <- if (!is.null(o_tri))  decision_curve_table(o_tri$obs,  o_tri$p)  else NULL
  
  t_full <- if (!is.null(df_full)) nb_zero_threshold(df_full) else NA_real_
  t_tri  <- if (!is.null(df_tri))  nb_zero_threshold(df_tri)  else NA_real_
  t_max  <- suppressWarnings(max(c(t_full, t_tri), na.rm = TRUE))
  if (!is.finite(t_max)) {
    t_max <- suppressWarnings(max(c(df_full$threshold, df_tri$threshold), na.rm = TRUE))
    if (!is.finite(t_max)) t_max <- 0.99
  }
  
  if (!is.null(df_full) && is.finite(t_full)) df_full <- dplyr::filter(df_full, threshold <= t_full)
  if (!is.null(df_tri)  && is.finite(t_tri))  df_tri  <- dplyr::filter(df_tri,  threshold <= t_tri)
  
  base <- if (!is.null(o_tri)) decision_curve_table(o_tri$obs, o_tri$p) else decision_curve_table(o_full$obs, o_full$p)
  base <- base %>% dplyr::filter(threshold <= t_max)
  
  y_max <- max(
    if (!is.null(df_full)) df_full$Net_Benefit else 0,
    if (!is.null(df_tri))  df_tri$Net_Benefit  else 0,
    base$NB_TreatAll,
    0, na.rm = TRUE
  )
  y_max <- if (is.finite(y_max) && y_max > 0) y_max else 0.05
  
  g <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60", show.legend = FALSE) +
    geom_line(data = base, aes(threshold, NB_TreatAll, linetype = "Treat-all"), color = "black", linewidth = 0.6, na.rm = TRUE) +
    geom_line(data = base, aes(threshold, NB_TreatNone, linetype = "Treat-none"), color = "grey40", linewidth = 0.6, na.rm = TRUE) +
    scale_linetype_manual(values = c("Treat-all" = "solid", "Treat-none" = "dotted"), name = NULL) +
    coord_cartesian(xlim = c(0, t_max), ylim = c(0, y_max)) +
    labs(title = NULL, x = "Threshold probability", y = "Net benefit") +
    theme_pub()
  
  if (!is.null(df_full) && nrow(df_full)) {
    g <- g + geom_line(data = df_full, aes(threshold, Net_Benefit, color = "Complete feature set"), linewidth = 0.9, show.legend = TRUE, na.rm = TRUE)
  }
  if (!is.null(df_tri) && nrow(df_tri)) {
    g <- g + geom_line(data = df_tri,  aes(threshold, Net_Benefit, color = "Triage feature set"), linewidth = 0.9, show.legend = TRUE, na.rm = TRUE)
  }
  
  g + scale_color_manual(values = c("Complete feature set" = pastel_red, "Triage feature set" = pastel_green), name = "Model")
}

# ===================== RUN AND DRAW ALL FIVE IN ONE GO =====================
df <- load_and_split(cfg$file_path, cfg$sheet, cfg$outcome, cfg$boundary, cfg$pct_test)

algos_avail <- names(model_specs(algos = cfg$algos_req))
if (!length(algos_avail)) stop("No requested algorithms are available. Please install: RWeka, kknn, kernlab, randomForest, glmnet")

hold_full <- lapply(algos_avail, function(a) run_holdout(df, cfg$outcome, features = feat_full,   algo_name = a))
names(hold_full) <- algos_avail; hold_full <- list(set_name = "Full",   objects = hold_full)

hold_tri  <- lapply(algos_avail, function(a) run_holdout(df, cfg$outcome, features = feat_triage, algo_name = a))
names(hold_tri)  <- algos_avail; hold_tri  <- list(set_name = "Triage", objects = hold_tri)

algos_all <- intersect(names(hold_full$objects), names(hold_tri$objects))
if (!length(algos_all)) stop("No algorithms available in both Full and Triage results")

make_two_panel_overlay <- function(a, show_legend = FALSE) {
  o_full_test <- hold_full$objects[[a]][["test"]]
  o_tri_test  <- hold_tri$objects[[a]][["test"]]
  o_full_ext  <- hold_full$objects[[a]][["external"]]
  o_tri_ext   <- hold_tri$objects[[a]][["external"]]
  p_test <- overlay_dca_until_zero(o_full_test, o_tri_test)
  p_ext  <- overlay_dca_until_zero(o_full_ext,  o_tri_ext)
  comb <- if (inherits(p_test, "ggplot") && inherits(p_ext, "ggplot")) {
    (p_test | p_ext) + plot_layout(guides = if (show_legend) "collect" else "keep")
  } else if (inherits(p_test, "ggplot")) {
    p_test + plot_layout(guides = if (show_legend) "collect" else "keep")
  } else if (inherits(p_ext, "ggplot")) {
    p_ext + plot_layout(guides = if (show_legend) "collect" else "keep")
  } else {
    ggplot() + theme_void()
  }
  if (!show_legend) comb <- comb & theme(legend.position = "none")
  comb & guides(
    linetype = guide_legend(
      override.aes = list(color = c("black","grey40"),
                          linetype = c("solid","dotted"))
    )
  )
}

plots_by_algo <- lapply(seq_along(algos_all), function(i) {
  make_two_panel_overlay(algos_all[[i]], show_legend = (i == 1))
})

panel <- wrap_plots(plots_by_algo, ncol = 1, guides = "collect") &
  theme(legend.position = "right")

# ---- draw, with a forced draw fallback in case RStudio swallows it ----
print(panel)
try({
  grid::grid.newpage()
  grid::grid.draw(patchwork::patchworkGrob(panel))
  grDevices::dev.flush()
}, silent = TRUE)

# Optional save (kept commented)
# ggsave("DCA_overlay_all_algos.png", panel,
#        width = 10, height = 5.5 * length(plots_by_algo), dpi = 300, bg = "white")
