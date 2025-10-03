# ============================================================
# Parsimony curves by algorithm (MCC) with ±1 SD ribbons
# Temporal split, identical CV folds across models and sets
# 2 rows layout: C4.5, k-NN, SVM on top, RF, LR bottom, legend bottom-right
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(glmnet)
  library(pROC)
  library(ranger)
  library(kknn)
  library(RWeka)
  library(kernlab)
  library(patchwork)
  library(cowplot)
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
  inner_k   = 3,
  tune_len  = 10
)

# algorithms in required facet order
algo_map <- c("C4.5" = "J48", "k-NN" = "kknn", "SVM" = "svmRadial", "RF" = "ranger", "LR" = "glmnet")

needs_dummies <- function(method) method %in% c("glmnet","svmRadial","kknn")
needs_scaling <- function(method) method %in% c("glmnet","svmRadial","kknn")

# predictors
feat_full   <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission",
                 "monocyte_abs_number","monocytes_perc",
                 "lymphocyte_abs_number","lymphocytes_perc",
                 "neutrophil_abs_number","neutrophils_perc")
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- date helpers ----------------
norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
try_parse_date <- function(x) {
  if (inherits(x, "Date"))   return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    d1 <- as.Date(x, origin = "1899-12-30")
    d2 <- as.Date(x, origin = "1970-01-01")
    r <- as.Date("2019-01-01"); R <- as.Date("2023-12-31")
    if (sum(d2>=r & d2<=R, na.rm=TRUE) > sum(d1>=r & d1<=R, na.rm=TRUE)) d2 else d1
  } else {
    xs <- trimws(as.character(x))
    fmts <- c("%d/%m/%Y","%Y-%m-%d","%m/%d/%Y","%d-%m-%Y","%d.%m.%Y","%Y/%m/%d")
    mats <- lapply(fmts, function(f) suppressWarnings(as.Date(xs, format=f)))
    r <- as.Date("2019-01-01"); R <- as.Date("2023-12-31")
    counts <- vapply(mats, function(d) sum(d>=r & d<=R, na.rm=TRUE), integer(1))
    if (all(counts==0)) suppressWarnings(as.Date(xs)) else mats[[which.max(counts)]]
  }
}
find_sampling_date <- function(df) {
  nms_raw <- names(df); nms <- norm_nm(nms_raw); has <- function(p) grepl(p, nms)
  cand <- which((has("sampling")|has("sample")|has("collection")|has("admission")|has("draw")|has("visit")) & has("date"))
  if (!length(cand)) cand <- which(has("samplingdate")|has("date"))
  score <- function(col){ d<-try_parse_date(df[[col]]); sum(d>=as.Date("2019-01-01") & d<=as.Date("2023-12-31"), na.rm=TRUE)}
  if (length(cand)) {
    sc <- vapply(cand, function(i) score(nms_raw[i]), integer(1))
    if (max(sc,na.rm=TRUE)>0) return(nms_raw[cand[which.max(sc)]])
  }
  sc_all <- vapply(nms_raw, score, integer(1))
  if (max(sc_all,na.rm=TRUE)>0) return(nms_raw[which.max(sc_all)])
  stop("A sampling date column was not found")
}

# ---------------- load + temporal split ----------------
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
      stopifnot(all(x %in% c(0,1), na.rm=TRUE))
      return(factor(ifelse(x==1,cfg$pos_label,cfg$neg_label), levels = c(cfg$pos_label, cfg$neg_label)))
    }
    if (is.character(x)) {
      s <- trimws(tolower(x))
      yes <- c("1","yes","y","true","pos","positive")
      no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% yes, cfg$pos_label, ifelse(s %in% no, cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome could not be coerced to Yes/No")
      factor(out, levels = c(cfg$pos_label, cfg$neg_label))
    } else stop("Unsupported outcome type")
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
  df
}

# ---------------- metrics ----------------
get_neg_level <- function(y) {
  levs <- levels(y); nl <- setdiff(levs, cfg$pos_label)
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
  aucv  <- {
    y01 <- as.integer(y == pos_label)
    if (length(unique(y01)) < 2) NA_real_
    else suppressMessages(as.numeric(pROC::auc(y01, as_num(p_pos), quiet = TRUE)))
  }
  c(MCC = mcc, AUC = aucv, `F1` = f1, Accuracy = acc, Precision = preci, Sensitivity = sens, Specificity = spec)
}

# caret summary for MCC, robust to degenerate folds
mccSummary <- function(data, lev = NULL, model = NULL) {
  pos <- cfg$pos_label
  neg <- setdiff(lev, pos)[1] %||% cfg$neg_label
  p   <- as_num(data[[pos]])
  one_class <- length(unique(as.character(data$obs))) < 2
  thr_obj <- try(best_threshold_mcc(data$obs, p, pos = pos, neg = neg), silent = TRUE)
  thr <- if (inherits(thr_obj, "try-error") || !is.finite(thr_obj$mcc)) 0.5 else thr_obj$t
  pred <- factor(ifelse(p >= thr, pos, neg), levels = c(neg, pos))
  met  <- compute_metrics_binary(data$obs, pred, p, pos_label = pos, neg_label = neg)
  m    <- unname(met["MCC"])
  if (!is.finite(m) || one_class) m <- 0
  c(MCC = m)
}

# ---------------- stratified folds, reused everywhere ----------------
safe_cv_index <- function(y, k_desired, seed = 2020) {
  set.seed(seed)
  k_safe <- max(2, min(k_desired, min(table(y))))
  caret::createFolds(y, k = k_safe, returnTrain = TRUE)
}

# ---------------- recipe per algorithm ----------------
make_recipe_for <- function(dat, yvar, caret_method) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c("data_split",".rid","group"), names(dat))
  if (length(ign)) rec <- rec |> recipes::update_role(all_of(ign), new_role = "ignore")
  rec <- rec |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors())  |>
    recipes::step_novel(recipes::all_nominal_predictors())        |>
    recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01)
  if (needs_dummies(caret_method)) rec <- rec |> recipes::step_dummy(recipes::all_nominal_predictors())
  rec <- rec |> recipes::step_zv(recipes::all_predictors())
  if (needs_scaling(caret_method)) {
    rec <- rec |>
      recipes::step_YeoJohnson(recipes::all_numeric_predictors()) |>
      recipes::step_normalize(recipes::all_numeric_predictors())
  }
  rec
}

# ---------------- permutation importance (ΔMCC) ----------------
perm_drop_mcc <- function(fit, data, yvar, features, threshold, pos_label = cfg$pos_label) {
  p0 <- as.numeric(predict(fit, newdata = data, type = "prob")[, pos_label])
  neg_label <- get_neg_level(data[[yvar]])
  pred0 <- factor(ifelse(p0 >= threshold, pos_label, neg_label),
                  levels = levels(data[[yvar]]))
  m0 <- as.numeric(compute_metrics_binary(data[[yvar]], pred0, p0, pos_label = pos_label, neg_label = neg_label)["MCC"])
  purrr::map_dfr(features, function(f) {
    dperm <- data
    dperm[[f]] <- sample(dperm[[f]])
    pp <- as.numeric(predict(fit, newdata = dperm, type = "prob")[, pos_label])
    pr <- factor(ifelse(pp >= threshold, pos_label, neg_label),
                 levels = levels(data[[yvar]]))
    m1 <- as.numeric(compute_metrics_binary(data[[yvar]], pr, pp, pos_label = pos_label, neg_label = neg_label)["MCC"])
    tibble(Feature = f, Drop_MCC = m0 - m1)
  }) %>% arrange(desc(Drop_MCC))
}

# ---- helper: minimum predictors across folds after recipe ----
min_predictors_across_folds <- function(rec, data, idx, yvar) {
  pvec <- vapply(idx, function(tr) {
    rec_i <- recipes::prep(rec, training = data[tr, , drop = FALSE], retain = TRUE)
    max(0L, ncol(recipes::juice(rec_i)) - 1L)
  }, integer(1))
  max(1L, min(pvec, na.rm = TRUE))
}

# ---- robust CV mean, SD for the best tuning row ----
cv_mean_sd <- function(fit) {
  rs <- fit$resample
  if (is.null(rs) || !nrow(rs) || !("MCC" %in% names(rs))) {
    m <- suppressWarnings(max(fit$results$MCC, na.rm = TRUE))
    m <- if (is.finite(m)) m else NA_real_
    return(c(mean = m, sd = NA_real_))
  }
  bt <- fit$bestTune
  if (!is.null(bt) && ncol(bt) > 0) {
    keep <- intersect(names(bt), names(rs))
    if (length(keep)) for (nm in keep) rs <- rs[rs[[nm]] == bt[[nm]], , drop = FALSE]
  }
  m <- mean(rs$MCC, na.rm = TRUE)
  s <- sd(rs$MCC, na.rm = TRUE)
  if (!is.finite(m)) {
    m <- suppressWarnings(max(fit$results$MCC, na.rm = TRUE))
    m <- if (is.finite(m)) m else NA_real_
    s <- NA_real_
  }
  c(mean = m, sd = s)
}

# ---- rank once on TRAIN with LR, reused for all algorithms ----
rank_by_perm_on_train_lr <- function(df, yvar, feats, idx) {
  feats <- intersect(feats, names(df))
  dtr <- df[df$data_split == "train", c(yvar, feats, intersect("group", names(df))), drop = FALSE]
  stopifnot(nrow(dtr) > 0)
  rec <- make_recipe_for(dtr, yvar, "glmnet")
  tr_ctrl <- caret::trainControl(method = "cv", index = idx,
                                 classProbs = TRUE, summaryFunction = twoClassSummary)
  p_min <- min_predictors_across_folds(rec, dtr, idx, yvar)
  if (p_min < 2) {
    fit <- caret::train(rec, data = dtr, method = "glm",
                        trControl = tr_ctrl, metric = "ROC", family = binomial())
  } else {
    fit <- caret::train(rec, data = dtr, method = "glmnet",
                        trControl = tr_ctrl, metric = "ROC", tuneLength = cfg$tune_len)
  }
  p_tr  <- as.numeric(predict(fit, newdata = dtr, type = "prob")[, cfg$pos_label])
  neg   <- get_neg_level(dtr[[yvar]])
  thr   <- best_threshold_mcc(dtr[[yvar]], p_tr, pos = cfg$pos_label, neg = neg)$t
  imp   <- perm_drop_mcc(fit, dtr, yvar, feats, thr)
  imp |> arrange(desc(Drop_MCC)) |> dplyr::pull(Feature)
}

# ---- parsimony curve for one algorithm and one set, returns mean and SD ----
parsimony_curve_MCC_algo <- function(df, feats, set_name, algo_key, order_feats, idx) {
  dtr <- df[df$data_split == "train", , drop = FALSE]
  feats <- intersect(feats, names(df))
  stopifnot(length(feats) > 0)
  
  caret_method <- unname(algo_map[algo_key])
  if (!length(order_feats)) order_feats <- feats
  
  tr_ctrl <- caret::trainControl(method = "cv", index = idx,
                                 classProbs = TRUE, summaryFunction = mccSummary)
  
  res <- purrr::map_dfr(seq_along(order_feats), function(k) {
    feats_k <- order_feats[seq_len(k)]
    rec_k <- recipes::recipe(
      stats::as.formula(paste(cfg$outcome, "~", paste(feats_k, collapse = "+"))),
      data = dtr
    ) |>
      recipes::step_impute_median(recipes::all_numeric_predictors()) |>
      recipes::step_impute_mode(recipes::all_nominal_predictors())  |>
      recipes::step_novel(recipes::all_nominal_predictors())        |>
      recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01)
    
    if (needs_dummies(caret_method)) rec_k <- rec_k |> recipes::step_dummy(recipes::all_nominal_predictors())
    rec_k <- rec_k |> recipes::step_zv(recipes::all_predictors())
    if (needs_scaling(caret_method)) {
      rec_k <- rec_k |>
        recipes::step_YeoJohnson(recipes::all_numeric_predictors()) |>
        recipes::step_normalize(recipes::all_numeric_predictors())
    }
    
    if (caret_method == "glmnet") {
      p_min_k <- min_predictors_across_folds(rec_k, dtr, idx, cfg$outcome)
      if (p_min_k < 2) {
        fit_k <- caret::train(rec_k, data = dtr, method = "glm",
                              trControl = tr_ctrl, metric = "MCC", family = binomial())
      } else {
        fit_k <- caret::train(rec_k, data = dtr, method = "glmnet",
                              trControl = tr_ctrl, metric = "MCC",
                              tuneLength = cfg$tune_len)
      }
    } else if (caret_method == "J48") {
      fit_k <- caret::train(rec_k, data = dtr, method = "J48",
                            trControl = tr_ctrl, metric = "MCC",
                            tuneLength = cfg$tune_len)
    } else if (caret_method == "svmRadial") {
      fit_k <- caret::train(rec_k, data = dtr, method = "svmRadial",
                            trControl = tr_ctrl, metric = "MCC",
                            tuneLength = cfg$tune_len)
    } else if (caret_method == "kknn") {
      fit_k <- caret::train(rec_k, data = dtr, method = "kknn",
                            trControl = tr_ctrl, metric = "MCC",
                            tuneLength = cfg$tune_len)
    } else if (caret_method == "ranger") {
      p_min_k <- min_predictors_across_folds(rec_k, dtr, idx, cfg$outcome)
      mtry_vals <- unique(pmax(1, pmin(p_min_k, round(seq(1, p_min_k, length.out = min(cfg$tune_len, p_min_k))))))
      fit_k <- caret::train(rec_k, data = dtr, method = "ranger",
                            trControl = tr_ctrl, metric = "MCC",
                            tuneGrid = expand.grid(mtry = mtry_vals, splitrule = "gini", min.node.size = 1),
                            num.trees = 1000, importance = "none")
    } else {
      fit_k <- caret::train(rec_k, data = dtr, method = caret_method,
                            trControl = tr_ctrl, metric = "MCC",
                            tuneLength = cfg$tune_len)
    }
    
    ms <- cv_mean_sd(fit_k)
    tibble(k = k, Score = unname(ms["mean"]), SD = unname(ms["sd"]))
  })
  res$Set <- set_name
  res$Algorithm <- algo_key
  res
}

# ---------------- plotting theme ----------------
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "grey30"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

# ============================================================
# RUN
# ============================================================
df <- load_and_split(cfg$file_path, cfg$sheet, cfg$outcome, cfg$boundary, cfg$pct_test)

# identical folds reused everywhere
y_tr <- df[[cfg$outcome]][df$data_split == "train"]
idx_once <- safe_cv_index(y_tr, cfg$inner_k, seed = 2020)

# rank once per set with LR
order_full <- rank_by_perm_on_train_lr(df, cfg$outcome, feat_full,   idx_once)
order_tri  <- rank_by_perm_on_train_lr(df, cfg$outcome, feat_triage, idx_once)

sets <- list(
  "Complete feature set" = list(feats = feat_full,   order = order_full),
  "Triage feature set"   = list(feats = feat_triage, order = order_tri)
)

res_all <- purrr::imap_dfr(sets, function(spec, set_name) {
  purrr::map_dfr(names(algo_map), function(ak) {
    parsimony_curve_MCC_algo(df, spec$feats, set_name, ak, order_feats = spec$order, idx = idx_once)
  })
}) |>
  dplyr::filter(is.finite(Score)) |>
  dplyr::mutate(Algorithm = factor(Algorithm, levels = names(algo_map)))

# ---------------- custom 2×3 layout with legend bottom-right ----------------
make_alg_plot <- function(alg, show_y = FALSE, show_x = FALSE) {
  dat <- dplyr::filter(res_all, Algorithm == alg)
  
  # in-panel axis lines
  x_min <- min(dat$k) - 0.5
  x_max <- max(dat$k) + 0.5
  y_min <- 0.25
  y_max <- 1.00
  
  ggplot(dat, aes(k, Score, color = Set)) +
    # black axes drawn inside the panel so they are not clipped
    annotate("segment", x = x_min, xend = x_max, y = y_min, yend = y_min,
             linewidth = 0.6, colour = "black") +
    annotate("segment", x = x_min, xend = x_min, y = y_min, yend = y_max,
             linewidth = 0.6, colour = "black") +
    # ±1 SD ribbons, no fill legend
    geom_ribbon(aes(ymin = pmax(0.25, Score - SD), ymax = pmin(1, Score + SD),
                    fill = Set, group = Set),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = function(x) seq_len(max(x)),
                       expand = ggplot2::expansion(mult = c(0.12, 0.03))) +
    scale_y_continuous(limits = c(0.25, 1),
                       expand = ggplot2::expansion(mult = c(0, 0))) +
    labs(x = if (show_x) "Number of attributes" else NULL,
         y = if (show_y) "MCC" else NULL,
         title = alg) +
    theme_pub() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.title.x = element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = element_text(margin = ggplot2::margin(r = 6))
    )
}


p1 <- make_alg_plot("C4.5", show_y = TRUE,  show_x = FALSE)
p2 <- make_alg_plot("k-NN", show_y = FALSE, show_x = FALSE)
p3 <- make_alg_plot("SVM",  show_y = FALSE, show_x = FALSE)
p4 <- make_alg_plot("RF",   show_y = TRUE,  show_x = TRUE)
p5 <- make_alg_plot("LR",   show_y = FALSE, show_x = TRUE)

# legend grob with no title
legend_src <- ggplot(res_all, aes(k, Score, color = Set)) +
  geom_line() +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")
pleg <- patchwork::wrap_elements(cowplot::get_legend(legend_src))

# arrange: 3 plots top row, 2 plots + legend bottom row
p_grid <- (p1 + p2 + p3) /
  (p4 + p5 + pleg) +
  plot_layout(widths = c(1,1,1), heights = c(1,1))

print(p_grid)



ggsave(filename = "20250924_parsimony.png", plot = p_grid,  width= 12, height = 6, units = "in", dpi = 300, bg = "white")
