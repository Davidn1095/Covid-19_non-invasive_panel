# ============================================================
# Feature importance figure (ΔMCC permutation) — complete script
# Bigger fonts, grey vertical lines every 0.125, black x=0 on top,
# one black horizontal line per facet, no y-axis label
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

# algorithms
algo_map <- c("C4.5" = "J48", "k-NN" = "kknn", "SVM" = "svmRadial", "RF" = "ranger", "LR" = "glmnet")
needs_dummies <- function(algo) algo %in% c("glmnet","svmRadial","kknn")
needs_scaling <- function(algo) algo %in% c("glmnet","svmRadial","kknn")

# predictors
feat_full   <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission",
                 "monocyte_abs_number","monocytes_perc",
                 "lymphocyte_abs_number","lymphocytes_perc",
                 "neutrophil_abs_number","neutrophils_perc")
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# labels
label_map <- c(
  "Diagnosis"              = "Diagnosis",
  "severity_admission"     = "Admission severity",
  "Age"                    = "Age",
  "Gender"                 = "Sex",
  "SpO2_admission"         = "Oxygen saturation",
  "monocyte_abs_number"    = "Monocyte count",
  "monocytes_perc"         = "Monocytes %",
  "lymphocyte_abs_number"  = "Lymphocyte count",
  "lymphocytes_perc"       = "Lymphocytes %",
  "neutrophil_abs_number"  = "Neutrophil count",
  "neutrophils_perc"       = "Neutrophils %"
)
f_label <- function(x) unname(label_map[x] %||% x)

# ---------------- date helpers ----------------
norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
try_parse_date <- function(x) {
  if (inherits(x,"Date"))   return(as.Date(x))
  if (inherits(x,"POSIXt")) return(as.Date(x))
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

# ---------------- modeling helpers ----------------
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

# caret recipe per algorithm
make_recipe_for <- function(dat, yvar, algo) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat)
  ign <- intersect(c("data_split",".rid"), names(dat))
  if (length(ign)) rec <- rec |> recipes::update_role(all_of(ign), new_role = "ignore")  # avoid warning
  rec <- rec |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors())  |>
    recipes::step_novel(recipes::all_nominal_predictors())        |>
    recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01)
  if (needs_dummies(algo)) rec <- rec |> recipes::step_dummy(recipes::all_nominal_predictors())
  rec <- rec |> recipes::step_zv(recipes::all_predictors())
  if (needs_scaling(algo)) rec <- rec |>
    recipes::step_YeoJohnson(recipes::all_numeric_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors())
  rec
}

# train and evaluate on holdouts
run_holdout_algo <- function(df, yvar, features, algo_key,
                             split_train = "train", split_test = c("test","external")) {
  caret_method <- unname(algo_map[algo_key])
  present_splits <- intersect(split_test, unique(as.character(df$data_split)))
  if (!length(present_splits)) return(NULL)
  base_cols <- unique(c(features, yvar, intersect("group", names(df))))
  base_cols <- intersect(base_cols, names(df))
  dtrain <- df[df$data_split == split_train, base_cols, drop = FALSE]
  if (!nrow(dtrain)) return(NULL)
  
  rec <- make_recipe_for(dtrain, yvar, caret_method)
  tr_ctrl <- caret::trainControl(method = "cv", number = cfg$inner_k,
                                 classProbs = TRUE, summaryFunction = twoClassSummary)
  
  fit <- caret::train(
    rec, data = dtrain, method = caret_method,
    trControl = tr_ctrl, metric = "ROC", tuneLength = cfg$tune_len
  )
  
  out <- list()
  for (sp in present_splits) {
    dtest <- df[df$data_split == sp, base_cols, drop = FALSE]
    if (!nrow(dtest)) next
    p_tr  <- as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label])
    neg_label <- get_neg_level(dtrain[[yvar]])
    thr   <- best_threshold_mcc(dtrain[[yvar]], p_tr, pos = cfg$pos_label, neg = neg_label)
    p_te  <- as.numeric(predict(fit, newdata = dtest, type = "prob")[, cfg$pos_label])
    pred  <- factor(ifelse(p_te >= thr$t, cfg$pos_label, neg_label),
                    levels = levels(dtrain[[yvar]]))
    out[[sp]] <- list(
      obs = factor(dtest[[yvar]], levels = levels(dtrain[[yvar]])),
      p   = p_te, pred = pred, threshold = thr$t,
      model = fit, holdout_df = dtest,
      algo  = algo_key
    )
  }
  out
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

# ---------------- plotting helpers ----------------
theme_pub <- function(base = 18) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title  = element_text(face = "bold"),
      axis.title  = element_text(color = "grey30"),
      axis.text   = element_text(color = "grey20"),
      strip.text  = element_text(face = "bold"),
      legend.position = "right",
      legend.title    = element_text(face = "bold")
    )
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  within <- interaction(within, drop = TRUE)
  x <- paste(x, within, sep = sep)
  stats::reorder(x, by, fun = fun)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(drop = TRUE,
                            labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# two panels with grey verticals under bars, black x=0 on top, black horizontal per facet, no y label
plot_perm_two_panel <- function(holds_full, holds_tri) {
  imp_full_all <- purrr::map_dfr(names(holds_full), function(ak) {
    H <- holds_full[[ak]][["test"]]
    if (is.null(H)) return(NULL)
    perm_drop_mcc(H$model, H$holdout_df, cfg$outcome,
                  intersect(feat_full, names(H$holdout_df)), H$threshold) |>
      dplyr::mutate(Algorithm = ak)
  }) |>
    dplyr::filter(is.finite(Drop_MCC)) |>
    dplyr::group_by(Algorithm) |>
    dplyr::arrange(dplyr::desc(Drop_MCC), .by_group = TRUE) |>
    dplyr::mutate(Feature_label = f_label(Feature),
                  y_key = reorder_within(Feature_label, Drop_MCC, within = Algorithm),
                  sign = ifelse(Drop_MCC >= 0, "Positive", "Negative")) |>
    dplyr::ungroup()
  
  imp_tri_all <- purrr::map_dfr(names(holds_tri), function(ak) {
    H <- holds_tri[[ak]][["test"]]
    if (is.null(H)) return(NULL)
    perm_drop_mcc(H$model, H$holdout_df, cfg$outcome,
                  intersect(feat_triage, names(H$holdout_df)), H$threshold) |>
      dplyr::mutate(Algorithm = ak)
  }) |>
    dplyr::filter(is.finite(Drop_MCC)) |>
    dplyr::group_by(Algorithm) |>
    dplyr::arrange(dplyr::desc(Drop_MCC), .by_group = TRUE) |>
    dplyr::mutate(Feature_label = f_label(Feature),
                  y_key = reorder_within(Feature_label, Drop_MCC, within = Algorithm),
                  sign = ifelse(Drop_MCC >= 0, "Positive", "Negative")) |>
    dplyr::ungroup()
  
  x_min <- min(imp_full_all$Drop_MCC, imp_tri_all$Drop_MCC, 0, na.rm = TRUE)
  x_max <- max(imp_full_all$Drop_MCC, imp_tri_all$Drop_MCC, 0, na.rm = TRUE)
  pad   <- 0.02 * (x_max - x_min + 1e-9)
  grid_breaks <- seq(floor((x_min - pad)/0.125)*0.125,
                     ceiling((x_max + pad)/0.125)*0.125,
                     by = 0.125)
  grid_breaks <- grid_breaks[grid_breaks >= (x_min - pad) & grid_breaks <= (x_max + pad)]  # avoid vline warnings
  
  grid_layer   <- geom_vline(xintercept = grid_breaks, linewidth = 0.3, color = "grey85")
  horiz_layer  <- geom_hline(yintercept = 0.5, linewidth = 0.6, color = "black")
  x0_top_layer <- geom_vline(xintercept = 0, linewidth = 0.8, color = "black")
  
  p_full <- ggplot(imp_full_all, aes(Drop_MCC, y_key, fill = sign)) +
    grid_layer + horiz_layer +
    geom_col(show.legend = FALSE) +
    x0_top_layer +
    facet_grid(Algorithm ~ ., scales = "free_y") +
    scale_y_reordered() +
    scale_x_continuous(limits = c(x_min - pad, x_max + pad), expand = c(0, 0)) +
    scale_fill_manual(values = c("Negative" = "#F4CCCC", "Positive" = "#D9EAD3")) +
    labs(title = "Complete feature set", x = "Drop in MCC when permuted", y = NULL) +
    theme_pub()
  
  p_tri <- ggplot(imp_tri_all, aes(Drop_MCC, y_key, fill = sign)) +
    grid_layer + horiz_layer +
    geom_col(show.legend = FALSE) +
    x0_top_layer +
    facet_grid(Algorithm ~ ., scales = "free_y") +
    scale_y_reordered() +
    scale_x_continuous(limits = c(x_min - pad, x_max + pad), expand = c(0, 0)) +
    scale_fill_manual(values = c("Negative" = "#F4CCCC", "Positive" = "#D9EAD3")) +
    labs(title = "Triage feature set", x = "Drop in MCC when permuted", y = NULL) +
    theme_pub()
  
  p_full + patchwork::plot_spacer() + p_tri +
    patchwork::plot_layout(widths = c(1, 0.12, 1)) &
    theme(plot.title = element_text(hjust = 0.5))
}

# ============================================================
# RUN
# ============================================================
df <- load_and_split(cfg$file_path, cfg$sheet, cfg$outcome, cfg$boundary, cfg$pct_test)

holds_full <- lapply(names(algo_map), function(ak) run_holdout_algo(df, cfg$outcome, features = feat_full,  algo_key = ak))
names(holds_full) <- names(algo_map)
holds_tri  <- lapply(names(algo_map), function(ak) run_holdout_algo(df, cfg$outcome, features = feat_triage, algo_key = ak))
names(holds_tri)  <- names(algo_map)

p_importance_combined <- plot_perm_two_panel(holds_full, holds_tri)
print(p_importance_combined)

# Optional save
ggsave("20250923_features.png", p_importance_combined, width = 12, height = 14, dpi = 300, bg = "white")
