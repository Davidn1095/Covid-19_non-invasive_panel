# ===============================================================
# PARSIMONY CURVES (MCC) — FIXED FOR LENGTH MISMATCH
# Key fix vs previous version:
#   - Align factor/character levels between df_train and df_test per k,
#     mapping unseen test levels -> "Unknown" BEFORE training/predicting.
# This prevents predict() from silently dropping rows (which caused
# "table(): all arguments must have the same length").
#
# Uses:
#   df_train, df_test, feat_full, feat_triage, cfg
# tuneLength = 5 for ALL models
# No files written (no ggsave)
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(pROC)
  library(kknn)
  library(randomForest)
  library(kernlab)
  library(RWeka)
  library(patchwork)
  library(cowplot)
  library(forcats)
  library(digest)
})

options(warn = -1)
options(na.action = "na.pass")

`%||%` <- function(a, b) if (!is.null(a)) a else b
seed_from_name <- function(name, offset = 0L) as.integer(offset + sum(utf8ToInt(as.character(name))))

# -----------------------------
# 0) Guards + config
# -----------------------------
stopifnot(exists("df_train", inherits = TRUE))
stopifnot(exists("df_test",  inherits = TRUE))
stopifnot(exists("feat_full",   inherits = TRUE))
stopifnot(exists("feat_triage", inherits = TRUE))
stopifnot(exists("cfg", inherits = TRUE))

df_train <- get("df_train", inherits = TRUE)
df_test  <- get("df_test",  inherits = TRUE)
feat_full   <- get("feat_full",   inherits = TRUE)
feat_triage <- get("feat_triage", inherits = TRUE)
cfg <- get("cfg", inherits = TRUE)

yvar <- "outcome"
stopifnot(yvar %in% names(df_train), yvar %in% names(df_test))

df_train[[yvar]] <- factor(df_train[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
df_test[[yvar]]  <- factor(df_test[[yvar]],  levels = c(cfg$neg_label, cfg$pos_label))

use_class_weights <- if (exists("use_class_weights", inherits = TRUE)) {
  isTRUE(get("use_class_weights", inherits = TRUE))
} else TRUE

max_n_train   <- 5000L
tune_len_all  <- 5L
algo_order    <- c("C4.5","k-NN","SVM","RF","LR")
cascade_label <- "Cascade model"

t_low_fixed  <- if (exists("t_low_fixed",  inherits = TRUE)) get("t_low_fixed",  inherits = TRUE) else (cfg$t_low_default %||% 0.33)
t_high_fixed <- if (exists("t_high_fixed", inherits = TRUE)) get("t_high_fixed", inherits = TRUE) else (cfg$t_high_default %||% 0.66)

# speed knobs (bootstrap SD)
speed_mode  <- TRUE
speed_ultra <- FALSE
B_parsimony <- if (speed_ultra) 0 else if (speed_mode) 200 else 1000

RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed_cv %||% 123)

# -----------------------------
# 1) Utilities: preprocessing + alignment (FIX), weights, calibration, threshold, metrics
# -----------------------------
prep_df_for_algo <- function(d, feats, yvar) {
  feats <- intersect(feats, names(d))
  d2 <- d[!is.na(d[[yvar]]), , drop = FALSE]
  d2[[yvar]] <- factor(d2[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  
  if (!length(feats)) return(list(df = d2, feats = feats))
  
  # characters -> keep as character for now; we will align later
  is_char <- vapply(d2[feats], is.character, logical(1))
  if (any(is_char)) {
    ch <- feats[is_char]
    d2[ch] <- lapply(d2[ch], function(x) replace(x, is.na(x), "Unknown"))
  }
  
  # factors -> add Unknown for NA, but alignment later will enforce training levels
  is_fac <- vapply(d2[feats], is.factor, logical(1))
  if (any(is_fac)) {
    fc <- feats[is_fac]
    d2[fc] <- lapply(d2[fc], function(x) forcats::fct_na_value_to_level(x, "Unknown"))
  }
  
  # numeric median impute
  is_num <- vapply(d2[feats], is.numeric, logical(1))
  for (cl in feats[is_num]) {
    if (all(is.na(d2[[cl]]))) d2[[cl]] <- 0 else {
      med <- stats::median(d2[[cl]], na.rm = TRUE)
      d2[[cl]][is.na(d2[[cl]])] <- med
    }
  }
  
  # drop constant predictors in THIS subset
  zv <- feats[vapply(d2[feats], function(x) {
    ux <- unique(x[!is.na(x)])
    length(ux) <= 1
  }, logical(1))]
  feats <- setdiff(feats, zv)
  
  list(df = d2, feats = feats)
}

# FIX: enforce consistent factor levels between train/test and map unseen test levels -> "Unknown"
align_train_test_levels <- function(dtr, dte, feats) {
  for (f in feats) {
    if (!(f %in% names(dtr)) || !(f %in% names(dte))) next
    
    is_cat <- is.factor(dtr[[f]]) || is.factor(dte[[f]]) || is.character(dtr[[f]]) || is.character(dte[[f]])
    if (!is_cat) next
    
    tr <- as.character(dtr[[f]])
    te <- as.character(dte[[f]])
    
    tr[is.na(tr)] <- "Unknown"
    te[is.na(te)] <- "Unknown"
    
    # training levels define the vocabulary (+ Unknown)
    lev_tr <- unique(tr)
    if (!("Unknown" %in% lev_tr)) lev_tr <- c(lev_tr, "Unknown")
    
    # map unseen test levels to Unknown
    te[!(te %in% lev_tr)] <- "Unknown"
    
    dtr[[f]] <- factor(tr, levels = lev_tr)
    dte[[f]] <- factor(te, levels = lev_tr)
  }
  list(dtr = dtr, dte = dte)
}

add_lr_shadow_if_needed <- function(df, feats, yvar, algo, seed) {
  if (identical(algo, "LR") && length(feats) == 1) {
    set.seed(seed)
    nm <- "shadow_LR"
    while (nm %in% names(df)) nm <- paste0("shadow_LR_", sample(1000:9999, 1))
    df[[nm]] <- rnorm(nrow(df), sd = 1e-8)
    list(df = df, feats = c(feats, nm))
  } else {
    list(df = df, feats = feats)
  }
}

compute_case_weights <- function(y) {
  if (!isTRUE(use_class_weights)) return(rep(1, length(y)))
  tab <- table(y)
  n_pos <- as.numeric(tab[cfg$pos_label])
  n_neg <- as.numeric(tab[cfg$neg_label])
  w_pos <- if (!is.na(n_pos) && n_pos > 0) n_neg / n_pos else 1
  ifelse(y == cfg$pos_label, w_pos, 1)
}

fit_calibrator <- function(obs, prob, pos_label) {
  y <- ifelse(obs == pos_label, 1, 0)
  p <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  suppressWarnings(stats::glm(y ~ qlogis(p), data = data.frame(y = y, p = p),
                              family = stats::binomial(link = "logit")))
}

apply_calibration <- function(prob, calib_model) {
  if (is.null(calib_model)) return(prob)
  p <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  lp <- try(predict(calib_model, newdata = data.frame(p = p), type = "link"), silent = TRUE)
  if (inherits(lp, "try-error")) return(p)
  out <- 1 / (1 + exp(-lp))
  out[!is.finite(out)] <- p[!is.finite(out)]
  out
}

best_thresh_roc <- function(obs, prob, pos = cfg$pos_label, neg = cfg$neg_label) {
  mask <- !is.na(prob)
  if (sum(mask) < 2L) return(0.5)
  roc_obj <- try(pROC::roc(response = obs[mask], predictor = as.numeric(prob[mask]),
                           levels = c(neg, pos), direction = "<"), silent = TRUE)
  if (inherits(roc_obj, "try-error")) return(0.5)
  thr <- try(pROC::coords(roc_obj, "best", best.method = "youden", ret = "threshold"), silent = TRUE)
  if (inherits(thr, "try-error") || length(thr) == 0L || is.na(thr[1])) return(0.5)
  as.numeric(thr[1])
}

mcc_from_tab <- function(obs, pred, pos, neg) {
  obs  <- factor(obs,  levels = c(neg, pos))
  pred <- factor(pred, levels = c(neg, pos))
  tab <- table(obs, pred)
  TN <- as.numeric(tab[neg, neg])
  FP <- as.numeric(tab[neg, pos])
  FN <- as.numeric(tab[pos, neg])
  TP <- as.numeric(tab[pos, pos])
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  ((TP * TN) - (FP * FN)) / denom
}

compute_metrics_binary <- function(obs, pred, p_pos, pos = cfg$pos_label, neg = cfg$neg_label) {
  # extra safety: enforce same length
  n <- min(length(obs), length(pred), length(p_pos))
  obs <- obs[seq_len(n)]
  pred <- pred[seq_len(n)]
  p_pos <- p_pos[seq_len(n)]
  
  mcc <- mcc_from_tab(obs, pred, pos, neg)
  p_pos <- as.numeric(p_pos); p_pos[is.na(p_pos)] <- 0.5
  auc <- NA_real_
  roc_obj <- try(pROC::roc(response = factor(obs, levels = c(neg, pos)),
                           predictor = p_pos, levels = c(neg, pos), direction = "<"),
                 silent = TRUE)
  if (!inherits(roc_obj, "try-error")) auc <- as.numeric(pROC::auc(roc_obj))
  c(MCC = mcc, AUC = auc)
}

boot_sd_mcc <- function(obs, pred, p, B = B_parsimony, seed = 123) {
  if (!is.numeric(B) || B <= 0) return(NA_real_)
  if (length(obs) < 2 || length(unique(obs)) < 2) return(NA_real_)
  n <- length(obs); set.seed(seed)
  tryCatch({
    sd(replicate(B, {
      idx <- sample.int(n, n, replace = TRUE)
      as.numeric(compute_metrics_binary(obs[idx], pred[idx], p[idx])["MCC"])
    }), na.rm = TRUE)
  }, error = function(e) NA_real_)
}

aucSummary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) lev <- levels(data$obs)
  if (length(lev) != 2L) stop("aucSummary requires a binary outcome.")
  pos <- lev[2]
  roc_obj <- try(pROC::roc(response = data$obs, predictor = data[[pos]],
                           levels = lev, direction = "<"), silent = TRUE)
  roc_auc <- if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj))
  c(ROC = roc_auc)
}

ctrl_auc <- caret::trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = aucSummary,
  savePredictions = "final",
  verboseIter     = FALSE,
  allowParallel   = TRUE
)

method_from_algo <- function(algo) {
  switch(
    algo,
    "LR"   = "glmnet",
    "RF"   = "rf",
    "SVM"  = "svmRadial",
    "k-NN" = "kknn",
    "C4.5" = "J48",
    stop("Unknown algorithm: ", algo)
  )
}

extract_oof_best <- function(fit) {
  oof <- fit$pred
  if (is.null(oof) || !nrow(oof)) stop("Out-of-fold predictions were not available from caret.")
  bt <- fit$bestTune
  if (!is.null(bt) && nrow(bt)) {
    for (nm in names(bt)) oof <- oof[oof[[nm]] == bt[[nm]], , drop = FALSE]
  }
  oof
}

# -----------------------------
# 2) Holdout training: fit on df_train, calibrate on OOF, threshold on dev, evaluate on df_test
# -----------------------------
run_holdout <- function(df_train, df_test, yvar, features, algo, set_label) {
  pp_tr <- prep_df_for_algo(df_train, features, yvar)
  pp_te <- prep_df_for_algo(df_test,  features, yvar)
  
  feats <- intersect(pp_tr$feats, pp_te$feats)
  if (!length(feats)) return(NULL)
  
  dtr <- pp_tr$df[, c(yvar, feats), drop = FALSE]
  dte <- pp_te$df[, c(yvar, feats), drop = FALSE]
  
  # FIX: align categorical levels (prevents predict() from dropping rows)
  al <- align_train_test_levels(dtr, dte, feats)
  dtr <- al$dtr
  dte <- al$dte
  
  if (nrow(dtr) > max_n_train) {
    set.seed(seed_from_name(paste(algo, set_label, "subsample", sep = "|")))
    dtr <- dplyr::slice_sample(dtr, n = max_n_train)
  }
  
  # LR single-feature guard
  sh_seed <- seed_from_name(paste(algo, set_label, "shadow", sep = "|"))
  sh_tr <- add_lr_shadow_if_needed(dtr, feats, yvar, algo, sh_seed)
  sh_te <- add_lr_shadow_if_needed(dte, feats, yvar, algo, sh_seed)
  dtr <- sh_tr$df; dte <- sh_te$df; feats2 <- sh_tr$feats
  
  w <- compute_case_weights(dtr[[yvar]])
  form <- stats::as.formula(paste(yvar, "~ ."))
  method <- method_from_algo(algo)
  
  set.seed(seed_from_name(paste(algo, set_label, "fit", sep = "|")))
  fit <- switch(
    algo,
    "LR" = caret::train(form, data = dtr, method = method,
                        metric = "ROC", trControl = ctrl_auc,
                        tuneLength = tune_len_all,
                        weights = w, family = "binomial"),
    "RF" = caret::train(form, data = dtr, method = method,
                        metric = "ROC", trControl = ctrl_auc,
                        tuneLength = tune_len_all,
                        weights = w, ntree = 500),
    caret::train(form, data = dtr, method = method,
                 metric = "ROC", trControl = ctrl_auc,
                 tuneLength = tune_len_all,
                 weights = w)
  )
  
  # OOF -> calibrate
  oof <- extract_oof_best(fit)
  pos_col <- cfg$pos_label
  if (!(pos_col %in% names(oof))) stop("Positive class column not found in OOF predictions.")
  
  p_oof_raw <- as.numeric(oof[[pos_col]]); p_oof_raw[is.na(p_oof_raw)] <- 0.5
  calib <- fit_calibrator(oof$obs, p_oof_raw, cfg$pos_label)
  p_oof_cal <- apply_calibration(p_oof_raw, calib)
  thr_dev <- best_thresh_roc(oof$obs, p_oof_cal)
  
  # external calibrated probabilities
  p_te_mat <- predict(fit, newdata = dte, type = "prob", na.action = na.pass)
  p_te_raw <- if (is.vector(p_te_mat)) as.numeric(p_te_mat) else as.numeric(p_te_mat[, cfg$pos_label])
  
  # hard check: prediction must match test rows (should hold after alignment fix)
  if (length(p_te_raw) != nrow(dte)) {
    stop(
      "Prediction length mismatch for algo=", algo, " set=", set_label,
      " | nrow(test)=", nrow(dte), " but length(prob)=", length(p_te_raw),
      ". This indicates rows were dropped during prediction."
    )
  }
  
  p_te_raw[is.na(p_te_raw)] <- 0.5
  p_te_cal <- apply_calibration(p_te_raw, calib)
  
  pred_te <- ifelse(p_te_cal >= thr_dev, cfg$pos_label, cfg$neg_label)
  pred_te <- factor(pred_te, levels = c(cfg$neg_label, cfg$pos_label))
  
  list(
    algo      = algo,
    set       = set_label,
    fit       = fit,
    threshold = thr_dev,
    oof       = data.frame(obs = oof$obs, p = p_oof_cal),
    obs       = dte[[yvar]],
    p         = p_te_cal,
    pred      = pred_te,
    feats     = feats2
  )
}

# memoize repeated k runs
.hold_cache <- new.env(parent = emptyenv())
run_holdout_memo <- function(df_train, df_test, yvar, features, algo, set_label) {
  key <- paste(algo, set_label, digest::digest(sort(features)), sep = "|")
  if (exists(key, envir = .hold_cache, inherits = FALSE)) return(get(key, envir = .hold_cache))
  res <- try(run_holdout(df_train, df_test, yvar, features, algo, set_label), silent = TRUE)
  if (inherits(res, "try-error")) res <- NULL
  assign(key, res, envir = .hold_cache)
  res
}

# -----------------------------
# 3) Cascade apply (fixed band)
# -----------------------------
apply_cascade_external_fixed <- function(tri_hold, com_hold, t_low = t_low_fixed, t_high = t_high_fixed) {
  stopifnot(length(tri_hold$obs) == length(com_hold$obs))
  
  y <- tri_hold$obs
  p_tri  <- as.numeric(tri_hold$p);  p_tri[is.na(p_tri)]   <- 0.5
  p_full <- as.numeric(com_hold$p);  p_full[is.na(p_full)] <- 0.5
  
  thr_full <- if (is.finite(com_hold$threshold)) com_hold$threshold else 0.5
  
  decision <- ifelse(
    p_tri < t_low, cfg$neg_label,
    ifelse(p_tri > t_high, cfg$pos_label, "Defer")
  )
  
  pred <- character(length(decision))
  nd <- decision != "Defer"
  pred[nd] <- decision[nd]
  pred[!nd] <- ifelse(p_full[!nd] >= thr_full, cfg$pos_label, cfg$neg_label)
  pred <- factor(pred, levels = c(cfg$neg_label, cfg$pos_label))
  
  list(obs = y, pred = pred, p = ifelse(nd, p_tri, p_full), defer_rate = mean(!nd))
}

# -----------------------------
# 4) Order features by dev importance (glmnet varImp), tuneLength=5
# -----------------------------
ref_feats <- intersect(feat_full, names(df_train))
stopifnot(length(ref_feats) > 0)

df_rank <- df_train[, unique(c(yvar, ref_feats)), drop = FALSE]

# light impute
num_cols <- intersect(ref_feats, names(df_rank)[vapply(df_rank, is.numeric, logical(1))])
for (cl in num_cols) {
  med <- stats::median(df_rank[[cl]], na.rm = TRUE)
  if (is.finite(med)) df_rank[[cl]][is.na(df_rank[[cl]])] <- med
}
fac_cols <- setdiff(ref_feats, num_cols)
if (length(fac_cols)) {
  df_rank[fac_cols] <- lapply(df_rank[fac_cols], function(x) {
    x <- as.factor(x)
    forcats::fct_na_value_to_level(x, "Unknown")
  })
}
df_rank[[yvar]] <- factor(df_rank[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))

form_ref <- stats::as.formula(paste(yvar, "~", paste(ref_feats, collapse = "+")))
ctrl_ref <- caret::trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = aucSummary,
  savePredictions = "final",
  allowParallel   = TRUE
)

set.seed(cfg$seed_cv %||% 123)
fit_ref <- caret::train(
  form_ref,
  data       = df_rank,
  method     = "glmnet",
  family     = "binomial",
  tuneLength = tune_len_all,
  metric     = "ROC",
  trControl  = ctrl_ref
)

imp_tbl <- caret::varImp(fit_ref)$importance
imp_tbl$Var <- rownames(imp_tbl)

map_to_base <- function(vn, base_feats) {
  hits <- base_feats[startsWith(vn, base_feats)]
  if (!length(hits)) return(NA_character_)
  hits[which.max(nchar(hits))]
}
imp_tbl$Base <- vapply(imp_tbl$Var, map_to_base, character(1), base_feats = ref_feats)

imp_agg <- imp_tbl %>%
  dplyr::filter(!is.na(Base)) %>%
  dplyr::group_by(Base) %>%
  dplyr::summarise(Overall = max(Overall, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Overall))

ref_ord <- unique(c(imp_agg$Base, ref_feats))

# triage-first ordering
tri_order <- intersect(ref_ord, feat_triage)
triage_k_cap <- 4
tri_order <- head(tri_order, triage_k_cap)

rest_pool  <- setdiff(intersect(feat_full, names(df_train)), tri_order)
rest_order <- unique(c(intersect(ref_ord, rest_pool),
                       sort(setdiff(rest_pool, ref_ord), decreasing = TRUE)))

order_by_set <- list(
  "Complete feature set" = c(tri_order, rest_order),
  "Triage feature set"   = tri_order
)
full_order <- c(tri_order, rest_order)

# -----------------------------
# 5) Parsimony curves
# -----------------------------
parsimony_external_curve <- function(set_name, feats_order, algo) {
  if (!length(feats_order)) {
    return(tibble(k = integer(), Score = numeric(), SD = numeric(), Set = character(), Algorithm = character()))
  }
  
  purrr::map_dfr(seq_along(feats_order), function(k) {
    feats_k <- feats_order[seq_len(k)]
    lab <- ifelse(set_name == "Complete feature set", "Complete", "Triage")
    
    hold <- run_holdout_memo(df_train, df_test, yvar, feats_k, algo, lab)
    
    if (is.null(hold)) {
      tibble(k = k, Score = NA_real_, SD = NA_real_, Set = set_name, Algorithm = algo)
    } else {
      mets <- compute_metrics_binary(hold$obs, hold$pred, hold$p)
      sd_mcc <- boot_sd_mcc(
        hold$obs, hold$pred, hold$p, B = B_parsimony,
        seed = seed_from_name(paste(algo, set_name, k, "boot", sep = "|"))
      )
      tibble(k = k, Score = as.numeric(mets["MCC"]), SD = sd_mcc, Set = set_name, Algorithm = algo)
    }
  })
}

cascade_external_curve <- function(tri_order, full_order, algo) {
  K_full <- length(full_order)
  if (!K_full) {
    return(tibble(k = integer(), Score = numeric(), SD = numeric(), Set = character(), Algorithm = character()))
  }
  
  purrr::map_dfr(seq_len(K_full), function(k) {
    tri_feats_k <- head(tri_order, min(k, length(tri_order)))
    com_feats_k <- head(full_order, k)
    
    tri <- run_holdout_memo(df_train, df_test, yvar, tri_feats_k, algo, "Triage")
    com <- run_holdout_memo(df_train, df_test, yvar, com_feats_k, algo, "Complete")
    
    if (is.null(tri) || is.null(com)) {
      return(tibble(k = k, Score = NA_real_, SD = NA_real_, Set = cascade_label, Algorithm = algo))
    }
    
    cas <- apply_cascade_external_fixed(tri, com, t_low = t_low_fixed, t_high = t_high_fixed)
    mets <- compute_metrics_binary(cas$obs, cas$pred, cas$p)
    sd_mcc <- boot_sd_mcc(
      cas$obs, cas$pred, cas$p, B = B_parsimony,
      seed = seed_from_name(paste(algo, "Cascade", k, "boot", sep = "|"))
    )
    
    tibble(k = k, Score = as.numeric(mets["MCC"]), SD = sd_mcc, Set = cascade_label, Algorithm = algo)
  })
}

res_core <- purrr::imap_dfr(order_by_set, function(ord, set_name) {
  purrr::map_dfr(algo_order, function(algo) parsimony_external_curve(set_name, ord, algo))
})
res_cas <- purrr::map_dfr(algo_order, function(algo) cascade_external_curve(tri_order, full_order, algo))

res_all <- dplyr::bind_rows(res_core, res_cas) %>%
  dplyr::mutate(
    Algorithm = factor(Algorithm, levels = algo_order),
    Set = factor(Set, levels = c("Complete feature set", "Triage feature set", cascade_label))
  )

# Optional sanity count
chk <- dplyr::count(res_all, Set, Algorithm)
print(chk %>% tidyr::pivot_wider(names_from = Set, values_from = n))

# -----------------------------
# 6) Plot
# -----------------------------
theme_pub <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(color = "grey30"),
      legend.position = "right",
      legend.title = ggplot2::element_blank()
    )
}

levels_set <- c("Complete feature set","Triage feature set", cascade_label)
pal <- c(
  "Complete feature set" = "#D55E00",
  "Triage feature set"   = "#0072B2",
  "Cascade model"        = "#009E73"
)

make_alg_plot <- function(alg, show_y = FALSE, show_x = FALSE) {
  dat <- dplyr::filter(res_all, Algorithm == alg)
  if (!nrow(dat) || all(is.na(dat$Score))) return(ggplot() + theme_void() + ggtitle(alg))
  
  tri_name <- "Triage feature set"
  dat_tri  <- dplyr::filter(dat, Set == tri_name)
  dat_oth  <- dplyr::filter(dat, Set != tri_name)
  x_max <- max(dat$k, na.rm = TRUE) + 0.5
  
  p <- ggplot(dat, aes(k, Score, color = Set, group = Set)) +
    annotate("segment", x = 0, xend = x_max, y = 0, yend = 0, linewidth = 0.6, colour = "black") +
    annotate("segment", x = 0, xend = 0, y = 0, yend = 1, linewidth = 0.6, colour = "black")
  
  if (!all(is.na(dat$SD))) {
    p <- p +
      geom_ribbon(data = dat_oth,
                  aes(ymin = pmax(0, Score - SD), ymax = pmin(1, Score + SD), fill = Set),
                  alpha = 0.15, colour = NA, show.legend = FALSE) +
      geom_ribbon(data = dat_tri,
                  aes(ymin = pmax(0, Score - SD), ymax = pmin(1, Score + SD), fill = Set),
                  alpha = 0.15, colour = NA, show.legend = FALSE) +
      scale_fill_manual(values = pal, breaks = levels_set)
  }
  
  p <- p +
    geom_line(data = dat_oth, linewidth = 1) +
    geom_point(data = dat_oth, size = 1.8, na.rm = TRUE) +
    geom_line(data = dat_tri, linewidth = 1.2) +
    geom_point(data = dat_tri, size = 2.1, na.rm = TRUE)
  
  p +
    scale_x_continuous(
      breaks = function(x) {
        hi <- floor(max(x, na.rm = TRUE) / 2) * 2
        seq(0, hi, by = 2)
      },
      expand = expansion(mult = c(0.12, 0.03))
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    scale_color_manual(values = pal, breaks = levels_set, limits = levels_set) +
    labs(x = if (show_x) "Number of attributes" else NULL,
         y = if (show_y) "MCC" else NULL,
         title = as.character(alg)) +
    theme_pub() +
    theme(
      legend.position = "none",
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

legend_src <- ggplot(res_all, aes(k, Score, color = Set)) +
  geom_line() +
  scale_color_manual(values = pal, breaks = levels_set, limits = levels_set) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")
pleg <- patchwork::wrap_elements(cowplot::get_legend(legend_src))

p_parsimony <- (p1 + p2 + p3) / (p4 + p5 + pleg) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1))

print(p_parsimony)

# Data behind curves
res_all
