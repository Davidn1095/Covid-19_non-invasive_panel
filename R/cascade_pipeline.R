suppressWarnings(options(warn = -1))
options(na.action = "na.pass")

`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(dplyr)
  library(pROC)
  library(randomForest)
  library(kernlab)
  library(kknn)
  library(rpart)
})

source(file.path("R", "common_utils.R"))
source(file.path("R", "utils_metrics.R"))
source(file.path("R", "utils_thresholds.R"))

prep_outcome_01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  
  if (is.numeric(x) || is.integer(x)) {
    out <- rep(NA_integer_, length(x))
    out[x == 1] <- 1L
    out[x == 0] <- 0L
    return(out)
  }
  
  x0 <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(x0))
  
  out[x0 %in% c("yes","y","1","true","t","positive","pos")] <- 1L
  out[x0 %in% c("no","n","0","false","f","negative","neg")] <- 0L
  out
}

harmonize_units <- function(df) {
  out <- df
  
  if ("temp_raw" %in% names(out)) {
    x <- coerce_numeric_robust(out$temp_raw)
    medx <- suppressWarnings(median(x, na.rm = TRUE))
    if (is.finite(medx) && medx > 60) x <- (x - 32) * (5/9)
    out$temp_raw <- x
  }
  
  if ("spo2_raw" %in% names(out)) {
    x <- coerce_numeric_robust(out$spo2_raw)
    medx <- suppressWarnings(median(x, na.rm = TRUE))
    if (is.finite(medx) && medx <= 1.2) x <- x * 100
    out$spo2_raw <- x
  }
  
  if ("sex_raw" %in% names(out)) {
    if (!(is.numeric(out$sex_raw) || is.integer(out$sex_raw))) {
      out$sex_raw <- coerce_sex_to_01(out$sex_raw)
    } else {
      out$sex_raw <- as.numeric(out$sex_raw)
    }
  }
  
  out
}

safe_div <- function(a, b, eps = 1e-6) {
  a <- coerce_numeric_robust(a)
  b <- coerce_numeric_robust(b)
  ok <- is.finite(a) & !is.na(a) & is.finite(b) & !is.na(b) & abs(b) > eps
  out <- rep(NA_real_, length(a))
  out[ok] <- a[ok] / b[ok]
  out
}

add_engineered_features <- function(df) {
  out <- df
  
  needed <- c(
    "neutrophils_raw","lymphocytes_raw","monocytes_raw","platelets_raw",
    "heart_rate_raw","sbp_raw","dbp_raw"
  )
  for (nm in needed) {
    if (!nm %in% names(out)) out[[nm]] <- NA_real_
  }
  
  neut  <- coerce_numeric_robust(out$neutrophils_raw)
  lymph <- coerce_numeric_robust(out$lymphocytes_raw)
  mono  <- coerce_numeric_robust(out$monocytes_raw)
  plat  <- coerce_numeric_robust(out$platelets_raw)
  hr    <- coerce_numeric_robust(out$heart_rate_raw)
  sbp   <- coerce_numeric_robust(out$sbp_raw)
  dbp   <- coerce_numeric_robust(out$dbp_raw)
  
  if (!"nlr" %in% names(out))  out$nlr  <- safe_div(neut, lymph)
  if (!"mlr" %in% names(out))  out$mlr  <- safe_div(mono, lymph)
  if (!"plr" %in% names(out))  out$plr  <- safe_div(plat, lymph)
  if (!"siri" %in% names(out)) out$siri <- safe_div(neut * mono, lymph)
  if (!"sii" %in% names(out))  out$sii  <- safe_div(plat * neut, lymph)
  
  if (!"shock_index" %in% names(out))    out$shock_index <- safe_div(hr, sbp)
  if (!"pulse_pressure" %in% names(out)) out$pulse_pressure <- sbp - dbp
  
  for (nm in c("nlr","mlr","plr","siri","sii","shock_index","pulse_pressure")) {
    if (nm %in% names(out)) out[[nm]] <- coerce_numeric_robust(out[[nm]])
  }
  
  out
}

prep_df_any <- function(df, label) {
  if (!("outcome" %in% names(df))) stop("Missing column outcome in ", label)
  
  names(df) <- trimws(names(df))
  out <- harmonize_units(df)
  
  out <- add_engineered_features(out)
  
  out$outcome <- prep_outcome_01(out$outcome)
  if (any(is.na(out$outcome))) stop("Outcome has NA after coercion in ", label)
  
  raw_cols <- setdiff(grep("_raw$", names(out), value = TRUE), "outcome")
  for (nm in raw_cols) {
    if (!(is.numeric(out[[nm]]) || is.integer(out[[nm]]))) {
      out[[nm]] <- coerce_numeric_robust(out[[nm]])
    } else {
      out[[nm]] <- as.numeric(out[[nm]])
    }
  }
  
  out
}

panel_noninv_base <- c(
  "age_raw","sex_raw","spo2_raw","heart_rate_raw","resp_rate_raw","temp_raw","sbp_raw","dbp_raw","map_raw"
)

panel_noninv_optional <- c(
  "triage_acuity_raw","triage_acuity","acuity_raw","acuity",
  "pain_raw","pain"
)

panel_labs <- c(
  "wbc_raw","neutrophils_raw","lymphocytes_raw","monocytes_raw","platelets_raw",
  "hemoglobin_raw",  "creatinine_raw","urea_raw",
  "sodium_raw","potassium_raw","chloride_raw","anion_gap_raw",
  "glucose_raw","crp_raw","d_dimer_raw",
  "ast_raw","alt_raw","lactate_raw","troponin_raw"
)

panel_noninv_engineered <- c(
  "shock_index","pulse_pressure"
)

panel_labs_engineered <- c(
  "nlr","mlr","plr","siri","sii"
)

add_missing_cols_as_na <- function(df, cols) {
  cols <- unique(cols)
  miss <- setdiff(cols, names(df))
  if (length(miss)) for (m in miss) df[[m]] <- NA_real_
  df
}

pick_predictors_with_signal <- function(df, cols) {
  cols <- unique(cols)
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(character(0))
  
  keep <- character(0)
  for (nm in cols) {
    x <- coerce_numeric_robust(df[[nm]])
    ok <- is.finite(x) & !is.na(x)
    if (!any(ok)) next
    if (stats::sd(x[ok]) <= 0) next
    keep <- c(keep, nm)
  }
  unique(keep)
}

get_panel_predictors <- function(df, panel) {
  if (panel == "Non-invasive") {
    cols <- c(panel_noninv_base, panel_noninv_optional, panel_noninv_engineered)
  } else if (panel == "Laboratory augmented") {
    cols <- c(panel_noninv_base, panel_noninv_optional, panel_noninv_engineered, panel_labs, panel_labs_engineered)
  } else {
    stop("Unknown panel: ", panel)
  }
  df2 <- add_missing_cols_as_na(df, cols)
  pick_predictors_with_signal(df2, cols)
}

row_complete_frac_global <- function(df) {
  cols <- unique(c(panel_noninv_base, panel_noninv_optional, panel_labs))
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(rep(0, nrow(df)))
  
  ok_mat <- vapply(cols, function(nm) {
    x <- coerce_numeric_robust(df[[nm]])
    is.finite(x) & !is.na(x)
  }, logical(nrow(df)))
  
  if (is.vector(ok_mat)) ok_mat <- matrix(ok_mat, nrow = nrow(df))
  rowMeans(ok_mat)
}

stratified_sample_idx <- function(y01, n_max, seed) {
  set.seed(seed)
  n <- length(y01)
  if (n <= n_max) return(seq_len(n))
  
  idx_pos <- which(y01 == 1L)
  idx_neg <- which(y01 == 0L)
  
  n_pos <- length(idx_pos)
  n_neg <- length(idx_neg)
  
  prop_pos <- n_pos / (n_pos + n_neg)
  take_pos <- max(1L, min(n_pos, as.integer(round(n_max * prop_pos))))
  take_neg <- max(1L, min(n_neg, n_max - take_pos))
  
  samp_pos <- sample(idx_pos, take_pos, replace = FALSE)
  samp_neg <- sample(idx_neg, take_neg, replace = FALSE)
  sort(c(samp_pos, samp_neg))
}

apply_global_filter_then_sample <- function(df, row_thr, n_max, seed) {
  rc <- row_complete_frac_global(df)
  keep <- which(is.finite(rc) & !is.na(rc) & rc >= row_thr & !is.na(df$outcome))
  df2 <- df[keep, , drop = FALSE]
  idx <- stratified_sample_idx(df2$outcome, n_max, seed)
  df2[idx, , drop = FALSE]
}

make_repeated_stratified_folds <- function(y01, k, repeats, seed) {
  set.seed(seed)
  y01 <- as.integer(y01)
  n <- length(y01)
  
  out <- list()
  
  for (r in seq_len(repeats)) {
    pos <- sample(which(y01 == 1L))
    neg <- sample(which(y01 == 0L))
    
    pos_folds <- split(pos, rep(seq_len(k), length.out = length(pos)))
    neg_folds <- split(neg, rep(seq_len(k), length.out = length(neg)))
    
    for (j in seq_len(k)) {
      te <- sort(c(pos_folds[[j]] %||% integer(0), neg_folds[[j]] %||% integer(0)))
      tr <- setdiff(seq_len(n), te)
      out[[paste0("Rep", r, "_Fold", j)]] <- list(train = tr, test = te)
    }
  }
  out
}

imputer_fit <- function(X) {
  X <- as.data.frame(X)
  med <- vapply(names(X), function(nm) {
    x <- coerce_numeric_robust(X[[nm]])
    m <- suppressWarnings(median(x, na.rm = TRUE))
    if (!is.finite(m) || is.na(m)) 0.0 else as.numeric(m)
  }, numeric(1))
  list(median = med)
}

imputer_apply <- function(X, imp) {
  out <- as.data.frame(X)
  for (nm in names(out)) {
    x <- coerce_numeric_robust(out[[nm]])
    x[!is.finite(x) | is.na(x)] <- imp$median[[nm]]
    out[[nm]] <- x
  }
  out
}

winsor_fit <- function(X, clip_q) {
  if (is.na(clip_q) || !is.finite(clip_q) || clip_q <= 0) {
    return(list(clip_q = NA_real_, lo = NULL, hi = NULL))
  }
  X <- as.data.frame(X)
  lo <- vapply(names(X), function(nm) {
    x <- X[[nm]]
    q <- suppressWarnings(stats::quantile(x, probs = clip_q, na.rm = TRUE, type = 7))
    if (!is.finite(q) || is.na(q)) -Inf else as.numeric(q)
  }, numeric(1))
  hi <- vapply(names(X), function(nm) {
    x <- X[[nm]]
    q <- suppressWarnings(stats::quantile(x, probs = 1 - clip_q, na.rm = TRUE, type = 7))
    if (!is.finite(q) || is.na(q)) Inf else as.numeric(q)
  }, numeric(1))
  list(clip_q = clip_q, lo = lo, hi = hi)
}

winsor_apply <- function(X, wf) {
  out <- as.data.frame(X)
  if (is.na(wf$clip_q)) return(out)
  for (nm in names(out)) {
    x <- out[[nm]]
    x <- pmin(pmax(x, wf$lo[[nm]]), wf$hi[[nm]])
    out[[nm]] <- x
  }
  out
}

scale_fit <- function(X) {
  X <- as.matrix(X)
  mu <- colMeans(X)
  sd <- apply(X, 2, stats::sd)
  sd[!is.finite(sd) | is.na(sd) | sd == 0] <- 1.0
  list(mu = mu, sd = sd)
}

scale_apply <- function(X, sf) {
  X <- as.matrix(X)
  out <- sweep(sweep(X, 2, sf$mu, "-"), 2, sf$sd, "/")
  out[!is.finite(out) | is.na(out)] <- 0.0
  as.data.frame(out)
}

prep_fit <- function(X_tr, clip_q, standardise) {
  X_tr <- as.data.frame(X_tr)
  imp <- imputer_fit(X_tr)
  X1  <- imputer_apply(X_tr, imp)
  wf  <- winsor_fit(X1, clip_q)
  X2  <- winsor_apply(X1, wf)
  sf  <- if (isTRUE(standardise)) scale_fit(X2) else NULL
  list(imp = imp, wf = wf, sf = sf, standardise = isTRUE(standardise))
}

prep_apply <- function(X, pp) {
  X <- as.data.frame(X)
  X1 <- imputer_apply(X, pp$imp)
  X2 <- winsor_apply(X1, pp$wf)
  if (isTRUE(pp$standardise)) {
    X3 <- scale_apply(X2, pp$sf)
  } else {
    X3 <- X2
    for (nm in names(X3)) {
      x <- X3[[nm]]
      x[!is.finite(x) | is.na(x)] <- 0.0
      X3[[nm]] <- x
    }
  }
  X3
}

make_case_weights <- function(y01) {
  y01 <- as.integer(y01)
  n_pos <- sum(y01 == 1L, na.rm = TRUE)
  n_neg <- sum(y01 == 0L, na.rm = TRUE)
  n <- length(y01)
  if (n_pos == 0L || n_neg == 0L) return(rep(1.0, n))
  w_pos <- n / (2 * n_pos)
  w_neg <- n / (2 * n_neg)
  ifelse(y01 == 1L, w_pos, w_neg)
}

fit_calibrator <- function(p_raw, y01, method = "none") {
  method <- tolower(as.character(method))
  p <- as.numeric(p_raw)
  y <- as.integer(y01)
  ok <- is.finite(p) & !is.na(p) & !is.na(y)
  
  if (method == "none" || sum(ok) < 30L) return(list(meth = "none", x = NULL, y = NULL))
  
  if (method == "isotonic") {
    p2 <- p[ok]; y2 <- y[ok]
    if (length(unique(p2)) < 5L) return(list(meth = "none", x = NULL, y = NULL))
    ord <- order(p2)
    iso <- stats::isoreg(p2[ord], y2[ord])
    list(meth = "isotonic", x = p2[ord], y = iso$yf)
  } else {
    list(meth = "none", x = NULL, y = NULL)
  }
}

apply_calibrator <- function(cal, p_raw) {
  p <- pmin(pmax(as.numeric(p_raw), 1e-6), 1 - 1e-6)
  if (is.null(cal$x) || is.null(cal$y) || cal$meth == "none") return(p)
  out <- stats::approx(x = cal$x, y = cal$y, xout = p, rule = 2, ties = "ordered")$y
  pmin(pmax(as.numeric(out), 1e-6), 1 - 1e-6)
}

to_factor_y <- function(y01) factor(ifelse(y01 == 1L, "Yes", "No"), levels = c("No","Yes"))

fit_model <- function(algo, X, y01, w_case = NULL) {
  algo <- as.character(algo)
  Xdf <- as.data.frame(X)
  for (nm in names(Xdf)) {
    Xdf[[nm]] <- as.numeric(coerce_numeric_robust(Xdf[[nm]]))
    Xdf[[nm]][!is.finite(Xdf[[nm]]) | is.na(Xdf[[nm]])] <- 0.0
  }
  yfac <- to_factor_y(y01)
  if (length(unique(yfac)) < 2L) return(list(algo = algo, obj = NULL))
  
  if (algo == "LR") {
    df <- cbind(outcome = yfac, Xdf)
    ww <- w_case %||% rep(1.0, nrow(df))
    obj <- suppressWarnings(stats::glm(outcome ~ ., data = df, family = stats::binomial(), weights = ww))
    return(list(algo = algo, obj = obj))
  }
  
  if (algo == "RF") {
    cls_w <- NULL
    if (!is.null(w_case)) {
      cls_w <- c(No = mean(w_case[yfac == "No"]), Yes = mean(w_case[yfac == "Yes"]))
      cls_w[!is.finite(cls_w) | is.na(cls_w)] <- 1.0
    }
    obj <- randomForest::randomForest(x = Xdf, y = yfac, ntree = 500, classwt = cls_w)
    return(list(algo = algo, obj = obj))
  }
  
  if (algo == "SVMRBF") {
    cls_w <- NULL
    if (!is.null(w_case)) {
      cls_w <- c(No = mean(w_case[yfac == "No"]), Yes = mean(w_case[yfac == "Yes"]))
      cls_w[!is.finite(cls_w) | is.na(cls_w)] <- 1.0
    }
    obj <- kernlab::ksvm(
      x = as.matrix(Xdf),
      y = yfac,
      type = "C-svc",
      kernel = "rbfdot",
      prob.model = FALSE,
      class.weights = cls_w
    )
    return(list(algo = algo, obj = obj))
  }
  
  if (algo == "kNN") {
    ntr <- nrow(Xdf)
    k0 <- as.integer(round(sqrt(ntr)))
    k0 <- max(5L, min(31L, k0))
    if (k0 %% 2L == 0L) k0 <- k0 + 1L
    
    train_df <- cbind(outcome = yfac, Xdf)
    
    if (!is.null(w_case)) {
      pr <- as.numeric(w_case)
      pr[!is.finite(pr) | is.na(pr) | pr < 0] <- 0
      if (sum(pr) > 0) {
        pr <- pr / sum(pr)
        idx <- sample(seq_len(nrow(train_df)), size = nrow(train_df), replace = TRUE, prob = pr)
        train_df <- train_df[idx, , drop = FALSE]
      }
    }
    
    return(list(algo = algo, obj = NULL, k = k0, train_df = train_df))
  }
  
  if (algo == "C4.5") {
    df <- cbind(outcome = yfac, Xdf)
    ww <- w_case %||% rep(1.0, nrow(df))
    obj <- rpart::rpart(
      outcome ~ .,
      data = df,
      weights = ww,
      method = "class",
      control = rpart::rpart.control(cp = 0.005, minsplit = 20)
    )
    return(list(algo = algo, obj = obj))
  }
  
  stop("Unknown algo: ", algo)
}

predict_prob <- function(model, X) {
  algo <- model$algo
  Xdf <- as.data.frame(X)
  for (nm in names(Xdf)) {
    Xdf[[nm]] <- as.numeric(coerce_numeric_robust(Xdf[[nm]]))
    Xdf[[nm]][!is.finite(Xdf[[nm]]) | is.na(Xdf[[nm]])] <- 0.0
  }
  
  if (algo == "kNN") {
    train_df <- model$train_df
    if (is.null(train_df)) return(rep(NA_real_, nrow(Xdf)))
    test_df <- cbind(outcome = factor(rep("No", nrow(Xdf)), levels = c("No","Yes")), Xdf)
    
    obj2 <- kknn::kknn(
      outcome ~ ., 
      train = train_df,
      test  = test_df,
      k = model$k %||% 15L,
      distance = 2,
      kernel = "rectangular"
    )
    
    pred_cls <- as.character(obj2$fitted.values)
    pr <- obj2$prob
    if (is.matrix(pr) || is.data.frame(pr)) {
      pr <- as.data.frame(pr)
      if ("Yes" %in% names(pr)) {
        p_yes <- as.numeric(pr[["Yes"]])
      } else {
        p_yes <- as.numeric(pr[[min(2, ncol(pr))]])
      }
    } else {
      p_win <- as.numeric(pr)
      p_yes <- ifelse(pred_cls == "Yes", p_win, 1 - p_win)
    }
    return(pmin(pmax(p_yes, 1e-6), 1 - 1e-6))
  }
  
  if (is.null(model$obj)) return(rep(NA_real_, nrow(Xdf)))
  
  if (algo == "LR") {
    p <- suppressWarnings(stats::predict(model$obj, newdata = Xdf, type = "response"))
    return(pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6))
  }
  
  if (algo == "RF") {
    pr <- suppressWarnings(stats::predict(model$obj, newdata = Xdf, type = "prob"))
    p <- pr[, "Yes"]
    return(pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6))
  }
  
  if (algo == "SVMRBF") {
    pr <- try(suppressWarnings(kernlab::predict(model$obj, newdata = as.matrix(Xdf), type = "probabilities")),
              silent = TRUE)
    if (!inherits(pr, "try-error") && !is.null(dim(pr)) && ("Yes" %in% colnames(pr))) {
      p <- pr[, "Yes"]
      return(pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6))
    }
    dv <- try(suppressWarnings(kernlab::predict(model$obj, newdata = as.matrix(Xdf), type = "decision")),
              silent = TRUE)
    if (inherits(dv, "try-error")) return(rep(NA_real_, nrow(Xdf)))
    dv <- as.numeric(dv)
    p <- 1 / (1 + exp(-dv))
    return(pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6))
  }
  
  if (algo == "C4.5") {
    pr <- suppressWarnings(stats::predict(model$obj, newdata = Xdf, type = "prob"))
    if (is.null(dim(pr)) || !("Yes" %in% colnames(pr))) return(rep(NA_real_, nrow(Xdf)))
    p <- pr[, "Yes"]
    return(pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6))
  }
  
  rep(NA_real_, nrow(Xdf))
}

apply_noise_filter <- function(X, y01, filter_rate, seed) {
  if (is.null(filter_rate) || is.na(filter_rate) || !is.finite(filter_rate) || filter_rate <= 0) {
    return(list(X = X, y01 = y01, keep = seq_len(nrow(X))))
  }
  
  n <- nrow(X)
  n_drop <- as.integer(floor(filter_rate * n))
  if (n_drop < 1L) return(list(X = X, y01 = y01, keep = seq_len(n)))
  
  set.seed(seed)
  
  df <- cbind(outcome = to_factor_y(y01), as.data.frame(X))
  fit <- try(stats::glm(outcome ~ ., data = df, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(X = X, y01 = y01, keep = seq_len(n)))
  
  p <- suppressWarnings(stats::predict(fit, type = "response"))
  p <- pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)
  
  y <- as.integer(y01)
  loss <- -(y * log(p) + (1 - y) * log(1 - p))
  loss[!is.finite(loss) | is.na(loss)] <- -Inf
  
  ord <- order(loss, decreasing = TRUE)
  drop_idx <- ord[seq_len(min(n_drop, length(ord)))]
  keep_idx <- setdiff(seq_len(n), drop_idx)
  
  list(X = X[keep_idx, , drop = FALSE], y01 = y01[keep_idx], keep = keep_idx)
}

make_Xy <- function(df, predictors) {
  predictors <- unique(predictors)
  df2 <- add_missing_cols_as_na(df, predictors)
  X <- df2[, predictors, drop = FALSE]
  for (nm in names(X)) X[[nm]] <- coerce_numeric_robust(X[[nm]])
  y01 <- as.integer(df2$outcome)
  list(X = X, y01 = y01, predictors = predictors)
}

run_oof_single_panel <- function(df_dev, algo, predictors,
                                 clip_q, standardise, use_weights, filter_rate, calibration,
                                 folds) {
  xy <- make_Xy(df_dev, predictors)
  X <- xy$X
  y01 <- xy$y01
  
  n <- nrow(X)
  p_oof <- rep(NA_real_, n)
  
  for (nm in names(folds)) {
    tr <- folds[[nm]]$train
    te <- folds[[nm]]$test
    
    y_tr <- y01[tr]
    if (length(unique(y_tr)) < 2L) {
      p_oof[te] <- NA_real_
      next
    }
    
    X_tr_raw <- X[tr, , drop = FALSE]
    X_te_raw <- X[te, , drop = FALSE]
    
    pp <- prep_fit(X_tr_raw, clip_q = clip_q, standardise = standardise)
    X_tr <- prep_apply(X_tr_raw, pp)
    X_te <- prep_apply(X_te_raw, pp)
    
    w_case <- if (isTRUE(use_weights)) make_case_weights(y_tr) else NULL
    
    nf <- apply_noise_filter(X_tr, y_tr, filter_rate = filter_rate, seed = seed_cv + nchar(nm))
    X_tr2 <- nf$X
    y_tr2 <- nf$y01
    w_case2 <- if (!is.null(w_case)) w_case[nf$keep] else NULL
    
    mod <- fit_model(algo, X_tr2, y_tr2, w_case = w_case2)
    p_te_raw <- predict_prob(mod, X_te)
    
    if (calibration == "isotonic") {
      p_tr_raw <- predict_prob(mod, X_tr2)
      cal <- fit_calibrator(p_tr_raw, y_tr2, method = calibration)
      p_te <- apply_calibrator(cal, p_te_raw)
    } else {
      p_te <- p_te_raw
    }
    
    p_oof[te] <- p_te
  }
  
  o <- orient_prob_global(y01, p_oof)
  list(p_oof = o$p, y01 = y01, flip = o$flip, auc = o$auc, predictors = predictors)
}

cascade_apply <- function(p1, p2, tau_low, tau_high, thr2) {
  p1 <- as.numeric(p1); p2 <- as.numeric(p2)
  defer <- (p1 > tau_low) & (p1 < tau_high)
  
  pred <- ifelse(p1 >= tau_high, 1L,
                 ifelse(p1 <= tau_low, 0L,
                        as.integer(p2 >= thr2)))
  
  score <- ifelse(p1 >= tau_high, 1,
                  ifelse(p1 <= tau_low, 0, p2))
  
  list(pred01 = as.integer(pred), defer = as.logical(defer), score = as.numeric(score))
}

build_cascade_context <- function(
  dev_csv = Sys.getenv("DEV_CSV", unset = "20260101_development_mimic_raw_rowcomplete60.csv"),
  ext_csv = Sys.getenv("EXT_CSV", unset = "20260101_external_mcmed_raw_rowcomplete60.csv"),
  seed_cv = 123,
  k_folds = 5L,
  repeats = 2L,
  dev_max_n = 10000L,
  ext_max_n = 10000L,
  row_thr_global = 0.80,
  clip_q = 0.01,
  standardise = TRUE,
  use_weights = TRUE,
  filter_rate = 0.00,
  calibration = "isotonic",
  thr_grid = seq(0.01, 0.99, by = 0.01),
  sens_min_target = -Inf,
  spec_min_target = -Inf,
  max_def_rate_target = 0.30,
  require_sens_ge_spec = TRUE
) {
  assign("seed_cv", seed_cv, envir = .GlobalEnv)
  set.seed(seed_cv)
  
  stopifnot(file.exists(dev_csv))
  stopifnot(file.exists(ext_csv))
  
  df_dev_full <- readr::read_csv(dev_csv, show_col_types = FALSE) |> prep_df_any("Development")
  df_ext_full <- readr::read_csv(ext_csv, show_col_types = FALSE) |> prep_df_any("External")
  
  df_dev <- apply_global_filter_then_sample(df_dev_full, row_thr = row_thr_global, n_max = dev_max_n, seed = seed_cv + 10001L)
  df_ext <- apply_global_filter_then_sample(df_ext_full, row_thr = row_thr_global, n_max = ext_max_n, seed = seed_cv + 20001L)
  
  if (nrow(df_dev) < 100L) stop("Too few development rows after global filter and sampling")
  if (nrow(df_ext) < 100L) stop("Too few external rows after global filter and sampling")
  
  folds <- make_repeated_stratified_folds(df_dev$outcome, k = k_folds, repeats = repeats, seed = seed_cv)
  
  algos <- c("LR", "RF", "SVMRBF", "kNN", "C4.5")
  
  dev_params_rows <- list()
  dev_cache <- list()
  
  for (algo in algos) {
    preds1 <- get_panel_predictors(df_dev, "Non-invasive")
    res1 <- run_oof_single_panel(
      df_dev = df_dev, algo = algo, predictors = preds1,
      clip_q = clip_q, standardise = standardise, use_weights = use_weights,
      filter_rate = filter_rate, calibration = calibration,
      folds = folds
    )
    
    preds2 <- get_panel_predictors(df_dev, "Laboratory augmented")
    res2 <- run_oof_single_panel(
      df_dev = df_dev, algo = algo, predictors = preds2,
      clip_q = clip_q, standardise = standardise, use_weights = use_weights,
      filter_rate = filter_rate, calibration = calibration,
      folds = folds
    )
    
    stopifnot(length(res1$y01) == length(res2$y01))
    
    thr1_sel <- best_thr_mcc_at_sens(
      y01 = res1$y01, p = res1$p_oof, thr_grid = thr_grid,
      sens_min = sens_min_target, spec_min = spec_min_target
    )
    
    thr2_sel <- best_thr_mcc_at_sens(
      y01 = res2$y01, p = res2$p_oof, thr_grid = thr_grid,
      sens_min = sens_min_target, spec_min = spec_min_target
    )
    
    sel <- best_cascade_params(
      y01 = res1$y01, p1 = res1$p_oof, p2 = res2$p_oof,
      tau_grid = thr_grid, thr2_grid = thr_grid,
      sens_min = sens_min_target, spec_min = spec_min_target,
      max_def_rate = max_def_rate_target,
      require_sens_ge_spec = require_sens_ge_spec
    )
    
    dev_params_rows[[length(dev_params_rows) + 1L]] <- tibble::tibble(
      Algorithm = algo,
      thr_stage1 = thr1_sel$thr,
      thr_stage2 = thr2_sel$thr,
      tau_low = sel$tau_low,
      tau_high = sel$tau_high,
      thr_cascade = sel$thr2
    )
    
    dev_cache[[algo]] <- list(
      preds1 = preds1, flip1 = res1$flip,
      preds2 = preds2, flip2 = res2$flip
    )
  }
  
  dev_params_tbl <- dplyr::bind_rows(dev_params_rows) |>
    dplyr::arrange(.data$Algorithm)
  
  list(
    df_dev = df_dev,
    df_ext = df_ext,
    folds = folds,
    algos = algos,
    dev_cache = dev_cache,
    dev_params_tbl = dev_params_tbl,
    clip_q = clip_q,
    standardise = standardise,
    use_weights = use_weights,
    filter_rate = filter_rate,
    calibration = calibration,
    thr_grid = thr_grid,
    sens_min_target = sens_min_target,
    spec_min_target = spec_min_target,
    max_def_rate_target = max_def_rate_target,
    require_sens_ge_spec = require_sens_ge_spec
  )
}
