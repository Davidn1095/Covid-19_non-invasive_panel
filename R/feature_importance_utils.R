mcc_safe_at_thr <- function(y01, p, thr) {
  met <- compute_metrics_at_thr(y01, p, thr)
  m <- met$MCC
  if (!is.finite(m) || is.na(m)) 0 else as.numeric(m)
}

mcc_safe_from_pred <- function(y01, pred01, score = NULL) {
  met <- compute_metrics_from_pred(y01, pred01, score = score)
  m <- met$MCC
  if (!is.finite(m) || is.na(m)) 0 else as.numeric(m)
}

fit_full_panel_obj <- function(df_dev, algo, predictors,
                               clip_q, standardise, use_weights, filter_rate, calibration,
                               flip_from_dev = FALSE, seed_fit = 1L) {
  set.seed(seed_fit)
  
  xy_dev <- make_Xy(df_dev, predictors)
  X_dev_raw <- xy_dev$X
  y_dev01 <- xy_dev$y01
  
  pp <- prep_fit(X_dev_raw, clip_q = clip_q, standardise = standardise)
  X_dev <- prep_apply(X_dev_raw, pp)
  
  w_case <- if (isTRUE(use_weights)) make_case_weights(y_dev01) else NULL
  
  nf <- apply_noise_filter(X_dev, y_dev01, filter_rate = filter_rate, seed = seed_fit + 777L)
  X_dev2 <- nf$X
  y_dev2 <- nf$y01
  w_case2 <- if (!is.null(w_case)) w_case[nf$keep] else NULL
  
  mod <- fit_model(algo, X_dev2, y_dev2, w_case = w_case2)
  
  cal <- list(meth = "none", x = NULL, y = NULL)
  if (tolower(as.character(calibration)) == "isotonic") {
    p_dev_raw <- predict_prob(mod, X_dev2)
    cal <- fit_calibrator(p_dev_raw, y_dev2, method = "isotonic")
  }
  
  list(
    algo = as.character(algo),
    predictors = unique(predictors),
    pp = pp,
    mod = mod,
    cal = cal,
    flip = isTRUE(flip_from_dev)
  )
}

predict_prob_from_fitobj_Xraw <- function(fitobj, X_raw) {
  X_raw <- as.data.frame(X_raw)
  X_new <- prep_apply(X_raw, fitobj$pp)
  
  p_raw <- predict_prob(fitobj$mod, X_new)
  
  if (!is.null(fitobj$cal) && !is.null(fitobj$cal$x) && !is.null(fitobj$cal$y) &&
      tolower(as.character(fitobj$cal$meth)) == "isotonic") {
    p <- apply_calibrator(fitobj$cal, p_raw)
  } else {
    p <- p_raw
  }
  
  if (isTRUE(fitobj$flip)) p <- 1 - p
  pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)
}

perm_drop_mcc_panel_fitobj <- function(fitobj, X_ext_raw, y_ext01, thr, seed_prefix) {
  X0 <- as.data.frame(X_ext_raw[, fitobj$predictors, drop = FALSE])
  
  p0 <- predict_prob_from_fitobj_Xraw(fitobj, X0)
  m0 <- mcc_safe_at_thr(y_ext01, p0, thr)
  
  out <- vector("list", length(fitobj$predictors))
  
  for (i in seq_along(fitobj$predictors)) {
    f <- fitobj$predictors[[i]]
    X1 <- X0
    
    set.seed(seed_from_name(paste0(seed_prefix, "|", f), offset = 99001L))
    perm_idx <- sample.int(nrow(X1), nrow(X1), replace = FALSE)
    X1[[f]] <- X1[[f]][perm_idx]
    
    p1 <- predict_prob_from_fitobj_Xraw(fitobj, X1)
    m1 <- mcc_safe_at_thr(y_ext01, p1, thr)
    
    out[[i]] <- tibble::tibble(
      Feature  = f,
      Drop_MCC = m0 - m1
    )
  }
  
  dplyr::bind_rows(out)
}

perm_drop_mcc_cascade_fitobjs <- function(fit1, fit2, X1_ext_raw, X2_ext_raw, y_ext01,
                                          tau_low, tau_high, thr2, seed_prefix) {
  X10 <- as.data.frame(X1_ext_raw[, fit1$predictors, drop = FALSE])
  X20 <- as.data.frame(X2_ext_raw[, fit2$predictors, drop = FALSE])
  
  p10 <- predict_prob_from_fitobj_Xraw(fit1, X10)
  p20 <- predict_prob_from_fitobj_Xraw(fit2, X20)
  
  cas0 <- cascade_apply(p10, p20, tau_low, tau_high, thr2)
  m0 <- mcc_safe_from_pred(y_ext01, cas0$pred01, score = cas0$score)
  
  feats <- sort(unique(c(names(X10), names(X20))))
  out <- vector("list", length(feats))
  
  for (i in seq_along(feats)) {
    f <- feats[[i]]
    X1 <- X10
    X2 <- X20
    
    set.seed(seed_from_name(paste0(seed_prefix, "|", f), offset = 99002L))
    perm_idx <- sample.int(length(y_ext01), length(y_ext01), replace = FALSE)
    
    if (f %in% names(X1)) X1[[f]] <- X1[[f]][perm_idx]
    if (f %in% names(X2)) X2[[f]] <- X2[[f]][perm_idx]
    
    p1 <- predict_prob_from_fitobj_Xraw(fit1, X1)
    p2 <- predict_prob_from_fitobj_Xraw(fit2, X2)
    
    cas1 <- cascade_apply(p1, p2, tau_low, tau_high, thr2)
    m1 <- mcc_safe_from_pred(y_ext01, cas1$pred01, score = cas1$score)
    
    out[[i]] <- tibble::tibble(
      Feature  = f,
      Drop_MCC = m0 - m1
    )
  }
  
  dplyr::bind_rows(out)
}

build_feature_importance_perm_dMCC <- function(df_dev, df_ext, dev_cache, dev_params_tbl, algos,
                                               clip_q, standardise, use_weights, filter_rate, calibration) {
  feats_noninv <- sort(unique(unlist(lapply(dev_cache, function(z) z$preds1))))
  feats_labaug <- sort(unique(unlist(lapply(dev_cache, function(z) z$preds2))))
  feats_cascade <- sort(unique(c(feats_noninv, feats_labaug)))
  
  y_ext01 <- as.integer(df_ext$outcome)
  
  fi_rows <- list()
  k_row <- 0L
  
  for (algo in algos) {
    cache <- dev_cache[[algo]]
    rowp <- dev_params_tbl[dev_params_tbl$Algorithm == algo, , drop = FALSE]
    stopifnot(nrow(rowp) == 1L)
    
    fit1 <- fit_full_panel_obj(
      df_dev = df_dev, algo = algo, predictors = cache$preds1,
      clip_q = clip_q, standardise = standardise, use_weights = use_weights,
      filter_rate = filter_rate, calibration = calibration,
      flip_from_dev = cache$flip1,
      seed_fit = seed_from_name(paste0(algo, "|fit|stage1"), offset = 3100L)
    )
    
    fit2 <- fit_full_panel_obj(
      df_dev = df_dev, algo = algo, predictors = cache$preds2,
      clip_q = clip_q, standardise = standardise, use_weights = use_weights,
      filter_rate = filter_rate, calibration = calibration,
      flip_from_dev = cache$flip2,
      seed_fit = seed_from_name(paste0(algo, "|fit|stage2"), offset = 3200L)
    )
    
    X1_ext_raw <- make_Xy(df_ext, fit1$predictors)$X
    X2_ext_raw <- make_Xy(df_ext, fit2$predictors)$X
    
    k_row <- k_row + 1L
    fi_rows[[k_row]] <- perm_drop_mcc_panel_fitobj(
      fitobj = fit1,
      X_ext_raw = X1_ext_raw,
      y_ext01 = y_ext01,
      thr = as.numeric(rowp$thr_stage1),
      seed_prefix = paste0("FI|", algo, "|Non-invasive")
    ) %>% dplyr::mutate(Algorithm = algo, Policy = "Non-invasive")
    
    k_row <- k_row + 1L
    fi_rows[[k_row]] <- perm_drop_mcc_panel_fitobj(
      fitobj = fit2,
      X_ext_raw = X2_ext_raw,
      y_ext01 = y_ext01,
      thr = as.numeric(rowp$thr_stage2),
      seed_prefix = paste0("FI|", algo, "|Laboratory augmented")
    ) %>% dplyr::mutate(Algorithm = algo, Policy = "Laboratory augmented")
    
    k_row <- k_row + 1L
    fi_rows[[k_row]] <- perm_drop_mcc_cascade_fitobjs(
      fit1 = fit1, fit2 = fit2,
      X1_ext_raw = X1_ext_raw,
      X2_ext_raw = X2_ext_raw,
      y_ext01 = y_ext01,
      tau_low = as.numeric(rowp$tau_low),
      tau_high = as.numeric(rowp$tau_high),
      thr2 = as.numeric(rowp$thr_cascade),
      seed_prefix = paste0("FI|", algo, "|Cascade")
    ) %>% dplyr::mutate(Algorithm = algo, Policy = "Cascade")
  }
  
  fi_raw <- dplyr::bind_rows(fi_rows) %>%
    dplyr::group_by(.data$Policy, .data$Algorithm, .data$Feature) %>%
    dplyr::summarise(Drop_MCC = mean(.data$Drop_MCC, na.rm = TRUE), .groups = "drop")
  
  grid_noninv <- tidyr::expand_grid(Policy = "Non-invasive", Algorithm = algos, Feature = feats_noninv)
  grid_labaug <- tidyr::expand_grid(Policy = "Laboratory augmented", Algorithm = algos, Feature = feats_labaug)
  grid_cascade <- tidyr::expand_grid(Policy = "Cascade", Algorithm = algos, Feature = feats_cascade)
  
  fi_grid <- dplyr::bind_rows(grid_noninv, grid_labaug, grid_cascade) %>%
    dplyr::left_join(fi_raw, by = c("Policy","Algorithm","Feature"))
  
  feature_importance_perm_dMCC <- fi_grid %>%
    dplyr::mutate(
      Policy = factor(.data$Policy, levels = c("Non-invasive","Laboratory augmented","Cascade"))
    ) %>%
    dplyr::group_by(.data$Policy, .data$Algorithm) %>%
    dplyr::mutate(
      Drop_MCC_rankval = ifelse(is.finite(.data$Drop_MCC), round(.data$Drop_MCC, 12), NA_real_),
      Rank = dplyr::if_else(
        is.na(.data$Drop_MCC_rankval),
        NA_integer_,
        dplyr::dense_rank(dplyr::desc(.data$Drop_MCC_rankval))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$Drop_MCC_rankval)
  
  consensus_tbl <- feature_importance_perm_dMCC %>%
    dplyr::group_by(.data$Policy, .data$Feature) %>%
    dplyr::summarise(
      ConsensusScore = mean(.data$Rank, na.rm = TRUE),
      n_algos        = sum(!is.na(.data$Rank)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ConsensusScore = ifelse(is.finite(.data$ConsensusScore), .data$ConsensusScore, Inf)
    ) %>%
    dplyr::group_by(.data$Policy) %>%
    dplyr::mutate(
      ConsensusRank = dplyr::dense_rank(.data$ConsensusScore)
    ) %>%
    dplyr::ungroup()
  
  feature_importance_perm_dMCC <- feature_importance_perm_dMCC %>%
    dplyr::left_join(consensus_tbl, by = c("Policy","Feature")) %>%
    dplyr::arrange(.data$Policy, .data$ConsensusRank, .data$ConsensusScore, dplyr::desc(.data$n_algos),
                   .data$Feature, .data$Algorithm)
  
  feature_importance_perm_dMCC_wide <- feature_importance_perm_dMCC %>%
    dplyr::select(.data$Policy, .data$ConsensusRank, .data$ConsensusScore, .data$n_algos,
                  .data$Feature, .data$Algorithm, .data$Drop_MCC, .data$Rank) %>%
    dplyr::arrange(.data$Policy, .data$ConsensusRank, .data$ConsensusScore, dplyr::desc(.data$n_algos),
                   .data$Feature, .data$Algorithm)
  
  list(
    feature_importance_perm_dMCC = feature_importance_perm_dMCC,
    feature_importance_perm_dMCC_wide = feature_importance_perm_dMCC_wide,
    consensus_tbl = consensus_tbl,
    feats_noninv = feats_noninv,
    feats_labaug = feats_labaug,
    feats_cascade = feats_cascade
  )
}
