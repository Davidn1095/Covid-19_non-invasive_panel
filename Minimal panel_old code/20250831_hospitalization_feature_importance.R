# ============================================================
# Permutation feature importance (ΔMCC) for hospitalization
# ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr); library(ggplot2) })

# 1) pick best algorithm by MCC from hosp_tbl
stopifnot(exists("hosp_tbl") && nrow(hosp_tbl) > 0)
algo_choice <- hosp_tbl %>% arrange(desc(MCC_Mean)) %>% slice(1) %>% pull(Algorithm)
cat(sprintf("\n[Permutation importance] Using algorithm: %s\n", algo_choice))

# 2) helpers specific to permutation
.make_pred_fun <- function(fit, positive_class) {
  force(positive_class)
  function(newdata) {
    pr <- stats::predict(fit, newdata = newdata, type = "prob")
    if (!positive_class %in% colnames(pr)) stop("positive_class column not found in predict(type='prob')")
    as.numeric(pr[[positive_class]])
  }
}
.mcc_at_threshold <- function(truth, prob, pos_label, threshold) {
  truth <- factor(truth)
  lev <- levels(truth)
  neg_label <- setdiff(lev, pos_label)[1]
  pred <- ifelse(prob >= threshold, pos_label, neg_label)
  mcc_binary_vec(truth, pred, pos_label, neg_label)  # uses your earlier helper
}
.permute_one_column <- function(df, col) {
  idx <- sample.int(nrow(df))
  out <- df
  out[[col]] <- out[[col]][idx]
  out
}

# 3) collect per-fold artifacts for the chosen algorithm, reusing your CV scheme
collect_fold_artifacts <- function(yvar, pos_label, algo_name, seed = 444) {
  dat <- df %>%
    filter(!is.na(.data[[yvar]])) %>%
    select(all_of(c(yvar, feature_vars, ".rid"))) %>%
    droplevels()
  lev <- levels(dat[[yvar]])
  splits <- build_cv_splits(dat[[yvar]], R = cv_R, seed = seed)
  spec   <- model_specs[[algo_name]]
  out    <- list()
  for (nm in names(splits)) {
    tr_idx <- splits[[nm]]$train_idx
    va_idx <- splits[[nm]]$test_idx
    dat_tr <- dat[tr_idx, , drop = FALSE]
    dat_va <- dat[va_idx, , drop = FALSE]
    
    inner <- train_with_inner(dat_tr, yvar, spec)
    if (is.null(inner)) next
    if (!all(lev %in% colnames(inner$oof))) next
    
    neg_label <- setdiff(lev, pos_label)[1]
    t_res <- best_threshold_mcc(inner$oof$obs,
                                as.numeric(inner$oof[[pos_label]]),
                                pos = pos_label, neg = neg_label)
    
    bestTune <- inner$fit$bestTune %||% NULL
    fit <- fit_fixed(dat_tr, yvar, spec, bestTune)
    if (is.null(fit)) next
    
    # carry .rid into x_test to satisfy recipe bake
    x_test <- dat_va[, c(feature_vars, ".rid"), drop = FALSE]
    y_test <- dat_va[[yvar]]
    out[[nm]] <- list(fit = fit, x_test = x_test, y_test = y_test, threshold = t_res$t)
  }
  out
}

arts <- collect_fold_artifacts(
  yvar = "Hosp_Bin",
  pos_label = "Yes",
  algo_name = algo_choice,
  seed = 444
)
stopifnot(length(arts) > 0)

# 4) per-fold permutation ΔMCC (no Brier), per-feature only
perm_imp_one_fold <- function(a, positive_class = "Yes", n_permute = 200) {
  x_test <- a$x_test
  y_test <- factor(a$y_test)
  thr    <- a$threshold
  pred_fun <- .make_pred_fun(a$fit, positive_class)
  p_base <- pred_fun(x_test)
  base_mcc <- .mcc_at_threshold(y_test, p_base, positive_class, thr)
  
  feats <- setdiff(colnames(x_test), ".rid")
  purrr::map_dfr(feats, function(fcol) {
    # repeat permutations, average the MCC drop within this fold
    drops <- replicate(n_permute, {
      p <- pred_fun(.permute_one_column(x_test, fcol))
      base_mcc - .mcc_at_threshold(y_test, p, positive_class, thr)
    })
    tibble(
      target = fcol,
      kind = "feature",
      delta_mcc_mean = mean(drops),
      delta_mcc_sd   = stats::sd(drops),
      base_mcc = base_mcc
    )
  })
}

list_of_imp_fold <- lapply(arts, perm_imp_one_fold)

# 5) aggregate across folds, compute 95% CIs across folds
raw_imp <- dplyr::bind_rows(list_of_imp_fold, .id = "fold_id")

perm_mcc_ci <- raw_imp %>%
  dplyr::filter(kind == "feature") %>%
  dplyr::group_by(target) %>%
  dplyr::summarise(
    n_folds = dplyr::n(),
    mean_delta_mcc = mean(delta_mcc_mean),
    se = stats::sd(delta_mcc_mean) / sqrt(n_folds),
    lo = mean_delta_mcc - 1.96 * se,
    hi = mean_delta_mcc + 1.96 * se,
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_delta_mcc))

# 6) plot top features with 95% CI, nothing is saved
top_n <- 12
p_perm <- ggplot2::ggplot(
  perm_mcc_ci %>% dplyr::slice_max(order_by = mean_delta_mcc, n = top_n),
  ggplot2::aes(x = stats::reorder(target, mean_delta_mcc), y = mean_delta_mcc)
) +
  ggplot2::geom_col() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Permutation feature importance",
    y = "ΔMCC after permutation", x = NULL
  ) +
  ggplot2::theme_minimal()

print(p_perm)

# 7) optional table in console, top 30 rows
perm_tbl <- perm_mcc_ci %>%
  dplyr::transmute(Feature = target,
                   `ΔMCC_mean` = round(mean_delta_mcc, 4),
                   `CI_low` = round(lo, 4),
                   `CI_high` = round(hi, 4))
perm_tbl_head <- dplyr::slice_head(perm_tbl, n = 30)
print(perm_tbl_head, n = nrow(perm_tbl_head))
