# ============================================================
# FEATURE IMPORTANCE (Permutation on MCC), no disk writes
# Requires in-memory: df, cfg, feat_full, feat_triage,
# and functions: model_specs, build_cv_splits, train_inner,
# train_final, compute_metrics_binary, best_threshold_mcc
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
})

# --- sanity checks ---
req_objs <- c("df","cfg","feat_full","feat_triage")
req_funs <- c("model_specs","build_cv_splits","train_inner",
              "train_final","compute_metrics_binary","best_threshold_mcc")
missing_objs <- req_objs[!vapply(req_objs, exists, logical(1), inherits = TRUE)]
missing_funs <- req_funs[!vapply(req_funs, function(f) exists(f, mode = "function"), logical(1))]
if (length(missing_objs) || length(missing_funs)) {
  stop(paste0("Missing in memory: ",
              if (length(missing_objs)) paste0("objects {", paste(missing_objs, collapse=", "), "} ") else "",
              if (length(missing_funs)) paste0("functions {", paste(missing_funs, collapse=", "), "}") else ""))
}

# --- utility: single-fold permutation importance ---
perm_importance_fold <- function(dat_tr, dat_va, spec, yvar, features, pos_label, inner_k, verbose = FALSE) {
  lev <- levels(dat_tr[[yvar]])
  neg_label <- setdiff(lev, pos_label)[1]
  
  inner <- train_inner(dat_tr, yvar, spec, inner_k)
  if (is.null(inner) || !all(lev %in% colnames(inner$oof))) return(NULL)
  
  t_res <- best_threshold_mcc(inner$oof$obs,
                              as.numeric(inner$oof[, pos_label]),
                              pos = pos_label, neg = neg_label)
  bestTune <- inner$fit$bestTune %||% NULL
  fit <- train_final(dat_tr, yvar, spec, bestTune, inner_k)
  if (is.null(fit)) return(NULL)
  
  # baseline on validation
  pr0 <- try(predict(fit, newdata = dat_va, type = "prob"), silent = TRUE)
  if (inherits(pr0, "try-error") || is.null(pr0) || !all(lev %in% colnames(pr0))) return(NULL)
  pred0 <- factor(ifelse(pr0[[pos_label]] >= t_res$t, pos_label, neg_label), levels = lev)
  base  <- compute_metrics_binary(dat_va[[yvar]], pred0, pr0[[pos_label]], pos_label, neg_label)[["MCC"]]
  
  # permute each feature in original space, then predict
  rows <- vector("list", length(features))
  for (i in seq_along(features)) {
    f <- features[i]
    dat_perm <- dat_va
    dat_perm[[f]] <- sample(dat_perm[[f]])
    prp <- try(predict(fit, newdata = dat_perm, type = "prob"), silent = TRUE)
    if (inherits(prp, "try-error") || is.null(prp) || !all(lev %in% colnames(prp))) next
    predp <- factor(ifelse(prp[[pos_label]] >= t_res$t, pos_label, neg_label), levels = lev)
    mcc_p <- compute_metrics_binary(dat_va[[yvar]], predp, prp[[pos_label]], pos_label, neg_label)[["MCC"]]
    rows[[i]] <- tibble(Feature = f, Drop_MCC = as.numeric(base - mcc_p))
  }
  dplyr::bind_rows(rows)
}

# --- CV-level runner for importance ---
feature_importance_cv <- function(df, yvar, features, algos, pos_label,
                                  cv_R, cv_k, inner_k, tune_len, seed, verbose = FALSE) {
  dat <- df |>
    dplyr::filter(!is.na(.data[[yvar]])) |>
    dplyr::select(dplyr::all_of(c(yvar, features, ".rid"))) |>
    droplevels()
  lev <- levels(dat[[yvar]])
  splits <- build_cv_splits(dat[[yvar]], R = cv_R, k_desired = cv_k, seed = seed)
  specs  <- model_specs(tune_len, algos)
  
  out <- list()
  for (algo in names(specs)) {
    spec <- specs[[algo]]
    for (nm in names(splits)) {
      tr_idx <- splits[[nm]]$train_idx
      va_idx <- splits[[nm]]$test_idx
      dat_tr <- dat[tr_idx, , drop = FALSE]
      dat_va <- dat[va_idx, , drop = FALSE]
      imp <- perm_importance_fold(dat_tr, dat_va, spec, yvar, features, pos_label, inner_k, verbose)
      if (is.null(imp) || !nrow(imp)) next
      out <- c(out, list(dplyr::mutate(imp, Algorithm = algo, Fold = nm)))
    }
  }
  if (!length(out)) return(tibble())
  dplyr::bind_rows(out)
}

# --- pretty printer ---
print_importance <- function(imp_tbl, title, top_n = 15) {
  cat(sprintf("\n=== %s — Permutation importance (MCC drop) ===\n", title))
  if (!nrow(imp_tbl)) return(cat("(no results)\n"))
  agg <- imp_tbl |>
    dplyr::group_by(Algorithm, Feature) |>
    dplyr::summarise(Mean_Drop = mean(Drop_MCC, na.rm = TRUE),
                     SD_Drop   = stats::sd(Drop_MCC, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::group_by(Algorithm) |>
    dplyr::mutate(Importance = Mean_Drop / pmax(sum(Mean_Drop, na.rm = TRUE), 1e-12)) |>
    dplyr::arrange(Algorithm, dplyr::desc(Mean_Drop)) |>
    dplyr::group_modify(~ dplyr::slice_head(.x, n = min(top_n, nrow(.x)))) |>
    dplyr::ungroup()
  print(agg, n = nrow(agg))
  invisible(agg)
}

# --- run for FULL and TRIAGE ---
yvar <- "Hosp_Bin"
imp_full <- feature_importance_cv(
  df = df, yvar = yvar, features = feat_full, algos = cfg$algos,
  pos_label = cfg$pos_label, cv_R = cfg$cv_R, cv_k = cfg$cv_k,
  inner_k = cfg$inner_k, tune_len = cfg$tune_len, seed = cfg$seed_cv,
  verbose = FALSE
)
imp_triage <- feature_importance_cv(
  df = df, yvar = yvar, features = feat_triage, algos = cfg$algos,
  pos_label = cfg$pos_label, cv_R = cfg$cv_R, cv_k = cfg$cv_k,
  inner_k = cfg$inner_k, tune_len = cfg$tune_len, seed = cfg$seed_cv,
  verbose = FALSE
)

# --- show top features per algorithm for each setup ---
top_full   <- print_importance(imp_full,   "FULL feature set")
top_triage <- print_importance(imp_triage, "TRIAGE feature set")

# Objects in memory: imp_full, imp_triage, top_full, top_triage




# ============================================================
# Top 5 feature importance per algorithm
# 2 rows: Full (row 1), Triage (row 2)
# Separate barplot per algorithm so each shows its own labels
# Requires: imp_full, imp_triage (from permutation-importance step)
# No disk writes
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(patchwork)
})

stopifnot(exists("imp_full", inherits = TRUE), exists("imp_triage", inherits = TRUE))

# mean MCC-drop then keep top 5 within each algorithm
compute_top5 <- function(imp_tbl) {
  imp_tbl %>%
    dplyr::group_by(Algorithm, Feature) %>%
    dplyr::summarise(Mean_Drop = mean(Drop_MCC, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Algorithm) %>%
    dplyr::slice_max(order_by = Mean_Drop, n = 5, with_ties = FALSE) %>%
    dplyr::ungroup()
}

# single horizontal barplot for one algorithm
plot_algo <- function(top_tbl, algo_name) {
  dd <- dplyr::filter(top_tbl, Algorithm == algo_name) %>%
    dplyr::arrange(Mean_Drop) %>%
    dplyr::mutate(Feature = factor(Feature, levels = Feature))
  if (!nrow(dd)) return(ggplot() + theme_void())
  ggplot(dd, aes(x = Mean_Drop, y = Feature)) +
    geom_col(width = 0.7) +
    labs(title = algo_name, x = "MCC drop on permutation", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(margin = ggplot2::margin(t = 6))
    )
}

# prepare lists of per-algorithm plots for full and triage
top_full   <- compute_top5(imp_full)
top_triage <- compute_top5(imp_triage)

algos <- if (exists("cfg", inherits = TRUE) && !is.null(cfg$algos)) {
  as.character(cfg$algos)
} else {
  sort(unique(c(top_full$Algorithm, top_triage$Algorithm)))
}

plots_full   <- lapply(algos, function(a) plot_algo(top_full, a))
plots_triage <- lapply(algos, function(a) plot_algo(top_triage, a))

# assemble 2-row figure
row1 <- patchwork::wrap_plots(plots_full,   nrow = 1)
row2 <- patchwork::wrap_plots(plots_triage, nrow = 1)

final_plot <- row1 / row2 +
  plot_annotation(
    title = "Top 5 features per algorithm",
    subtitle = "Row 1 Full, row 2 Triage"
  )

print(final_plot)






