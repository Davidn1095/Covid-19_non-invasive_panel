# Compact ΔMCC heatmap with triage-first ordering, reverse alphabetical inside blocks,
# color-blind friendly palette, contrast stretched around zero

# ==== Build holdout objects for permutation ΔMCC importance ====

# Train once per algorithm and return fit + external data + threshold
build_holds_for_set <- function(df, yvar, features, set_label) {
  specs  <- model_specs(cfg$tune_len)                # uses your algo_order & grids
  algos  <- names(specs)
  out    <- setNames(vector("list", length(algos)), algos)
  
  row_train <- which(df$data_split == "train")
  row_ext   <- which(df$data_split == "external")
  
  # data for CV fold construction (time-blocked if enabled)
  dtrain_with_date <- df[row_train, , drop = FALSE]
  
  for (algo in algos) {
    spec <- specs[[algo]]
    
    # columns needed for modeling
    model_cols <- unique(c(features, yvar))
    model_cols <- intersect(model_cols, names(df))
    
    dtrain <- df[row_train, model_cols, drop = FALSE]
    dext   <- df[row_ext,   model_cols, drop = FALSE]
    if (!nrow(dtrain) || !nrow(dext)) next
    
    # same CV/time-blocking strategy as the main pipeline
    seed_base <- seed_for(algo, set_label)
    cv_idx <- if (isTRUE(cfg$cv_blocked_by_time)) {
      make_timeblock_indices(dtrain_with_date, k = cfg$inner_k)
    } else {
      make_cv_indices(dtrain[[yvar]], k = cfg$inner_k, seed = seed_base,
                      repeats = effective_repeats())
    }
    
    # keep ENN logging consistent
    old_key <- getOption("enn_log_key", NULL)
    options(enn_log_key = paste(algo, set_label, sep = "|"))
    on.exit(options(enn_log_key = old_key), add = TRUE)
    
    # fit the final model with your train_inner()
    fit <- try(train_inner(dtrain, yvar, spec, cfg$inner_k, cv_idx, seed_base), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit)) next
    
    # derive a threshold on TRAIN OOF (uncalibrated, to match most perm routines)
    oof <- get_oof_bestTune(fit, pos_label = cfg$pos_label)
    thr <- if (!is.null(oof) && nrow(oof)) {
      best_threshold_mcc_bacc_med(oof$obs, as_num(oof$p))$t
    } else {
      # fallback: compute on raw train predictions
      p_tr <- try(as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label]), silent = TRUE)
      if (inherits(p_tr, "try-error")) 0.5 else best_threshold_mcc_bacc_med(dtrain[[yvar]], p_tr)$t
    }
    
    out[[algo]] <- list(
      fit      = fit,
      threshold= thr,
      data     = dext[, unique(c(features, yvar)), drop = FALSE],
      features = intersect(features, names(dext)),
      yvar     = yvar
    )
  }
  out
}

# Build the two objects your feature-importance code expects
holds_full   <- build_holds_for_set(df, yvar, feat_full,   "Complete")
holds_triage <- build_holds_for_set(df, yvar, feat_triage, "Triage")


suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
})

# safe permutation ΔMCC
perm_drop_mcc_safe <- function(fit, data, yvar, features, threshold,
                               pos_label = cfg$pos_label, neg_label = cfg$neg_label) {
  p0 <- as.numeric(predict(fit, newdata = data, type = "prob")[, pos_label])
  pr0 <- factor(ifelse(p0 >= threshold, pos_label, neg_label),
                levels = c(neg_label, pos_label))
  m0 <- as.numeric(compute_metrics_binary(data[[yvar]], pr0, p0,
                                          pos_label = pos_label, neg_label = neg_label)["MCC"])
  purrr::map_dfr(features, function(f) {
    dperm <- data
    dperm[[f]] <- sample(dperm[[f]])
    pp <- as.numeric(predict(fit, newdata = dperm, type = "prob")[, pos_label])
    pr <- factor(ifelse(pp >= threshold, pos_label, neg_label),
                 levels = c(neg_label, pos_label))
    m1 <- as.numeric(compute_metrics_binary(data[[yvar]], pr, pp,
                                            pos_label = pos_label, neg_label = neg_label)["MCC"])
    tibble(Feature = f, Drop_MCC = m0 - m1)
  })
}

# human readable labels
pretty_feature_label <- function(x) {
  map <- c(
    "WHO_score_admission_mod" = "WHO score at admission",
    "SpO2_admission"          = "SpO\u2082 at admission",
    "monocyte_abs_number"     = "Monocytes (abs #)",
    "lymphocyte_abs_number"   = "Lymphocytes (abs #)",
    "neutrophil_abs_number"   = "Neutrophils (abs #)",
    "monocytes_perc"          = "Monocytes (%)",
    "lymphocytes_perc"        = "Lymphocytes (%)",
    "neutrophils_perc"        = "Neutrophils (%)",
    "D_Dimer" = "D-dimer", "CRP" = "CRP", "albumin" = "Albumin",
    "Diagnosis" = "Diagnosis", "Gender" = "Gender", "Age" = "Age",
    "NLR"="NLR","MLR"="MLR","SIRI"="SIRI","NMR"="NMR","CAR"="CAR","CLR"="CLR","PNI"="PNI","DLR"="DLR"
  )
  y <- unname(map[match(x, names(map))])
  idx <- is.na(y)
  if (any(idx)) {
    y[idx] <- x[idx]
    y[idx] <- gsub("_perc$", " (%)", y[idx])
    y[idx] <- gsub("_abs_number$", " (abs #)", y[idx])
    y[idx] <- gsub("_", " ", y[idx], fixed = TRUE)
    y[idx] <- tools::toTitleCase(y[idx])
  }
  y
}

# build ΔMCC from existing holdouts
importance_from_holds <- function(holds, feat_vec, set_label) {
  purrr::map_dfr(names(holds), function(ak) {
    H <- holds[[ak]]
    # skip if this algo failed to train or returned nothing
    if (is.null(H) || is.null(H$fit) || is.null(H$data)) return(NULL)
    
    perm_drop_mcc_safe(
      fit       = H$fit,                          # was H$model
      data      = H$data,                         # was H$holdout_df
      yvar      = H$yvar %||% yvar,
      features  = intersect(feat_vec, H$features %||% names(H$data)),
      threshold = H$threshold
    ) %>%
      dplyr::mutate(
        Algorithm     = factor(ak, levels = algo_order),
        Feature_label = pretty_feature_label(Feature),
        Set           = set_label
      )
  })
}

# importance for both sets
imp_full <- importance_from_holds(holds_full, feat_full,  "Complete")
imp_tri <- importance_from_holds(holds_triage, feat_triage, "Triage")


# ----- Correct feature order -----
# Triage panel = reverse alphabetical
# Complete panel = triage features first (reverse alphabetical), then the rest (reverse alphabetical)

# split features
tri_feats_raw   <- intersect(feat_triage, feat_full)
other_feats_raw <- setdiff(feat_full, tri_feats_raw)

# human labels, then reverse alphabetical
tri_rev   <- sort(pretty_feature_label(tri_feats_raw),   decreasing = TRUE)
other_rev <- sort(pretty_feature_label(other_feats_raw), decreasing = TRUE)

# panel-specific factor levels (note: ggplot places the FIRST level at the BOTTOM)
# so to have triage at the TOP of Complete, others are put first in the factor levels
levels_full <- c(other_rev, tri_rev)   # bottom→top: others then triage  → TOP shows triage
levels_tri  <- rev(tri_rev)            # bottom→top: A→W              → TOP shows W→A (reverse alpha)

# complete per-panel grids, join ΔMCC
grid_full <- tidyr::expand_grid(
  Set = factor("Complete", levels = c("Complete","Triage")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_full, levels = levels_full)
)
grid_tri <- tidyr::expand_grid(
  Set = factor("Triage", levels = c("Complete","Triage")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_tri, levels = levels_tri)
)

df_full <- dplyr::left_join(
  grid_full,
  imp_full %>% dplyr::mutate(Feature_label = factor(Feature_label, levels = levels_full)) %>%
    dplyr::select(Set, Algorithm, Feature_label, Drop_MCC),
  by = c("Set","Algorithm","Feature_label")
)
df_tri <- dplyr::left_join(
  grid_tri,
  imp_tri %>% dplyr::mutate(Feature_label = factor(Feature_label, levels = levels_tri)) %>%
    dplyr::select(Set, Algorithm, Feature_label, Drop_MCC),
  by = c("Set","Algorithm","Feature_label")
)

df_plot <- dplyr::bind_rows(df_full, df_tri)

# use the original orange–blue palette with stretch around zero
lim <- stats::quantile(abs(df_plot$Drop_MCC), 0.95, na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- max(abs(df_plot$Drop_MCC), na.rm = TRUE)

p_heat <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Algorithm, y = Feature_label, fill = Drop_MCC)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::scale_fill_gradient2(
    low = "#D55E00",  # orange for negative
    mid = "white",
    high = "#0072B2", # blue for positive
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish,
    name = expression(Delta*"MCC"),
    na.value = "white"
  ) +
  ggplot2::facet_grid(Set ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x.top = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8),
    strip.text.y.right = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold")
  )+ ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "grey60", fill = NA, linewidth = 0.6)
  )


print(p_heat)
