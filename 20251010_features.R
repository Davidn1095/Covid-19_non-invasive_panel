# Compact ΔMCC heatmap with triage-first ordering, now including the Cascade set
# Reverse alphabetical inside blocks, color-blind friendly palette, contrast stretched around zero
# Expects: cfg, df, yvar, feat_full, feat_triage, algo_order, model_specs, seed_for,
#          effective_repeats, make_timeblock_indices, make_cv_indices, train_inner,
#          get_oof_bestTune, compute_metrics_binary, best_threshold_mcc_bacc_med

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
  library(purrr); library(tibble); library(scales); library(grid)
})

# ---------- helper for NULL defaults ----------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- deterministic seed helper ----------
seed_from_name <- function(name, offset = 0L) {
  as.integer(offset + sum(utf8ToInt(as.character(name))))
}

# ---------- Cascade helpers (band tuning and combination) ----------
tune_band_from_triage_oof <- function(oof_df, thr_tri,
                                      grid_tl = cfg$cascade$grid_tl %||% seq(0.05, 0.45, 0.05),
                                      grid_th = cfg$cascade$grid_th %||% seq(0.55, 0.95, 0.05),
                                      max_defer = cfg$cascade$max_defer %||% 1) {
  stopifnot(all(c("obs","p") %in% names(oof_df)))
  pos <- cfg$pos_label; neg <- cfg$neg_label
  best <- NULL
  for (tl in grid_tl) for (th in grid_th) if (tl < th) {
    defer <- (oof_df$p > tl & oof_df$p < th)
    if (mean(defer) > max_defer) next
    keep <- !defer; if (!any(keep)) next
    pred_keep <- factor(ifelse(oof_df$p[keep] >= thr_tri, pos, neg), levels = c(neg, pos))
    obs_keep  <- factor(oof_df$obs[keep], levels = c(neg, pos))
    m <- compute_metrics_binary(obs_keep, pred_keep, p_pos = oof_df$p[keep], pos_label = pos, neg_label = neg)
    score <- as.numeric(m["MCC"])
    cand <- list(t_low = tl, t_high = th, mcc = score, defer = mean(defer), n_keep = sum(keep))
    if (is.null(best) ||
        cand$mcc > best$mcc + 1e-12 ||
        (abs(cand$mcc - best$mcc) < 1e-12 &&
         (cand$defer < best$defer - 1e-12 ||
          (abs(cand$defer - best$defer) < 1e-12 && cand$n_keep > best$n_keep)))) best <- cand
  }
  if (is.null(best)) list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high, mcc = NA, defer = NA, n_keep = NA) else best
}

# ---------- Build holdout objects for permutation ΔMCC importance ----------
build_holds_for_set <- function(df, yvar, features, set_label) {
  specs  <- model_specs(cfg$tune_len)
  algos  <- names(specs)
  out    <- setNames(vector("list", length(algos)), algos)
  
  row_train <- which(df$data_split == "train")
  row_ext   <- which(df$data_split == "external")
  dtrain_with_date <- df[row_train, , drop = FALSE]
  
  for (algo in algos) {
    spec <- specs[[algo]]
    model_cols <- unique(c(features, yvar))
    model_cols <- intersect(model_cols, names(df))
    dtrain <- df[row_train, model_cols, drop = FALSE]
    dext   <- df[row_ext,   model_cols, drop = FALSE]
    if (!nrow(dtrain) || !nrow(dext)) next
    
    seed_base <- seed_for(algo, set_label)
    cv_idx <- if (isTRUE(cfg$cv_blocked_by_time)) {
      make_timeblock_indices(dtrain_with_date, k = cfg$inner_k)
    } else {
      make_cv_indices(dtrain[[yvar]], k = cfg$inner_k, seed = seed_base, repeats = effective_repeats())
    }
    
    old_key <- getOption("enn_log_key", NULL)
    options(enn_log_key = paste(algo, set_label, sep = "|"))
    on.exit(options(enn_log_key = old_key), add = TRUE)
    
    fit <- try(train_inner(dtrain, yvar, spec, cfg$inner_k, cv_idx, seed_base), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit)) next
    
    oof <- get_oof_bestTune(fit, pos_label = cfg$pos_label)
    thr <- if (!is.null(oof) && nrow(oof)) {
      best_threshold_mcc_bacc_med(oof$obs, as.numeric(oof$p))$t
    } else {
      p_tr <- try(as.numeric(predict(fit, newdata = dtrain, type = "prob")[, cfg$pos_label]), silent = TRUE)
      if (inherits(p_tr, "try-error")) 0.5 else best_threshold_mcc_bacc_med(dtrain[[yvar]], p_tr)$t
    }
    
    out[[algo]] <- list(
      fit       = fit,
      threshold = thr,
      oof       = oof,
      data      = dext[, unique(c(features, yvar)), drop = FALSE],
      features  = intersect(features, names(dext)),
      yvar      = yvar
    )
  }
  out
}

# ---------- permutation ΔMCC for single-model sets ----------
perm_drop_mcc_safe <- function(fit, data, yvar, features, threshold,
                               pos_label = cfg$pos_label, neg_label = cfg$neg_label) {
  p0  <- as.numeric(predict(fit, newdata = data, type = "prob")[, pos_label])
  pr0 <- factor(ifelse(p0 >= threshold, pos_label, neg_label), levels = c(neg_label, pos_label))
  m0  <- as.numeric(compute_metrics_binary(data[[yvar]], pr0, p0,
                                           pos_label = pos_label, neg_label = neg_label)["MCC"])
  
  feats <- intersect(unique(as.character(features)), names(data))
  purrr::map_dfr(feats, function(f) {
    dperm <- data
    set.seed(seed_from_name(f, 1007L))
    idx <- sample.int(nrow(dperm), nrow(dperm))
    dperm[[f]] <- dperm[[f]][idx]
    
    pp <- as.numeric(predict(fit, newdata = dperm, type = "prob")[, pos_label])
    pr <- factor(ifelse(pp >= threshold, pos_label, neg_label), levels = c(neg_label, pos_label))
    m1 <- as.numeric(compute_metrics_binary(dperm[[yvar]], pr, pp,
                                            pos_label = pos_label, neg_label = neg_label)["MCC"])
    tibble(Feature = f, Drop_MCC = m0 - m1)
  })
}

# ---------- permutation ΔMCC for Cascade (uses union of triage+complete features) ----------
perm_drop_mcc_cascade_safe <- function(fit_tri, fit_full,
                                       data_tri, data_full, yvar,
                                       features, thr_tri, thr_full,
                                       t_low, t_high,
                                       pos_label = cfg$pos_label, neg_label = cfg$neg_label) {
  stopifnot(nrow(data_tri) == nrow(data_full))
  
  # baseline decision: triage outside band, complete inside band
  p_tri0  <- as.numeric(predict(fit_tri,  newdata = data_tri,  type = "prob")[, pos_label])
  p_full0 <- as.numeric(predict(fit_full, newdata = data_full, type = "prob")[, pos_label])
  defer0  <- (p_tri0 > t_low & p_tri0 < t_high)
  pr0 <- factor(ifelse(!defer0,
                       ifelse(p_tri0  >= thr_tri,  pos_label, neg_label),
                       ifelse(p_full0 >= thr_full, pos_label, neg_label)),
                levels = c(neg_label, pos_label))
  p_used0 <- ifelse(!defer0, p_tri0, p_full0)
  m0 <- as.numeric(compute_metrics_binary(data_full[[yvar]], pr0, p_used0,
                                          pos_label = pos_label, neg_label = neg_label)["MCC"])
  
  # use union of features across stages, permute where available
  feats <- intersect(unique(as.character(features)), union(names(data_tri), names(data_full)))
  
  purrr::map_dfr(feats, function(f) {
    dtri  <- data_tri
    dfull <- data_full
    
    set.seed(seed_from_name(f, 20251010L))
    idx <- sample.int(nrow(dtri), nrow(dtri))
    if (f %in% names(dtri))  dtri[[f]]  <- dtri[[f]][idx]
    if (f %in% names(dfull)) dfull[[f]] <- dfull[[f]][idx]
    
    p_tri  <- as.numeric(predict(fit_tri,  newdata = dtri,  type = "prob")[, pos_label])
    p_full <- as.numeric(predict(fit_full, newdata = dfull, type = "prob")[, pos_label])
    
    defer <- (p_tri > t_low & p_tri < t_high)
    pr <- factor(ifelse(!defer,
                        ifelse(p_tri  >= thr_tri,  pos_label, neg_label),
                        ifelse(p_full >= thr_full, pos_label, neg_label)),
                 levels = c(neg_label, pos_label))
    p_used <- ifelse(!defer, p_tri, p_full)
    
    m1 <- as.numeric(compute_metrics_binary(dfull[[yvar]], pr, p_used,
                                            pos_label = pos_label, neg_label = neg_label)["MCC"])
    tibble(Feature = f, Drop_MCC = m0 - m1)
  })
}

# ---------- build per-set permutation ΔMCC ----------
importance_from_holds <- function(holds, features, set_label) {
  algos <- names(holds)
  purrr::map_dfr(algos, function(algo) {
    h <- holds[[algo]]
    if (is.null(h) || is.null(h$fit) || !is.data.frame(h$data) || !nrow(h$data)) return(tibble())
    imp <- try(
      perm_drop_mcc_safe(
        fit       = h$fit,
        data      = h$data,
        yvar      = h$yvar,
        features  = h$features,
        threshold = h$threshold,
        pos_label = cfg$pos_label,
        neg_label = cfg$neg_label
      ),
      silent = TRUE
    )
    if (inherits(imp, "try-error") || !nrow(imp)) return(tibble())
    imp %>%
      mutate(
        Set           = set_label,
        Algorithm     = factor(algo, levels = algo_order),
        Feature_label = pretty_feature_label(Feature)
      ) %>%
      select(Set, Algorithm, Feature_label, Drop_MCC)
  })
}

# ---------- Cascade permutation ΔMCC using triage band tuning when available ----------
importance_cascade_from_holds <- function(holds_tri, holds_full, features_union) {
  algos <- intersect(names(holds_tri), names(holds_full))
  purrr::map_dfr(algos, function(algo) {
    ht <- holds_tri[[algo]]
    hf <- holds_full[[algo]]
    if (is.null(ht) || is.null(hf) || is.null(ht$fit) || is.null(hf$fit)) return(tibble())
    
    tuned <- try(
      if (!is.null(ht$oof) && is.data.frame(ht$oof) && nrow(ht$oof)) {
        tune_band_from_triage_oof(ht$oof, thr_tri = ht$threshold)
      } else {
        list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high)
      },
      silent = TRUE
    )
    if (inherits(tuned, "try-error") || is.null(tuned$t_low) || is.null(tuned$t_high)) {
      tuned <- list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high)
    }
    
    cols_tri  <- unique(c(ht$features, features_union, ht$yvar))
    cols_full <- unique(c(hf$features, features_union, hf$yvar))
    dtri  <- ht$data[,  intersect(cols_tri,  names(ht$data)),  drop = FALSE]
    dfull <- hf$data[, intersect(cols_full, names(hf$data)), drop = FALSE]
    
    imp <- try(
      perm_drop_mcc_cascade_safe(
        fit_tri  = ht$fit,
        fit_full = hf$fit,
        data_tri = dtri,
        data_full = dfull,
        yvar = hf$yvar,
        features = features_union,
        thr_tri = ht$threshold,
        thr_full = hf$threshold,
        t_low = tuned$t_low,
        t_high = tuned$t_high,
        pos_label = cfg$pos_label,
        neg_label = cfg$neg_label
      ),
      silent = TRUE
    )
    if (inherits(imp, "try-error") || !nrow(imp)) return(tibble())
    imp %>%
      mutate(
        Set           = "Cascade",
        Algorithm     = factor(algo, levels = algo_order),
        Feature_label = pretty_feature_label(Feature)
      ) %>%
      select(Set, Algorithm, Feature_label, Drop_MCC)
  })
}

# ---------- human readable labels ----------
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

# ---------- build holdouts ----------
holds_full   <- build_holds_for_set(df, yvar, feat_full,   "Complete")
holds_triage <- build_holds_for_set(df, yvar, feat_triage, "Triage")

# ---------- importance for all three sets ----------
imp_full <- purrr::quietly(importance_from_holds)(holds_full,   feat_full,   "Complete")$result
imp_tri  <- purrr::quietly(importance_from_holds)(holds_triage, feat_triage, "Triage")$result
feat_union <- union(feat_full, feat_triage)
imp_cas  <- purrr::quietly(importance_cascade_from_holds)(holds_triage, holds_full, feat_union)$result

# ---------- feature order ----------
tri_feats_raw   <- intersect(feat_triage, feat_full)
other_feats_raw <- setdiff(feat_full, tri_feats_raw)
tri_rev   <- sort(pretty_feature_label(tri_feats_raw),   decreasing = TRUE)
other_rev <- sort(pretty_feature_label(other_feats_raw), decreasing = TRUE)
levels_full_like <- c(other_rev, tri_rev)
levels_tri       <- rev(tri_rev)

# ---------- panel grids ----------
grid_full <- tidyr::expand_grid(
  Set = factor("Complete", levels = c("Complete","Triage","Cascade")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_full_like, levels = levels_full_like)
)
grid_tri <- tidyr::expand_grid(
  Set = factor("Triage", levels = c("Complete","Triage","Cascade")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_tri, levels = levels_tri)
)
grid_cas <- tidyr::expand_grid(
  Set = factor("Cascade", levels = c("Complete","Triage","Cascade")),
  Algorithm = factor(algo_order, levels = algo_order),
  Feature_label = factor(levels_full_like, levels = levels_full_like)
)

# ---------- join ΔMCC onto grids ----------
df_full <- dplyr::left_join(
  grid_full,
  imp_full %>% dplyr::mutate(Feature_label = factor(Feature_label, levels = levels_full_like)) %>%
    dplyr::select(Set, Algorithm, Feature_label, Drop_MCC),
  by = c("Set","Algorithm","Feature_label")
)
df_tri <- dplyr::left_join(
  grid_tri,
  imp_tri %>% dplyr::mutate(Feature_label = factor(Feature_label, levels = levels_tri)) %>%
    dplyr::select(Set, Algorithm, Feature_label, Drop_MCC),
  by = c("Set","Algorithm","Feature_label")
)
df_cas <- dplyr::left_join(
  grid_cas,
  imp_cas %>% dplyr::mutate(Feature_label = factor(Feature_label, levels = levels_full_like)) %>%
    dplyr::select(Set, Algorithm, Feature_label, Drop_MCC),
  by = c("Set","Algorithm","Feature_label")
)

df_plot <- dplyr::bind_rows(df_full, df_tri, df_cas)

df_plot <- df_plot %>%
  mutate(Set = forcats::fct_relevel(Set, "Complete","Triage","Cascade"))

# ---------- palette stretch around zero ----------
lim <- stats::quantile(abs(df_plot$Drop_MCC), 0.95, na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- max(abs(df_plot$Drop_MCC), na.rm = TRUE)

# ---------- superscript dagger for triage attributes ----------
tri_labels <- pretty_feature_label(feat_triage)
super_y_labels <- function(labs, tri = tri_labels, mark = "\u2020") {
  lapply(labs, function(lab) if (lab %in% tri) bquote(.(lab)^.(mark)) else bquote(.(lab)))
}

# ---------- plot ----------
p_heat <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Algorithm, y = Feature_label, fill = Drop_MCC)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::scale_fill_gradient2(
    low = "#D55E00", mid = "white", high = "#0072B2",
    midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
    name = expression(Delta*"MCC"), na.value = "white"
  ) +
  ggplot2::facet_grid(rows = vars(Set), scales = "free_y", space = "free_y") +
  ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
  ggplot2::scale_y_discrete(labels = super_y_labels) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 8) +
  ggplot2::theme(
    panel.grid         = ggplot2::element_blank(),
    panel.border       = ggplot2::element_rect(color = "grey60", fill = NA, linewidth = 0.6),
    axis.text.x.top    = ggplot2::element_text(size = 7, face = "plain"),
    axis.text.y        = ggplot2::element_text(size = 7),
    strip.text.y.right = ggplot2::element_text(size = 8, face = "plain"),
    legend.title       = ggplot2::element_text(size = 7, face = "plain"),
    legend.text        = ggplot2::element_text(size = 6)
  ) +
  ggplot2::guides(fill = ggplot2::guide_colorbar(
    barheight = grid::unit(18, "mm"),
    barwidth  = grid::unit(2.5, "mm"),
    ticks = TRUE
  ))

print(p_heat)

ggplot2::ggsave(
  filename = "20251015_features.png",
  plot = p_heat,
  width = 100, height = 140, units = "mm",
  dpi = 300, bg = "white"
)
