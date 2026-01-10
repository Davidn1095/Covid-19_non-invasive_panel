# ============================================================
# 16b) PARSIMONY CURVES (DEV OOF MCC), consensus-ordered features
# - Uses consensus feature order from feature_importance_perm_dMCC
# - Evaluates on DEVELOPMENT OOF only (no repeated EXTERNAL use)
# - Produces:
#   parsimony_tbl
#   parsimony_plot
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --------
# Pipeline context
# --------
source(file.path("R", "common_utils.R"))
source(file.path("R", "cascade_pipeline.R"))
source(file.path("R", "feature_importance_utils.R"))

ctx <- build_cascade_context()

df_dev <- ctx$df_dev
folds <- ctx$folds
algos <- ctx$algos
dev_cache <- ctx$dev_cache
clip_q <- ctx$clip_q
standardise <- ctx$standardise
use_weights <- ctx$use_weights
filter_rate <- ctx$filter_rate
calibration <- ctx$calibration
thr_grid <- ctx$thr_grid
sens_min_target <- ctx$sens_min_target
spec_min_target <- ctx$spec_min_target
max_def_rate_target <- ctx$max_def_rate_target
require_sens_ge_spec <- ctx$require_sens_ge_spec

fi_out <- build_feature_importance_perm_dMCC(
  df_dev = df_dev,
  df_ext = ctx$df_ext,
  dev_cache = dev_cache,
  dev_params_tbl = ctx$dev_params_tbl,
  algos = algos,
  clip_q = clip_q,
  standardise = standardise,
  use_weights = use_weights,
  filter_rate = filter_rate,
  calibration = calibration
)

feature_importance_perm_dMCC <- fi_out$feature_importance_perm_dMCC

fi_src <- feature_importance_perm_dMCC

policy_col <- if ("Policy" %in% names(fi_src)) "Policy" else if ("Panel" %in% names(fi_src)) "Panel" else NA_character_
stopifnot(!is.na(policy_col))
stopifnot(all(c("Algorithm","Feature") %in% names(fi_src)))

# --------
# Build consensus order per policy, then a global order for additions
# - Non-invasive first, then the remaining lab-augmented features
# --------
get_consensus_order <- function(fi, policy_name) {
  fi2 <- fi %>%
    dplyr::filter(.data[[policy_col]] == policy_name) %>%
    dplyr::mutate(
      Algorithm = as.character(.data$Algorithm),
      Feature = as.character(.data$Feature)
    )
  
  if ("ConsensusRank" %in% names(fi2)) {
    cons <- fi2 %>%
      dplyr::group_by(.data$Feature) %>%
      dplyr::summarise(
        ConsensusRank = suppressWarnings(min(.data$ConsensusRank, na.rm = TRUE)),
        ConsensusScore = if ("ConsensusScore" %in% names(fi2)) suppressWarnings(min(.data$ConsensusScore, na.rm = TRUE)) else NA_real_,
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$ConsensusRank, .data$ConsensusScore, .data$Feature)
    return(cons$Feature)
  }
  
  # fallback if ConsensusRank not present: mean per-algorithm Rank
  stopifnot("Rank" %in% names(fi2))
  cons <- fi2 %>%
    dplyr::group_by(.data$Feature) %>%
    dplyr::summarise(ConsensusScore = mean(.data$Rank, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(ConsensusRank = dplyr::dense_rank(.data$ConsensusScore)) %>%
    dplyr::arrange(.data$ConsensusRank, .data$ConsensusScore, .data$Feature)
  cons$Feature
}

order_noninv <- get_consensus_order(fi_src, "Non-invasive")
order_labaug <- get_consensus_order(fi_src, "Laboratory augmented")

# Global order used for k additions: all non-invasive first, then remaining lab-augmented
order_full <- c(order_noninv, setdiff(order_labaug, order_noninv))

K_noninv <- length(order_noninv)
K_full   <- length(order_full)
stopifnot(K_noninv >= 1L, K_full >= 1L)

# --------
# Pretty x labels (feature-addition axis)
# --------
labels_full <- vapply(order_full, pretty_feature_label, character(1), style = "parsimony")
labels_full <- if (length(labels_full) >= 1L) c(labels_full[1], paste0("+", labels_full[-1])) else character(0)

# --------
# OOF memoisation to avoid recomputing the same subset repeatedly
# --------
.oof_cache <- new.env(parent = emptyenv())

oof_key <- function(algo, stage, feats) {
  paste0(as.character(algo), "|", as.character(stage), "|", paste(feats, collapse = ","))
}

get_oof <- function(algo, stage, feats) {
  feats <- unique(as.character(feats))
  key <- oof_key(algo, stage, feats)
  if (exists(key, envir = .oof_cache, inherits = FALSE)) return(get(key, envir = .oof_cache))
  
  res <- try(
    run_oof_single_panel(
      df_dev = df_dev, algo = algo, predictors = feats,
      clip_q = clip_q, standardise = standardise, use_weights = use_weights,
      filter_rate = filter_rate, calibration = calibration,
      folds = folds
    ),
    silent = TRUE
  )
  if (inherits(res, "try-error")) res <- NULL
  assign(key, res, envir = .oof_cache)
  res
}

# --------
# Parsimony evaluation on DEV OOF only
# --------
subset_feats <- function(order_vec, k, allowed) {
  intersect(order_vec[seq_len(k)], allowed)
}

par_rows <- list()
kk <- 0L

for (algo in algos) {
  allowed1 <- dev_cache[[algo]]$preds1
  allowed2 <- dev_cache[[algo]]$preds2
  
  # Non-invasive: k = 1..K_noninv
  for (k in seq_len(K_noninv)) {
    feats1 <- subset_feats(order_noninv, k, allowed1)
    if (!length(feats1)) next
    
    res1 <- get_oof(algo, "stage1", feats1)
    if (is.null(res1)) next
    
    thr1_sel <- best_thr_mcc_at_sens(
      y01 = res1$y01, p = res1$p_oof, thr_grid = thr_grid,
      sens_min = sens_min_target, spec_min = spec_min_target
    )
    met1 <- compute_metrics_at_thr(res1$y01, res1$p_oof, thr1_sel$thr)
    
    kk <- kk + 1L
    par_rows[[kk]] <- tibble::tibble(
      Algorithm = algo,
      Policy = "Non-invasive",
      k = k,
      n_used_stage1 = length(feats1),
      n_used_stage2 = NA_integer_,
      tau_low = NA_real_, tau_high = NA_real_,
      thr = as.numeric(thr1_sel$thr),
      Deferral = 0.0,
      MCC = as.numeric(met1$MCC),
      AUC = as.numeric(met1$AUC)
    )
  }
  
  # Laboratory augmented + Cascade: k = 1..K_full
  for (k in seq_len(K_full)) {
    feats2 <- subset_feats(order_full, k, allowed2)
    if (!length(feats2)) next
    
    # Lab-augmented only
    res2 <- get_oof(algo, "stage2", feats2)
    if (!is.null(res2)) {
      thr2_sel <- best_thr_mcc_at_sens(
        y01 = res2$y01, p = res2$p_oof, thr_grid = thr_grid,
        sens_min = sens_min_target, spec_min = spec_min_target
      )
      met2 <- compute_metrics_at_thr(res2$y01, res2$p_oof, thr2_sel$thr)
      
      kk <- kk + 1L
      par_rows[[kk]] <- tibble::tibble(
        Algorithm = algo,
        Policy = "Laboratory augmented",
        k = k,
        n_used_stage1 = NA_integer_,
        n_used_stage2 = length(feats2),
        tau_low = NA_real_, tau_high = NA_real_,
        thr = as.numeric(thr2_sel$thr),
        Deferral = 0.0,
        MCC = as.numeric(met2$MCC),
        AUC = as.numeric(met2$AUC)
      )
    }
    
    # Cascade at this k: stage1 uses up to min(k, K_noninv) from the same global prefix
    k1 <- min(k, K_noninv)
    feats1k <- subset_feats(order_noninv, k1, allowed1)
    if (!length(feats1k) || is.null(res2)) next
    
    res1k <- get_oof(algo, "stage1", feats1k)
    if (is.null(res1k)) next
    
    sel <- best_cascade_params(
      y01 = res1k$y01, p1 = res1k$p_oof, p2 = res2$p_oof,
      tau_grid = thr_grid, thr2_grid = thr_grid,
      sens_min = sens_min_target, spec_min = spec_min_target,
      max_def_rate = max_def_rate_target,
      require_sens_ge_spec = require_sens_ge_spec
    )
    
    cas <- cascade_apply(res1k$p_oof, res2$p_oof, sel$tau_low, sel$tau_high, sel$thr2)
    metc <- compute_metrics_from_pred(res1k$y01, cas$pred01, score = cas$score)
    def_rate <- mean(cas$defer, na.rm = TRUE)
    
    kk <- kk + 1L
    par_rows[[kk]] <- tibble::tibble(
      Algorithm = algo,
      Policy = "Cascade",
      k = k,
      n_used_stage1 = length(feats1k),
      n_used_stage2 = length(feats2),
      tau_low = as.numeric(sel$tau_low),
      tau_high = as.numeric(sel$tau_high),
      thr = as.numeric(sel$thr2),
      Deferral = as.numeric(def_rate),
      MCC = as.numeric(metc$MCC),
      AUC = as.numeric(metc$AUC)
    )
  }
}

parsimony_tbl <- dplyr::bind_rows(par_rows) %>%
  dplyr::mutate(
    Policy = factor(.data$Policy, levels = c("Non-invasive","Laboratory augmented","Cascade")),
    Algorithm = factor(.data$Algorithm, levels = algos)
  ) %>%
  dplyr::arrange(.data$Algorithm, .data$Policy, .data$k)

par_plot_df <- parsimony_tbl

print(parsimony_tbl, n = Inf, width = Inf)

# ============================================================
# PLOT (FINAL ITERATION ONLY)
# - algo order: C4.5 ... LR
# - x tick labels in ALL panels, x title only on LR
# - y title only on C4.5, SVMRBF, LR
# - y range: 0 to 0.5
# - black x/y axis lines
# - legend order: Non-invasive, Laboratory augmented, Cascade
# - draw order (top): Non-invasive over Cascade over Laboratory augmented
# - printed only, no file writing
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

algo_levels_plot <- rev(algos)  # C4.5 ... LR

policy_levels_plot <- c("Non-invasive", "Laboratory augmented", "Cascade")
pal <- c(
  "Non-invasive"         = "#0072B2",
  "Laboratory augmented" = "#D55E00",
  "Cascade"              = "#009E73"
)

# ensure correct factor levels for facet order + legend order
par_plot_df2 <- par_plot_df %>%
  dplyr::mutate(
    Algorithm = factor(as.character(.data$Algorithm), levels = algo_levels_plot),
    Policy    = factor(as.character(.data$Policy), levels = policy_levels_plot)
  )

theme_axes_black <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      
      axis.line  = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = grid::unit(2, "pt"),
      
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7, lineheight = 0.9),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      
      legend.title = ggplot2::element_blank()
    )
}

make_alg_plot_final <- function(alg) {
  dat <- dplyr::filter(par_plot_df2, .data$Algorithm == alg) %>%
    dplyr::mutate(
      .draw_order = dplyr::case_when(
        as.character(.data$Policy) == "Laboratory augmented" ~ 1L,
        as.character(.data$Policy) == "Cascade"              ~ 2L,
        as.character(.data$Policy) == "Non-invasive"         ~ 3L,
        TRUE ~ 99L
      )
    ) %>%
    dplyr::arrange(.data$.draw_order, .data$k)
  
  show_x_title <- identical(as.character(alg), "LR")
  show_y_title <- as.character(alg) %in% c("C4.5", "SVMRBF", "LR")
  
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$k, y = .data$MCC, color = .data$Policy, group = .data$Policy)) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(K_full),
      labels = labels_full[seq_len(K_full)],
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 0.5),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_color_manual(
      values = pal,
      breaks = policy_levels_plot,
      limits = policy_levels_plot
    ) +
    ggplot2::labs(
      title = as.character(alg),
      x = if (show_x_title) "Added attribute" else NULL,
      y = if (show_y_title) "MCC" else NULL,
      color = NULL
    ) +
    theme_axes_black()
}

plots <- lapply(algo_levels_plot, make_alg_plot_final)

parsimony_plot <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") &
  ggplot2::theme(legend.position = "right")

print(parsimony_plot)

# Objects returned:
# parsimony_tbl
# parsimony_plot
# order_noninv, order_full
