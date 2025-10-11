# ===============================================================
# DCA curves (Complete, Triage, Cascade) per algorithm
# Self-contained block: safely defines any missing helpers
# Requires in the environment: df, yvar, cfg, feat_full, feat_triage,
#   run_holdout(), model_specs(), algo_order (optional)
# ===============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(tibble)
})

# ---- defaults for cfg$cascade if missing ----
if (is.null(cfg$cascade)) {
  cfg$cascade <- list(
    enabled   = TRUE, tune_band = TRUE,
    t_low     = 0.20, t_high = 0.80,
    grid_tl   = seq(0.05, 0.45, by = 0.05),
    grid_th   = seq(0.55, 0.95, by = 0.05),
    max_defer = 1.0
  )
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- helper: decision curve table (define if missing) ----
if (!exists("decision_curve_table", mode = "function")) {
  decision_curve_table <- function(obs, p, thresholds = seq(0.01, 0.99, by = 0.01)) {
    y_raw <- as.integer(obs == cfg$pos_label)
    keep  <- is.finite(as.numeric(p)) & !is.na(y_raw)
    y <- y_raw[keep]; p <- as.numeric(p)[keep]
    if (!length(y)) {
      return(tibble::tibble(threshold=thresholds, Net_Benefit=NA_real_,
                            NB_TreatAll=NA_real_, NB_TreatNone=0, prevalence=NA_real_))
    }
    N <- length(y); prev <- mean(y)
    purrr::map_dfr(thresholds, function(pt){
      pred <- as.integer(p >= pt)
      TP <- sum(pred==1 & y==1)
      FP <- sum(pred==1 & y==0)
      NB <- TP/N - FP/N * (pt/(1-pt))
      NB_all <- prev - (1 - prev) * (pt/(1-pt))
      tibble::tibble(threshold=pt, Net_Benefit=NB,
                     NB_TreatAll=NB_all, NB_TreatNone=0, prevalence=prev)
    })
  }
}

# ---- tighter ribbons: STRATIFIED bootstrap (fix prevalence) ----
dca_boot_band_stratified <- function(obj,
                                     B = cfg$B_boot %||% 1000,
                                     thresholds = seq(0.01, 0.99, by = 0.01)) {
  if (is.null(obj)) return(NULL)
  y <- as.integer(obj$obs == cfg$pos_label)
  p <- as.numeric(obj$p)
  
  idx_pos <- which(y == 1L)
  idx_neg <- which(y == 0L)
  P  <- length(idx_pos)
  Nn <- length(idx_neg)
  N  <- P + Nn
  if (P == 0L || Nn == 0L) {
    warning("Stratified bootstrap: one class missing; skipping ribbons for ", obj$set)
    return(NULL)
  }
  
  kvec   <- thresholds / (1 - thresholds)
  nb_mat <- matrix(NA_real_, nrow = B, ncol = length(thresholds))
  
  for (b in seq_len(B)) {
    ip   <- sample(idx_pos, P,  replace = TRUE)   # resample positives
    ineg <- sample(idx_neg, Nn, replace = TRUE)   # resample negatives
    idx  <- c(ip, ineg)
    yb   <- y[idx]
    pb   <- p[idx]
    
    for (j in seq_along(thresholds)) {
      pt <- thresholds[j]
      pred1 <- (pb >= pt)
      TP <- sum(pred1 & (yb == 1L))
      FP <- sum(pred1 & (yb == 0L))
      nb_mat[b, j] <- TP / N - FP / N * kvec[j]
    }
  }
  
  qs <- apply(nb_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  tibble::tibble(
    threshold = thresholds,
    NB_L = qs[1, ],
    NB_U = qs[2, ],
    set  = obj$set,
    algo = obj$algo
  )
}

# ---- helper: tune uncertainty band on TRIAGE OOF (define if missing) ----
if (!exists("tune_band_from_triage_oof", mode = "function")) {
  tune_band_from_triage_oof <- function(oof_df, thr_tri,
                                        grid_tl   = cfg$cascade$grid_tl %||% seq(0.05, 0.45, 0.05),
                                        grid_th   = cfg$cascade$grid_th %||% seq(0.55, 0.95, 0.05),
                                        max_defer = cfg$cascade$max_defer %||% 1) {
    stopifnot(all(c("obs","p") %in% names(oof_df)))
    pos <- cfg$pos_label; neg <- cfg$neg_label
    best <- NULL
    for (tl in grid_tl) for (th in grid_th) if (tl < th) {
      defer <- (oof_df$p > tl & oof_df$p < th)
      defer_rate <- mean(defer)
      if (defer_rate > max_defer) next
      keep <- !defer
      if (!any(keep)) next
      pred_keep <- factor(ifelse(oof_df$p[keep] >= thr_tri, pos, neg), levels = c(neg, pos))
      obs_keep  <- factor(oof_df$obs[keep], levels = c(neg, pos))
      N  <- length(obs_keep)
      y  <- as.integer(obs_keep == pos)
      pr <- as.integer(pred_keep == pos)
      TP <- sum(pr==1 & y==1); FP <- sum(pr==1 & y==0)
      NB <- TP/N - FP/N * (thr_tri/(1-thr_tri))
      cand <- list(t_low = tl, t_high = th, score = NB, defer = defer_rate, n_keep = sum(keep))
      if (is.null(best) ||
          cand$score > best$score + 1e-12 ||
          (abs(cand$score - best$score) < 1e-12 &&
           (cand$defer < best$defer - 1e-12 ||
            (abs(cand$defer - best$defer) < 1e-12 && cand$n_keep > best$n_keep)))) best <- cand
    }
    if (is.null(best)) list(t_low = cfg$cascade$t_low, t_high = cfg$cascade$t_high) else best
  }
}

# ---------- get external predictions for Complete + Triage ----------
algos_avail <- names(model_specs(cfg$tune_len))
get_hold <- function(a, set_feats, set_label) {
  h <- run_holdout(df, yvar, set_feats, a, set_label)
  if (is.null(h)) return(NULL)
  h$algo <- a
  h$set  <- set_label
  h
}

holds_external <- purrr::compact(unlist(lapply(algos_avail, function(a) {
  list(
    full = get_hold(a, feat_full,   "Complete feature set"),
    tri  = get_hold(a, feat_triage, "Triage feature set")
  )
}), recursive = FALSE))

stopifnot(length(holds_external) > 0)

# ---------- build Cascade predictions ----------
by_algo <- split(holds_external, vapply(holds_external, function(x) x$algo, character(1)))
cascade_holds <- list()

for (a in names(by_algo)) {
  hx   <- by_algo[[a]]
  full <- purrr::detect(hx, ~ .x$set == "Complete feature set")
  tri  <- purrr::detect(hx, ~ .x$set == "Triage feature set")
  if (is.null(full) || is.null(tri)) next
  
  if (isTRUE(cfg$cascade$enabled) &&
      isTRUE(cfg$cascade$tune_band) &&
      !is.null(tri$oof) && nrow(tri$oof) > 0 && is.finite(tri$threshold)) {
    tuned <- tune_band_from_triage_oof(tri$oof, thr_tri = tri$threshold)
    t_low  <- tuned$t_low  %||% cfg$cascade$t_low
    t_high <- tuned$t_high %||% cfg$cascade$t_high
  } else {
    t_low  <- cfg$cascade$t_low
    t_high <- cfg$cascade$t_high
  }
  
  p_tri    <- as.numeric(tri$p)
  p_full   <- as.numeric(full$p)
  stopifnot(length(p_tri) == length(p_full))
  use_full <- (p_tri > t_low & p_tri < t_high)
  p_cas    <- ifelse(use_full, p_full, p_tri)
  
  cascade_holds[[length(cascade_holds) + 1]] <- list(
    obs  = tri$obs,
    p    = p_cas,
    algo = a,
    set  = "Cascade"
  )
}

holds_external <- c(holds_external, cascade_holds)

# ---------- build DCA tables + stratified bootstrap ribbons ----------
by_algo <- split(holds_external, vapply(holds_external, function(x) x$algo, character(1)))
combo <- lapply(names(by_algo), function(a) {
  hx <- by_algo[[a]]
  full <- purrr::detect(hx, ~ .x$set == "Complete feature set")
  tri  <- purrr::detect(hx, ~ .x$set == "Triage feature set")
  cas  <- purrr::detect(hx, ~ .x$set == "Cascade")
  
  d_full <- if (!is.null(full)) decision_curve_table(full$obs, full$p) %>% dplyr::mutate(set = full$set, algo = a) else NULL
  d_tri  <- if (!is.null(tri))  decision_curve_table(tri$obs,  tri$p)  %>% dplyr::mutate(set = tri$set,  algo = a) else NULL
  d_cas  <- if (!is.null(cas))  decision_curve_table(cas$obs,  cas$p)  %>% dplyr::mutate(set = cas$set,  algo = a) else NULL
  
  ref <- if (!is.null(full)) full else if (!is.null(tri)) tri else cas
  t_df <- decision_curve_table(ref$obs, ref$p) %>% dplyr::transmute(algo = a, threshold, NB_TreatAll, NB_TreatNone)
  
  list(dca = dplyr::bind_rows(d_full, d_tri, d_cas), treat = t_df)
})

dca_df   <- dplyr::bind_rows(purrr::map(combo, "dca"))
treat_df <- dplyr::bind_rows(purrr::map(combo, "treat"))
bands_df <- holds_external %>%
  purrr::map(dca_boot_band_stratified) %>%
  purrr::compact() %>%
  { if (!length(.)) NULL else dplyr::bind_rows(.) }

# ---------- factor levels and palettes ----------
if (exists("algo_order") && length(algo_order) > 0) {
  dca_df$algo   <- factor(as.character(dca_df$algo),   levels = algo_order, ordered = TRUE)
  treat_df$algo <- factor(as.character(treat_df$algo), levels = algo_order, ordered = TRUE)
  if (!is.null(bands_df)) bands_df$algo <- factor(as.character(bands_df$algo), levels = algo_order, ordered = TRUE)
}
dca_df$set <- factor(dca_df$set, levels = c("Complete feature set","Triage feature set","Cascade"))
if (!is.null(bands_df)) bands_df$set <- factor(bands_df$set, levels = levels(dca_df$set))

# color-blind friendly
cb_orange <- "#D55E00"  # Complete
cb_blue   <- "#0072B2"  # Triage
cb_green  <- "#009E73"  # Cascade

# ---------- panel builder ----------
y_cap <- max(c(dca_df$Net_Benefit, treat_df$NB_TreatAll %||% 0, 0), na.rm = TRUE)
if (!is.finite(y_cap) || y_cap <= 0) y_cap <- 0.05
y0 <- -0.02

axis_breaks_y <- function(y_cap) if (y_cap <= 0.06) seq(0, y_cap, by = 0.01) else scales::breaks_extended()(c(0, y_cap))
theme_ref <- function(base_size = 15) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.line          = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "plain", color = "grey20"),
      legend.title       = ggplot2::element_blank()
    )
}
safe_ribbon <- function(df_b) {
  if (is.null(df_b) || !nrow(df_b)) return(list())
  df_b <- df_b %>% dplyr::filter(is.finite(NB_L), is.finite(NB_U), !is.na(set))
  layers <- list()
  for (s in c("Complete feature set","Triage feature set","Cascade")) {
    d <- dplyr::filter(df_b, set == s)
    if (!nrow(d)) next
    fill_col <- if (s == "Complete feature set") cb_orange else if (s == "Triage feature set") cb_blue else cb_green
    layers[[length(layers) + 1]] <- ggplot2::geom_ribbon(
      data = d, ggplot2::aes(threshold, ymin = NB_L, ymax = NB_U),
      fill = fill_col, alpha = 0.18, inherit.aes = FALSE
    )
  }
  layers
}

dca_panel <- function(a, show_y = FALSE) {
  df_d <- dplyr::filter(dca_df, algo == a, is.finite(Net_Benefit))
  df_t <- dplyr::filter(treat_df, algo == a, is.finite(NB_TreatAll), is.finite(NB_TreatNone))
  df_b <- if (!is.null(bands_df)) dplyr::filter(bands_df, algo == a) else NULL
  
  ggplot2::ggplot() +
    safe_ribbon(df_b) +
    ggplot2::geom_line(data = df_t, ggplot2::aes(threshold, NB_TreatAll, linetype = "Treat-all"),
                       linewidth = 0.6, color = "black", na.rm = TRUE) +
    ggplot2::geom_line(data = df_t, ggplot2::aes(threshold, NB_TreatNone, linetype = "Treat-none"),
                       linewidth = 0.6, color = "grey40", na.rm = TRUE) +
    ggplot2::geom_line(data = df_d, ggplot2::aes(threshold, Net_Benefit, color = set),
                       linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0, y = y0, yend = y_cap),
                          inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = y0, yend = y0),
                          inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    ggplot2::scale_color_manual(values = c("Complete feature set" = cb_orange,
                                           "Triage feature set"   = cb_blue,
                                           "Cascade"              = cb_green),
                                name = "") +
    ggplot2::scale_linetype_manual(values = c("Treat-all" = "solid", "Treat-none" = "dashed"), name = "") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    ggplot2::scale_y_continuous(breaks = axis_breaks_y(y_cap), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::coord_cartesian(ylim = c(y0, y_cap)) +
    ggplot2::labs(x = "Threshold probability", y = if (show_y) "Net benefit" else NULL, title = a) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA)) +
    theme_ref()
}

# ---------- build panels & legend ----------
algos_in_data <- levels(dca_df$algo) %||% unique(as.character(dca_df$algo))
panels <- lapply(seq_along(algos_in_data), function(i) dca_panel(algos_in_data[i], show_y = i %% 3 == 1))
panels <- lapply(panels, function(p) p + theme(
  legend.position = "none",
  axis.title.x = element_blank(),
  plot.title = element_text(size = 11),
  axis.title = element_text(size = 9),
  axis.text  = element_text(size = 8)
))

legend_levels <- c("Treat-all","Treat-none","Complete feature set","Triage feature set","Cascade")
leg_df <- data.frame(legend = factor(legend_levels, levels = legend_levels), y = seq_along(legend_levels))
cb_orange <- "#D55E00"; cb_blue <- "#0072B2"; cb_green <- "#009E73"

cb_leg <- ggplot(leg_df) +
  geom_segment(aes(x = 0, xend = 1, y = y, yend = y, color = legend, linetype = legend), linewidth = 1.1) +
  scale_color_manual(values = c("Treat-all" = "black",
                                "Treat-none" = "grey40",
                                "Complete feature set" = cb_orange,
                                "Triage feature set"   = cb_blue,
                                "Cascade"              = cb_green),
                     breaks = legend_levels, name = NULL) +
  scale_linetype_manual(values = c("Treat-all" = "solid",
                                   "Treat-none" = "dashed",
                                   "Complete feature set" = "solid",
                                   "Triage feature set"   = "solid",
                                   "Cascade"              = "solid"),
                        breaks = legend_levels, name = NULL) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.height = grid::unit(3, "mm"),
        legend.key.width  = grid::unit(10, "mm"))
leg_plot <- cowplot::ggdraw(cowplot::get_legend(cb_leg))

blank <- ggplot() + theme_void()
fill_panel <- function(i) if (i <= length(panels)) panels[[i]] else blank
p1 <- fill_panel(1)
p2 <- fill_panel(2)
p3 <- fill_panel(3)
p4 <- fill_panel(4) + theme(axis.title.x = element_text())
p5 <- fill_panel(5) + theme(axis.title.x = element_text())

p_dca_grid <- (p1 | p2 | p3) / (p4 | p5 | leg_plot)

print(p_dca_grid)
ggsave("external_dca_complete_triage_cascade.png", p_dca_grid,
       width = 250, height = 120, units = "mm", dpi = 300, bg = "white")
