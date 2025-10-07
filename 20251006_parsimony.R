# ===============================================================
# PARSIMONY CURVES (MCC) — evaluated on the EXTERNAL set
# Reuses: df, yvar, feat_full, feat_triage, algo_order,
#         run_holdout(), compute_metrics_binary(), bootstrap_ci_all_metrics()
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
})

# --- deterministic RNG (same as your pipeline) ---
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# --- feature ordering: triage-first, reverse alphabetical inside blocks ---
tri_in_df  <- intersect(feat_triage, names(df))
tri_order  <- sort(tri_in_df, decreasing = TRUE)
full_in_df <- intersect(feat_full, names(df))
rest_order <- sort(setdiff(full_in_df, tri_order), decreasing = TRUE)

order_by_set <- list(
  "Complete feature set" = c(tri_order, rest_order),
  "Triage feature set"   = tri_order
)

# --- one external parsimony curve (k = 1..p), using your run_holdout() ---
parsimony_external_curve <- function(set_name, feats_order, algo) {
  if (!length(feats_order)) {
    return(tibble(k = integer(), Score = numeric(), SD = numeric(),
                  Set = character(), Algorithm = character()))
  }
  purrr::map_dfr(seq_along(feats_order), function(k) {
    feats_k <- feats_order[seq_len(k)]
    # run_holdout() trains with time-blocked CV, ENN+ADASYN, calibrates, picks τ, and predicts on external
    hold <- run_holdout(df, yvar, feats_k, algo,
                        ifelse(set_name == "Complete feature set", "Complete", "Triage"))
    if (is.null(hold)) {
      tibble(k = k, Score = NA_real_, SD = NA_real_,
             Set = set_name, Algorithm = algo)
    } else {
      mets <- compute_metrics_binary(hold$obs, hold$pred, hold$p)
      ci   <- bootstrap_ci_all_metrics(hold$obs, hold$pred, hold$p)  # sd for ribbon
      tibble(k = k,
             Score = as.numeric(mets["MCC"]),
             SD    = as.numeric(ci$sd["MCC"]),
             Set = set_name,
             Algorithm = algo)
    }
  })
}

# --- build all curves (algorithms in your declared order) ---
res_all <- purrr::imap_dfr(order_by_set, function(ord, set_name) {
  purrr::map_dfr(algo_order, function(algo) parsimony_external_curve(set_name, ord, algo))
}) %>%
  mutate(Algorithm = factor(Algorithm, levels = algo_order))

# --- simple publication theme ---
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "grey30"),
      legend.position = "right",
      legend.title = element_blank()
    )
}

# --- panel builder (one algorithm) ---
make_alg_plot <- function(alg, show_y = FALSE, show_x = FALSE) {
  dat <- dplyr::filter(res_all, Algorithm == alg)
  if (!nrow(dat) || all(is.na(dat$Score))) {
    return(ggplot() + theme_void() + ggtitle(alg))
  }
  # in-panel axes
  x_min <- 0.5
  x_max <- max(dat$k, na.rm = TRUE) + 0.5
  y_min <- 0
  y_max <- 1
  
  ggplot(dat, aes(k, Score, color = Set)) +
    annotate("segment", x = x_min, xend = x_max, y = y_min, yend = y_min,
             linewidth = 0.6, colour = "black") +
    annotate("segment", x = x_min, xend = x_min, y = y_min, yend = y_max,
             linewidth = 0.6, colour = "black") +
    geom_ribbon(aes(ymin = pmax(0, Score - SD), ymax = pmin(1, Score + SD),
                    fill = Set, group = Set),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = function(x) seq_len(max(x, na.rm = TRUE)),
                       expand = expansion(mult = c(0.12, 0.03))) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    labs(x = if (show_x) "Number of attributes" else NULL,
         y = if (show_y) "MCC" else NULL,
         title = as.character(alg)) +
    theme_pub() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      
    )
}

# --- draw 2×3 grid (legend bottom-right) ---
p1 <- make_alg_plot("C4.5", show_y = TRUE,  show_x = FALSE)
p2 <- make_alg_plot("k-NN", show_y = FALSE, show_x = FALSE)
p3 <- make_alg_plot("SVM",  show_y = FALSE, show_x = FALSE)
p4 <- make_alg_plot("RF",   show_y = TRUE,  show_x = TRUE)
p5 <- make_alg_plot("LR",   show_y = FALSE, show_x = TRUE)

legend_src <- ggplot(res_all, aes(k, Score, color = Set)) +
  geom_line() + theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")
pleg <- patchwork::wrap_elements(cowplot::get_legend(legend_src))

p_grid <- (p1 + p2 + p3) /
  (p4 + p5 + pleg) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1))

print(p_grid)
