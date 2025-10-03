# ===== DCA, panels by Algorithm, x from 0..1, minimal < 0 on y, print only =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("margin","ggplot2", quiet = TRUE)
}

# pred_all must have: Algorithm, FeatureSet, Fold, obs (factor No/Yes), prob_pos (0..1)

# thresholds computed on a stable range, axis later shown as 0..1
ths <- seq(0.05, 0.80, by = 0.01)

nb_curve <- function(obs, prob, thresholds){
  y <- as.integer(obs == levels(obs)[2])
  N <- length(y)
  lapply(thresholds, function(t){
    pred_pos <- prob >= t
    TP <- sum(pred_pos & y == 1)
    FP <- sum(pred_pos & y == 0)
    NB <- (TP / N) - (FP / N) * t / (1 - t)
    c(t = t, NB = NB)
  }) |>
    do.call(rbind, args = _) |>
    as.data.frame() |>
    dplyr::mutate(t = as.numeric(t), NB = as.numeric(NB))
}

# model curves, averaged across outer folds
dca_model <- pred_all |>
  dplyr::group_by(Algorithm, FeatureSet, Fold) |>
  dplyr::group_modify(function(d, key) nb_curve(d$obs, d$prob_pos, ths)) |>
  dplyr::ungroup() |>
  dplyr::group_by(Algorithm, FeatureSet, t) |>
  dplyr::summarise(NB_mean = mean(NB, na.rm = TRUE), .groups = "drop")

# references, averaged across folds
dca_ref <- pred_all |>
  dplyr::group_by(Algorithm, Fold) |>
  dplyr::group_modify(function(d, key){
    y <- as.integer(d$obs == levels(d$obs)[2]); prev <- mean(y)
    tibble::tibble(t = ths,
                   NB_all  = prev - (1 - prev) * ths / (1 - ths),
                   NB_none = 0)
  }) |>
  dplyr::ungroup() |>
  dplyr::group_by(Algorithm, t) |>
  dplyr::summarise(NB_all = mean(NB_all, na.rm = TRUE),
                   NB_none = 0, .groups = "drop")

# order of panels
algo_order_inv <- c("C4.5","k-NN","SVM","RF","LR")
dca_model$Algorithm <- factor(dca_model$Algorithm, levels = algo_order_inv)
dca_ref$Algorithm   <- factor(dca_ref$Algorithm,   levels = algo_order_inv)

# dynamic top so the headroom is kept, floor near zero to avoid deep negatives
y_min <- -0.02
y_max <- max(dca_model$NB_mean, dca_ref$NB_all, na.rm = TRUE) + 0.02

p_dca <- ggplot() +
  geom_line(data = dca_ref, aes(t, NB_all),  color = "grey60", linetype = "dashed") +
  geom_line(data = dca_ref, aes(t, NB_none), color = "grey70") +
  geom_line(data = dca_model, aes(t, NB_mean, color = FeatureSet), linewidth = 0.9) +
  facet_wrap(~ Algorithm, nrow = 1) +
  scale_color_manual(values = c("Full" = "#1f78b4", "Triage" = "#33a02c"), name = NULL) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(x = "Threshold probability", y = "Net benefit") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text = element_text(color = "black")
  )

print(p_dca)

# save PNG for Overleaf
ggplot2::ggsave(
  filename = "fig_dca_panels.png",
  plot     = p_dca,
  width    = 190,   # mm, two-column friendly
  height   = 70,    # mm
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)
