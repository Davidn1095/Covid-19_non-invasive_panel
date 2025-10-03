# ===== Calibration curves (binned, per-fold weighted), print only =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("margin","ggplot2", quiet = TRUE)
}

# pred_all must exist with: Algorithm, FeatureSet, Fold, obs (factor "No","Yes"), prob_pos (0..1)

# --- Build calibration data: equal-frequency bins per Algorithm×FeatureSet,
#     compute per-fold bin stats, then average weighting by bin size
calibration_df <- function(pred_all, nbins = 10){
  stopifnot(all(c("Algorithm","FeatureSet","Fold","obs","prob_pos") %in% names(pred_all)))
  pred_all %>%
    mutate(y = as.integer(obs == levels(obs)[2])) %>%
    group_by(Algorithm, FeatureSet) %>%
    group_modify(function(d, key){
      brk <- quantile(d$prob_pos, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
      brk <- unique(brk)
      if (length(brk) < 2) brk <- seq(0, 1, length.out = nbins + 1)
      d %>%
        mutate(bin = cut(prob_pos, breaks = brk, include.lowest = TRUE, labels = FALSE)) %>%
        group_by(Fold, bin, .add = TRUE) %>%
        summarise(
          p_hat = mean(prob_pos, na.rm = TRUE),
          y_bar = mean(y, na.rm = TRUE),
          n     = dplyr::n(),
          .groups = "drop"
        ) %>%
        group_by(bin) %>%
        summarise(
          p_hat = weighted.mean(p_hat, w = n, na.rm = TRUE),
          y_bar = weighted.mean(y_bar, w = n, na.rm = TRUE),
          n     = sum(n),
          .groups = "drop"
        ) %>%
        arrange(p_hat)
    }) %>%
    ungroup()
}

cal_df <- calibration_df(pred_all, nbins = 10)

# panel order
algo_order_inv <- c("C4.5","k-NN","SVM","RF","LR")
cal_df$Algorithm <- factor(cal_df$Algorithm, levels = algo_order_inv)

# --- Plot ---
p_cal <- ggplot(cal_df, aes(x = p_hat, y = y_bar, color = FeatureSet, group = FeatureSet)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8, stroke = 0.2) +
  facet_wrap(~ Algorithm, nrow = 1) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_color_manual(values = c("Full" = "#1f78b4", "Triage" = "#33a02c"), name = NULL) +
  labs(x = "Predicted probability", y = "Observed frequency") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text = element_text(color = "black")
  )

print(p_cal)

# If later you want to save:
ggplot2::ggsave("fig_calibration_panels.png", p_cal, width = 190, height = 70, units = "mm", dpi = 300, bg = "white")
