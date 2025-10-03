# ==== Two-row-per-metric heatmap, white=0.5 diverging, bold max per row ====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("margin","ggplot2", quiet = TRUE)
}

metric_order <- c("MCC","AUC","F1","Accuracy","Precision","Sensitivity","Specificity")
metric_label <- c(MCC="MCC", AUC="AUC-ROC", F1="F1-score",
                  Accuracy="Accuracy", Precision="Precision",
                  Sensitivity="Sensitivity", Specificity="Specificity")
algo_order_inv <- c("C4.5","k-NN","SVM","RF","LR")

prep_two_row_df <- function(triage_agg, full_agg){
  tri <- triage_agg %>% mutate(FeatureSet = "Triage features")
  ful <- full_agg   %>% mutate(FeatureSet = "Full features")
  
  df <- bind_rows(ful, tri) %>%
    mutate(
      Metric    = factor(Metric, levels = metric_order,
                         labels = metric_label[metric_order]),
      Algorithm = factor(Algorithm, levels = algo_order_inv),
      Label     = sprintf("%.4f \u00B1 %.4f", Mean, SD),
      Row       = paste0(Metric, " – ",
                         ifelse(FeatureSet == "Full features","Full","Triage"))
    )
  
  # order rows as MCC–Full, MCC–Triage, AUC–Full, AUC–Triage, ... then flipped so MCC at top
  metric_lvls <- levels(df$Metric)
  row_levels  <- unlist(lapply(metric_lvls, function(m) c(paste0(m," – Full"),
                                                          paste0(m," – Triage"))))
  df$Row <- factor(df$Row, levels = rev(row_levels))
  
  # mark the max within each Row (one bold per row, ties allowed)
  df %>%
    group_by(Row) %>%
    mutate(IsMax = Mean == max(Mean, na.rm = TRUE) & !is.na(Mean)) %>%
    ungroup()
}

plot_perf_heatmap_midwhite <- function(triage_agg, full_agg,
                                       ratio = 0.35,
                                       limits = c(0,1),
                                       breaks = seq(0, 1, by = 0.2)) {
  dfp <- prep_two_row_df(triage_agg, full_agg) %>%
    group_by(Row) %>% filter(!all(is.na(Mean))) %>% ungroup()
  
  ggplot(dfp, aes(x = Algorithm, y = Row, fill = Mean)) +
    geom_tile(color = "white", width = 0.98, height = 0.98) +
    geom_text(aes(label = Label, fontface = ifelse(IsMax, "bold", "plain")),
              size = 2.6, lineheight = 0.95) +
    scale_fill_gradient2(
      name      = "Mean",
      low       = "#f46d43",   # 0 orange
      mid       = "white",     # 0.5 white
      high      = "#1a9850",   # 1 green
      midpoint  = 0.5,
      limits    = limits,
      breaks    = breaks,
      labels    = number_format(accuracy = 0.1),
      oob       = squish,
      na.value  = "grey90",
      guide     = guide_colorbar(ticks = TRUE)
    ) +
    scale_x_discrete(position = "top") +
    coord_fixed(ratio = ratio, clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_blank(),
      axis.text.x = element_text(hjust = 0.5),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
}

# ---- Usage (active) ----
p <- plot_perf_heatmap_midwhite(triage_agg, full_agg)
print(p)

# assume the functions/objects from before exist:
# prep_two_row_df(), plot_perf_heatmap_midwhite(), triage_agg, full_agg

p <- plot_perf_heatmap_midwhite(triage_agg, full_agg, ratio = 0.35)

# save for Overleaf, two-column friendly size
ggplot2::ggsave(
  filename = "fig_perf_heatmap.png",
  plot     = p,
  width    = 300,   # mm
  height   = 120,   # mm
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)

