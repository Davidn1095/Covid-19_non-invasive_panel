# ==== Top-5 permutation importance (panel-specific), 2 rows × 5 columns ====
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(cowplot)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("margin","ggplot2", quiet = TRUE)
}

# expect: imp_full, imp_triage with cols: Algorithm, Feature, MeanDrop, SD
pretty_names <- c(
  Diagnosis                = "Diagnosis",
  severity_admission       = "Severity at admission",
  Age                      = "Age",
  Gender                   = "Sex",
  SpO2_admission           = "SpO2 at admission",
  albumin                  = "Albumin",
  CRP                      = "C-reactive protein",
  D_Dimer                  = "D-dimer",
  monocyte_abs_number      = "Monocyte count",
  neutrophil_abs_number    = "Neutrophil count",
  lymphocyte_abs_number    = "Lymphocyte count",
  monocytes_perc           = "Monocytes (%)",
  neutrophils_perc         = "Neutrophils (%)",
  lymphocytes_perc         = "Lymphocytes (%)"
)

imp_full$FeatureSet   <- "Full"
imp_triage$FeatureSet <- "Triage"

imp_all <- bind_rows(imp_full, imp_triage) %>%
  mutate(
    Feature    = dplyr::recode(Feature, !!!pretty_names),
    Algorithm  = factor(Algorithm,  levels = c("C4.5","k-NN","SVM","RF","LR")),
    FeatureSet = factor(FeatureSet, levels = c("Full","Triage"))
  )

# one panel: local top-5, panel-specific x-limits, titles only for top row
panel_plot <- function(data, feature_set, algo, show_title = TRUE){
  d <- data %>%
    filter(FeatureSet == feature_set, Algorithm == algo) %>%
    filter(is.finite(MeanDrop), is.finite(SD)) %>%
    arrange(desc(MeanDrop)) %>%
    slice_head(n = 5) %>%                                # <-- constant 5 fixes your error
    mutate(Feature = factor(Feature, levels = rev(Feature)))  # top on top
  
  if (nrow(d) == 0) return(ggplot() + theme_void())
  
  xmax <- max(d$MeanDrop + d$SD, na.rm = TRUE) * 1.12
  
  ggplot(d, aes(x = MeanDrop, y = Feature)) +
    geom_col(fill = "grey65") +
    geom_errorbarh(aes(xmin = pmax(0, MeanDrop - SD), xmax = MeanDrop + SD), height = 0.18) +
    geom_text(aes(label = sprintf("%.3f", MeanDrop)), hjust = -0.10, size = 3) +
    scale_x_continuous(limits = c(0, xmax), breaks = scales::pretty_breaks()) +
    labs(x = "Mean AUC drop when permuted", y = NULL,
         title = if (show_title) as.character(algo) else NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 2)),
      axis.text.y = element_text(size = 9)
    )
}

algos <- levels(imp_all$Algorithm)

# build two rows with the same column order
top_row    <- lapply(algos, function(a) panel_plot(imp_all, "Full",   a, show_title = TRUE))
bottom_row <- lapply(algos, function(a) panel_plot(imp_all, "Triage", a, show_title = FALSE))

# combine strictly as 2 rows × 5 columns
p_imp_2x5 <- cowplot::plot_grid(
  plotlist = c(top_row, bottom_row),
  nrow = 2, ncol = 5, align = "hv"
)

print(p_imp_2x5)

# save PNG for Overleaf
ggplot2::ggsave(
  filename = "fig_importance_top5.png",
  plot     = p_imp_2x5,
  width    = 330,   # mm, two-column friendly
  height   = 120,   # mm
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)

