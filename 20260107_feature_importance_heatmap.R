# ============================================================
# 16) FEATURE IMPORTANCE: permutation Drop in MCC on EXTERNAL
# Paste after you print external_policy_tbl
#
# Creates:
#   feature_importance_perm_dMCC
#   feature_importance_perm_dMCC_wide
#   feature_importance_perm_dMCC_plotdata
#   p_fi_heatmap
#
# Adds:
#   ConsensusScore, ConsensusRank across algorithms per Policy, Feature
# Plot:
#   columns: C4.5, kNN, SVMRBF, RF, LR, Consensus
#   non-invasive features marked with † in row labels
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
})

source(file.path("R", "common_utils.R"))
source(file.path("R", "cascade_pipeline.R"))
source(file.path("R", "feature_importance_utils.R"))

ctx <- build_cascade_context()

df_dev <- ctx$df_dev
df_ext <- ctx$df_ext
dev_cache <- ctx$dev_cache
dev_params_tbl <- ctx$dev_params_tbl
algos <- ctx$algos
clip_q <- ctx$clip_q
standardise <- ctx$standardise
use_weights <- ctx$use_weights
filter_rate <- ctx$filter_rate
calibration <- ctx$calibration

fi_out <- build_feature_importance_perm_dMCC(
  df_dev = df_dev,
  df_ext = df_ext,
  dev_cache = dev_cache,
  dev_params_tbl = dev_params_tbl,
  algos = algos,
  clip_q = clip_q,
  standardise = standardise,
  use_weights = use_weights,
  filter_rate = filter_rate,
  calibration = calibration
)

feature_importance_perm_dMCC <- fi_out$feature_importance_perm_dMCC
feature_importance_perm_dMCC_wide <- fi_out$feature_importance_perm_dMCC_wide
consensus_tbl <- fi_out$consensus_tbl
feats_noninv <- fi_out$feats_noninv
feats_labaug <- fi_out$feats_labaug
feats_cascade <- fi_out$feats_cascade

# --------
# Plot data with Consensus as 6th column
# --------
lab_tbl <- consensus_tbl %>%
  dplyr::mutate(
    Feature_label_base = pretty_feature_label(.data$Feature),
    Feature_label = dplyr::if_else(.data$Feature %in% feats_noninv,
                                   paste0(.data$Feature_label_base, " \u2020"),
                                   .data$Feature_label_base)
  ) %>%
  dplyr::arrange(.data$Policy, .data$ConsensusRank, .data$ConsensusScore, dplyr::desc(.data$n_algos), .data$Feature_label)

# Feature_key levels, top is rank 1
mk_levels <- function(policy_name) {
  keys_top_to_bottom <- lab_tbl %>%
    dplyr::filter(as.character(.data$Policy) == policy_name) %>%
    dplyr::mutate(Feature_key = paste(policy_name, .data$Feature, sep = "||")) %>%
    dplyr::pull(.data$Feature_key)
  rev(keys_top_to_bottom)
}

levels_all <- c(
  mk_levels("Non-invasive"),
  mk_levels("Laboratory augmented"),
  mk_levels("Cascade")
)

lab_vec <- lab_tbl %>%
  dplyr::mutate(Feature_key = paste(as.character(.data$Policy), .data$Feature, sep = "||")) %>%
  dplyr::select(.data$Feature_key, .data$Feature_label) %>%
  tibble::deframe()

# Algorithm columns, reversed, then Consensus
pref <- c("C4.5","kNN","SVMRBF","RF","LR")
algo_levels <- c(pref[pref %in% as.character(algos)], setdiff(as.character(algos), pref), "Consensus")

plot_alg <- feature_importance_perm_dMCC %>%
  dplyr::mutate(
    Policy = as.character(.data$Policy),
    Feature_key = paste(.data$Policy, .data$Feature, sep = "||"),
    Algorithm = as.character(.data$Algorithm),
    CellLabel = as.character(.data$Rank)
  ) %>%
  dplyr::select(.data$Policy, .data$Algorithm, .data$Feature_key, .data$Drop_MCC, .data$CellLabel)

plot_cons <- consensus_tbl %>%
  dplyr::mutate(
    Policy = as.character(.data$Policy),
    Feature_key = paste(.data$Policy, .data$Feature, sep = "||"),
    Algorithm = "Consensus",
    Drop_MCC = NA_real_,
    CellLabel = as.character(.data$ConsensusRank)
  ) %>%
  dplyr::select(.data$Policy, .data$Algorithm, .data$Feature_key, .data$Drop_MCC, .data$CellLabel)

feature_importance_perm_dMCC_plotdata <- dplyr::bind_rows(plot_alg, plot_cons) %>%
  dplyr::mutate(
    Policy = factor(.data$Policy, levels = c("Non-invasive","Laboratory augmented","Cascade")),
    Algorithm = factor(.data$Algorithm, levels = algo_levels),
    Feature_key = factor(.data$Feature_key, levels = levels_all)
  )

lim <- stats::quantile(abs(feature_importance_perm_dMCC_plotdata$Drop_MCC), 0.95, na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- max(abs(feature_importance_perm_dMCC_plotdata$Drop_MCC), na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- 1e-6

p_fi_heatmap <- ggplot2::ggplot(
  feature_importance_perm_dMCC_plotdata,
  ggplot2::aes(x = .data$Algorithm, y = .data$Feature_key, fill = .data$Drop_MCC)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = .data$CellLabel), size = 2.2, na.rm = TRUE) +
  ggplot2::scale_fill_gradient2(
    low = "#D55E00",
    mid = "white",
    high = "#0072B2",
    midpoint = 0,
    limits = c(-lim, lim),
    oob = scales::squish,
    name = expression(Delta * "MCC"),
    na.value = "white"
  ) +
  ggplot2::facet_grid(rows = vars(Policy), scales = "free_y", space = "free_y") +
  ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
  ggplot2::scale_y_discrete(labels = function(x) unname(lab_vec[x])) +
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

print(p_fi_heatmap)

