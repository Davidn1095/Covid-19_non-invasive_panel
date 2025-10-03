# ============================================================
# 2) ABLATION PIPELINE FOR FIG 7 (parsimony), robust version
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2)
})

stopifnot(exists("cfg"), exists("feat_triage"), exists("feat_full"), exists("df"))
yvar <- "Hosp_Bin"

safe_auc_lr <- function(feat) {
  agg <- try(
    run_cv(df, yvar,
           features = feat, algos = "LR",
           type = "binary", pos_label = cfg$pos_label,
           cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
           tune_len = cfg$tune_len, seed = cfg$seed_cv,
           verbose = FALSE, return = "agg"),
    silent = TRUE
  )
  if (inherits(agg, "try-error") || !is.data.frame(agg) ||
      !"Metric" %in% names(agg) || !nrow(agg)) return(NA_real_)
  val <- agg %>% dplyr::filter(Metric == "AUC") %>% dplyr::pull(Mean)
  if (length(val) == 0) NA_real_ else as.numeric(val[1])
}

rank_tbl <- tibble(Feature = feat_triage) %>%
  rowwise() %>% mutate(AUC = safe_auc_lr(Feature)) %>% ungroup() %>%
  arrange(dplyr::desc(AUC))

feat_order <- rank_tbl$Feature
cat("\nOrden de features para ablation (mejor primero):\n"); print(rank_tbl, n = nrow(rank_tbl))

m <- length(feat_order)
abl_rows <- vector("list", m)
for (k in seq_len(m)) {
  feats_k <- feat_order[seq_len(k)]
  agg <- try(
    run_cv(df, yvar,
           features = feats_k, algos = cfg$algos,
           type = "binary", pos_label = cfg$pos_label,
           cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
           tune_len = cfg$tune_len, seed = cfg$seed_cv,
           verbose = FALSE, return = "agg"),
    silent = TRUE
  )
  if (!(inherits(agg, "try-error")) && is.data.frame(agg) && nrow(agg)) {
    abl_rows[[k]] <- agg %>% mutate(K = k) %>% dplyr::filter(Metric %in% c("MCC","AUC"))
  }
}
abl <- bind_rows(abl_rows)
stopifnot(nrow(abl) > 0)

abl_mcc <- abl %>% dplyr::filter(Metric == "MCC") %>% select(K, Algorithm, Mean, SD) %>% arrange(K, Algorithm)
cat("\n=== Fig 7 datos — MCC por K (media ± SD) ===\n")
pretty_tbl <- abl_mcc %>% mutate(Value = sprintf("%.4f ± %.4f", Mean, SD)) %>%
  select(K, Algorithm, Value) %>% tidyr::pivot_wider(names_from = Algorithm, values_from = Value)
print(pretty_tbl, n = nrow(pretty_tbl))

p_fig7 <- abl %>%
  dplyr::filter(Metric == "MCC") %>%
  ggplot(aes(x = K, y = Mean, group = Algorithm)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)), width = 0.15) +
  scale_x_continuous(breaks = seq_len(m)) +
  labs(title = "Parsimony curve, triage", subtitle = "MCC vs número de features, refit en cada paso",
       x = "Número de features (K)", y = "MCC (media ± SD)") +
  theme_minimal(base_size = 12)
print(p_fig7)

incl_tbl <- tibble(K = seq_len(m),
                   Included = purrr::map_chr(seq_len(m), ~ paste(feat_order[seq_len(.x)], collapse = ", ")))
cat("\n=== Inclusión por K ===\n"); print(incl_tbl, n = m)

# Objects: rank_tbl, abl, abl_mcc, p_fig7, incl_tbl


ggsave(filename = "fig_ablation_triage.png",
       plot = p_fig7, width = 3.5, height = 2.6, units = "in", dpi = 300)
