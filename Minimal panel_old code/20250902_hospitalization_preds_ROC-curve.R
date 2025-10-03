# ===== ROC only, 5 panels (one per algorithm), color = FeatureSet, print only =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(pROC)
})
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("margin","ggplot2", quiet = TRUE)
}

# helper: outer-CV predictions for one feature set (uses your existing build_cv_splits, model_specs, train_inner, train_final)
cv_oof_predictions <- function(df, yvar, features, algos,
                               pos_label = "Yes", cv_R = 1, cv_k = 3, inner_k = 3,
                               tune_len = 3, seed = 444, verbose = FALSE) {
  dat <- df %>%
    dplyr::filter(!is.na(.data[[yvar]])) %>%
    dplyr::select(dplyr::all_of(c(yvar, features, ".rid")))
  splits <- build_cv_splits(dat[[yvar]], R = cv_R, k_desired = cv_k, seed = seed)
  specs  <- model_specs(tune_len, algos)
  
  rows <- list()
  for (algo in names(specs)) {
    spec <- specs[[algo]]
    for (nm in names(splits)) {
      tr_idx <- splits[[nm]]$train_idx
      va_idx <- splits[[nm]]$test_idx
      dat_tr <- dat[tr_idx, , drop = FALSE]
      dat_va <- dat[va_idx, , drop = FALSE]
      
      inner <- train_inner(dat_tr, yvar, spec, inner_k); if (is.null(inner)) next
      bestTune <- inner$fit$bestTune %||% NULL
      
      fit <- train_final(dat_tr, yvar, spec, bestTune, inner_k); if (is.null(fit)) next
      
      lev <- levels(dat_tr[[yvar]])
      pr  <- try(predict(fit, newdata = dat_va, type = "prob"), silent = TRUE)
      if (inherits(pr, "try-error") || is.null(pr) || !(pos_label %in% colnames(pr))) next
      
      rows[[length(rows) + 1]] <- tibble::tibble(
        Algorithm = algo,
        Fold      = nm,
        obs       = factor(dat_va[[yvar]], levels = lev),
        prob_pos  = as.numeric(pr[[pos_label]])
      )
      if (isTRUE(verbose)) cat(sprintf("[Fold %s] %s done\n", nm, algo))
    }
  }
  if (!length(rows)) return(tibble::tibble())
  dplyr::bind_rows(rows)
}

# --- Build outer OOF predictions for both feature sets ---
triage_pred <- cv_oof_predictions(
  df, yvar, feat_triage, cfg$algos,
  pos_label = cfg$pos_label, cv_R = cfg$cv_R, cv_k = cfg$cv_k,
  inner_k = cfg$inner_k, tune_len = cfg$tune_len, seed = cfg$seed_cv
) %>% dplyr::mutate(FeatureSet = "Triage")

full_pred <- cv_oof_predictions(
  df, yvar, feat_full, cfg$algos,
  pos_label = cfg$pos_label, cv_R = cfg$cv_R, cv_k = cfg$cv_k,
  inner_k = cfg$inner_k, tune_len = cfg$tune_len, seed = cfg$seed_cv
) %>% dplyr::mutate(FeatureSet = "Full")

pred_all <- dplyr::bind_rows(full_pred, triage_pred)

# --- ROC curve points and AUC per Algorithm × FeatureSet ---
roc_df <- pred_all %>%
  dplyr::group_by(Algorithm, FeatureSet) %>%
  dplyr::group_modify(function(d, key){
    r <- pROC::roc(d$obs, d$prob_pos, quiet = TRUE, levels = levels(d$obs), direction = "<")
    tibble::tibble(
      fpr = 1 - r$specificities,
      tpr = r$sensitivities,
      AUC = as.numeric(pROC::auc(r))
    )
  }) %>% dplyr::ungroup()

# --- Invert algorithm order in panels, then plot ---
algo_order_inv <- c("C4.5","k-NN","SVM","RF","LR")

roc_df$Algorithm <- factor(roc_df$Algorithm, levels = algo_order_inv)

auc_lab <- roc_df %>%
  dplyr::group_by(Algorithm, FeatureSet) %>%
  dplyr::summarise(AUC = max(AUC, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(x = 0.60,
                y = ifelse(FeatureSet == "Full", 0.15, 0.05),
                lab = sprintf("%s AUC=%.3f", FeatureSet, AUC))

p_roc_only <- ggplot(roc_df, aes(x = fpr, y = tpr, color = FeatureSet)) +
  geom_path(linewidth = 0.9, alpha = 0.95) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_text(data = auc_lab, aes(x = x, y = y, label = lab, color = FeatureSet),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ Algorithm, nrow = 1) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("Full" = "#1f78b4", "Triage" = "#33a02c"), name = NULL) +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text = element_text(color = "black")
  )

print(p_roc_only)

# save PNG for Overleaf
ggplot2::ggsave(
  filename = "fig_roc_panels.png",
  plot     = p_roc_only,
  width    = 190,   # mm
  height   = 70,    # mm
  units    = "mm",
  dpi      = 300,
  bg       = "white"
)
