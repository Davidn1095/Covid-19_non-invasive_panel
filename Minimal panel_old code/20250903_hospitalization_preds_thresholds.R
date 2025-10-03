# ==== Operating thresholds for MCC, F1, Accuracy, Precision, Sensitivity, Specificity
# ==== plus AUC-ROC operating point, with AUC_ROC column appended ====

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(pROC)
})

# pred_all must exist with columns: Algorithm, FeatureSet, obs, prob_pos
stopifnot(all(c("Algorithm","FeatureSet","obs","prob_pos") %in% names(pred_all)))

# -------- helpers --------
.as_factor01 <- function(y){
  y <- factor(y)
  if (length(levels(y)) != 2) stop("Outcome must have 2 levels")
  y
}
.mcc <- function(tp, tn, fp, fn){
  num <- tp*tn - fp*fn
  den <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if (den == 0) return(NA_real_)
  num / den
}
.metrics_at <- function(y, prob, thr){
  lev <- levels(y); pos <- lev[2]; neg <- lev[1]
  pred <- factor(ifelse(prob >= thr, pos, neg), levels = lev)
  cm <- table(y, pred)
  TP <- as.numeric(cm[pos, pos]); TN <- as.numeric(cm[neg, neg])
  FP <- as.numeric(cm[neg, pos]);  FN <- as.numeric(cm[pos, neg])
  eps <- 1e-9
  sens <- TP / pmax(TP + FN, eps)
  spec <- TN / pmax(TN + FP, eps)
  tibble::tibble(
    Threshold   = thr,
    MCC         = .mcc(TP, TN, FP, FN),
    F1          = 2*TP / pmax(2*TP + FP + FN, eps),
    Accuracy    = (TP + TN) / pmax(sum(cm), eps),
    Precision   = TP / pmax(TP + FP, eps),
    Sensitivity = sens,
    Specificity = spec
  )
}
.pick_best <- function(df, target){
  stopifnot(target %in% names(df))
  ord <- order(-df[[target]], -df$MCC, -df$Precision, -df$Specificity, na.last = TRUE)
  df[ord, , drop = FALSE] |> dplyr::slice(1)
}
.find_selected_cuts <- function(y, prob){
  p <- sort(unique(prob)); grid <- c(0, p, 1)
  mets <- purrr::map_dfr(grid, ~ .metrics_at(y, prob, .x))
  
  # AUC-ROC operating point: closest to (0,1)
  roc_closest <- mets %>%
    mutate(Dist2 = (1 - Sensitivity)^2 + (1 - Specificity)^2) %>%
    { .[order(.$Dist2, -.$MCC, -.$Precision, -.$Specificity), , drop = FALSE] } %>%
    slice(1) %>%
    mutate(Selection = "AUC-ROC closest (0,1)") %>%
    select(-Dist2)
  
  # maximize core thresholded metrics
  targets <- c("MCC","F1","Accuracy","Precision","Sensitivity","Specificity")
  per_metric <- purrr::map_dfr(targets, function(targ){
    .pick_best(mets, targ) %>% mutate(Selection = paste0("Max ", targ))
  })
  
  bind_rows(per_metric, roc_closest)
}

# -------- build operating-point table (op_core) --------
op_core <- pred_all %>%
  mutate(obs = .as_factor01(obs)) %>%
  group_by(FeatureSet, Algorithm) %>%
  group_modify(~ .find_selected_cuts(.x$obs, .x$prob_pos)) %>%
  ungroup() %>%
  relocate(Selection, .after = Algorithm)

# -------- append AUC_ROC per Algorithm × FeatureSet --------
auc_tbl <- pred_all %>%
  mutate(y = as.integer(obs == levels(obs)[2])) %>%
  group_by(FeatureSet, Algorithm) %>%
  summarise(AUC_ROC = as.numeric(pROC::auc(y, prob_pos, quiet = TRUE)), .groups = "drop")

op_core_auc <- op_core %>%
  left_join(auc_tbl, by = c("FeatureSet","Algorithm")) %>%
  relocate(AUC_ROC, .after = Threshold) %>%
  mutate(across(c(Threshold, AUC_ROC, MCC, F1, Accuracy, Precision, Sensitivity, Specificity),
                ~ round(.x, 3))) %>%
  arrange(FeatureSet, Algorithm, Selection)

print(op_core_auc, n = nrow(op_core_auc))
