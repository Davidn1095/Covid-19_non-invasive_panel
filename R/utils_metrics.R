source(file.path("R", "common_utils.R"))

confusion_counts <- function(y01, pred01) {
  y01 <- as.integer(y01)
  pred01 <- as.integer(pred01)
  TP <- sum(y01 == 1L & pred01 == 1L, na.rm = TRUE)
  TN <- sum(y01 == 0L & pred01 == 0L, na.rm = TRUE)
  FP <- sum(y01 == 0L & pred01 == 1L, na.rm = TRUE)
  FN <- sum(y01 == 1L & pred01 == 0L, na.rm = TRUE)
  c(TP = TP, TN = TN, FP = FP, FN = FN)
}

auc_from_prob_fixed <- function(y01, p) {
  ok <- is.finite(p) & !is.na(p) & !is.na(y01)
  if (sum(ok, na.rm = TRUE) < 10L) return(NA_real_)
  y <- factor(y01[ok], levels = c(0, 1))
  roc_obj <- suppressWarnings(pROC::roc(response = y, predictor = p[ok], levels = c(0, 1), direction = "<", quiet = TRUE))
  as.numeric(pROC::auc(roc_obj))
}

compute_metrics_at_thr <- function(y01, p, thr) {
  ok <- is.finite(p) & !is.na(p) & !is.na(y01)
  if (sum(ok, na.rm = TRUE) == 0L) {
    return(list(TP = NA, TN = NA, FP = NA, FN = NA, MCC = NA, AUC = NA, F1 = NA,
                Accuracy = NA, Precision = NA, Sensitivity = NA, Specificity = NA))
  }

  y <- as.integer(y01[ok])
  pr <- as.numeric(p[ok])
  pred <- as.integer(pr >= thr)

  cc <- confusion_counts(y, pred)
  mcc <- mcc_from_counts(cc)

  TP <- cc["TP"]; TN <- cc["TN"]; FP <- cc["FP"]; FN <- cc["FN"]
  acc <- (TP + TN) / (TP + TN + FP + FN)
  prec <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
  f1 <- if (!is.finite(prec) || !is.finite(sens) || (prec + sens) == 0) NA_real_ else 2 * prec * sens / (prec + sens)
  auc <- auc_from_prob_fixed(y, pr)

  list(
    TP = as.integer(TP), TN = as.integer(TN), FP = as.integer(FP), FN = as.integer(FN),
    MCC = mcc, AUC = auc, F1 = f1,
    Accuracy = acc, Precision = prec, Sensitivity = sens, Specificity = spec
  )
}

compute_metrics_from_pred <- function(y01, pred01, score = NULL) {
  ok <- is.finite(pred01) & !is.na(pred01) & !is.na(y01)
  if (sum(ok) == 0L) {
    return(list(TP=NA, TN=NA, FP=NA, FN=NA, MCC=NA, AUC=NA, F1=NA,
                Accuracy=NA, Precision=NA, Sensitivity=NA, Specificity=NA))
  }

  y <- as.integer(y01[ok])
  pr <- as.integer(pred01[ok])

  cc <- confusion_counts(y, pr)
  mcc <- mcc_from_counts(cc)

  TP <- cc["TP"]; TN <- cc["TN"]; FP <- cc["FP"]; FN <- cc["FN"]
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  prec <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
  f1   <- if (!is.finite(prec) || !is.finite(sens) || (prec + sens) == 0) NA_real_ else 2 * prec * sens / (prec + sens)

  auc <- NA_real_
  if (!is.null(score)) {
    sc <- as.numeric(score)[ok]
    auc <- auc_from_prob_fixed(y, sc)
  }

  list(
    TP = as.integer(TP), TN = as.integer(TN), FP = as.integer(FP), FN = as.integer(FN),
    MCC = mcc, AUC = auc, F1 = f1,
    Accuracy = acc, Precision = prec, Sensitivity = sens, Specificity = spec
  )
}
