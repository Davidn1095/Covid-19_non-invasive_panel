source(file.path("R", "utils_metrics.R"))

best_thr_by_mcc <- function(y01, p, thr_grid) {
  best_mcc <- -Inf
  best_thr <- thr_grid[1]
  best_sens <- -Inf

  for (thr in thr_grid) {
    met <- compute_metrics_at_thr(y01, p, thr)
    mcc <- met$MCC
    sens <- met$Sensitivity
    if (!is.finite(mcc) || is.na(mcc)) next

    if (mcc > best_mcc) {
      best_mcc <- mcc
      best_thr <- thr
      best_sens <- sens
    } else if (mcc == best_mcc && is.finite(sens) && !is.na(sens) && sens > best_sens) {
      best_thr <- thr
      best_sens <- sens
    }
  }

  list(thr = best_thr, mcc = best_mcc, sens = best_sens)
}

best_thr_mcc_at_sens <- function(y01, p, thr_grid, sens_min = 0.70, spec_min = -Inf) {
  best_thr  <- thr_grid[1]
  best_mcc  <- -Inf
  best_sens <- -Inf
  best_spec <- -Inf

  for (thr in thr_grid) {
    met  <- compute_metrics_at_thr(y01, p, thr)
    mcc  <- met$MCC
    sens <- met$Sensitivity
    spec <- met$Specificity

    if (!is.finite(mcc)  || is.na(mcc))  next
    if (!is.finite(sens) || is.na(sens)) next
    if (!is.finite(spec) || is.na(spec)) next

    if (sens < sens_min) next
    if (is.finite(spec_min) && spec < spec_min) next

    if (sens < spec) next

    if (mcc > best_mcc ||
        (mcc == best_mcc && sens > best_sens) ||
        (mcc == best_mcc && sens == best_sens && spec > best_spec)) {
      best_mcc  <- mcc
      best_thr  <- thr
      best_sens <- sens
      best_spec <- spec
    }
  }

  if (!is.finite(best_mcc) || is.na(best_mcc)) {
    out <- best_thr_by_mcc(y01, p, thr_grid)
    return(list(thr = out$thr, mcc = out$mcc, sens = out$sens, spec = NA_real_, hit_constraint = FALSE))
  }

  list(thr = best_thr, mcc = best_mcc, sens = best_sens, spec = best_spec, hit_constraint = TRUE)
}

orient_prob_global <- function(y01, p) {
  auc <- auc_from_prob_fixed(y01, p)
  flip <- is.finite(auc) && !is.na(auc) && auc < 0.5
  list(p = if (flip) 1 - p else p, flip = flip, auc = auc)
}

best_cascade_params <- function(y01, p1, p2,
                                tau_grid, thr2_grid,
                                sens_min = -Inf, spec_min = -Inf,
                                max_def_rate = Inf,
                                require_sens_ge_spec = FALSE) {
  y <- as.integer(y01)
  p1 <- as.numeric(p1)
  p2 <- as.numeric(p2)

  ok <- is.finite(p1) & !is.na(p1) & is.finite(p2) & !is.na(p2) & !is.na(y)
  y  <- y[ok]; p1 <- p1[ok]; p2 <- p2[ok]

  if (length(y) < 50L) {
    return(list(tau_low = 0.49, tau_high = 0.51, thr2 = 0.50, mcc = NA_real_,
                sens = NA_real_, spec = NA_real_, def_rate = NA_real_, hit_constraint = FALSE))
  }

  tau_low_grid  <- tau_grid
  tau_high_grid <- tau_grid

  best <- list(mcc = -Inf, sens = -Inf, spec = -Inf, def_rate = Inf,
               tau_low = tau_low_grid[1], tau_high = tau_high_grid[1], thr2 = thr2_grid[1],
               hit_constraint = TRUE)

  found <- FALSE

  for (tau_low in tau_low_grid) {
    for (tau_high in tau_high_grid) {
      if (!(tau_low < tau_high)) next

      defer <- (p1 > tau_low) & (p1 < tau_high)
      def_rate <- mean(defer)
      if (is.finite(max_def_rate) && def_rate > max_def_rate) next

      for (thr2 in thr2_grid) {
        pred <- ifelse(p1 >= tau_high, 1L,
                       ifelse(p1 <= tau_low, 0L,
                              as.integer(p2 >= thr2)))

        cc <- confusion_counts(y, pred)
        mcc <- mcc_from_counts(cc)
        if (!is.finite(mcc) || is.na(mcc)) next

        TP <- as.numeric(cc["TP"]); TN <- as.numeric(cc["TN"])
        FP <- as.numeric(cc["FP"]); FN <- as.numeric(cc["FN"])
        sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
        spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)

        if (!is.finite(sens) || is.na(sens) || !is.finite(spec) || is.na(spec)) next
        if (is.finite(sens_min) && sens < sens_min) next
        if (is.finite(spec_min) && spec < spec_min) next
        if (isTRUE(require_sens_ge_spec) && sens < spec) next

        found <- TRUE

        if (mcc > best$mcc ||
            (mcc == best$mcc && sens > best$sens) ||
            (mcc == best$mcc && sens == best$sens && def_rate < best$def_rate) ||
            (mcc == best$mcc && sens == best$sens && def_rate == best$def_rate && spec > best$spec)) {
          best <- list(
            mcc = mcc, sens = sens, spec = spec, def_rate = def_rate,
            tau_low = tau_low, tau_high = tau_high, thr2 = thr2,
            hit_constraint = TRUE
          )
        }
      }
    }
  }

  if (!found) {
    best2 <- list(mcc = -Inf, sens = -Inf, spec = -Inf, def_rate = Inf,
                  tau_low = tau_low_grid[1], tau_high = tau_high_grid[1], thr2 = thr2_grid[1],
                  hit_constraint = FALSE)

    for (tau_low in tau_low_grid) {
      for (tau_high in tau_high_grid) {
        if (!(tau_low < tau_high)) next

        defer <- (p1 > tau_low) & (p1 < tau_high)
        def_rate <- mean(defer)
        if (is.finite(max_def_rate) && def_rate > max_def_rate) next

        for (thr2 in thr2_grid) {
          pred <- ifelse(p1 >= tau_high, 1L,
                         ifelse(p1 <= tau_low, 0L,
                                as.integer(p2 >= thr2)))

          cc <- confusion_counts(y, pred)
          mcc <- mcc_from_counts(cc)
          if (!is.finite(mcc) || is.na(mcc)) next

          TP <- as.numeric(cc["TP"]); TN <- as.numeric(cc["TN"])
          FP <- as.numeric(cc["FP"]); FN <- as.numeric(cc["FN"])
          sens <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
          spec <- if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
          if (!is.finite(sens) || is.na(sens) || !is.finite(spec) || is.na(spec)) next

          if (mcc > best2$mcc ||
              (mcc == best2$mcc && sens > best2$sens) ||
              (mcc == best2$mcc && sens == best2$sens && def_rate < best2$def_rate) ||
              (mcc == best2$mcc && sens == best2$sens && def_rate == best2$def_rate && spec > best2$spec)) {
            best2 <- list(
              mcc = mcc, sens = sens, spec = spec, def_rate = def_rate,
              tau_low = tau_low, tau_high = tau_high, thr2 = thr2,
              hit_constraint = FALSE
            )
          }
        }
      }
    }
    return(best2)
  }

  best
}
