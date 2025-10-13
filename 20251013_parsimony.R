# ===============================================================
# PARSIMONY CURVES (MCC), importance ordered, LR point at k = 1
# Needs: df, yvar, feat_full, feat_triage, algo_order already defined
#        run_holdout(), compute_metrics_binary(), tune_band_from_triage_oof(), apply_cascade_external()
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(patchwork)
  library(cowplot)
  library(forcats)
  library(digest)
})

# switches
speed_mode       <- TRUE
speed_ultra      <- TRUE
B_parsimony      <- if (speed_ultra) 0 else if (speed_mode) 200 else 1000
rank_tune_length <- if (speed_ultra) 8 else 20

# memoized holdout
.hold_cache <- new.env(parent = emptyenv())
run_holdout_memo <- function(df, yvar, features, algo, set_label) {
  key <- paste(algo, set_label, digest::digest(sort(features)), sep = "|")
  if (exists(key, envir = .hold_cache, inherits = FALSE))
    return(get(key, envir = .hold_cache))
  res <- try(run_holdout(df, yvar, features, algo, set_label), silent = TRUE)
  if (inherits(res, "try-error")) res <- NULL
  assign(key, res, envir = .hold_cache)
  res
}

# robust preprocessor
prep_df_for_algo <- function(d, feats, algo) {
  feats <- intersect(feats, names(d))
  d2 <- d
  if (!is.null(d2[[yvar]])) {
    d2 <- d2[!is.na(d2[[yvar]]), , drop = FALSE]
    d2[[yvar]] <- factor(d2[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))
  }
  if (!length(feats)) return(list(df = d2, feats = feats))
  is_char <- vapply(d2[feats], is.character, logical(1))
  if (any(is_char)) {
    ch <- feats[is_char]
    d2[ch] <- lapply(d2[ch], function(x) factor(replace(x, is.na(x), "Unknown")))
  }
  is_fac <- vapply(d2[feats], is.factor, logical(1))
  if (any(is_fac)) {
    fc <- feats[is_fac]
    d2[fc] <- lapply(d2[fc], function(x) forcats::fct_na_value_to_level(x, "Unknown"))
  }
  is_num <- vapply(d2[feats], is.numeric, logical(1))
  for (cl in feats[is_num]) {
    if (all(is.na(d2[[cl]]))) d2[[cl]] <- 0 else {
      med <- median(d2[[cl]], na.rm = TRUE)
      d2[[cl]][is.na(d2[[cl]])] <- med
    }
  }
  zv <- feats[vapply(d2[feats], function(x) { ux <- unique(x[!is.na(x)]); length(ux) <= 1 }, logical(1))]
  feats <- setdiff(feats, zv)
  list(df = d2, feats = feats)
}

# cascade label
cascade_label <- "Cascade model"

# reproducibility
RNGkind("L'Ecuyer-CMRG"); set.seed(123)

# cohort filter guard
stopifnot("group" %in% names(df))
df <- dplyr::filter(df, group %in% c("CTRL_noCOVID","COVID"))

# TRAIN rows for ranking
get_train_mask <- function(d) {
  if ("data_split" %in% names(d)) {
    tolower(as.character(d$data_split)) %in% c("train","training","cv","inner","dev")
  } else if ("split" %in% names(d)) {
    tolower(as.character(d$split)) %in% c("train","training","cv","inner","dev")
  } else TRUE
}
train_mask <- get_train_mask(df)
df_train   <- df[train_mask, , drop = FALSE]

# reference importance for ordering
ref_feats <- intersect(feat_full, names(df_train)); stopifnot(length(ref_feats) > 0)
cols_use  <- unique(c(yvar, ref_feats))
df_rank   <- df_train[, cols_use, drop = FALSE]

num_cols <- intersect(ref_feats, names(df_rank)[vapply(df_rank, is.numeric, logical(1))])
for (cl in num_cols) {
  med <- median(df_rank[[cl]], na.rm = TRUE)
  if (!is.na(med)) df_rank[[cl]][is.na(df_rank[[cl]])] <- med
}
chr_fac <- setdiff(ref_feats, num_cols)
if (length(chr_fac)) {
  df_rank[chr_fac] <- lapply(df_rank[chr_fac], function(x) {
    x <- as.factor(x); forcats::fct_na_value_to_level(x, level = "Unknown")
  })
}
df_rank[[yvar]] <- factor(df_rank[[yvar]], levels = c(cfg$neg_label, cfg$pos_label))

form_ref <- as.formula(paste(yvar, "~", paste(ref_feats, collapse = "+")))
ctrl_ref <- caret::trainControl(method = "cv",
                                number = 5,
                                classProbs = TRUE,
                                summaryFunction = twoClassSummary,
                                savePredictions = "final",
                                allowParallel = TRUE)
options(na.action = "na.pass")
set.seed(123)
fit_ref <- caret::train(form_ref,
                        data = df_rank,
                        method = "glmnet",
                        family = "binomial",
                        tuneLength = rank_tune_length,
                        metric = "ROC",
                        trControl = ctrl_ref)

imp_tbl <- caret::varImp(fit_ref)$importance
ref_ord <- rownames(imp_tbl)[order(imp_tbl$Overall, decreasing = TRUE)]
ref_ord <- unique(c(ref_ord, ref_feats))

# build orders
tri_order    <- intersect(ref_ord, feat_triage)
triage_k_cap <- 4
tri_order    <- head(tri_order, triage_k_cap)

rest_pool  <- setdiff(intersect(feat_full, names(df)), tri_order)
rest_order <- unique(c(intersect(ref_ord, rest_pool),
                       sort(setdiff(rest_pool, ref_ord), decreasing = TRUE)))

order_by_set <- list(
  "Complete feature set" = c(tri_order, rest_order),
  "Triage feature set"   = tri_order
)
full_order <- c(tri_order, rest_order)

# τ rule safety wrapper
best_threshold_mcc_bacc_med <- local({
  f_old <- get("best_threshold_mcc_bacc_med", inherits = TRUE)
  function(obs, p, grid = seq(0.01, 0.99, by = 0.001), tol = 1e-12) {
    pos <- cfg$pos_label; neg <- cfg$neg_label
    y <- factor(obs, levels = c(neg, pos))
    if (length(unique(stats::na.omit(y))) < 2 || sd(p, na.rm = TRUE) < .Machine$double.eps)
      return(list(t = 0.5, mcc = NA_real_))
    out <- try(f_old(obs, p, grid, tol), silent = TRUE)
    if (inherits(out, "try-error") || !is.finite(out$t)) list(t = 0.5, mcc = NA_real_) else out
  }
})

# bootstrap SD for ribbons
boot_sd_mcc <- function(obs, pred, p, B = B_parsimony, seed = 123) {
  if (!is.numeric(B) || B <= 0) return(NA_real_)
  if (length(obs) < 2 || length(unique(obs)) < 2) return(NA_real_)
  n <- length(obs); set.seed(seed)
  tryCatch({
    sd(replicate(B, {
      idx  <- sample.int(n, n, replace = TRUE)
      as.numeric(compute_metrics_binary(obs[idx], pred[idx], p[idx])["MCC"])
    }), na.rm = TRUE)
  }, error = function(e) NA_real_)
}

# helper to guarantee LR works with a single feature
add_lr_shadow_if_needed <- function(pp, algo, seed = 123) {
  if (identical(algo, "LR") && length(pp$feats) == 1) {
    df_aug <- pp$df
    set.seed(seed)
    nm <- "shadow_LR"
    while (nm %in% names(df_aug)) nm <- paste0("shadow_LR_", sample(1000:9999, 1))
    df_aug[[nm]] <- rnorm(nrow(df_aug), sd = 1e-8)
    list(df = df_aug, feats = c(pp$feats, nm))
  } else pp
}

# curves
parsimony_external_curve <- function(set_name, feats_order, algo) {
  if (!length(feats_order)) {
    return(tibble(k = integer(), Score = numeric(), SD = numeric(),
                  Set = character(), Algorithm = character()))
  }
  purrr::map_dfr(seq_along(feats_order), function(k) {
    feats_k <- feats_order[seq_len(k)]
    lab     <- ifelse(set_name == "Complete feature set", "Complete", "Triage")
    pp <- prep_df_for_algo(df, feats_k, algo)
    pp <- add_lr_shadow_if_needed(pp, algo)
    hold <- if (algo == "RF")
      suppressWarnings(run_holdout_memo(pp$df, yvar, pp$feats, algo, lab))
    else
      run_holdout_memo(pp$df, yvar, pp$feats, algo, lab)
    if (is.null(hold)) {
      tibble(k = k, Score = NA_real_, SD = NA_real_,
             Set = set_name, Algorithm = algo)
    } else {
      mets   <- compute_metrics_binary(hold$obs, hold$pred, hold$p)
      sd_mcc <- boot_sd_mcc(hold$obs, hold$pred, hold$p, B = B_parsimony)
      tibble(k = k,
             Score = as.numeric(mets["MCC"]),
             SD    = sd_mcc,
             Set = set_name,
             Algorithm = algo)
    }
  })
}

# cascade over FULL order with TRIAGE gate
cascade_external_curve <- function(feats_tri_order, full_order, algo) {
  K_full <- length(full_order)
  if (!K_full) {
    return(tibble(k = integer(), Score = numeric(), SD = numeric(),
                  Set = character(), Algorithm = character()))
  }
  purrr::map_dfr(seq_len(K_full), function(k) {
    tri_feats_k <- head(feats_tri_order, min(k, length(feats_tri_order)))
    com_feats_k <- head(full_order, k)
    pp_tri <- prep_df_for_algo(df, tri_feats_k, algo)
    pp_com <- prep_df_for_algo(df, com_feats_k,  algo)
    pp_tri <- add_lr_shadow_if_needed(pp_tri, algo)
    pp_com <- add_lr_shadow_if_needed(pp_com, algo)
    tri <- if (algo == "RF")
      suppressWarnings(run_holdout_memo(pp_tri$df, yvar, pp_tri$feats, algo, "Triage"))
    else
      run_holdout_memo(pp_tri$df, yvar, pp_tri$feats, algo, "Triage")
    com <- if (algo == "RF")
      suppressWarnings(run_holdout_memo(pp_com$df, yvar, pp_com$feats, algo, "Complete"))
    else
      run_holdout_memo(pp_com$df, yvar, pp_com$feats, algo, "Complete")
    if (is.null(tri) || is.null(com)) {
      return(tibble(k = k, Score = NA_real_, SD = NA_real_,
                    Set = cascade_label, Algorithm = algo))
    }
    tl <- cfg$cascade$t_low; th <- cfg$cascade$t_high
    if (!is.null(tri$oof) && nrow(tri$oof) && is.finite(tri$threshold)) {
      tuned <- try(tune_band_from_triage_oof(tri$oof, thr_tri = tri$threshold), silent = TRUE)
      if (!inherits(tuned, "try-error") && is.finite(tuned$t_low) && is.finite(tuned$t_high)) {
        tl <- tuned$t_low; th <- tuned$t_high
      }
    }
    cas <- try(apply_cascade_external(tri, com, t_low = tl, t_high = th), silent = TRUE)
    if (inherits(cas, "try-error")) {
      return(tibble(k = k, Score = NA_real_, SD = NA_real_,
                    Set = cascade_label, Algorithm = algo))
    }
    mets   <- compute_metrics_binary(cas$obs, cas$pred, cas$p)
    sd_mcc <- boot_sd_mcc(cas$obs, cas$pred, cas$p, B = B_parsimony)
    tibble(k = k, Score = as.numeric(mets["MCC"]), SD = sd_mcc,
           Set = cascade_label, Algorithm = algo)
  })
}

# compute panels
if (!exists("algo_order")) algo_order <- c("C4.5","k-NN","SVM","RF","LR")

res_core <- purrr::imap_dfr(order_by_set, function(ord, set_name) {
  purrr::map_dfr(algo_order, function(algo) parsimony_external_curve(set_name, ord, algo))
})
res_cas  <- purrr::map_dfr(algo_order, function(algo)
  cascade_external_curve(tri_order, full_order, algo))

res_all <- dplyr::bind_rows(res_core, res_cas) %>%
  dplyr::mutate(
    Algorithm = factor(Algorithm, levels = algo_order),
    Set = factor(Set, levels = c("Complete feature set","Triage feature set", cascade_label))
  )

# assert expected counts
chk <- dplyr::count(res_all, Set, Algorithm)
print(chk %>% tidyr::pivot_wider(names_from = Set, values_from = n))
stopifnot(all(chk$n[chk$Set == "Triage feature set"] == length(tri_order)))
stopifnot(all(chk$n[chk$Set == cascade_label]          == length(full_order)))

# plotting
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "grey30"),
          legend.position = "right",
          legend.title = element_blank())
}

levels_set <- c("Complete feature set","Triage feature set", cascade_label)
pal <- c("Complete feature set" = "#D55E00",  # orange
         "Triage feature set"   = "#0072B2",  # blue
         "Cascade model"        = "#009E73")  # green

make_alg_plot <- function(alg, show_y = FALSE, show_x = FALSE) {
  dat <- dplyr::filter(res_all, Algorithm == alg)
  if (!nrow(dat) || all(is.na(dat$Score)))
    return(ggplot() + theme_void() + ggtitle(alg))
  tri_name <- "Triage feature set"
  dat_tri  <- dplyr::filter(dat, Set == tri_name)
  dat_oth  <- dplyr::filter(dat, Set != tri_name)
  x_max <- max(dat$k, na.rm = TRUE) + 0.5
  
  p <- ggplot(dat, aes(k, Score, color = Set, group = Set)) +
    annotate("segment", x = 0.5, xend = x_max, y = 0, yend = 0, linewidth = 0.6, colour = "black") +
    annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 1, linewidth = 0.6, colour = "black")
  
  # ribbons under everything
  if (!all(is.na(dat$SD))) {
    p <- p +
      geom_ribbon(data = dat_oth, aes(ymin = pmax(0, Score - SD), ymax = pmin(1, Score + SD), fill = Set),
                  alpha = 0.15, colour = NA, show.legend = FALSE) +
      geom_ribbon(data = dat_tri, aes(ymin = pmax(0, Score - SD), ymax = pmin(1, Score + SD), fill = Set),
                  alpha = 0.15, colour = NA, show.legend = FALSE) +
      scale_fill_manual(values = pal, breaks = levels_set)
  }
  
  # other sets first
  p <- p +
    geom_line(data = dat_oth, linewidth = 1) +
    geom_point(data = dat_oth, size = 1.8, na.rm = TRUE)
  
  # triage on top
  p <- p +
    geom_line(data = dat_tri, linewidth = 1.2) +
    geom_point(data = dat_tri, size = 2.1, na.rm = TRUE)
  
  p +
    scale_x_continuous(breaks = function(x) seq_len(max(x, na.rm = TRUE)),
                       expand = expansion(mult = c(0.12, 0.03))) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    scale_color_manual(values = pal, breaks = levels_set, limits = levels_set) +
    labs(x = if (show_x) "Number of attributes" else NULL,
         y = if (show_y) "MCC" else NULL,
         title = as.character(alg)) +
    theme_pub() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)))
}

p1 <- make_alg_plot("C4.5", show_y = TRUE,  show_x = FALSE)
p2 <- make_alg_plot("k-NN", show_y = FALSE, show_x = FALSE)
p3 <- make_alg_plot("SVM",  show_y = FALSE, show_x = FALSE)
p4 <- make_alg_plot("RF",   show_y = TRUE,  show_x = TRUE)
p5 <- make_alg_plot("LR",   show_y = FALSE, show_x = TRUE)

legend_src <- ggplot(res_all, aes(k, Score, color = Set)) +
  geom_line() +
  scale_color_manual(values = pal, breaks = levels_set, limits = levels_set) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")
pleg <- patchwork::wrap_elements(cowplot::get_legend(legend_src))

p_parsimony <- (p1 + p2 + p3) / (p4 + p5 + pleg) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1))

print(p_parsimony)

# data behind the curves
res_all
