# ============================================================
# Push se_NB down further (increase k) + lighten multiplicity (per-metric BH)
# Treat all algorithms & metrics equally
# ============================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(stringr) })

# ---------- 1) Try more folds (k) to shrink NB variance ----------
# Note: your helper choose_k_for_task() will cap k at the max feasible
cfg$cv_k <- 20          # target 20 folds (capped automatically if class counts are smaller)
cfg$cv_R <- 10          # repeats (kept same; increasing R helps a little vs increasing k)
cfg$fast_mode <- FALSE

# Rebuild CV results with the new k
full_folds <- run_cv(
  df, yvar, features = feat_full, algos = cfg$algos,
  type = "binary", pos_label = cfg$pos_label,
  cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
  tune_len = cfg$tune_len, seed = cfg$seed_cv,
  verbose = FALSE, return = "folds"
)
triage_folds <- run_cv(
  df, yvar, features = feat_triage, algos = cfg$algos,
  type = "binary", pos_label = cfg$pos_label,
  cv_R = cfg$cv_R, cv_k = cfg$cv_k, inner_k = cfg$inner_k,
  tune_len = cfg$tune_len, seed = cfg$seed_cv,
  verbose = FALSE, return = "folds"
)

# ---------- helpers (define if missing) ----------
extract_num_first <- function(x, patterns) {
  x <- as.character(x)
  for (p in patterns) {
    m <- stringr::str_match(x, p)
    if (!is.null(m)) {
      v <- suppressWarnings(as.integer(m[, 2]))
      if (any(!is.na(v))) return(v)
    }
  }
  rep(NA_integer_, length(x))
}
if (!exists("add_ids", mode = "function")) {
  add_ids <- function(df) {
    df <- dplyr::mutate(df, Fold = as.character(Fold))
    x <- df$Fold
    comb   <- stringr::str_match(x, "(?i)r\\s*([0-9]+)\\s*[_-]?\\s*f\\s*([0-9]+)")
    rep_c  <- suppressWarnings(as.integer(comb[, 2]))
    fold_c <- suppressWarnings(as.integer(comb[, 3]))
    rep_fb  <- extract_num_first(x, c("(?i)(?:^|[_-])r(?:ep|epeat)?\\D*([0-9]+)", "(?i)\\brep(?:eat)?\\D*([0-9]+)\\b"))
    fold_fb <- extract_num_first(x, c("(?i)(?:^|[_-])f(?:old)?\\D*([0-9]+)", "(?i)\\bfold\\D*([0-9]+)\\b", "(?i).*?_([0-9]+)$"))
    dplyr::mutate(df,
                  rep_id  = dplyr::coalesce(rep_c,  rep_fb,  1L),
                  fold_id = dplyr::coalesce(fold_c, fold_fb)
    )
  }
}
if (!exists("nadeau_bengio_summarise", mode = "function")) {
  nadeau_bengio_summarise <- function(full_folds_df, triage_folds_df,
                                      alpha_two_sided = 0.05, alpha_one_sided_NI = 0.025) {
    paired <- dplyr::inner_join(
      dplyr::rename(add_ids(full_folds_df),   Value_full   = Value),
      dplyr::rename(add_ids(triage_folds_df), Value_triage = Value),
      by = c("Algorithm","Fold","Metric","fold_id","rep_id")
    ) |>
      dplyr::mutate(Diff = Value_triage - Value_full)
    
    out <- paired |>
      dplyr::group_by(Algorithm, Metric) |>
      dplyr::summarise(
        mean_full   = mean(Value_full,   na.rm = TRUE),
        mean_triage = mean(Value_triage, na.rm = TRUE),
        diff_mean   = mean(Diff,         na.rm = TRUE),
        sd_diff     = stats::sd(Diff,    na.rm = TRUE),
        N           = dplyr::n(),
        k_est       = dplyr::n_distinct(fold_id),
        R_est       = dplyr::n_distinct(rep_id),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        c_NB   = 1/N + 1/pmax(k_est - 1, 1),
        se_NB  = sd_diff * sqrt(c_NB),
        df     = pmax(N - 1, 1),
        tcrit2 = stats::qt(1 - alpha_two_sided/2, df),
        ci95_lo = diff_mean - tcrit2 * se_NB,
        ci95_hi = diff_mean + tcrit2 * se_NB,
        tcrit_NI = stats::qt(1 - alpha_one_sided_NI, df),
        rel_diff_pct = 100 * diff_mean / pmax(abs(mean_full), 1e-12)
      )
    if (any(out$k_est < 2)) warning("Some k_est < 2; NB correction becomes conservative.")
    out
  }
}

alpha_one_sided_NI <- 0.025
nb_tbl <- nadeau_bengio_summarise(full_folds, triage_folds,
                                  alpha_two_sided = 0.05,
                                  alpha_one_sided_NI = alpha_one_sided_NI)

cat("\n=== NB diagnostics (k, R, N, c_NB, se_NB) with target k=20 ===\n")
print(nb_tbl |> select(Metric, Algorithm, N, k_est, R_est, c_NB, se_NB) |> arrange(Metric, Algorithm), n = Inf)

# ---------- 2) Lighten multiplicity: BH within each metric (treats algos equally) ----------
adjust_scope  <- "by_metric"   # adjust per metric family (lighter than global)
adjust_method <- "BH"          # Benjamini–Hochberg FDR

all_metrics <- c("MCC","F1","Accuracy","Precision","Sensitivity","Specificity","AUC")
ni_margin_by_metric <- setNames(rep(0.03, length(all_metrics)), all_metrics)  # same δ for everyone

nb_ni <- nb_tbl |>
  dplyr::mutate(
    NI_margin   = unname(ni_margin_by_metric[Metric]),
    t_NI        = (diff_mean + NI_margin) / pmax(se_NB, 1e-12),
    p_NI        = 1 - stats::pt(t_NI, df),
    lower_97_5  = diff_mean - tcrit_NI * se_NB,
    NI_pass_raw = lower_97_5 > -NI_margin
  ) |>
  dplyr::group_by(Metric) |>
  dplyr::mutate(
    p_NI_adj    = {
      p <- p_NI; idx <- which(is.finite(p))
      out <- rep(NA_real_, length(p))
      if (length(idx)) out[idx] <- p.adjust(p[idx], method = adjust_method)
      out
    },
    NI_pass_adj = p_NI_adj < alpha_one_sided_NI
  ) |>
  dplyr::ungroup()

cat("\n=== NI results | δ = 0.03 | α = 0.025 | BH within metric ===\n")
pretty_ni <- nb_ni |>
  mutate(
    mean_full   = round(mean_full, 4),
    mean_triage = round(mean_triage, 4),
    diff_mean   = round(diff_mean, 4),
    se_NB       = round(se_NB, 4),
    t_NI        = round(t_NI, 3),
    p_NI        = round(p_NI, 4),
    p_NI_adj    = round(p_NI_adj, 4),
    lower_97_5  = round(lower_97_5, 4)
  ) |>
  select(Metric, Algorithm, mean_full, mean_triage, diff_mean, se_NB, df,
         lower_97_5, NI_margin, t_NI, p_NI, p_NI_adj, NI_pass_raw, NI_pass_adj) |>
  arrange(Metric, Algorithm)
print(pretty_ni, n = nrow(pretty_ni))

cat("\n=== NI pass counts after BH within metric ===\n")
print(dplyr::count(nb_ni, Metric, NI_pass_adj), n = nrow(dplyr::count(nb_ni, Metric, NI_pass_adj)))
