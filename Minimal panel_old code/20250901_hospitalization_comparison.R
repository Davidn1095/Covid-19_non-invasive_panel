# ============================================================
# STATISTICAL ANALYSIS: paired NB comparison, no disk writes
# Requires in-memory full_folds and triage_folds from predictions
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
})

if (!exists("full_folds", inherits = TRUE) || !exists("triage_folds", inherits = TRUE)) {
  stop("full_folds and triage_folds must exist in memory. Run the predictions block first.")
}

alpha_two_sided    <- 0.05
alpha_one_sided_NI <- 0.025

nadeau_bengio_summarise <- function(full_folds, triage_folds,
                                    alpha_two_sided = 0.05, alpha_one_sided_NI = 0.025) {
  paired <- dplyr::inner_join(
    dplyr::rename(full_folds,  Value_full   = Value),
    dplyr::rename(triage_folds,Value_triage = Value),
    by = c("Algorithm","Fold","Metric")
  ) |>
    dplyr::mutate(Diff = Value_triage - Value_full)
  
  paired |>
    dplyr::group_by(Algorithm, Metric) |>
    dplyr::summarise(
      mean_full   = mean(Value_full,   na.rm = TRUE),
      mean_triage = mean(Value_triage, na.rm = TRUE),
      diff_mean   = mean(Diff,         na.rm = TRUE),
      sd_diff     = stats::sd(Diff,    na.rm = TRUE),
      N           = dplyr::n(),
      k_est       = { f_ids <- stringr::str_extract(Fold, "(?<=_f)\\d+"); length(unique(f_ids)) },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      c_NB   = 1/N + 1/pmax(k_est - 1, 1),
      se_NB  = sd_diff * sqrt(c_NB),
      df     = pmax(N - 1, 1),
      t_val  = diff_mean / pmax(se_NB, 1e-12),
      p_two  = 2 * (1 - stats::pt(abs(t_val), df)),
      tcrit95 = stats::qt(1 - alpha_two_sided/2, df),
      ci95_lo = diff_mean - tcrit95 * se_NB,
      ci95_hi = diff_mean + tcrit95 * se_NB,
      p_one_worse = ifelse(diff_mean < 0, 1 - stats::pt(t_val, df), stats::pt(t_val, df)),
      tcrit_NI = stats::qt(1 - alpha_one_sided_NI, df),
      delta_min_abs      = pmax(0, tcrit_NI * se_NB - diff_mean),
      delta_min_rel_pct  = 100 * delta_min_abs / pmax(abs(mean_full), 1e-12),
      rel_diff_pct       = 100 * diff_mean / pmax(abs(mean_full), 1e-12),
      sig_two_sided_0.05       = p_two < 0.05,
      sig_worse_one_sided_0.05 = (diff_mean < 0) & (p_one_worse < 0.05)
    )
}

nb_tbl <- nadeau_bengio_summarise(full_folds, triage_folds,
                                  alpha_two_sided = alpha_two_sided,
                                  alpha_one_sided_NI = alpha_one_sided_NI)

cat("\n=== How inferior or improved is triage vs full, triage − full ===\n")
pretty_tbl <- nb_tbl |>
  dplyr::mutate(
    mean_full    = round(mean_full, 4),
    mean_triage  = round(mean_triage, 4),
    diff_mean    = sprintf("%+.4f", diff_mean),
    rel_diff_pct = sprintf("%+.1f%%", rel_diff_pct),
    ci95         = sprintf("[%+.4f, %+.4f]", ci95_lo, ci95_hi),
    p_two        = sprintf("%.4f", p_two),
    p_one_worse  = sprintf("%.4f", p_one_worse),
    se_NB        = round(se_NB, 4),
    delta_min_abs     = round(delta_min_abs, 4),
    delta_min_rel_pct = sprintf("%.1f%%", delta_min_rel_pct)
  ) |>
  dplyr::select(
    Metric, Algorithm,
    mean_full, mean_triage,
    diff_mean, rel_diff_pct, ci95,
    p_two, sig_two_sided_0.05,
    p_one_worse, sig_worse_one_sided_0.05,
    se_NB, df,
    delta_min_abs, delta_min_rel_pct
  ) |>
  dplyr::arrange(Metric, Algorithm)

print(pretty_tbl, n = nrow(pretty_tbl))

cat("\n=== Wide view, triage − full mean differences ===\n")
pretty_wide <- nb_tbl |>
  dplyr::transmute(Metric, Algorithm, Value = sprintf("%+.4f", diff_mean)) |>
  tidyr::pivot_wider(names_from = Algorithm, values_from = Value) |>
  dplyr::arrange(Metric)
print(pretty_wide, n = nrow(pretty_wide))
