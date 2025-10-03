# ============================================================
# 3) BUILDERS FOR TABLE 2 AND TABLE 8 (overall)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(tibble); library(purrr); library(stringr)
})

stopifnot(exists("cfg"), exists("feat_full"), exists("feat_triage"))

df0 <- readxl::read_excel(cfg$file_path, sheet = cfg$sheet) %>% as_tibble()

all_features <- unique(c(feat_triage, setdiff(feat_full, feat_triage)))
missing_cols <- setdiff(all_features, names(df0))
if (length(missing_cols)) message("Warning: not in data -> ", paste(missing_cols, collapse = ", "))
present_features <- intersect(all_features, names(df0))
stopifnot(length(present_features) > 0)

is_categorical <- function(x) {
  is.factor(x) || is.logical(x) || (is.character(x) && dplyr::n_distinct(x, na.rm = TRUE) <= 25)
}
derive_type <- function(v) {
  if (is.numeric(v) || inherits(v, "Date") || inherits(v, "POSIXct")) return("numeric")
  if (is_categorical(v)) return("categorical")
  "other"
}
preproc_for_type <- function(type) {
  if (type == "numeric")     return("median impute, Yeo–Johnson, center-scale, zero-variance removal")
  if (type == "categorical") return("mode impute, novel level, rare merge (1%), one-hot, zero-variance removal")
  "inspect"
}

tbl2 <- tibble(
  Predictor = present_features,
  Set = case_when(
    Predictor %in% feat_triage ~ "Triage",
    Predictor %in% setdiff(feat_full, feat_triage) ~ "Full-only",
    TRUE ~ "Unspecified"
  ),
  Type = purrr::map_chr(present_features, ~ derive_type(df0[[.x]]))
) %>%
  mutate(Preprocessing = vapply(Type, preproc_for_type, character(1))) %>%
  arrange(factor(Set, levels = c("Triage","Full-only","Unspecified")), Predictor)

miss_summary <- function(v) {
  n <- length(v); n_mis <- sum(is.na(v)); n_avl <- n - n_mis
  tibble(
    Missing_n = n_mis,
    Missing_pct = round(100 * n_mis / n, 1),
    Available_n = n_avl,
    Available_pct = round(100 * n_avl / n, 1)
  )
}

tbl8_overall <- purrr::map_dfr(present_features, function(f) {
  bind_cols(
    tibble(Predictor = f,
           Set = if (f %in% feat_triage) "Triage" else "Full-only",
           Type = derive_type(df0[[f]])),
    miss_summary(df0[[f]])
  )
}) %>%
  arrange(factor(Set, levels = c("Triage","Full-only")), Predictor)

cat("\n=== Table 2 — Predictors & preprocessing ===\n"); print(tbl2, n = nrow(tbl2))
cat("\n=== Table 8 — Missingness & availability (overall) ===\n"); print(tbl8_overall, n = nrow(tbl8_overall))

# Objects: tbl2, tbl8_overall
