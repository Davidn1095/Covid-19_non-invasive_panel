# ======================= TABLE 1 + TABLE (Missingness & Availability) =======================
# Builds exactly:
#   • tbl1_baseline            — Baseline characteristics by data_split (numeric summaries + factor %)
#   • tbl_missingness          — Missingness & availability per predictor (overall + by split)
# Assumes your data uses Group ∈ {CTRL_noCOVID, COVID}. Outcome coerced to {Yes, No}.
# Optional: enforce a temporal split at 21-Apr-2020 if you provide the admission date column name.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(writexl)   # only if you want to export at the end
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))

# ---------------- config ----------------
cfg <- list(
  file_path = "biomarkers_acuteCOVID_meta.xlsx",
  sheet     = "meta",
  outcome   = "Hospital_ID",   # factor {Yes, No} after coercion
  pos_label = "Yes",
  neg_label = "No"
)

# ---------------- feature sets used for missingness table ----------------
feat_full <- c(
  "Diagnosis","severity_admission","Age","Gender","SpO2_admission",
  "monocyte_abs_number","monocytes_perc",
  "lymphocyte_abs_number","lymphocytes_perc",
  "neutrophil_abs_number","neutrophils_perc"
)
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- loader (as in your pipeline, minimal) ----------------
load_data <- function(path, sheet, outcome, keep) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet = sheet) |> as.data.frame()
  
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found, filtering cannot be applied")
  names(df)[names(df) == grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  df$.rid <- seq_len(nrow(df))
  
  guess_factor <- intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))
  for (v in guess_factor) if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' not in data", outcome))
  y <- df[[outcome]]
  to_YN <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) x <- ifelse(x, cfg$pos_label, cfg$neg_label)
    if (is.numeric(x)) {
      if (!all(x %in% c(0, 1), na.rm = TRUE)) stop("Outcome numeric but not 0/1")
      out <- ifelse(x == 1, cfg$pos_label, cfg$neg_label)
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    if (is.character(x)) {
      s <- trimws(tolower(x))
      map_yes <- c("1","yes","y","true","pos","positive")
      map_no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% map_yes, cfg$pos_label,
                    ifelse(s %in% map_no,  cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome must be Yes/No, Y/N, TRUE/FALSE, or 0/1 after trimming")
      return(factor(out, levels = c(cfg$pos_label, cfg$neg_label)))
    }
    stop("Outcome type unsupported for coercion")
  }
  df[[outcome]] <- to_YN(y)
  
  keepx <- unique(c(keep, outcome, "data_split", "group", ".rid"))
  keepx <- intersect(keepx, names(df))
  df <- df[, keepx, drop = FALSE]
  df
}

# ---------------- optional: temporal split helper (set admit_col or skip) ----------------
enforce_temporal_split <- function(df, admit_col = NULL, boundary = as.Date("2020-04-21")) {
  if (is.null(admit_col)) return(df)  # no change if not provided
  if (!(admit_col %in% names(df))) stop(sprintf("Column '%s' not found in df.", admit_col))
  x <- df[[admit_col]]
  to_date <- function(z) {
    if (inherits(z, "Date")) return(z)
    if (inherits(z, "POSIXt")) return(as.Date(z))
    zc <- as.character(z)
    d  <- suppressWarnings(as.Date(zc))
    if (all(is.na(d))) d <- suppressWarnings(lubridate::ymd(zc))
    if (all(is.na(d))) d <- suppressWarnings(lubridate::dmy(zc))
    if (all(is.na(d))) d <- suppressWarnings(lubridate::mdy(zc))
    as.Date(d)
  }
  dts <- to_date(x)
  if (anyNA(dts)) stop("Admission dates could not be parsed cleanly.")
  df$data_split <- factor(ifelse(dts <= boundary, "train", "external"), levels = c("train","external"))
  df
}

# ---------------- Table builders ----------------
baseline_table <- function(df, yvar, split_col = "data_split") {
  df2 <- df
  drop_cols <- c(".rid", yvar, split_col)
  drop_cols <- intersect(drop_cols, names(df2))
  df2 <- df2[, setdiff(names(df2), drop_cols), drop = FALSE]
  has_split <- split_col %in% names(df)
  
  num_cols <- names(df2)[sapply(df2, is.numeric)]
  fac_cols <- names(df2)[sapply(df2, is.factor)]
  
  num_tbl <- if (length(num_cols)) {
    agg <- function(z) c(
      N = sum(is.finite(z)),
      Mean = mean(z, na.rm = TRUE),
      SD = stats::sd(z, na.rm = TRUE),
      Median = stats::median(z, na.rm = TRUE),
      Q1 = stats::quantile(z, 0.25, na.rm = TRUE),
      Q3 = stats::quantile(z, 0.75, na.rm = TRUE)
    )
    if (has_split) {
      sp_names <- unique(as.character(df[[split_col]]))
      parts <- lapply(sp_names, function(s) {
        idx <- df[[split_col]] == s
        as.data.frame(t(vapply(df2[idx, num_cols, drop=FALSE], agg, numeric(6))))
      })
      names(parts) <- sp_names
      parts <- Map(function(x, nm) { colnames(x) <- paste0(colnames(x), "_", nm); x }, parts, names(parts))
      tbl <- Reduce(function(a,b) cbind(a,b), parts)
      tibble::tibble(Variable = rownames(tbl)) |> cbind(tbl, row.names = NULL)
    } else {
      tbl <- as.data.frame(t(vapply(df2[, num_cols, drop=FALSE], agg, numeric(6))))
      tibble::tibble(Variable = rownames(tbl)) |> cbind(tbl, row.names = NULL)
    }
  } else tibble::tibble()
  
  fac_tbl <- if (length(fac_cols)) {
    count_one <- function(var) {
      if (has_split) {
        by(df[[var]], df[[split_col]], function(x) round(100*prop.table(table(x)), 1)) |>
          lapply(function(pct) tibble::tibble(Level = names(pct), Pct = as.numeric(pct))) |>
          purrr::imap(~dplyr::rename(.x, !!paste0("pct_", .y) := Pct)) |>
          purrr::reduce(dplyr::left_join, by = "Level") |>
          dplyr::mutate(Variable = var, .before = 1)
      } else {
        pct <- round(100*prop.table(table(df[[var]])), 1)
        tibble::tibble(Variable = var, Level = names(pct), pct_overall = as.numeric(pct))
      }
    }
    purrr::map_dfr(fac_cols, count_one)
  } else tibble::tibble()
  
  dplyr::bind_rows(num_tbl, fac_tbl)
}

missingness_table <- function(df, predictors, split_col = "data_split") {
  present <- predictors[predictors %in% names(df)]
  base <- tibble::tibble(
    Predictor = present,
    Overall_N = vapply(present, function(v) length(df[[v]]), integer(1)),
    Overall_Missing_pct = vapply(present, function(v) 100*mean(is.na(df[[v]])), numeric(1))
  )
  if (split_col %in% names(df)) {
    spl <- levels(factor(df[[split_col]]))
    for (sp in spl) {
      idx <- df[[split_col]] == sp
      base[[paste0("Missing_", sp)]] <- vapply(present, function(v) 100*mean(is.na(df[[v]][idx])), numeric(1))
    }
  }
  base
}

# ============================ RUN: build exactly 2 tables ============================
df <- load_data(cfg$file_path, cfg$sheet, cfg$outcome, unique(c(feat_full, feat_triage)))

# Optional: uncomment and set your admission date column to enforce the temporal split
# df <- enforce_temporal_split(df, admit_col = "AdmissionDate")  # <-- set the exact column name

tbl1_baseline   <- baseline_table(df, cfg$outcome, split_col = "data_split")
tbl_missingness <- missingness_table(df, unique(c(feat_full, feat_triage)), split_col = "data_split")

# Quick sanity prints
cat(sprintf("\n[Table 1] Baseline: %d x %d\n", nrow(tbl1_baseline), ncol(tbl1_baseline)))
print(utils::head(tbl1_baseline, 12))
cat(sprintf("\n[Missingness] Predictors: %d x %d\n", nrow(tbl_missingness), ncol(tbl_missingness)))
print(utils::head(tbl_missingness, 12))

# Optional: export to a single Excel file
# writexl::write_xlsx(list(Table1_Baseline = tbl1_baseline,
#                          Table_Missingness = tbl_missingness),
#                     "tables_baseline_missingness.xlsx")
