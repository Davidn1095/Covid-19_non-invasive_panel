suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})

root <- "~/20250610_conchi"
WINDOW_HOURS <- 24

# =========================
# Helpers
# =========================
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fread_any <- function(path, ...) {
  stopifnot(file.exists(path))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    if (requireNamespace("R.utils", quietly = TRUE)) {
      return(data.table::fread(path, ...))
    }
    if (nzchar(Sys.which("zcat")))  return(data.table::fread(cmd = paste("zcat",  shQuote(path)), ...))
    if (nzchar(Sys.which("gzcat"))) return(data.table::fread(cmd = paste("gzcat", shQuote(path)), ...))
    stop("Cannot read .gz: install R.utils or ensure zcat/gzcat is available.")
  }
  data.table::fread(path, ...)
}

normalize_time_string <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaT", "NULL", "null", "nan", "NaN")] <- NA_character_
  x <- gsub("Z$", "", x)
  x <- gsub("T", " ", x, fixed = TRUE)
  x
}

parse_time <- function(x, tz = "UTC") {
  x0 <- normalize_time_string(x)
  out <- rep(as.POSIXct(NA, tz = tz), length(x0))

  is_numish <- !is.na(x0) & grepl("^[0-9]+(\\.[0-9]+)?$", x0)
  if (any(is_numish)) {
    xn <- suppressWarnings(as.numeric(x0[is_numish]))
    xs <- ifelse(xn > 1e12, xn / 1000, xn)
    out[is_numish] <- as.POSIXct(xs, origin = "1970-01-01", tz = tz)
  }

  fmts <- c(
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M:%OS",
    "%Y-%m-%d %H:%M",
    "%Y/%m/%d %H:%M:%S",
    "%Y/%m/%d %H:%M:%OS",
    "%m/%d/%Y %H:%M:%S",
    "%m/%d/%Y %H:%M:%OS"
  )

  for (f in fmts) {
    idx <- is.na(out) & !is.na(x0)
    if (!any(idx)) next
    tmp <- suppressWarnings(as.POSIXct(x0[idx], tz = tz, format = f))
    out[idx] <- tmp
  }

  out
}

harmonize_temp_c <- function(x) {
  x <- to_num(x)
  ok <- !is.na(x)
  if (sum(ok) >= 100) {
    frac_gt60 <- mean(x[ok] > 60)
    if (is.finite(frac_gt60) && frac_gt60 > 0.80) x <- (x - 32) * 5/9
  }
  x
}

disposition_to_outcome <- function(x) {
  z <- stringr::str_to_lower(trimws(as.character(x)))
  z[is.na(z) | !nzchar(z)] <- NA_character_

  dplyr::case_when(
    is.na(z) ~ NA_character_,
    stringr::str_detect(z, "admit|admitted|observation|\\bobs\\b|transfer|icu") ~ "Yes",
    stringr::str_detect(z, "discharg|home|left|elop|ama|against\\s+medical") ~ "No",
    TRUE ~ NA_character_
  )
}

resolve_itemid <- function(d_labitems, prefer_itemid = NA_integer_, label_pat) {
  d_labitems <- dplyr::mutate(d_labitems, itemid = as.integer(.data$itemid))
  if (!is.na(prefer_itemid) && prefer_itemid %in% d_labitems$itemid) return(as.integer(prefer_itemid))

  hit <- d_labitems %>%
    dplyr::mutate(label_l = stringr::str_to_lower(.data$label)) %>%
    dplyr::filter(stringr::str_detect(.data$label_l, label_pat)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(.data$itemid)

  if (length(hit) == 0L) NA_integer_ else as.integer(hit[[1]])
}

missing_rate <- function(df) {
  tibble::tibble(
    variable = names(df),
    missing  = vapply(df, function(z) sum(is.na(z)), numeric(1)),
    total    = nrow(df)
  ) %>%
    dplyr::mutate(missing_rate = .data$missing / .data$total) %>%
    dplyr::arrange(dplyr::desc(.data$missing_rate))
}

pick_col <- function(nms, candidates, label) {
  hit <- candidates[candidates %in% nms][1]
  if (length(hit) == 0L || is.na(hit)) {
    stop(
      "Missing required column for ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n- ", paste(nms, collapse = "\n- ")
    )
  }
  hit
}

opt_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms][1]
  if (length(hit) == 0L || is.na(hit)) NA_character_ else hit
}

# =========================
# Build MIMIC development cohort
# =========================
build_mimic_development <- function(root, window_hours = 24) {
  f_edstays  <- file.path(root, "mimic-iv-ed-2.2/ed/edstays.csv.gz")
  f_triage   <- file.path(root, "mimic-iv-ed-2.2/ed/triage.csv.gz")
  f_patients <- file.path(root, "mimic-iv-hosp-2.2/patients.csv.gz")
  f_labs     <- file.path(root, "mimic-iv-hosp-2.2/labevents.csv.gz")
  f_items    <- file.path(root, "mimic-iv-hosp-2.2/d_labitems.csv.gz")
  stopifnot(file.exists(f_edstays), file.exists(f_triage), file.exists(f_patients), file.exists(f_labs), file.exists(f_items))

  ed <- fread_any(
    f_edstays,
    select = c("subject_id", "hadm_id", "stay_id", "intime", "gender", "disposition"),
    na.strings = c("", "NA")
  )
  ed[, intime := parse_time(intime)]

  tri <- fread_any(
    f_triage,
    select = c("subject_id", "stay_id", "temperature", "heartrate", "resprate", "o2sat", "sbp", "dbp"),
    na.strings = c("", "NA")
  )

  pat <- fread_any(
    f_patients,
    select = c("subject_id", "anchor_age"),
    na.strings = c("", "NA")
  )

  dt <- merge(ed, tri, by = c("subject_id", "stay_id"), all.x = TRUE)
  dt <- merge(dt, pat, by = "subject_id", all.x = TRUE)

  dt[, `:=`(
    dataset        = "MIMIC",
    encounter_id   = as.character(stay_id),
    outcome        = disposition_to_outcome(disposition),
    age_raw        = to_num(anchor_age),
    sex_raw        = fifelse(gender %in% c("M", "m"), 1, fifelse(gender %in% c("F", "f"), 0, NA_real_)),
    temp_raw       = harmonize_temp_c(temperature),
    heart_rate_raw = to_num(heartrate),
    resp_rate_raw  = to_num(resprate),
    spo2_raw       = to_num(o2sat),
    sbp_raw        = to_num(sbp),
    dbp_raw        = to_num(dbp)
  )]
  dt[, map_raw := fifelse(!is.na(sbp_raw) & !is.na(dbp_raw), (2 * dbp_raw + sbp_raw) / 3, NA_real_)]
  dt <- dt[!is.na(intime) & !is.na(outcome)]

  d_labitems <- fread_any(
    f_items,
    select = c("itemid", "label"),
    na.strings = c("", "NA")
  ) %>% dplyr::as_tibble()

  item_ids <- tibble::tibble(
    feature = c(
      "neutrophils_raw","lymphocytes_raw","hemoglobin_raw",
      "creatinine_raw","glucose_raw","urea_raw",
      "sodium_raw","potassium_raw","anion_gap_raw","chloride_raw",
      "bilirubin_raw","albumin_raw","alk_phos_raw"
    ),
    prefer = c(
      51256L,51244L,51222L,
      50912L,50931L,51006L,
      50983L,50971L,50868L,50902L,
      50885L,50862L,50863L
    ),
    label_pat = c(
      "neutrophils",
      "lymphocytes",
      "hemoglobin",
      "creatinine",
      "glucose",
      "urea|bun",
      "\\bsodium\\b",
      "\\bpotassium\\b",
      "anion\\s*gap",
      "\\bchloride\\b",
      "bilirubin",
      "\\balbumin\\b",
      "alkaline\\s*phosphatase|alk\\s*phos"
    )
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(itemid = resolve_itemid(d_labitems, .data$prefer, .data$label_pat)) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$itemid, .data$feature) %>%
    dplyr::filter(!is.na(.data$itemid))

  keep_itemids <- unique(as.integer(item_ids$itemid))
  if (length(keep_itemids) == 0L) stop("No lab itemids resolved from d_labitems.")

  key_stays <- dt[, .(subject_id, stay_id, intime)]
  key_stays[, t1 := intime + window_hours * 3600]
  key_stays[, `:=`(start = intime, end = t1)]
  key_stays <- key_stays[!is.na(start) & !is.na(end)]
  setkey(key_stays, subject_id, start, end)

  labs <- fread_any(
    f_labs,
    select = c("subject_id", "itemid", "charttime", "valuenum"),
    na.strings = c("", "NA")
  )
  labs <- labs[itemid %in% keep_itemids]
  labs[, charttime := parse_time(charttime)]
  labs[, valuenum  := to_num(valuenum)]
  labs <- labs[!is.na(charttime) & !is.na(valuenum)]
  labs[, `:=`(start = charttime, end = charttime)]
  setkey(labs, subject_id, start, end)

  lab_join <- data.table::foverlaps(
    x = labs,
    y = key_stays,
    by.x = c("subject_id", "start", "end"),
    by.y = c("subject_id", "start", "end"),
    type = "within",
    nomatch = 0L
  )

  item_ids_dt <- as.data.table(item_ids)
  setkey(item_ids_dt, itemid)
  lab_join <- merge(lab_join, item_ids_dt, by = "itemid", all.x = TRUE)
  lab_join <- lab_join[!is.na(feature)]

  setorder(lab_join, stay_id, feature, charttime)
  labs_first <- lab_join[, .SD[1], by = .(stay_id, feature)]
  labs_wide  <- dcast(labs_first, stay_id ~ feature, value.var = "valuenum")

  dt2 <- merge(dt, labs_wide, by = "stay_id", all.x = TRUE)

  out <- dt2 %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(
      dataset,
      encounter_id,
      outcome,
      age_raw,
      sex_raw,
      spo2_raw,
      heart_rate_raw,
      resp_rate_raw,
      temp_raw,
      sbp_raw,
      dbp_raw,
      map_raw,
      neutrophils_raw,
      lymphocytes_raw,
      hemoglobin_raw,
      creatinine_raw,
      glucose_raw,
      urea_raw,
      sodium_raw,
      potassium_raw,
      anion_gap_raw,
      chloride_raw,
      bilirubin_raw,
      albumin_raw,
      alk_phos_raw
    )

  for (nm in item_ids$feature) if (!nm %in% names(out)) out[[nm]] <- NA_real_
  out
}

# =========================
# Build MC-MED external cohort
# Key change: strict component selection to avoid unit mixing
# =========================
build_mcmed_external <- function(root, window_hours = 24) {
  f_visits <- file.path(root, "mc-med-1.0.1/data/visits.csv")
  f_labs   <- file.path(root, "mc-med-1.0.1/data/labs.csv")
  stopifnot(file.exists(f_visits), file.exists(f_labs))

  v0 <- data.table::fread(f_visits, na.strings = c("", "NA"))
  l0 <- data.table::fread(f_labs,   na.strings = c("", "NA"))

  # required visit columns
  csn_col   <- pick_col(names(v0), c("CSN","csn","Encounter_ID","encounter_id"), "visit encounter id")
  arr_col   <- pick_col(names(v0), c("Arrival_time","ArrivalTime","arrival_time"), "Arrival_time")
  admit_col <- pick_col(names(v0), c("Admit_time","AdmitTime","admit_time"), "Admit_time")
  age_col   <- pick_col(names(v0), c("Age","age"), "Age")
  sex_col   <- pick_col(names(v0), c("Gender","gender","Sex","sex"), "Gender")
  spo2_col  <- pick_col(names(v0), c("Triage_SpO2","Triage_SPO2","SpO2","spo2"), "Triage SpO2")
  hr_col    <- pick_col(names(v0), c("Triage_HR","HR","HeartRate","heartrate"), "Triage HR")
  rr_col    <- pick_col(names(v0), c("Triage_RR","RR","RespRate","resprate"), "Triage RR")
  temp_col  <- pick_col(names(v0), c("Triage_Temp","Temp","Temperature","temperature"), "Triage Temp")
  sbp_col   <- pick_col(names(v0), c("Triage_SBP","SBP","sbp"), "Triage SBP")
  dbp_col   <- pick_col(names(v0), c("Triage_DBP","DBP","dbp"), "Triage DBP")

  # optional visit columns
  hospid_col <- opt_col(names(v0), c("Hospital_ID","HospitalID","Hosp_ID","Admission_ID","Admit_ID"))
  los_col    <- opt_col(names(v0), c("Hosp_LOS","Hospital_LOS","LOS","Length_of_stay"))
  dispo_col  <- opt_col(names(v0), c("ED_dispo","ED_Dispo","Disposition","ED_disposition"))

  nV <- nrow(v0)

  hospid_present <- if (!is.na(hospid_col)) {
    hv <- v0[[hospid_col]]
    if (is.numeric(hv) || is.integer(hv)) !is.na(hv) & hv > 0
    else {
      hs <- trimws(as.character(hv))
      !is.na(hs) & nzchar(hs)
    }
  } else rep(FALSE, nV)

  los_num   <- if (!is.na(los_col)) to_num(v0[[los_col]]) else rep(NA_real_, nV)
  dispo_str <- if (!is.na(dispo_col)) str_to_lower(as.character(v0[[dispo_col]])) else rep(NA_character_, nV)

  v <- v0 %>%
    dplyr::mutate(
      Arrival_time = parse_time(.data[[arr_col]]),
      Admit_time   = parse_time(.data[[admit_col]]),
      age_raw      = to_num(.data[[age_col]]),
      sex_raw      = dplyr::case_when(
        as.character(.data[[sex_col]]) %in% c("M","Male","m","male") ~ 1,
        as.character(.data[[sex_col]]) %in% c("F","Female","f","female") ~ 0,
        TRUE ~ NA_real_
      ),
      spo2_raw       = to_num(.data[[spo2_col]]),
      heart_rate_raw = to_num(.data[[hr_col]]),
      resp_rate_raw  = to_num(.data[[rr_col]]),
      temp_raw       = harmonize_temp_c(.data[[temp_col]]),
      sbp_raw        = to_num(.data[[sbp_col]]),
      dbp_raw        = to_num(.data[[dbp_col]]),
      map_raw        = dplyr::if_else(
        !is.na(.data$sbp_raw) & !is.na(.data$dbp_raw),
        (2 * .data$dbp_raw + .data$sbp_raw) / 3,
        NA_real_
      ),
      outcome = dplyr::case_when(
        hospid_present ~ "Yes",
        !is.na(.data$Admit_time) ~ "Yes",
        !is.na(los_num) & los_num > 0 ~ "Yes",
        !is.na(dispo_str) & stringr::str_detect(dispo_str, "admit|observation|\\bobs\\b") ~ "Yes",
        TRUE ~ "No"
      )
    ) %>%
    dplyr::transmute(
      dataset      = "MC-MED",
      encounter_id = as.character(.data[[csn_col]]),
      outcome,
      Arrival_time,
      age_raw,
      sex_raw,
      spo2_raw,
      heart_rate_raw,
      resp_rate_raw,
      temp_raw,
      sbp_raw,
      dbp_raw,
      map_raw
    ) %>%
    dplyr::filter(!is.na(.data$encounter_id), !is.na(.data$Arrival_time), !is.na(.data$outcome))

  # required lab columns
  l_csn   <- pick_col(names(l0), c("CSN","csn","Encounter_ID","encounter_id"), "lab encounter id")
  l_time  <- pick_col(names(l0), c("Result_time","ResultTime","result_time"), "Result_time")
  l_name  <- pick_col(names(l0), c("Component_name","ComponentName","component_name"), "Component_name")
  l_value <- pick_col(names(l0), c("Component_value","ComponentValue","component_value","Value","value"), "Component_value")
  l_unit  <- opt_col(names(l0), c("Component_units","component_units","Units","units","Unit","unit"))

  l2 <- l0 %>%
    dplyr::mutate(
      encounter_id = as.character(.data[[l_csn]]),
      Result_time  = parse_time(.data[[l_time]]),
      Component_u  = stringr::str_to_upper(as.character(.data[[l_name]])),
      unit_u       = if (!is.na(l_unit)) stringr::str_to_upper(as.character(.data[[l_unit]])) else NA_character_,
      value_num    = to_num(.data[[l_value]])
    ) %>%
    dplyr::filter(!is.na(.data$encounter_id), !is.na(.data$Result_time), !is.na(.data$value_num), !is.na(.data$Component_u))

  # Priority-based component selection per feature (avoid wrong analytes / units)
  spec <- tibble::tibble(
    feature = c(
      "neutrophils_raw","lymphocytes_raw","hemoglobin_raw",
      "creatinine_raw","glucose_raw","urea_raw",
      "sodium_raw","potassium_raw","anion_gap_raw","chloride_raw",
      "bilirubin_raw","albumin_raw","alk_phos_raw"
    ),
    unit_expected = c(
      "%","%","G/DL",
      "MG/DL","MG/DL","MG/DL",
      "MMOL/L","MMOL/L","MMOL/L","MMOL/L",
      "MG/DL","G/DL","U/L"
    ),
    # ordered patterns, first match is preferred
    pat1 = c(
      "^NEUTROPHIL % \\(AUTO DIFF\\)$",
      "^LYMPHOCYTE % \\(AUTO DIFF\\)$",
      "^HEMOGLOBIN \\(HGB\\)$",
      "^CREATININE$",
      "^GLUCOSE$",
      "^BLOOD UREA NITROGEN \\(BUN\\)$",
      "^SODIUM$",
      "^POTASSIUM$",
      "^ANION GAP$",
      "^CHLORIDE$",
      "^BILIRUBIN, TOTAL$",
      "^ALBUMIN$",
      "^ALKALINE PHOSPHATASE$"
    ),
    pat2 = c(
      "^NEUTROPHILS, SEGMENTED % \\(MANUAL DIFF\\)$",
      "^LYMPHOCYTES % \\(MANUAL DIFF\\)$",
      NA_character_,
      "^POC:CREATININE",
      "^POC:GLUCOSE",
      NA_character_,
      "^POC:SODIUM",
      "^POC:POTASSIUM",
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_
    )
  )

  assign_rank <- function(comp, p1, p2) {
    r <- rep(Inf, length(comp))
    if (!is.na(p1) && nzchar(p1)) r[str_detect(comp, p1)] <- 1
    if (!is.na(p2) && nzchar(p2)) r[str_detect(comp, p2)] <- pmin(r[str_detect(comp, p2)], 2)
    r
  }

  extract_feature <- function(feat, unit_expected, p1, p2) {
    df <- l2

    # reduce obvious wrong contexts first
    df <- df %>%
      dplyr::filter(
        !stringr::str_detect(.data$Component_u, "URINE|BODY FLU|BODY FLD|CSF|SYNOVIAL|RATIO|STOOL|ISOENZYMES|BONE|LIVER"),
        !stringr::str_detect(.data$Component_u, "RAW COUNT"),
        !stringr::str_detect(.data$Component_u, "ABSOLUTE") | feat %in% c("creatinine_raw","glucose_raw") # keep POC for those via rank
      )

    # candidate patterns
    keep <- stringr::str_detect(df$Component_u, p1)
    if (!is.na(p2) && nzchar(p2)) keep <- keep | stringr::str_detect(df$Component_u, p2)
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) return(tibble::tibble())

    # enforce unit when present
    if (!is.na(l_unit) && !is.na(unit_expected) && nzchar(unit_expected)) {
      df <- df %>% dplyr::filter(!is.na(.data$unit_u), .data$unit_u == unit_expected)
    }

    # value sanity for percentages
    if (feat %in% c("neutrophils_raw","lymphocytes_raw")) {
      df <- df %>% dplyr::filter(.data$value_num >= 0, .data$value_num <= 100)
    }

    df$rank <- assign_rank(df$Component_u, p1, p2)

    # keep within arrival window, pick preferred component then earliest time
    out <- df %>%
      dplyr::inner_join(v %>% dplyr::select(.data$encounter_id, .data$Arrival_time), by = "encounter_id") %>%
      dplyr::filter(
        .data$Result_time >= .data$Arrival_time,
        .data$Result_time <= .data$Arrival_time + window_hours * 3600
      ) %>%
      dplyr::arrange(.data$encounter_id, .data$rank, .data$Result_time) %>%
      dplyr::group_by(.data$encounter_id) %>%
      dplyr::summarise(value = dplyr::first(.data$value_num), .groups = "drop") %>%
      dplyr::mutate(feature = feat)

    out
  }

  keepers <- dplyr::bind_rows(lapply(seq_len(nrow(spec)), function(i) {
    extract_feature(spec$feature[i], spec$unit_expected[i], spec$pat1[i], spec$pat2[i])
  }))

  lab_base <- v %>% dplyr::distinct(.data$encounter_id)

  if (nrow(keepers) > 0) {
    lab_vals <- keepers %>%
      tidyr::pivot_wider(names_from = .data$feature, values_from = .data$value)
  } else {
    lab_vals <- lab_base
  }

  # ensure all expected columns exist
  lab_wide <- lab_base %>% dplyr::left_join(lab_vals, by = "encounter_id")
  for (nm in spec$feature) if (!nm %in% names(lab_wide)) lab_wide[[nm]] <- NA_real_

  v %>%
    dplyr::left_join(lab_wide, by = "encounter_id") %>%
    dplyr::select(
      .data$dataset,
      .data$encounter_id,
      .data$outcome,
      .data$age_raw,
      .data$sex_raw,
      .data$spo2_raw,
      .data$heart_rate_raw,
      .data$resp_rate_raw,
      .data$temp_raw,
      .data$sbp_raw,
      .data$dbp_raw,
      .data$map_raw,
      dplyr::all_of(spec$feature)
    )
}

# =========================
# Build datasets
# =========================
development_mimic_raw <- build_mimic_development(root = root, window_hours = WINDOW_HOURS)
external_mcmed_raw    <- build_mcmed_external(root = root, window_hours = WINDOW_HOURS)

cat("\n=== Outcome counts: Development (MIMIC) ===\n")
print(dplyr::count(development_mimic_raw, .data$outcome))

cat("\n=== Outcome counts: External (MC-MED) ===\n")
print(dplyr::count(external_mcmed_raw, .data$outcome))

cat("\n=== Missingness: Development (MIMIC), top 15 ===\n")
print(missing_rate(development_mimic_raw) %>% dplyr::slice(1:15))

cat("\n=== Missingness: External (MC-MED), top 15 ===\n")
print(missing_rate(external_mcmed_raw) %>% dplyr::slice(1:15))

# =========================
# Row-completeness filter (>= 60% non-missing across predictors)
# =========================
row_complete_frac <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0L) return(rep(1, nrow(df)))
  X <- as.data.frame(df)[, cols, drop = FALSE]
  mean_non_missing <- function(v) mean(!is.na(v))
  as.numeric(apply(X, 1, mean_non_missing))
}

filter_by_row_completeness <- function(df, frac = 0.60, exclude_cols = c("dataset","encounter_id","outcome","Arrival_time")) {
  cols <- setdiff(names(df), exclude_cols)
  rc <- row_complete_frac(df, cols)
  keep <- rc >= frac
  list(df = df[keep, , drop = FALSE], keep = keep, row_complete = rc)
}

dev_f60 <- filter_by_row_completeness(development_mimic_raw, frac = 0.60)
ext_f60 <- filter_by_row_completeness(external_mcmed_raw,    frac = 0.60)

df_dev_60 <- dev_f60$df
df_ext_60 <- ext_f60$df

cat("\n=== Kept rows (>=60% complete predictors) ===\n")
cat("DEV:", nrow(df_dev_60), "/", nrow(development_mimic_raw), "\n")
cat("EXT:", nrow(df_ext_60), "/", nrow(external_mcmed_raw), "\n")

cat("\n=== Quick quantiles check (DEV vs EXT) ===\n")
quick_q <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, 7))
  as.numeric(stats::quantile(x, probs = c(0, .01, .05, .5, .95, .99, 1), na.rm = TRUE))
}
feat_check <- c("hemoglobin_raw","neutrophils_raw","lymphocytes_raw","anion_gap_raw","albumin_raw","bilirubin_raw","alk_phos_raw")
for (f in feat_check) {
  if (!f %in% names(df_dev_60) || !f %in% names(df_ext_60)) next
  cat("\n---", f, "---\n")
  cat("DEV:", paste(round(quick_q(df_dev_60[[f]]), 3), collapse = " "), "\n")
  cat("EXT:", paste(round(quick_q(df_ext_60[[f]]), 3), collapse = " "), "\n")
}

# =========================
# Save outputs only if WRITE_OUT=1
# =========================
if (identical(Sys.getenv("WRITE_OUT"), "1")) {
  dev_raw_out <- "20260106_development_mimic_raw.csv"
  ext_raw_out <- "20260106_external_mcmed_raw.csv"
  dev_60_out  <- "20260106_development_mimic_raw_rowcomplete60.csv"
  ext_60_out  <- "20260106_external_mcmed_raw_rowcomplete60.csv"

  readr::write_csv(development_mimic_raw, dev_raw_out, na = "")
  readr::write_csv(external_mcmed_raw,    ext_raw_out, na = "")
  readr::write_csv(df_dev_60, dev_60_out, na = "")
  readr::write_csv(df_ext_60, ext_60_out, na = "")

  cat("\nSaved:\n- ", dev_raw_out, "\n- ", ext_raw_out, "\n- ", dev_60_out, "\n- ", ext_60_out, "\n", sep = "")
}
