coerce_numeric_robust <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x1 <- trimws(as.character(x))
  x1 <- gsub(",", ".", x1, fixed = TRUE)
  x1 <- gsub("[^0-9eE+\\-\\.]", "", x1)
  suppressWarnings(as.numeric(x1))
}

coerce_sex_to_01 <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x0 <- trimws(tolower(as.character(x)))
  out <- rep(NA_real_, length(x0))
  out[x0 %in% c("m","male","man","1")] <- 1
  out[x0 %in% c("f","female","woman","0")] <- 0
  out
}

mcc_from_counts <- function(cc) {
  TP <- as.numeric(cc["TP"]); TN <- as.numeric(cc["TN"])
  FP <- as.numeric(cc["FP"]); FN <- as.numeric(cc["FN"])
  num <- TP * TN - FP * FN
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (!is.finite(den) || is.na(den) || den == 0) return(NA_real_)
  num / den
}

seed_from_name <- function(name, offset = 0L) {
  as.integer(offset + sum(utf8ToInt(as.character(name))))
}

pretty_feature_label <- function(x, style = c("heatmap", "parsimony")) {
  style <- match.arg(style)
  x <- as.character(x)
  
  if (style == "parsimony") {
    out <- dplyr::case_when(
      x == "age_raw"        ~ "Age",
      x == "sex_raw"        ~ "Sex",
      x == "spo2_raw"       ~ "SpO2",
      x == "resp_rate_raw"  ~ "Respiratory rate",
      x == "heart_rate_raw" ~ "Heart rate",
      x == "temp_raw"       ~ "Temperature",
      x == "sbp_raw"        ~ "Systolic BP",
      x == "dbp_raw"        ~ "Diastolic BP",
      x == "map_raw"        ~ "Mean arterial pressure",
      
      x == "shock_index"    ~ "Shock index",
      x == "pulse_pressure" ~ "Pulse pressure",
      x == "nlr"            ~ "Neutrophil:lymphocyte ratio",
      x == "mlr"            ~ "Monocyte:lymphocyte ratio",
      x == "plr"            ~ "Platelet:lymphocyte ratio",
      x == "siri"           ~ "Systemic inflammation response index",
      x == "sii"            ~ "Systemic immune inflammation index",
      
      x == "wbc_raw"        ~ "White blood cells",
      x == "neutrophils_raw"~ "Neutrophils",
      x == "lymphocytes_raw"~ "Lymphocytes",
      x == "monocytes_raw"  ~ "Monocytes",
      x == "platelets_raw"  ~ "Platelets",
      x == "hemoglobin_raw" ~ "Hemoglobin",
      x == "creatinine_raw" ~ "Creatinine",
      x == "urea_raw"       ~ "Urea",
      x == "sodium_raw"     ~ "Sodium",
      x == "potassium_raw"  ~ "Potassium",
      x == "chloride_raw"   ~ "Chloride",
      x == "anion_gap_raw"  ~ "Anion gap",
      x == "glucose_raw"    ~ "Glucose",
      x == "crp_raw"        ~ "CRP",
      x == "d_dimer_raw"    ~ "D-dimer",
      x == "ast_raw"        ~ "AST",
      x == "alt_raw"        ~ "ALT",
      x == "lactate_raw"    ~ "Lactate",
      x == "troponin_raw"   ~ "Troponin",
      
      TRUE ~ {
        z <- gsub("_raw$", "", x)
        z <- gsub("_", " ", z)
        tools::toTitleCase(z)
      }
    )
    return(out)
  }
  
  base <- x
  base <- gsub("_raw_rank$", "", base)
  base <- gsub("_rank$", "", base)
  base <- gsub("_raw$", "", base)
  base <- gsub("\\.", "_", base)
  base <- gsub("__+", "_", base)
  
  out <- base
  
  is_log1p <- grepl("^log1p", out, ignore.case = TRUE)
  if (any(is_log1p)) {
    b2 <- out[is_log1p]
    b2 <- sub("^log1p[_]*", "", b2, ignore.case = TRUE)
    b2 <- gsub("([a-z])([A-Z])", "\\1 \\2", b2)
    b2 <- gsub("_", " ", b2, fixed = TRUE)
    b2 <- trimws(b2)
    b2 <- tools::toTitleCase(tolower(b2))
    out[is_log1p] <- paste0("log(1 + ", b2, ")")
  }
  
  map_special <- c(
    "age"            = "Age",
    "sex"            = "Sex, M=1",
    "spo2"           = "SpO\u2082",
    "rr_spo2"        = "RR/SpO\u2082",
    "rrspo2"         = "RR/SpO\u2082",
    "heart_rate"     = "Heart rate",
    "resp_rate"      = "Respiratory rate",
    "temperature"    = "Temperature",
    "temp"           = "Temperature",
    "sbp"            = "Systolic BP",
    "dbp"            = "Diastolic BP",
    "map"            = "MAP",
    "pulsepressure"  = "Pulse pressure",
    "pulse_pressure" = "Pulse pressure",
    "shockindex"     = "Shock index",
    "shock_index"    = "Shock index",
    "aniongap"       = "Anion gap",
    "chloride"       = "Chloride",
    "creatinine"     = "Creatinine",
    "glucose"        = "Glucose",
    "hemoglobin"     = "Hemoglobin",
    "lymphocytes"    = "Lymphocytes",
    "neutrophils"    = "Neutrophils",
    "potassium"      = "Potassium",
    "sodium"         = "Sodium",
    "urea"           = "Urea",
    "nlr"            = "NLR"
  )
  
  key <- tolower(base)
  key <- gsub("[^a-z0-9_]", "", key)
  key <- gsub("_+", "_", key)
  
  mapped <- unname(map_special[match(key, names(map_special))])
  has_map <- !is.na(mapped)
  out[has_map] <- mapped[has_map]
  
  rest <- !has_map & !is_log1p
  if (any(rest)) {
    tmp <- base[rest]
    tmp <- gsub("([a-z])([A-Z])", "\\1 \\2", tmp)
    tmp <- gsub("_", " ", tmp, fixed = TRUE)
    tmp <- trimws(tmp)
    tmp <- tools::toTitleCase(tolower(tmp))
    tmp <- gsub("\\bBp\\b", "BP", tmp)
    tmp <- gsub("\\bMap\\b", "MAP", tmp)
    tmp <- gsub("\\bSpo2\\b", "SpO\u2082", tmp)
    tmp <- gsub("\\bRr\\b", "RR", tmp)
    out[rest] <- tmp
  }
  
  out
}
