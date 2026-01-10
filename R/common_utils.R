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
