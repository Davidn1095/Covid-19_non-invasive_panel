## Run this AFTER your benchmark script has created:
## df_train, df_test, feat_triage, feat_full (and optionally eng_created / eng_skipped)

stopifnot(exists("df_train", inherits = TRUE))
stopifnot(exists("df_test",  inherits = TRUE))
stopifnot(exists("feat_triage", inherits = TRUE))
stopifnot(exists("feat_full",   inherits = TRUE))

df_train    <- get("df_train", inherits = TRUE)
df_test     <- get("df_test",  inherits = TRUE)
feat_triage <- get("feat_triage", inherits = TRUE)
feat_full   <- get("feat_full",   inherits = TRUE)

`%||%` <- function(a,b) if (!is.null(a)) a else b
eng_created <- (if (exists("eng_created", inherits = TRUE)) get("eng_created", inherits = TRUE) else NULL) %||% character(0)
eng_skipped <- (if (exists("eng_skipped", inherits = TRUE)) get("eng_skipped", inherits = TRUE) else NULL) %||% character(0)

# ----------------------------
# 1) What exactly are the attributes used?
# ----------------------------
cat("\n=== TRIAGE FEATURE SET ===\n")
cat("n =", length(feat_triage), "\n")
print(sort(feat_triage))

cat("\n=== FULL FEATURE SET ===\n")
cat("n =", length(feat_full), "\n")
print(sort(feat_full))

full_only <- setdiff(feat_full, feat_triage)
cat("\n=== FULL-ONLY (beyond triage) ===\n")
cat("n =", length(full_only), "\n")
print(sort(full_only))

# ----------------------------
# 2) Verify they are numeric in train/test
# ----------------------------
is_num_train <- vapply(df_train[, feat_full, drop = FALSE], is.numeric, logical(1))
is_num_test  <- vapply(df_test[,  feat_full, drop = FALSE], is.numeric, logical(1))

bad_train <- names(is_num_train)[!is_num_train]
bad_test  <- names(is_num_test)[!is_num_test]

cat("\n=== NON-NUMERIC FEATURES (should be empty) ===\n")
cat("train:", if (length(bad_train)) paste(bad_train, collapse = ", ") else "NONE", "\n")
cat("test :", if (length(bad_test))  paste(bad_test,  collapse = ", ") else "NONE", "\n")

# ----------------------------
# 3) Show engineered features actually created/skipped
# ----------------------------
cat("\n=== ENGINEERED FEATURES ===\n")
cat("Created:", if (length(eng_created)) paste(eng_created, collapse = ", ") else "None", "\n")
cat("Skipped:", if (length(eng_skipped)) paste(eng_skipped, collapse = ", ") else "None", "\n")

# ----------------------------
# 4) Produce a copy/paste table (column name + human-ish label)
# ----------------------------
strip_suffix <- function(x) sub("(_raw_rank|_rank|_raw)$", "", x, perl = TRUE)

pretty_label <- function(x) {
  b <- strip_suffix(x)
  b <- gsub("_", " ", b)
  # common clinical abbreviations
  b <- gsub("\\bspo2\\b", "SpO2", b, ignore.case = TRUE)
  b <- gsub("\\bsbp\\b",  "SBP",  b, ignore.case = TRUE)
  b <- gsub("\\bdbp\\b",  "DBP",  b, ignore.case = TRUE)
  b <- gsub("\\bmap\\b",  "MAP",  b, ignore.case = TRUE)
  b <- gsub("\\brr\\b",   "RR",   b, ignore.case = TRUE)
  b <- gsub("\\bhr\\b",   "HR",   b, ignore.case = TRUE)
  b <- gsub("\\bnlr\\b",  "NLR",  b, ignore.case = TRUE)
  b <- gsub("\\bpni\\b",  "PNI",  b, ignore.case = TRUE)
  b <- trimws(b)
  b
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

feat_tbl <- tibble(
  feature = sort(unique(feat_full)),
  base    = strip_suffix(feature),
  panel   = ifelse(feature %in% feat_triage, "Triage", "Augmented"),
  label   = vapply(feature, pretty_label, character(1)),
  engineered = feature %in% eng_created
) %>%
  dplyr::arrange(dplyr::desc(panel == "Triage"), feature)

cat("\n=== FEATURE TABLE (copy/paste into Methods) ===\n")
print(feat_tbl, n = Inf, width = Inf)

# ----------------------------
# 5) LaTeX-friendly strings (column names)
# ----------------------------
to_tex_inline <- function(x) paste0("\\texttt{", x, "}", collapse = ", ")

cat("\n=== LaTeX (TRIAGE, column names) ===\n")
cat(to_tex_inline(sort(feat_triage)), "\n")

cat("\n=== LaTeX (AUGMENTED-ONLY, column names) ===\n")
cat(to_tex_inline(sort(full_only)), "\n")


