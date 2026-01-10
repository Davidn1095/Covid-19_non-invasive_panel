#!/usr/bin/env Rscript

set_env_if_exists <- function(var, path) {
  if (file.exists(path)) {
    Sys.setenv(var = path)
  }
}

compute_checksum <- function(lines) {
  text <- paste(lines, collapse = "\n")
  sum(utf8ToInt(text)) %% 1e9
}

print_key_lines <- function(lines, patterns, max_matches = 10) {
  regex <- paste(patterns, collapse = "|")
  matches <- grep(regex, lines, value = TRUE)
  if (length(matches) == 0) {
    cat("  (none)\n")
    return(invisible(NULL))
  }
  matches <- head(matches, max_matches)
  for (line in matches) {
    cat("  ", line, "\n", sep = "")
  }
  invisible(NULL)
}

scripts <- c(
  "20260107_performance_cascade.R",
  "20260107_questionnaire.R",
  "20260109_parsimony.R",
  "20260107_feature_importance_heatmap.R"
)

set_env_if_exists("DEV_CSV", "20260106_development_raw_rowcomplete60.csv")
set_env_if_exists("EXT_CSV", "20260106_external_raw_rowcomplete60.csv")

patterns <- c(
  "DEV",
  "OOF",
  "EXT",
  "AUC",
  "MCC",
  "Deferral",
  "tau",
  "thr",
  "Cascade",
  "Non-invasive",
  "Laboratory augmented"
)

failures <- 0L

for (script in scripts) {
  cat("=== ", script, " ===\n", sep = "")
  output <- system2(
    "Rscript",
    c("--vanilla", script),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  line_count <- length(output)
  char_count <- nchar(paste(output, collapse = "\n"), type = "chars")
  checksum <- compute_checksum(output)

  cat("Status: ", status, "\n", sep = "")
  cat("Lines: ", line_count, "\n", sep = "")
  cat("Chars: ", char_count, "\n", sep = "")
  cat("Checksum: ", checksum, "\n", sep = "")
  cat("Key lines:\n")
  print_key_lines(output, patterns, max_matches = 10)
  cat("\n")

  if (status != 0) {
    failures <- failures + 1L
  }
}

if (failures == 0) {
  cat("OVERALL: OK (0 failed)\n")
  quit(status = 0)
}

cat("OVERALL: FAIL (", failures, " failed)\n", sep = "")
quit(status = 1)
