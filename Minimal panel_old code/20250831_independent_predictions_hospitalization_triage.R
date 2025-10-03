# ============================================================
# Triage model: bedside features only
# ============================================================

# define triage feature set
triage_features <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")
stopifnot(all(triage_features %in% colnames(df)))

# keep a copy of the full set, then switch to triage
feature_vars_full <- feature_vars
feature_vars      <- triage_features

# run exactly the same training and evaluation
hosp_tbl_triage <- run_task(
  task_name = "Hospitalization_Triage_bedside",
  yvar      = "Hosp_Bin",
  type      = "binary",
  pos_label = "Yes",
  seed      = 444
)
print_task_table(hosp_tbl_triage, "Hospitalization triage, bedside only")

# optional, show full vs triage side by side
if (exists("hosp_tbl") && nrow(hosp_tbl)) {
  comp_tbl <- dplyr::bind_rows(hosp_tbl, hosp_tbl_triage)
  print_task_table(comp_tbl, "Full vs triage comparison")
}

# restore the full feature set for any subsequent steps
feature_vars <- feature_vars_full
