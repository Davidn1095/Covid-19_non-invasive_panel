suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# ---------------- config ----------------
cfg <- list(
  file_path = "biomarkers_acuteCOVID_meta.xlsx",
  sheet     = "meta",
  outcome   = "Hospital_ID",
  pos_label = "Yes",
  neg_label = "No",
  boundary  = as.Date("2020-04-21"),
  pct_test  = 0.20,
  seed      = 444
)

# ---------------- helpers ----------------
`%||%` <- function(a,b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))
norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

try_parse_date <- function(x) {
  if (inherits(x, "Date"))   return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    # Excel vs Unix origin, pick the one with plausible 2019–2023 coverage
    d1 <- as.Date(x, origin = "1899-12-30")
    d2 <- as.Date(x, origin = "1970-01-01")
    r <- as.Date("2019-01-01"); R <- as.Date("2023-12-31")
    if (sum(d2>=r & d2<=R, na.rm=TRUE) > sum(d1>=r & d1<=R, na.rm=TRUE)) d2 else d1
  } else {
    xs <- trimws(as.character(x))
    fmts <- c("%d/%m/%Y","%Y-%m-%d","%m/%d/%Y","%d-%m-%Y","%d.%m.%Y","%Y/%m/%d")
    mats <- lapply(fmts, function(f) suppressWarnings(as.Date(xs, format=f)))
    r <- as.Date("2019-01-01"); R <- as.Date("2023-12-31")
    counts <- vapply(mats, function(d) sum(d>=r & d<=R, na.rm=TRUE), integer(1))
    if (all(counts==0)) suppressWarnings(as.Date(xs)) else mats[[which.max(counts)]]
  }
}

find_sampling_date <- function(df) {
  nms_raw <- names(df); nms <- norm_nm(nms_raw); has <- function(p) grepl(p, nms)
  cand <- which((has("sampling")|has("sample")|has("collection")|has("admission")|has("draw")|has("visit")) & has("date"))
  if (!length(cand)) cand <- which(has("samplingdate")|has("date"))
  score <- function(col){
    d <- try_parse_date(df[[col]])
    sum(d >= as.Date("2019-01-01") & d <= as.Date("2023-12-31"), na.rm=TRUE)
  }
  if (length(cand)) {
    sc <- vapply(cand, function(i) score(nms_raw[i]), integer(1))
    if (max(sc, na.rm=TRUE) > 0) return(nms_raw[cand[which.max(sc)]])
  }
  sc_all <- vapply(nms_raw, score, integer(1))
  if (max(sc_all, na.rm=TRUE) > 0) return(nms_raw[which.max(sc_all)])
  stop("A sampling date column was not found")
}

to_YN <- function(x, pos = cfg$pos_label, neg = cfg$neg_label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(factor(ifelse(x, pos, neg), levels = c(pos, neg)))
  if (is.numeric(x)) {
    stopifnot(all(x %in% c(0,1), na.rm=TRUE))
    return(factor(ifelse(x==1, pos, neg), levels = c(pos, neg)))
  }
  if (is.character(x)) {
    s <- trimws(tolower(x))
    yes <- c("1","yes","y","true","pos","positive")
    no  <- c("0","no","n","false","neg","negative")
    out <- ifelse(s %in% yes, pos, ifelse(s %in% no, neg, NA_character_))
    if (any(is.na(out))) stop("Outcome could not be coerced to Yes/No")
    return(factor(out, levels = c(pos, neg)))
  }
  stop("Unsupported outcome type for outcome conversion")
}

# ---------------- main maker ----------------
produce_analysis_df <- function(cfg) {
  if (!file.exists(cfg$file_path)) stop(sprintf("File not found: %s", cfg$file_path))
  df <- readxl::read_excel(cfg$file_path, sheet = cfg$sheet) |> as.data.frame()
  
  # normalize Group column to 'group' and enforce cohort filter
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df)==grp_col[1]] <- "group"
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
  df$group <- droplevels(factor(df$group))
  
  # coerce key categoricals
  for (v in intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  # coerce outcome to Yes/No with positive level first
  if (!(cfg$outcome %in% names(df))) stop(sprintf("Outcome '%s' not found", cfg$outcome))
  df[[cfg$outcome]] <- to_YN(df[[cfg$outcome]], cfg$pos_label, cfg$neg_label)
  
  # detect and parse sampling date
  dcol <- find_sampling_date(df)
  samp_dt <- try_parse_date(df[[dcol]])
  if (all(is.na(samp_dt))) stop(sprintf("Dates in '%s' could not be parsed", dcol))
  df$sampling_date <- samp_dt
  
  # temporal split: pre-boundary -> train+test, post-boundary -> external
  pre_idx  <- which(df$sampling_date <= cfg$boundary)
  post_idx <- which(df$sampling_date >  cfg$boundary)
  
  set.seed(cfg$seed)
  n_pre <- length(pre_idx)
  n_te  <- max(1L, floor(cfg$pct_test * n_pre))
  te_ids <- if (n_pre > 1) sample(pre_idx, n_te) else integer(0)
  tr_ids <- setdiff(pre_idx, te_ids)
  
  sp <- rep(NA_character_, nrow(df))
  sp[tr_ids]   <- "train"
  sp[te_ids]   <- "test"
  sp[post_idx] <- "external"
  df$data_split <- factor(sp, levels = c("train","test","external"))
  
  # ensure numeric where expected, keep only relevant columns
  triage_feats  <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")
  full_feats    <- c(triage_feats,
                     "monocyte_abs_number","monocytes_perc",
                     "lymphocyte_abs_number","lymphocytes_perc",
                     "neutrophil_abs_number","neutrophils_perc")
  
  keep_feats <- intersect(unique(c(triage_feats, full_feats)), names(df))
  
  # numeric coercion for lab-vitals, leave categoricals untouched
  numeric_guess <- setdiff(keep_feats, c("Diagnosis","severity_admission","Gender"))
  for (v in intersect(numeric_guess, names(df))) {
    if (!is.numeric(df[[v]])) df[[v]] <- as_num(df[[v]])
  }
  
  df$.rid <- seq_len(nrow(df))
  
  df |>
    dplyr::select(
      .rid, sampling_date, group, all_of(cfg$outcome),
      data_split, all_of(keep_feats)
    )
}

# ---------------- run ----------------
df <- produce_analysis_df(cfg)

# Optional persistence for reuse across scripts
# saveRDS(df, "analysis_df.rds")
# readr::write_csv(df, "analysis_df.csv")

# df is now ready for subsequent analysis
