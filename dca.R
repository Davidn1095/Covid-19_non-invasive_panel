# ============================================================
# DCA — Complete vs Triage with temporal split
# Rows = algorithms, Cols = Test | External
# Each model curve is truncated where its own net benefit reaches 0
# Treat-all is solid black up to prevalence, Treat-none is dashed at 0
# y-axis is clipped at >= 0
# Colors: Complete #F4CCCC, Triage #D9EAD3
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(recipes)
  library(caret)
  library(pROC)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
as_num <- function(x) suppressWarnings(as.numeric(x))
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}
set.seed(123)

# ---------------- config ----------------
cfg <- list(
  file_path = "biomarkers_acuteCOVID_meta.xlsx",
  sheet     = "meta",
  outcome   = "Hospital_ID",
  pos_label = "Yes",
  neg_label = "No",
  boundary  = as.Date("2020-04-21"),
  pct_test  = 0.20,
  algos_req = c("C4.5","k-NN","SVM","RF","LR")
)

# predictors
feat_full   <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission",
                 "monocyte_abs_number","monocytes_perc","lymphocyte_abs_number",
                 "lymphocytes_perc","neutrophil_abs_number","neutrophils_perc")
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- date helpers ----------------
norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

try_parse_date <- function(x) {
  if (inherits(x,"Date"))   return(as.Date(x))
  if (inherits(x,"POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    d1 <- as.Date(x, origin = "1899-12-30")
    d2 <- as.Date(x, origin = "1970-01-01")
    inrng1 <- sum(d1 >= as.Date("2019-01-01") & d1 <= as.Date("2023-12-31"), na.rm=TRUE)
    inrng2 <- sum(d2 >= as.Date("2019-01-01") & d2 <= as.Date("2023-12-31"), na.rm=TRUE)
    return(if (inrng2 > inrng1) d2 else d1)
  }
  xs <- trimws(as.character(x))
  fmts <- c("%d/%m/%Y","%Y-%m-%d","%m/%d/%Y","%d-%m-%Y","%d.%m.%Y","%Y/%m/%d")
  mats <- lapply(fmts, function(f) suppressWarnings(as.Date(xs, format=f)))
  counts <- vapply(mats, function(d) sum(d >= as.Date("2019-01-01") & d <= as.Date("2023-12-31"), na.rm=TRUE), integer(1))
  if (all(counts==0)) return(suppressWarnings(as.Date(xs)))
  mats[[which.max(counts)]]
}

find_sampling_date <- function(df) {
  nms_raw <- names(df)
  nms     <- norm_nm(nms_raw)
  has <- function(p) grepl(p, nms)
  cand_idx <- which(
    (has("sampling") | has("sample") | has("collection") | has("admission") | has("draw") | has("visit")) & has("date")
  )
  if (!length(cand_idx)) cand_idx <- which(has("samplingdate") | has("date"))
  score <- function(col) {
    d <- try_parse_date(df[[col]])
    sum(d >= as.Date("2019-01-01") & d <= as.Date("2023-12-31"), na.rm=TRUE)
  }
  if (length(cand_idx)) {
    sc <- vapply(cand_idx, score, integer(1))
    if (max(sc, na.rm=TRUE) > 0) return(nms_raw[cand_idx[which.max(sc)]])
  }
  sc_all <- vapply(nms_raw, score, integer(1))
  if (max(sc_all, na.rm=TRUE) > 0) return(nms_raw[which.max(sc_all)])
  stop("A sampling date column was not found")
}

# ---------------- load, filter, temporal split ----------------
load_and_split <- function(path, sheet, outcome, boundary, pct_test=0.2){
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  df <- readxl::read_excel(path, sheet=sheet) |> as.data.frame()
  
  grp_col <- intersect(c("Group","group","GROUP"), names(df))
  if (!length(grp_col)) stop("Column 'Group' not found")
  names(df)[names(df)==grp_col[1]] <- "group"
  
  df <- df[df$group %in% c("CTRL_noCOVID","COVID"), , drop=FALSE]
  df$group <- droplevels(factor(df$group))
  
  for (v in intersect(c("group","Gender","Diagnosis","severity_admission","data_split"), names(df))) {
    if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }
  
  if (!(outcome %in% names(df))) stop(sprintf("Outcome '%s' not found", outcome))
  
  to_YN <- function(x){
    if (is.factor(x)) x <- as.character(x)
    if (is.logical(x)) return(factor(ifelse(x, cfg$pos_label, cfg$neg_label), levels=c(cfg$pos_label,cfg$neg_label)))
    if (is.numeric(x)){
      stopifnot(all(x %in% c(0,1), na.rm=TRUE))
      return(factor(ifelse(x==1, cfg$pos_label, cfg$neg_label), levels=c(cfg$pos_label,cfg$neg_label)))
    }
    if (is.character(x)){
      s <- trimws(tolower(x))
      yes <- c("1","yes","y","true","pos","positive")
      no  <- c("0","no","n","false","neg","negative")
      out <- ifelse(s %in% yes, cfg$pos_label, ifelse(s %in% no, cfg$neg_label, NA_character_))
      if (any(is.na(out))) stop("Outcome could not be coerced to Yes/No")
      return(factor(out, levels=c(cfg$pos_label,cfg$neg_label)))
    }
    stop("Unsupported outcome type")
  }
  df[[outcome]] <- to_YN(df[[outcome]])
  
  dcol <- find_sampling_date(df)
  dts  <- try_parse_date(df[[dcol]])
  if (all(is.na(dts))) stop(sprintf("Dates in '%s' could not be parsed", dcol))
  
  pre_idx  <- which(dts <= boundary)
  post_idx <- which(dts >  boundary)
  
  set.seed(444)
  n_pre <- length(pre_idx)
  n_te  <- max(1L, floor(pct_test*n_pre))
  te_ids <- if (n_pre > 1) sample(pre_idx, n_te) else integer(0)
  tr_ids <- setdiff(pre_idx, te_ids)
  
  sp <- rep(NA_character_, nrow(df))
  sp[tr_ids]   <- "train"
  sp[te_ids]   <- "test"
  sp[post_idx] <- "external"
  df$data_split <- factor(sp, levels=c("train","test","external"))
  
  df
}

# ---------------- modeling ----------------
make_recipe <- function(dat, yvar){
  rec <- recipes::recipe(stats::as.formula(paste(yvar,"~ .")), data=dat)
  ign <- intersect(c("data_split"), names(dat))
  if (length(ign)) rec <- rec |> update_role(all_of(ign), new_role="ignore")
  rec |>
    step_impute_median(all_numeric_predictors()) |>
    step_impute_mode(all_nominal_predictors())  |>
    step_novel(all_nominal_predictors())        |>
    step_other(all_nominal_predictors(), threshold=0.01) |>
    step_dummy(all_nominal_predictors())        |>
    step_zv(all_predictors())                   |>
    step_normalize(all_numeric_predictors())
}

algo_requires <- function(method){
  switch(method,
         "glmnet"    = "glmnet",
         "rf"        = "randomForest",
         "svmRadial" = "kernlab",
         "kknn"      = "kknn",
         "J48"       = "RWeka",
         NULL)
}
method_available <- function(method){
  pkg <- algo_requires(method); if (is.null(pkg)) return(TRUE); requireNamespace(pkg, quietly=TRUE)
}

model_specs <- function(tune_len=5, algos=cfg$algos_req){
  all <- list(
    "LR"   = list(method="glmnet",    tuneLength=tune_len),
    "RF"   = list(method="rf",        tuneLength=tune_len),
    "SVM"  = list(method="svmRadial", tuneLength=tune_len),
    "k-NN" = list(method="kknn",      tuneLength=tune_len),
    "C4.5" = list(method="J48",       tuneLength=tune_len)
  )
  keep <- names(all) %in% algos & vapply(all, \(s) method_available(s$method), logical(1))
  all[keep]
}

train_inner <- function(dat_train, yvar, spec, inner_k=3){
  rec <- make_recipe(dat_train, yvar)
  tr_ctrl <- caret::trainControl(method="cv", number=inner_k, classProbs=TRUE,
                                 summaryFunction=twoClassSummary, savePredictions="final",
                                 allowParallel=TRUE)
  caret::train(rec, data=dat_train, method=spec$method, trControl=tr_ctrl, metric="ROC",
               tuneLength=spec$tuneLength)
}

# Platt scaling applied to both Test and External
safe_logit  <- function(p, eps=1e-6) qlogis(pmin(pmax(as.numeric(p), eps), 1 - eps))
platt_fit   <- function(obs, p_train){ y <- as.integer(obs == cfg$pos_label); s <- safe_logit(p_train); suppressWarnings(glm(y ~ s, family=binomial())) }
platt_apply <- function(fit, p_new){ s <- safe_logit(p_new); as.numeric(plogis(drop(cbind(1, s) %*% coef(fit)))) }

run_holdout <- function(df, yvar, features, algo_name, split_train="train", split_test=c("test","external")){
  present_splits <- intersect(split_test, unique(as.character(df$data_split)))
  if (!length(present_splits)) return(NULL)
  specs <- model_specs(algos=algo_name); if (!length(specs)) return(NULL)
  spec  <- specs[[algo_name]]
  
  base_cols <- unique(c(features, yvar))   # do not include 'group' as predictor
  base_cols <- intersect(base_cols, names(df))
  
  dtrain <- df[df$data_split==split_train, base_cols, drop=FALSE]
  if (!nrow(dtrain)) return(NULL)
  
  fit  <- train_inner(dtrain, yvar, spec, inner_k=3)
  p_tr <- as.numeric(predict(fit, newdata=dtrain, type="prob")[, cfg$pos_label])
  
  # MCC threshold on train, stored but not used for DCA
  neg  <- setdiff(levels(dtrain[[yvar]]), cfg$pos_label)[1] %||% cfg$neg_label
  best_t <- 0.5; best_m <- -Inf
  for (t in seq(0.01,0.99,by=0.001)){
    pred <- factor(ifelse(p_tr >= t, cfg$pos_label, neg), levels=c(neg,cfg$pos_label))
    tab <- table(factor(dtrain[[yvar]], levels=c(neg,cfg$pos_label)), pred)
    TP <- as_num(tab[cfg$pos_label,cfg$pos_label] %||% 0); TN <- as_num(tab[neg,neg] %||% 0)
    FP <- as_num(tab[neg,cfg$pos_label] %||% 0); FN <- as_num(tab[cfg$pos_label,neg] %||% 0)
    den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); m <- if (den==0) NA_real_ else (TP*TN - FP*FN)/den
    if (is.finite(m) && m > best_m) { best_m <- m; best_t <- t }
  }
  thr <- best_t
  
  platt <- tryCatch(platt_fit(dtrain[[yvar]], p_tr), error=function(e) NULL)
  
  out <- list()
  for (sp in present_splits){
    dtest <- df[df$data_split==sp, base_cols, drop=FALSE]
    if (!nrow(dtest)) next
    p_raw <- as.numeric(predict(fit, newdata=dtest, type="prob")[, cfg$pos_label])
    p_use <- if (!is.null(platt)) platt_apply(platt, p_raw) else p_raw
    out[[sp]] <- list(obs=factor(dtest[[yvar]]), p=p_use, p_raw=p_raw, thr=thr)
  }
  out
}

# ---------------- decision curve helpers ----------------
decision_curve_table <- function(obs, p, thresholds=seq(0.01,0.99,by=0.01)){
  y_raw <- as.integer(obs == cfg$pos_label)
  keep  <- is.finite(as.numeric(p)) & !is.na(y_raw)
  y <- y_raw[keep]; p <- as.numeric(p)[keep]
  if (!length(y)) return(tibble(threshold=thresholds, Net_Benefit=NA_real_, NB_TreatAll=NA_real_, NB_TreatNone=0, prevalence=NA_real_))
  N <- length(y); prev <- mean(y)
  purrr::map_dfr(thresholds, function(pt){
    pred <- as.integer(p >= pt)
    TP <- sum(pred==1 & y==1)
    FP <- sum(pred==1 & y==0)
    NB <- TP/N - FP/N * (pt/(1-pt))
    NB_all <- prev - (1 - prev) * (pt/(1-pt))
    tibble(threshold=pt, Net_Benefit=NB, NB_TreatAll=NB_all, NB_TreatNone=0, prevalence=prev)
  })
}

nb_zero_threshold <- function(dca_df){
  d <- dca_df[order(dca_df$threshold),]
  nb <- d$Net_Benefit; th <- d$threshold
  bad <- !is.finite(nb) | !is.finite(th); nb <- nb[!bad]; th <- th[!bad]
  if (!length(nb)) return(NA_real_)
  i_neg <- which(nb < 0)
  if (!length(i_neg)) return(max(th, na.rm=TRUE))
  i <- i_neg[1]; if (i==1) return(th[1])
  t0 <- th[i-1]; n0 <- nb[i-1]; t1 <- th[i]; n1 <- nb[i]
  if (!is.finite(n0) || !is.finite(n1) || t1==t0) return(th[i])
  t0 - n0 * (t1 - t0) / (n1 - n0)
}

# ---------------- RUN: train, score, DCA data ----------------
df <- load_and_split(cfg$file_path, cfg$sheet, cfg$outcome, cfg$boundary, cfg$pct_test)

algos_avail <- names(model_specs(algos=cfg$algos_req))
if (!length(algos_avail)) stop("No requested algorithms are available, please install RWeka, kknn, kernlab, randomForest, glmnet")

hold_full <- lapply(algos_avail, function(a) run_holdout(df, cfg$outcome, feat_full,   a))
names(hold_full) <- algos_avail
hold_tri  <- lapply(algos_avail, function(a) run_holdout(df, cfg$outcome, feat_triage, a))
names(hold_tri)  <- algos_avail

# build per-algo DCA data, each model truncated at its own zero crossing
mk_dca_with_t0 <- function(obj, set_name, split_name){
  if (is.null(obj)) return(list(d=NULL, t0=NA_real_, prev=NA_real_))
  d <- decision_curve_table(obj$obs, obj$p)
  t0 <- nb_zero_threshold(d)
  prev <- suppressWarnings(unique(d$prevalence)[1])
  d <- d %>% filter(threshold <= t0)
  if (!nrow(d)) return(list(d=NULL, t0=t0, prev=prev))
  d$set <- set_name
  d$split <- split_name
  list(d=d, t0=t0, prev=prev)
}

build_algo_df <- function(a){
  o_te_full <- hold_full[[a]][["test"]]
  o_ex_full <- hold_full[[a]][["external"]]
  o_te_tri  <- hold_tri [[a]][["test"]]
  o_ex_tri  <- hold_tri [[a]][["external"]]
  
  te_full <- mk_dca_with_t0(o_te_full, "Complete feature set", "Test")
  te_tri  <- mk_dca_with_t0(o_te_tri,  "Triage feature set",   "Test")
  ex_full <- mk_dca_with_t0(o_ex_full, "Complete feature set", "External")
  ex_tri  <- mk_dca_with_t0(o_ex_tri,  "Triage feature set",   "External")
  
  tmax_te <- safe_max(c(te_full$t0, te_tri$t0))
  prev_te <- safe_max(c(te_full$prev, te_tri$prev))
  tmax_ex <- safe_max(c(ex_full$t0, ex_tri$t0))
  prev_ex <- safe_max(c(ex_full$prev, ex_tri$prev))
  
  d_all <- bind_rows(te_full$d, te_tri$d, ex_full$d, ex_tri$d)
  
  mk_treat <- function(ref_obj, split_label, tmax, prev){
    if (is.null(ref_obj) || !is.finite(tmax) || !is.finite(prev)) return(NULL)
    base <- decision_curve_table(ref_obj$obs, ref_obj$p)
    t_lim <- min(tmax, prev, na.rm=TRUE)
    base %>% filter(threshold <= t_lim) %>%
      transmute(split = split_label, threshold, NB_TreatAll, NB_TreatNone)
  }
  
  treat_te <- mk_treat(if (!is.null(o_te_tri)) o_te_tri else o_te_full, "Test",     tmax_te, prev_te)
  treat_ex <- mk_treat(if (!is.null(o_ex_tri)) o_ex_tri else o_ex_full, "External", tmax_ex, prev_ex)
  
  list(
    dca   = d_all %>% mutate(algo = a),
    treat = bind_rows(treat_te, treat_ex) %>% mutate(algo = a)
  )
}

combo    <- lapply(algos_avail, build_algo_df)
dca_df   <- bind_rows(lapply(combo, `[[`, "dca"))
treat_df <- bind_rows(lapply(combo, `[[`, "treat"))

if (is.null(dca_df) || !nrow(dca_df)) stop("No decision-curve data were produced")

# factor order
alg_order <- intersect(cfg$algos_req, unique(dca_df$algo))
dca_df$algo   <- factor(dca_df$algo,   levels = alg_order)
treat_df$algo <- factor(treat_df$algo, levels = alg_order)
dca_df$split   <- factor(dca_df$split,   levels=c("Test","External"))
treat_df$split <- factor(treat_df$split, levels=c("Test","External"))

# y cap at positive values only
y_cap <- max(c(dca_df$Net_Benefit, treat_df$NB_TreatAll, 0), na.rm=TRUE)
if (!is.finite(y_cap) || y_cap <= 0) y_cap <- 0.05

# ---------------- plot ----------------
pastel_red   <- "#F4CCCC"
pastel_green <- "#D9EAD3"

theme_pub <- function(){
  theme_minimal(base_size=15) +
    theme(
      panel.grid.major = element_line(color="grey85", linewidth=0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color="grey30"),
      legend.position = "right",
      legend.title = element_text(face="bold")
    )
}

p <- ggplot() +
  geom_line(data=treat_df, aes(threshold, NB_TreatAll, linetype="Treat-all"), color="black") +
  geom_line(data=treat_df, aes(threshold, NB_TreatNone, linetype="Treat-none"), color="grey40") +
  geom_line(data=dca_df, aes(threshold, Net_Benefit, color=set), linewidth=0.9) +
  scale_color_manual(values=c("Complete feature set"=pastel_red, "Triage feature set"=pastel_green), name="Model") +
  scale_linetype_manual(values=c("Treat-all"="solid","Treat-none"="dashed"), name=NULL) +
  labs(x="Threshold probability", y="Net benefit") +
  facet_grid(algo ~ split) +
  coord_cartesian(ylim=c(0, y_cap)) +
  theme_pub()

print(p)

# Optional save
ggsave("20250922_DCA.png", p, width=12, height=10, dpi=300, bg="white")

