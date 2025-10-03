# ======================= ROC + PR with temporal boundary (print only, pastel red/green, single caption) =======================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(recipes); library(caret); library(glmnet)
  library(pROC); library(ggplot2); library(lubridate)
  suppressWarnings(requireNamespace("patchwork", quietly = TRUE))
  suppressWarnings(requireNamespace("gridExtra", quietly = TRUE))
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---------------- config ----------------
cfg <- list(
  file_path   = "biomarkers_acuteCOVID_meta.xlsx",
  sheet       = "meta",
  admit_col   = "Sampling _date",
  boundary    = as.Date("2020-04-21"),
  outcome     = "Hospital_ID",
  pos_label   = "Yes",
  neg_label   = "No",
  feature_set = "triage",   # "triage" or "full"
  seed        = 444
)

# ---------------- feature sets ----------------
feat_full  <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission",
                "monocyte_abs_number","monocytes_perc",
                "lymphocyte_abs_number","lymphocytes_perc",
                "neutrophil_abs_number","neutrophils_perc")
feat_triage <- c("Diagnosis","severity_admission","Age","Gender","SpO2_admission")

# ---------------- helpers ----------------
coerce_outcome_YN <- function(y, pos = cfg$pos_label, neg = cfg$neg_label) {
  if (is.factor(y)) y <- as.character(y)
  if (is.logical(y)) y <- ifelse(y, pos, neg)
  if (is.numeric(y)) {
    if (!all(y %in% c(0,1), na.rm = TRUE)) stop("Outcome numeric but not 0 or 1")
    return(factor(ifelse(y == 1, pos, neg), levels = c(pos, neg)))
  }
  if (is.character(y)) {
    s <- trimws(tolower(y))
    map_yes <- c("1","yes","y","true","pos","positive","si","sí")
    map_no  <- c("0","no","n","false","neg","negative")
    out <- ifelse(s %in% map_yes, pos, ifelse(s %in% map_no, neg, NA_character_))
    if (any(is.na(out))) stop("Outcome not coercible to Yes or No")
    return(factor(out, levels = c(pos, neg)))
  }
  stop("Outcome type unsupported")
}

parse_sampling_date <- function(z) {
  if (inherits(z, "Date"))   return(z)
  if (inherits(z, "POSIXt")) return(as.Date(z))
  n <- length(z)
  out <- rep(as.Date(NA), n)
  zn <- suppressWarnings(as.numeric(z))
  if (any(!is.na(zn))) {
    out_num <- as.Date(zn, origin = "1899-12-30")
    out[!is.na(zn)] <- out_num[!is.na(zn)]
  }
  zc <- as.character(z)
  d <- suppressWarnings(lubridate::dmy(zc)); out[is.na(out) & !is.na(d)] <- d[is.na(out) & !is.na(d)]
  d <- suppressWarnings(lubridate::ymd(zc)); out[is.na(out) & !is.na(d)] <- d[is.na(out) & !is.na(d)]
  d <- suppressWarnings(lubridate::mdy(zc)); out[is.na(out) & !is.na(d)] <- d[is.na(out) & !is.na(d)]
  out
}

enforce_temporal_split <- function(df, admit_col, boundary) {
  if (!(admit_col %in% names(df))) stop(sprintf("Admission column '%s' not found", admit_col))
  dts <- parse_sampling_date(df[[admit_col]])
  if (all(is.na(dts))) stop("Admission dates could not be parsed")
  keep <- !is.na(dts)
  if (any(!keep)) { df <- df[keep, , drop = FALSE]; dts <- dts[keep] }
  df$.temporal_split <- factor(ifelse(dts <= boundary, "pre", "post"), levels = c("pre","post"))
  df
}

carve_internal_test <- function(df, yvar, pre_label = "pre", test_frac = 0.2, seed = cfg$seed) {
  set.seed(seed)
  df$data_split <- NA_character_
  pre_idx  <- which(df$.temporal_split == pre_label)
  post_idx <- which(df$.temporal_split != pre_label)
  if (!length(pre_idx)) stop("No pre-boundary rows to create train and test")
  y_pre <- droplevels(df[[yvar]][pre_idx])
  take_test <- unlist(tapply(pre_idx, y_pre, function(ix) {
    n <- length(ix); size <- max(1, floor(test_frac * n)); sample(ix, size = size)
  }))
  df$data_split[pre_idx]   <- "train"
  df$data_split[take_test] <- "test"
  df$data_split[post_idx]  <- "external"
  df$data_split <- factor(df$data_split, levels = c("train","test","external"))
  df
}

make_recipe <- function(dat, yvar) {
  rec <- recipes::recipe(stats::as.formula(paste(yvar, "~ .")), data = dat) |>
    recipes::update_role(dplyr::any_of(c("group",".temporal_split","data_split")), new_role = "ignore") |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_impute_mode(recipes::all_nominal_predictors()) |>
    recipes::step_novel(recipes::all_nominal_predictors()) |>
    recipes::step_other(recipes::all_nominal_predictors(), threshold = 0.01) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_zv(recipes::all_predictors())
  if (exists("step_YeoJohnson", where = asNamespace("recipes"), inherits = FALSE)) {
    rec <- rec |> recipes::step_YeoJohnson(recipes::all_numeric_predictors())
  } else if (exists("step_yeojohnson", where = asNamespace("recipes"), inherits = FALSE)) {
    rec <- rec |> recipes::step_yeojohnson(recipes::all_numeric_predictors())
  }
  rec |> recipes::step_normalize(recipes::all_numeric_predictors())
}

train_lr_glmnet <- function(dat, yvar, seed = cfg$seed) {
  set.seed(seed)
  tr_ctrl <- caret::trainControl(
    method = "cv", number = 5,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", allowParallel = TRUE
  )
  grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 0, length.out = 10))
  caret::train(
    make_recipe(dat, yvar),
    data = dat,
    method = "glmnet",
    trControl = tr_ctrl,
    metric = "ROC",
    tuneGrid = grid
  )
}

# ===== curves without yardstick dependencies =====
roc_curve_pROC <- function(truth_factor, probs) {
  y <- as.integer(truth_factor == cfg$pos_label)
  r <- pROC::roc(y, probs, quiet = TRUE, direction = "<")
  tibble::tibble(FPR = 1 - as.numeric(r$specificities),
                 TPR = as.numeric(r$sensitivities))
}

pr_curve_basic <- function(truth_factor, probs) {
  y <- as.integer(truth_factor == cfg$pos_label)
  ord <- order(probs, decreasing = TRUE, na.last = NA)
  y <- y[ord]
  tp <- cumsum(y)
  fp <- cumsum(1 - y)
  P  <- sum(y)
  if (P == 0) return(tibble::tibble(Recall = numeric(0), Precision = numeric(0)))
  tibble::tibble(Recall = tp / P, Precision = tp / (tp + fp))
}

average_precision_basic <- function(truth_factor, probs) {
  y <- as.integer(truth_factor == cfg$pos_label)
  ord <- order(probs, decreasing = TRUE, na.last = NA)
  y <- y[ord]
  tp <- cumsum(y)
  fp <- cumsum(1 - y)
  P  <- sum(y)
  if (P == 0) return(NA_real_)
  recall    <- tp / P
  precision <- tp / (tp + fp)
  precision_env <- rev(cummax(rev(precision)))
  recall_ext    <- c(0, recall)
  precision_ext <- c(precision_env[1], precision_env)
  sum((recall_ext[-1] - recall_ext[-length(recall_ext)]) * precision_ext[-1], na.rm = TRUE)
}

# ---------------- run ----------------
set.seed(cfg$seed)

# load and harmonize
raw <- readxl::read_excel(cfg$file_path, sheet = cfg$sheet) |> as.data.frame()
grp_col <- intersect(c("Group","group","GROUP"), names(raw))
if (!length(grp_col)) stop("Column 'Group' not found")
names(raw)[names(raw) == grp_col[1]] <- "group"
raw <- raw[raw$group %in% c("CTRL_noCOVID","COVID"), , drop = FALSE]
if (!(cfg$outcome %in% names(raw))) stop(sprintf("Outcome '%s' not found", cfg$outcome))
raw[[cfg$outcome]] <- coerce_outcome_YN(raw[[cfg$outcome]])

# temporal split at 21-Apr-2020
raw_ts <- enforce_temporal_split(raw, admit_col = cfg$admit_col, boundary = cfg$boundary)
cat("[Counts] pre: ", sum(raw_ts$.temporal_split == "pre", na.rm = TRUE),
    " post: ", sum(raw_ts$.temporal_split == "post", na.rm = TRUE), "\n", sep = "")

# keep predictors that truly exist
all_feats <- unique(c(feat_full, feat_triage))
present_feats <- intersect(all_feats, names(raw_ts))
if (!length(present_feats)) stop("None of the expected predictors are present")
dat0 <- raw_ts[, unique(c(present_feats, cfg$outcome, "group", ".temporal_split")), drop = FALSE]

# carve internal test from pre-boundary
dat1 <- carve_internal_test(dat0, yvar = cfg$outcome, pre_label = "pre", test_frac = 0.2, seed = cfg$seed)

# train on pre-boundary train
features <- if (tolower(cfg$feature_set) == "triage") feat_triage else feat_full
features <- intersect(features, names(dat1))
train_set <- dplyr::filter(dat1, data_split == "train") |>
  dplyr::select(dplyr::all_of(c(features, cfg$outcome, "group", ".temporal_split", "data_split")))
fit <- train_lr_glmnet(train_set, yvar = cfg$outcome)

# predict helper
predict_split <- function(d, split_name) {
  d2 <- dplyr::filter(d, data_split == split_name)
  if (!nrow(d2)) return(NULL)
  p <- as.numeric(predict(fit, newdata = d2, type = "prob")[, cfg$pos_label])
  tibble::tibble(truth = factor(d2[[cfg$outcome]], levels = c(cfg$pos_label, cfg$neg_label)),
                 .pred_Yes = p)
}
pred_test <- predict_split(dat1, "test")
pred_ext  <- predict_split(dat1, "external")

# curves and areas with requested labels
roc_df <- pr_df <- NULL
lbls <- character(0)

add_split <- function(pred_df, split_lab) {
  if (is.null(pred_df) || !nrow(pred_df)) return(invisible(NULL))
  r <- roc_curve_pROC(pred_df$truth, pred_df$.pred_Yes) %>% mutate(Split = split_lab)
  p <- pr_curve_basic(pred_df$truth, pred_df$.pred_Yes) %>% mutate(Split = split_lab)
  y_bin <- as.integer(pred_df$truth == cfg$pos_label)
  aucv <- as.numeric(pROC::auc(y_bin, pred_df$.pred_Yes, quiet = TRUE))
  ap   <- average_precision_basic(pred_df$truth, pred_df$.pred_Yes)
  assign("roc_df", dplyr::bind_rows(get("roc_df", inherits = TRUE), r), inherits = TRUE)
  assign("pr_df",  dplyr::bind_rows(get("pr_df",  inherits = TRUE), p), inherits = TRUE)
  assign("lbls", c(get("lbls", inherits = TRUE), sprintf("%s  AUC %.4f, AP %.4f", split_lab, aucv, ap)), inherits = TRUE)
  invisible(NULL)
}
add_split(pred_test, "Internal test")
add_split(pred_ext,  "External validation")

if (is.null(roc_df) || is.null(pr_df) || !nrow(roc_df) || !nrow(pr_df)) stop("Could not compute curves")

# enforce split order and palette
roc_df$Split <- factor(roc_df$Split, levels = c("Internal test","External validation"))
pr_df$Split  <- factor(pr_df$Split,  levels = c("Internal test","External validation"))
pastel_red   <- "#FF9999"
pastel_green <- "#99CC99"
cap <- paste(lbls, collapse = "  |  ")

# plots, no per-panel title or caption
g_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Split, linetype = Split)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_blank()) +
  scale_color_manual(values = c("Internal test" = pastel_red, "External validation" = pastel_green)) +
  scale_linetype_manual(values = c("Internal test" = "solid", "External validation" = "dashed")) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))

g_pr <- ggplot(pr_df, aes(x = Recall, y = Precision, color = Split, linetype = Split)) +
  geom_line(linewidth = 0.9) +
  labs(x = "Recall", y = "Precision", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_blank()) +
  scale_color_manual(values = c("Internal test" = pastel_red, "External validation" = pastel_green)) +
  scale_linetype_manual(values = c("Internal test" = "solid", "External validation" = "dashed")) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))

# print both panels with a single shared caption, save nothing
if (requireNamespace("patchwork", quietly = TRUE)) {
  p <- g_roc + g_pr + patchwork::plot_layout(ncol = 2)
  p <- p + patchwork::plot_annotation(
    caption = cap,
    theme = theme(plot.caption = element_text(hjust = 0.5))
  )
  print(p)
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  gridExtra::grid.arrange(
    g_roc, g_pr, ncol = 2,
    bottom = grid::textGrob(cap, x = 0.5, hjust = 0.5, gp = grid::gpar(cex = 0.9))
  )
} else {
  # fallback without caption on figure
  print(g_roc); print(g_pr)
  cat(cap, "\n")  # caption printed once
}
