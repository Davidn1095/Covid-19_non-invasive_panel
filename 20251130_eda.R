## =====================================================================
## EDA plots: outcome prevalence, missingness, correlation and PCA
## PCA is run separately in each dataset
## =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
})

## ---------------------------------------------------------------
## 0. Data preparation
## ---------------------------------------------------------------
if (!exists("dev_raw") || !exists("val_raw")) {
  dev_raw <- readr::read_csv("mimic_final_22_features_60pct.csv", show_col_types = FALSE)
  val_raw <- readr::read_csv("mcmed_final_22_features_60pct.csv", show_col_types = FALSE)
}

prep_df <- function(df) {
  df %>%
    dplyr::mutate(outcome = factor(outcome, levels = c("No", "Yes"))) %>%
    dplyr::select(outcome, dplyr::everything())
}

df_dev <- prep_df(dev_raw) %>% dplyr::mutate(dataset = "dev")
df_val <- prep_df(val_raw) %>% dplyr::mutate(dataset = "val")
df_all <- dplyr::bind_rows(df_dev, df_val)

## ---------------------------------------------------------------
## 1. Bar plot of outcome prevalence by dataset
## ---------------------------------------------------------------
df_all %>%
  dplyr::count(dataset, outcome) %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  ggplot2::ggplot(ggplot2::aes(x = dataset, y = prop, fill = outcome)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    x = "Dataset",
    y = "Outcome prevalence",
    fill = "Outcome"
  ) +
  ggplot2::theme_bw()

## ---------------------------------------------------------------
## 2. Bar plot of percent missing per variable, dev vs val
## ---------------------------------------------------------------
missing_summary <- function(df, dataset_name) {
  df %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(is.na(.)) * 100)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "pct_missing") %>%
    dplyr::mutate(dataset = dataset_name)
}

miss_dev <- missing_summary(df_dev, "dev")
miss_val <- missing_summary(df_val, "val")
miss_all <- dplyr::bind_rows(miss_dev, miss_val)

top_miss <- miss_all %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(max_missing = max(pct_missing)) %>%
  dplyr::arrange(dplyr::desc(max_missing)) %>%
  dplyr::slice(1:20)

miss_all %>%
  dplyr::semi_join(top_miss, by = "variable") %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = stats::reorder(variable, pct_missing),
      y = pct_missing,
      fill = dataset
    )
  ) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = "Variable",
    y = "Percent missing",
    fill = "Dataset"
  ) +
  ggplot2::theme_bw()

## ---------------------------------------------------------------
## 3. Correlation heatmap of numeric features in dev
## ---------------------------------------------------------------
num_dev_mat <- df_dev %>%
  dplyr::select(where(is.numeric), -outcome, -dataset) %>%
  dplyr::select(where(~ !all(is.na(.))))

if (ncol(num_dev_mat) > 1) {
  cor_mat <- stats::cor(num_dev_mat, use = "pairwise.complete.obs")
  
  corrplot::corrplot(
    cor_mat,
    method = "color",
    type   = "lower",
    tl.cex = 0.7
  )
}

## ---------------------------------------------------------------
## 4. PCA scatter plots, one per dataset
##    color = outcome, shape = dataset
## ---------------------------------------------------------------
pca_plot_dataset <- function(df, dataset_name) {
  num_mat <- df %>%
    dplyr::select(where(is.numeric), -outcome, -dataset) %>%
    dplyr::select(where(~ !all(is.na(.))))
  
  cc_idx  <- stats::complete.cases(num_mat)
  num_cc  <- num_mat[cc_idx, , drop = FALSE]
  out_cc  <- df$outcome[cc_idx]
  
  if (nrow(num_cc) > 1 && ncol(num_cc) > 1) {
    pca_fit <- stats::prcomp(num_cc, center = TRUE, scale. = TRUE)
    
    scores <- as_tibble(pca_fit$x[, 1:2]) %>%
      dplyr::mutate(
        outcome = out_cc,
        dataset = dataset_name
      )
    
    ggplot2::ggplot(
      scores,
      ggplot2::aes(x = PC1, y = PC2, color = outcome, shape = dataset)
    ) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::labs(
        x = "PC1",
        y = "PC2",
        color = "Outcome",
        shape = "Dataset",
        title = paste("PCA of numeric features -", dataset_name)
      ) +
      ggplot2::theme_bw()
  } else {
    message("PCA not run for ", dataset_name, " because of insufficient complete cases or variables")
    invisible(NULL)
  }
}

## run PCA separately for dev and val
pca_plot_dataset(df_dev, "dev")
pca_plot_dataset(df_val, "val")
