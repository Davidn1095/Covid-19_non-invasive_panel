# conchi_colab

A compact triage-time model for hospital admission risk in acute COVID-19. Code is organized for reproducible figures and tables from RStudio.

## What is here
- Analysis scripts for triage and complete feature sets
- Temporal split with boundary on 2020-04-15
- Training on development split, evaluation on temporal external split
- Models: C4.5, k-NN, SVM, Random Forest and Logistic Regression
- Tuning and threshold selection optimized by MCC
- Calibration was disabled in current experiments

## Quick start
```r
# R >= 4.3 is recommended
install.packages(c("tidyverse","readxl","recipes","caret","pROC","kknn","RWeka","randomForest","glmnet"))

# open the Rproj in RStudio, run the target script, inspect outputs in the plotting pane

