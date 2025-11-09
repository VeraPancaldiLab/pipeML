# pipeML

A robust R machine learning pipeline for classification tasks and
survival analysis

## Installation

You can install the development version of `pipeML` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pkg_install("VeraPancaldiLab/pipeML")
```

## Description

`pipeML` is a robust R-based pipeline for building, evaluating, and
interpreting machine learning models in classification tasks. It is
designed for fast, user-friendly deployment while maintaining the
flexibility and rigor required for research-grade analyses. The pipeline
integrates all essential steps — from preprocessing to feature
selection, cross-validation, hyperparameter tuning, and model
interpretation — into a single, consistent framework (Figure 1).

![](reference/figures/MLpipeline.png?raw=true)

*Figure 1. General structure of the `pipeML` machine learning pipeline.*

## Key Features

- Stratified data partitioning to preserve class balance in both
  training and testing sets.

- Iterative Boruta algorithm for robust feature selection.

- Customizable cross-validation with repeated k-fold CV or stratified

  105. 

- Hyperparameter tuning driven by AUROC, AUPRC, or Accuracy.

- Parallelization support for faster cross-validation and model training
  across multiple cores.

- Custom fold construction functions: users can inject their own
  fold-building logic.

These functions can also accept a bestune argument, which is
automatically passed after hyperparameter optimization, allowing
seamless retraining on the full training set with the best parameters.

- Preprocessing utilities for feature filtering (e.g., correlation
  pruning).

- Model interpretation via SHAP values for feature importance.

- Model stacking based on GLM for ensemble learning.

- Visualization functions for ROC and PR curves, plus performance
  summaries.

- Support for 13 machine learning methods, including:

  - Bagged CART
  - Random Forest (RF)
  - C50
  - Logistic regression (LG)
  - CART
  - Naive Bayes (NB)
  - Regularized Lasso
  - Ridge regression
  - Linear Discriminant Analysis (LDA)
  - Regularized Logistic Regression (Elastic net)
  - K-nearest neighbors (KNN)
  - Support vector machine with radial kernel (SVMr)
  - Support vector machine with linear kernel (SVMl)
  - Extreme Gradient Boosting (XGboost)

## General usage

These are basic examples which shows you how to use `pipeML` for
different tasks. For a detailed tutorial, see [Get
started](https://VeraPancaldiLab.github.io/pipeML/articles/pipeML.html)

``` r
library(pipeML)
```

[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md):
This function is designed for training machine learning models on a
single dataset using repeated k-fold cross-validation. It supports
feature selection via Boruta, optional model stacking, and flexible
hyperparameter tuning and the construction of k-folds stratified by
cohorts when this information is available. It can be used when the user
do not account with a prediction dataset, in order to train different
folds on the same dataset and evaluate performance.

``` r
res_ml = compute_features.training.ML(features_train, clinical$Response, "CR", 
                                      metric = "AUROC", stack = F, k_folds = 5, 
                                      n_rep = 10, file_name = "Test", ncores = 2, return = T)
```

After training, predictions on new data can be computed using the
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
function. You can specify which metric to maximize when determining the
optimal classification threshold. Supported values for maximize include:
“Accuracy”, “Precision”, “Recall”, “Specificity”, “Sensitivity”, “F1”,
and “MCC”.

``` r
pred = compute_prediction(res_ml, features_test, traitData_test$Response, 
                          "CR", stack = F, file.name = "Test", 
                          maximize = "Accuracy", return = T)
```

[`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md):
This function is intended for training on a dataset and evaluating on a
separate test dataset when is available. It automatically computes the
prediction using the trained model in the testing set provided. It
includes both previous functions.

``` r
res = compute_features.ML(tme_features_train[[i]], tme_features_test[[i]], 
                          clinical = traitData, trait = "Response", 
                          trait.positive = "R", metric = "AUROC", stack = F, 
                          k_folds = 2, n_rep = 1, LODO = T, batch_id = "Cohort", 
                          ncores = 2, maximize = "Accuracy", return = F)
```

## Issues

If you encounter any problems or have questions about the package, we
encourage you to open an issue
[here](https://github.com/VeraPancaldiLab/pipeML/issues). We’ll do our
best to assist you!

## Authors

`pipeML` was developed by [Marcelo
Hurtado](https://github.com/mhurtado13) in supervision of [Vera
Pancaldi](https://github.com/VeraPancaldi) and is part of the
[Pancaldi](https://github.com/VeraPancaldiLab) team. Currently, Marcelo
is the primary maintainer of this package.

## Citing pipeML

If you use `pipeML` in a scientific publication, we would appreciate
citation to the :

Hurtado M, Pancaldi V (2025). pipeML: A robust R machine learning
pipeline for classification tasks and survival analysis. R package
version 0.0.1, <https://github.com/VeraPancaldiLab/pipeML>.
