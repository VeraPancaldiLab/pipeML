# pipeML

A flexible and modular machine learning framework designed to support
leakage-free model training through custom cross-validation fold
construction

## Installation

You can install the development version of `pipeML` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pkg_install("VeraPancaldiLab/pipeML")
```

## Description

`pipeML` is a flexible and leakage-aware machine learning framework for
R designed for predictive modeling in high-dimensional biological data.
The package integrates all key steps of the machine learning workflow —
feature selection, model training, validation, prediction, and
interpretation — into a single reproducible pipeline.

A key design goal of `pipeML` is to support fold-aware feature
construction, allowing features that depend on the dataset
(e.g. enrichment scores, correlation-based features, or network-derived
features) to be recomputed within each cross-validation fold. This
prevents information leakage and ensures reliable performance
estimation.

The framework is designed to integrate naturally with R/Bioconductor
workflows, making it particularly suitable for omics and biomedical
machine learning applications.

![](reference/figures/pipeML.svg?raw=true)

*Figure 1. General structure of the `pipeML` machine learning pipeline.*

## Key Features

### End-to-end ML workflow

- Integrated pipeline for feature selection, model training, validation,
  prediction, and interpretation

### Leakage-aware validation

- Custom cross-validation fold construction
- Support for fold-aware feature recomputation
- Prevents information leakage when using dataset-dependent features

### Flexible model evaluation

- Repeated and stratified k-fold cross-validation
- Leave-one-dataset-out (LODO) evaluation for cross-cohort
  generalization

### Feature selection

- Boruta-based feature selection
- Optional correlation-based feature filtering

### Hyperparameter tuning

- Automatic optimization based on:

  - AUROC
  - AUPRC
  - Accuracy

### Model interpretation

- SHAP-based feature importance
- Variable importance summaries
- Performance visualization (ROC and PR curves)

### Ensemble learning

- Model stacking

### Parallel computing

- Multi-core support for faster model training and cross-validation

### Custom workflows

- Users can define custom fold construction functions
- These functions can receive a bestTune argument after hyperparameter
  optimization to retrain models on the full training dataset.

## Supported Machine Learning Methods

### Classification algorithms:

For classification tasks, we implemented a diverse set of classification
algorithms that are benchmarked on the fly making extensive use of the R
package `caret`.

- Bagged classification trees
- Random forests
- C5.0 decision trees
- Regularized logistic regression (elastic net)
- k-nearest neighbors (KNN)
- Classification and regression trees (CART)
- Lasso regression
- Ridge regression
- Support vector machines with linear and radial kernels
- Extreme Gradient Boosting (XGBoost)

### Survival algorithms:

For time-to-event outcomes, `pipeML` implements a unified survival
modeling framework based on the `parsnip` and `workflows` ecosystems,
enabling consistent training, hyperparameter tuning, and evaluation
across multiple survival model families.

- Cox proportional hazards model
- Elastic net–regularized Cox regression
- Parametric accelerated failure time (AFT) models
- Conditional inference survival trees
- Bagged CART survival models
- Random survival forests
- Gradient boosting for censored outcomes

## General usage

Below are basic examples showing how to use `pipeML`

For a detailed tutorial, see [Get
started](https://VeraPancaldiLab.github.io/pipeML/articles/pipeML.html)

``` r
library(pipeML)
```

### Training models

``` r
res_ml = compute_features.training.ML(features_train, clinical$Response, "CR", 
                                      metric = "AUROC", stack = F, k_folds = 5, 
                                      n_rep = 10, file_name = "Test", ncores = 2, return = T)
```

### Predicting on new data

``` r
pred = compute_prediction(res_ml, features_test, traitData_test$Response, 
                          "CR", stack = F, file.name = "Test", 
                          maximize = "Accuracy", return = T)
```

### Training and Testing Workflow

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
