# Train and evaluate machine learning models on custom cross-validation folds

Internal function that trains and evaluates machine learning models
using pre-constructed k-folds. This function is intended for
**cohort-aware or custom fold strategies** (see package vignette for
details). It supports hyperparameter tuning over a grid and returns a
model object that mimics the structure of caret's `train()` output,
including performance metrics and predictions.

## Usage

``` r
compute_custom_k_fold_CV(processed_folds, ml_method, tuneGrid)
```

## Arguments

- processed_folds:

  A list of folds. Each fold should contain processed training and test
  datasets with features.

- ml_method:

  Character string specifying the machine learning model to use, as
  supported by the `caret` package (e.g., `"rf"`, `"svmRadial"`,
  `"glmnet"`).

- tuneGrid:

  Optional. A data frame specifying the grid of hyperparameters to
  evaluate. If `NULL`, a default grid of length 3 is generated using
  [`caret::getModelInfo()`](https://rdrr.io/pkg/caret/man/modelLookup.html).

## Value

A list containing:

- `Results_folds`: Data frame summarizing average cross-validated
  Accuracy, Kappa, and standard deviations for each hyperparameter
  combination.

- `Prediction_folds`: Data frame of predictions from each fold,
  including class probabilities, observed and predicted labels, and
  hyperparameter values.

- `Resample_matrix`: Data frame summarizing Accuracy and Kappa per fold
  for the best-tuned model.

- `Besttune`: List of optimized hyperparameters.

## Details

The function performs the following steps:

1.  Train models for each fold and hyperparameter combination.

2.  Predict on the held-out test data for each fold.

3.  Aggregate predictions and evaluate Accuracy and Kappa for each fold
    and hyperparameter set.

4.  Select the best-performing hyperparameter set based on mean Accuracy
    across folds.

5.  Train the final model on the full dataset using the selected
    hyperparameters.
