# Train one machine learning model on one custom cross-validation fold

Internal function used with **custom fold construction functions** (see
package vignette for details). For one fold built by the custom function
and one machine learning model, it trains the model with each
hyperparameter combination of the grid on the training part of the fold
and predicts the test part.

## Usage

``` r
compute_custom_k_fold_CV(processed_folds, ml_method, tuneGrid)
```

## Arguments

- processed_folds:

  A list with the data of one fold, as saved by the custom fold
  construction function: `train_data` (features and `target`),
  `test_data` (features), `rowIndex` (row indices of the test samples),
  `fold_name`, `obs_test` (observed labels of the test samples, in the
  same order as the rows of `test_data`) and, when the function has
  tunable arguments, `params` (the parameter combination used to build
  the features).

- ml_method:

  Character string specifying the machine learning model to use, as
  supported by the `caret` package (e.g., `"rf"`, `"svmRadial"`,
  `"glmnet"`).

- tuneGrid:

  A data frame specifying the grid of hyperparameters to evaluate (one
  row per combination), as returned by
  [`get_tune_grid()`](https://verapancaldilab.github.io/pipeML/reference/get_tune_grid.md).

## Value

A list with two elements:

1.  Data frame of predictions on the test part of the fold, with one row
    per test sample and hyperparameter combination: `rowIndex`,
    `Resample` (fold name), `obs` and `pred` (observed and predicted
    labels), the hyperparameter values, the class probabilities (`no`,
    `yes`) and the columns of `params` (if any).

2.  Character vector with the names of the hyperparameters in
    `tuneGrid`.

## Details

The function performs the following steps for each row of `tuneGrid`:

1.  Train the model on the training part of the fold with that
    hyperparameter combination (a single fit, without resampling).

2.  Predict the class and the class probabilities of the test part of
    the fold.

3.  Store the predictions together with the hyperparameter values.

The predictions of all folds are combined by
[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md),
and the hyperparameters are selected afterwards
([`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md)).
