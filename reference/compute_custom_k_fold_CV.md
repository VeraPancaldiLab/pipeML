# Train and evaluate machine learning models on previously constructed k folds

This function performs k-fold cross-validation using custom folds
created from custom functions to be used for cohort-dependent algorithms
(see vignette for more information about this). It supports
hyperparameter tuning over a grid and returns a model object that
mimicks the caret's training output, including performance metrics and
predictions.

## Usage

``` r
compute_custom_k_fold_CV(processed_folds, ml_method, tuneGrid)
```

## Arguments

- processed_folds:

  A list of folds. Each fold contains processed training and test data
  with features.

- ml_method:

  A character string indicating the machine learning model to use, as
  supported by the `caret` package (e.g., `"rf"`, `"svmRadial"`,
  `"glmnet"`).

- tuneGrid:

  Optional. A data frame specifying the grid of hyperparameters to
  evaluate. If `NULL`, a default grid of length 3 is generated using
  caret's `getModelInfo()`.

## Value

A list with the following components:

- `Results_folds`: A data frame summarizing average cross-validated
  Accuracy, Kappa, and their standard deviations for each hyperparameter
  combination.

- `Prediction_folds`: A data frame of predictions from each fold,
  including class probabilities, observed and predicted labels, and
  hyperparameter values.

- `Resample_matrix`: A data frame summarizing Accuracy and Kappa per
  fold for the best-tuned model.

- `Besttune`: List of optimized hyperparameters.

## Details

This function performs the following:

1.  Trains models for each fold and hyperparameter combination.

2.  Predicts on the held-out test data of each fold.

3.  Aggregates prediction results and evaluates Accuracy and Kappa for
    each fold and hyperparameter set.

4.  Selects the best-performing hyperparameter set based on mean
    Accuracy across folds.

5.  Trains the final model on the full dataset using the selected
    hyperparameters.
