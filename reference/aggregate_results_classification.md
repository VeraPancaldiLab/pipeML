# Aggregate Cross-Validation Results for Classification Models

This function aggregates the results of cross-validation runs across
folds, parameter combinations, and machine learning models for
classification tasks. It summarizes per-fold performance metrics
(Accuracy and Kappa), computes median and variability measures (MAD),
and identifies the best hyperparameter configuration for each model.

## Usage

``` r
aggregate_results_classification(all_loaded)
```

## Arguments

- all_loaded:

  A nested list containing model evaluation results from
  cross-validation:

  - `all_loaded[[fold]][[param]][[model]][[1]]` — Data frame of
    predictions, including columns `obs`, `pred`, and resample metadata.

  - `all_loaded[[fold]][[param]][[model]][[2]]` — Character vector of
    hyperparameter column names.

  Each prediction data frame must include columns such as: `"obs"`,
  `"pred"`, `"Resample"`, and any hyperparameter columns.

## Value

A list (one element per model), where each element includes:

- `Prediction_folds`:

  Data frame containing all per-fold predictions and metrics.

- `Results_folds`:

  Summary table of median Accuracy and Kappa (and their MADs).

- `Besttune`:

  Data frame with the best-performing hyperparameter configuration.

- `Resample_matrix`:

  Per-fold Accuracy and Kappa for the best configuration.

## Details

It assumes that the input `all_loaded` structure follows the hierarchy:

1.  Folds

2.  Parameter combinations

3.  Models (list per fold and parameter set)

Each model entry should contain:

- A data frame with predictions and observed labels.

- A character vector of hyperparameter names.

The function:

1.  Concatenates predictions across folds and parameter sets.

2.  Computes Accuracy and Kappa per resample using
    `calculate_accuracy_kappa_resample()`.

3.  Aggregates median and MAD (robust variability) of Accuracy and Kappa
    across folds.

4.  Identifies the hyperparameter configuration with the highest median
    Accuracy.

5.  Extracts per-fold metrics for the best configuration.

This allows ranking models by performance stability and selecting the
most reliable hyperparameter configuration.

## See also

`calculate_accuracy_kappa_resample()`,
`compute_k_fold_CV_classification()`

## Examples

``` r
if (FALSE) { # \dontrun{
# Aggregate results from nested CV
aggregated_results <- aggregate_results_classification(all_loaded)

# View best hyperparameters for model 2
aggregated_results[[2]]$Besttune

# Inspect Accuracy and Kappa per fold for best configuration
aggregated_results[[2]]$Resample_matrix
} # }
```
