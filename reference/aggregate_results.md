# Aggregate Nested Cross-Validation Results for Classification or Survival Tasks

Aggregates performance metrics from nested cross-validation experiments
for either **classification** or **survival** models. This function
reads the resampled predictions (often saved per fold and hyperparameter
combination) and computes overall performance summaries, identifies the
best hyperparameter configuration, and collates per-resample metrics for
detailed inspection.

## Usage

``` r
aggregate_results(all_loaded, task = c("classification", "survival"))
```

## Arguments

- all_loaded:

  A nested list containing the results of cross-validation runs. Each
  element typically corresponds to a fold and may include one or more
  hyperparameter configurations and model predictions:

  - For classification: `all_loaded[[fold]][[param]][[model]][[1]]`
    contains predictions, and `[[2]]` the corresponding hyperparameters.

  - For survival: `all_loaded[[fold]][[param]][[model]]` contains the
    C-index and hyperparameter columns.

- task:

  Character string specifying the task type. Must be one of:
  `"classification"` or `"survival"`.

## Value

A list of length equal to the number of models evaluated. Each element
contains:

- `Prediction_folds`:

  All predictions or C-index values per fold.

- `Results_folds`:

  Aggregated performance summaries across folds.

- `bestTune`:

  Best-performing hyperparameter combination.

- `Resample_matrix`:

  Fold-level metrics for the best configuration.

## Details

The function automatically detects whether the input structure includes
tunable hyperparameters (`has_params`) by inspecting the nested list
depth.

For **classification** tasks:

- Aggregates fold-level predictions and computes Accuracy and Kappa.

- Computes the median and MAD (robust SD) across resamples.

- Selects the hyperparameter set with the highest median Accuracy.

For **survival** tasks:

- Aggregates per-fold C-index values across hyperparameter
  configurations.

- Computes the median and MAD of the C-index per configuration.

- Identifies and returns the best-performing parameter combination.

The function is compatible with results produced by
[`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)
and analogous classification CV pipelines.

## See also

[`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md),
`calculate_accuracy_kappa_resample()`
