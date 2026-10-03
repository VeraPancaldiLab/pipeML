# Aggregate Cross-Validation Results for Classification or Survival Tasks

Aggregates performance metrics from cross-validation experiments for
either **classification** or **survival** models. This function takes
the resampled predictions (per fold and hyperparameter combination) and
computes overall performance summaries, identifies the best
hyperparameter configuration, and collates per-resample metrics for
detailed inspection.

## Usage

``` r
aggregate_results(all_loaded, task = c("classification", "survival"))
```

## Arguments

- all_loaded:

  A nested list containing the results of cross-validation runs. Each
  element corresponds to a fold (resample) and holds the results of
  every model:

  - For classification (output of
    [`compute_custom_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_custom_k_fold_CV.md)):
    `all_loaded[[fold]][[param]][[model]][[1]]` is a data frame with one
    row per test sample and hyperparameter configuration (`rowIndex`,
    `Resample`, `obs`, `pred`, the class probabilities, the
    hyperparameter values and, with tunable arguments, the feature
    parameters), and `[[2]]` the **names** of the model hyperparameters.

  - For survival: `all_loaded[[fold]][[param]][[model]]` is a data frame
    with one row per test sample and hyperparameter configuration:
    `predictions`, the C-index of the fold repeated on each row
    (`c_index`), `Resample`, `model`, the hyperparameter values and,
    with tunable arguments, the feature parameters. A failed fit gives a
    single row with its error message in `fit_error`.

  The `[[param]]` level only exists when the custom fold construction
  function has tunable arguments; otherwise the structure is
  `all_loaded[[fold]][[model]]`.

- task:

  Character string specifying the task type. Must be one of:
  `"classification"` or `"survival"`.

## Value

A list of length equal to the number of models evaluated. Each element
(`NULL` for a survival model without any configuration left) contains:

- `Prediction_folds`:

  The rows of all folds and configurations bound together (see
  `all_loaded`).

- `Results_folds`:

  One row per configuration: classification `Accuracy`, `Kappa`,
  `AccuracySD`, `KappaSD` (median and MAD across resamples); survival
  `c_index_median` and `c_index_mad`.

- `bestTune`:

  The selected configuration: model hyperparameters plus, with tunable
  arguments, the feature parameters.

- `Resample_matrix`:

  Results of the selected configuration: classification, `Accuracy` and
  `Kappa` per resample; survival, its rows of `Prediction_folds` (one
  per test sample, with the C-index of the fold).

## Details

The function detects from the nesting of `all_loaded` whether the custom
fold construction function had tunable arguments (`has_params`). The
model hyperparameters are always there; the feature parameters (tunable
arguments) are then handled as extra hyperparameter columns.

For **classification** tasks:

- Binds the predictions of all folds and computes Accuracy and Kappa per
  resample and configuration.

- Computes their median and MAD (scaled, robust SD) across resamples.

- Selects the configuration with the highest median Accuracy. This
  `bestTune` is only provisional:
  [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
  selects it again by AUROC or AUPRC with
  [`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md).

For **survival** tasks:

- Binds the C-index rows of all folds and configurations.

- Reports the fits that failed, and the resamples where the C-index
  could not be computed (`NA`, e.g. a held-out fold without events), and
  excludes every configuration with a failed fit or an `NA` C-index in
  at least one resample, so all configurations are compared over the
  same resamples. A model without any configuration left is `NULL`.

- Computes the median and MAD (scaled, robust SD, as for classification)
  of the C-index per configuration.

- Selects the configuration with the highest median C-index.

The function is compatible with results produced by
[`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)
and analogous classification CV pipelines.
