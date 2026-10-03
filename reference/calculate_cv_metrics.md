# Calculate Cross-Validated Metrics for Machine Learning Models

This function computes cross-validated performance metrics for a trained
machine learning model, including AUROC, AUPRC, Accuracy, and other
threshold-based metrics derived from out-of-fold predictions. It also
handles hyperparameter tuning by selecting the optimal parameter set
based on a specified metric.

## Usage

``` r
calculate_cv_metrics(ml_model, metric, hyperparameters = NULL)
```

## Arguments

- ml_model:

  A trained machine learning model object containing `Prediction_folds`
  (out-of-fold predictions), `Results_folds` and `Resample_matrix`.

- metric:

  Character string specifying the metric to optimize when selecting
  hyperparameters. Either `"AUROC"` or `"AUPRC"`.

- hyperparameters:

  Character vector of hyperparameter column names of `Prediction_folds`
  to evaluate. For a model without hyperparameters (treebag),
  `"parameter"` is used: the column that holds `"none"`. With `NULL`, no
  hyperparameter tuning is performed and all resamples get the same
  aggregated value.

## Value

A list containing:

- `Prediction_folds`:

  Data frame of out-of-fold predictions with computed metrics for each
  resample and hyperparameter combination, sorted by resample,
  hyperparameters and decreasing predicted probability.

- `Resample_matrix`:

  Data frame of per-resample performance metrics (AUROC, AUPRC,
  Accuracy) for the tuned model.

- `Results_folds`:

  Median AUROC, AUPRC and Accuracy of each hyperparameter combination.

- `bestTune`:

  Optimal hyperparameter configuration based on the selected metric. Set
  to `"none"` if no hyperparameters are provided.

## Details

The function performs the following steps:

1.  Calculates AUROC and AUPRC for each resample in `Prediction_folds`
    using
    [`calculate_auc_roc_resample()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auc_roc_resample.md)
    and
    [`calculate_auc_prc_resample()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auc_prc_resample.md).

2.  Computes threshold-based metrics (Accuracy, Sensitivity,
    Specificity, Precision, F1-score, MCC) for each resample using
    [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md).

3.  If hyperparameters are provided:

    - Aggregates metrics across all resamples per hyperparameter
      combination.

    - Selects the best hyperparameters based on the `metric` argument.

    - Filters predictions and recomputes metrics for the tuned
      configuration.

4.  If no hyperparameters are provided, aggregates metrics across all
    resamples.

5.  Updates `Resample_matrix` and returns all relevant components.

This workflow ensures robust cross-validation evaluation and
hyperparameter tuning in a reproducible way.
