# Compute SHAP Values for Machine Learning Models

This function calculates SHAP (SHapley Additive exPlanations) values to
assess feature importance for a trained machine learning model. It
supports both classification and survival tasks, and performs
calculations on cross-validation resamples in parallel. The results can
be summarized and optionally saved with a stability plot.

## Usage

``` r
compute_shap_values(
  model_trained,
  data_train,
  task_type = "classification",
  target_col = NULL,
  trait.positive,
  time_col = NULL,
  event_col = NULL,
  n_cores = 2,
  file.name = NULL
)
```

## Arguments

- model_trained:

  A trained machine learning model object (e.g., output from caret or
  custom ML pipeline), which includes cross-validation resamples.

- data_train:

  A data frame containing the training data used for the model.

- task_type:

  Character. Either `"classification"` (default) or `"survival"`.

- target_col:

  Character. Name of the target column for classification tasks.
  Required if `task_type = "classification"`.

- trait.positive:

  Value representing the positive class in classification tasks.

- time_col:

  Character. Column name representing survival time. Required if
  `task_type = "survival"`.

- event_col:

  Character. Column name representing survival event indicator. Required
  if `task_type = "survival"`.

- n_cores:

  Integer. Number of cores for parallel computation. Default is 2.

- file.name:

  Character. Optional filename prefix for saving SHAP stability plots.
  If `NULL`, plots are not saved.

## Value

A data frame containing SHAP values for all features, averaged across
resamples, with rows corresponding to training samples and columns to
features.

## Details

The function performs the following steps:

1.  Sets up classification or survival prediction functions based on the
    task type.

2.  Loops over all cross-validation resamples in parallel, refitting
    models on training folds.

3.  Computes SHAP values using
    [`fastshap::explain()`](https://bgreenwell.github.io/fastshap/reference/explain.html)
    for each resample, skipping trivial predictions.

4.  Combines SHAP values across resamples and summarizes them (median
    per sample).

5.  Generates and optionally saves a SHAP stability plot if `file.name`
    is provided.

Trivial predictions (constant probability for all samples) are skipped,
and a warning is issued if SHAP values cannot be computed.
