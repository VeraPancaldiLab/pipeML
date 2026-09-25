# Compute SHAP Values for Machine Learning Models

This function calculates SHAP (SHapley Additive exPlanations) values to
assess feature importance for a trained machine learning model. It
supports both classification and survival tasks, and performs
calculations on cross-validation resamples in parallel. Each sample is
explained only by the fold model(s) that held it out, and its SHAP
values are summarized across repeats by the median.

## Usage

``` r
compute_shap_values(
  model_trained,
  data_train = NULL,
  task_type = "classification",
  time_col = NULL,
  event_col = NULL,
  n_cores = 2,
  file.name = NULL,
  fold_models_dir = NULL
)
```

## Arguments

- model_trained:

  A trained machine learning model object, which includes
  cross-validation resamples. For classification, this must be the caret
  `train` object returned as `$Model` by
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  (e.g. `res$Model`): the training data (`$trainingData`), the outcome
  (`.outcome`, coded `"no"`/`"yes"`) and the positive class (`"yes"`)
  are all taken from it.

- data_train:

  Survival only. A data frame containing the training data used for the
  model. Ignored for classification, where the training data is taken
  from `model_trained$trainingData`.

- task_type:

  Character. Either `"classification"` (default) or `"survival"`.

- time_col:

  Character. Column name representing survival time. Required if
  `task_type = "survival"`.

- event_col:

  Character. Column name representing survival event indicator. Required
  if `task_type = "survival"`.

- n_cores:

  Integer. Number of cores for parallel computation. Default is 2.

- file.name:

  Character. Currently unused (the SHAP stability plot is disabled);
  kept for backward compatibility.

- fold_models_dir:

  Character. Directory where per-fold models saved during training (by
  `compute_custom_k_fold_CV` or `compute_k_fold_CV_survival`) are read
  from, to avoid retraining each resample. If `NULL` (default), uses
  `"Results/fold_models/<task_type>"`, the same default used by
  `compute_features.training.ML`.

## Value

A data frame containing SHAP values for all features, summarized
(median) across resamples, with rows corresponding to training samples
and columns to features. Values are in the units of the model output
(probability of the positive class for classification, risk score for
survival).

## Details

The function performs the following steps:

1.  Sets up classification or survival prediction functions based on the
    task type.

2.  Loops over all cross-validation resamples in parallel, loading the
    saved fold model (or refitting it on the training fold if no
    matching file is found).

3.  Computes SHAP values using
    [`fastshap::explain()`](https://bgreenwell.github.io/fastshap/reference/explain.html)
    on the held-out samples of each resample, skipping resamples with
    trivial predictions.

4.  Combines SHAP values across resamples and summarizes them (median
    per sample).

A summary of how many resamples were loaded, retrained or skipped is
printed. Samples that could not be explained in any resample are
reported with a warning. If no resample could be explained, `NULL` is
returned.
