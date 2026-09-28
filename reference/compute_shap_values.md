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
  task_type = "classification",
  n_cores = 2,
  file.name = NULL,
  fold_models_dir = NULL,
  return_resamples = FALSE
)
```

## Arguments

- model_trained:

  The trained model returned as `$Model` by
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  (e.g. `res$Model`). Everything else is taken from it:

  - Classification: the caret `train` object; training data from
    `$trainingData`, outcome `.outcome` (coded `"no"`/`"yes"`), positive
    class `"yes"`.

  - Survival: training data from `$trainingData` (features plus `time`
    and `event` columns, as standardized by
    [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)).

- task_type:

  Character. Either `"classification"` (default) or `"survival"`.

- n_cores:

  Integer. Number of cores for parallel computation. Default is 2.

- file.name:

  Character. Optional filename suffix. If provided, a SHAP stability
  plot (see
  [`plot_shap_stability()`](https://verapancaldilab.github.io/pipeML/reference/plot_shap_stability.md))
  is saved as `"Results/SHAP_stability_resample_<file.name>.pdf"`. If
  `NULL` (default), no plot is saved.

- fold_models_dir:

  Character. Directory where per-fold models saved during training (by
  `compute_custom_k_fold_CV` or `compute_k_fold_CV_survival`) are read
  from, to avoid retraining each resample. If `NULL` (default), uses
  `"Results/fold_models/<task_type>"`, the same default used by
  `compute_features.training.ML`.

- return_resamples:

  Logical. If `TRUE`, also return the per-resample SHAP values. Default
  `FALSE`.

## Value

If `return_resamples = FALSE` (default), a data frame containing SHAP
values for all features, summarized (median) across resamples, with rows
corresponding to training samples (sample IDs as rownames) and columns
to features. Values are in the units of the model output (probability of
the positive class for classification, risk score for survival).

If `return_resamples = TRUE`, a list with:

- `shap`: the summarized data frame described above.

- `shap_resamples`: a long data frame with one row per held-out sample
  and resample, with columns `Resample`, `Samples` and one column per
  feature. Can be passed to
  [`plot_shap_stability()`](https://verapancaldilab.github.io/pipeML/reference/plot_shap_stability.md).

## Details

The function performs the following steps:

1.  Sets up classification or survival prediction functions based on the
    task type.

2.  Loops over all cross-validation resamples in parallel, loading the
    saved fold model. For models trained with `fold_construction_fun`
    (custom folds), every fold model must be found in `fold_models_dir`,
    otherwise the function stops: their fold features were built per
    fold, so a fold cannot be rebuilt from the final training data. For
    models trained on the standard CV path (no saved fold models), each
    fold is refitted on its training samples.

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
