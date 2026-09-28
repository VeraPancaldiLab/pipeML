# Plot SHAP Feature Importance Stability Across Resamples

This function visualizes how stable SHAP feature importance is across
cross-validation resamples. For each resample it computes the global
importance of each feature (mean absolute SHAP value over the samples
held out in that resample), then summarizes these per-resample
importances across resamples as mean +/- standard deviation.

## Usage

``` r
plot_shap_stability(shap_resamples, file.name = NULL, top_n = 20)
```

## Arguments

- shap_resamples:

  A long data frame of per-resample SHAP values, with one row per sample
  and resample, a `Resample` column, a `Samples` column and one numeric
  column per feature. This is the `$shap_resamples` element returned by
  `compute_shap_values(..., return_resamples = TRUE)`.

- file.name:

  Character. Optional filename suffix. If provided, the plot is saved as
  `"Results/SHAP_stability_resample_<file.name>.pdf"`. If `NULL`
  (default), nothing is saved.

- top_n:

  Integer. Number of most important features to show (by mean importance
  across resamples). Default 20. Use `NULL` to show all features.

## Value

A ggplot object (horizontal bar plot of the mean per-resample
importance, with error bars showing the standard deviation across
resamples, truncated at 0).
