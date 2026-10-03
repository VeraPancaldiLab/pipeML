# Bar plot of the cross-validation performance of each model

Internal helper used by
[`compute_cv_AUC()`](https://verapancaldilab.github.io/pipeML/reference/compute_cv_AUC.md)
(AUROC, AUPRC) and
[`compute_cv_CINDEX()`](https://verapancaldilab.github.io/pipeML/reference/compute_cv_CINDEX.md)
(C-index), so the cross-validation plots of classification and survival
look the same: models sorted from best to worst, the selected model in
blue, the median of each model written above its bar, error bars of one
MAD (cut at 0 and 1, the limits of the metrics) and a dashed line for
the performance of a random prediction.

## Usage

``` r
plot_cv_metric(
  res,
  metric,
  chance,
  selected_model,
  selected_by,
  n_resamples,
  label = metric,
  wrap_names = FALSE,
  chance_label = "random classifier"
)
```

## Arguments

- res:

  Data frame with one row per model, sorted from best to worst: `model`,
  `Median_<metric>` and `MAD_<metric>`.

- metric:

  Character. Suffix of the median and MAD columns of `res` (`"AUROC"`,
  `"AUPRC"` or `"CINDEX"`).

- chance:

  Numeric. Value of a random prediction (dashed line).

- selected_model:

  Character. Name of the selected model (blue bar).

- selected_by:

  Character. Metric used to select the model, shown in the subtitle.

- n_resamples:

  Integer. Number of cross-validation resamples, shown in the subtitle.

- label:

  Character. Name of the metric in the title and axis. Default:
  `metric`.

- wrap_names:

  Logical. If `TRUE`, the model names are split at `"_"` over several
  lines (long survival model names).

- chance_label:

  Character. Name of the random prediction in the caption. Default:
  `"random classifier"`.

## Value

A ggplot object.
