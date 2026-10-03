# Summarize and Visualize C-index from Survival Model Cross-Validation (Internal)

Aggregates and visualizes C-index (concordance index) results from
cross-validation of multiple survival models. Computes median and MAD
(median absolute deviation) per model, identifies the top-performing
model, and optionally generates a bar plot summarizing performance.

## Usage

``` r
compute_cv_CINDEX(models, file_name = NULL, plot_results = TRUE)
```

## Arguments

- models:

  Named list of survival model objects (from
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)).
  Each element contains a `Resample_matrix` data frame with one row per
  test sample of each resample (for the selected configuration), with
  columns:

  `c_index`

  :   C-index of the resample, repeated on each of its rows.

  `Resample`

  :   Resample identifier (e.g., "Fold1.Rep1").

  Excluded models are `NULL` and are ignored.

- file_name:

  Optional character string to name the output PDF saved under
  `"Results/CINDEX_CV_methods_<file_name>.pdf"`. If `NULL`, the file is
  `"Results/CINDEX_CV_methods_.pdf"`.

- plot_results:

  Logical (default = TRUE). If `TRUE`, saves a bar plot of the median
  C-index +/- MAD of each model, with the same layout as the
  classification AUROC/AUPRC plots
  ([`plot_cv_metric()`](https://verapancaldilab.github.io/pipeML/reference/plot_cv_metric.md)):
  models sorted from best to worst, selected model in blue, values above
  the bars and a dashed line at 0.5 (random prediction).

## Value

A list with:

- `CINDEX_summary`:

  Tibble with one row per model (`model`, `Median_CINDEX`,
  `MAD_CINDEX`), sorted from best to worst.

- `All_folds`:

  Tibble with the rows of all models (`model`, `c_index`, `Resample`):
  one row per test sample of each resample.

- `Top_model`:

  Character string of the model with highest median C-index.

## Details

- The median and MAD are computed over the rows of `Resample_matrix`,
  i.e. each resample counts in proportion to its number of test samples
  (the same as per resample when the folds have the same size).

- The MAD is scaled ([`stats::mad()`](https://rdrr.io/r/stats/mad.html)
  default, comparable to a standard deviation), as for the
  classification AUROC and AUPRC.

- Models without cross-validation results (excluded because they could
  not be fitted) are ignored.

- The optional plot displays model performance with error bars +/- MAD
  (see `plot_results`).

## See also

[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md),
[`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)
