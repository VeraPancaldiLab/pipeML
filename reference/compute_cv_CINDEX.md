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

  Named list of survival model objects. Each element must contain a
  `Resample_matrix` data frame with columns:

  `c_index`

  :   Numeric C-index per fold/resample.

  `Resample`

  :   Fold or resample identifier (e.g., "Fold1").

- file_name:

  Optional character string to name the output PDF saved under
  `"Results/CINDEX_CV_methods_<file_name>.pdf"`. If `NULL`, the file
  uses a default naming convention.

- plot_results:

  Logical (default = TRUE). If `TRUE`, generates a PDF bar plot showing
  median C-index ± MAD per model.

## Value

A list with:

- `CINDEX_summary`:

  Tibble summarizing median and MAD per model.

- `All_folds`:

  Tibble of raw C-index values for all models and folds.

- `Top_model`:

  Character string of the model with highest median C-index.

## Details

- Median C-index represents typical discrimination performance across
  folds.

- MAD provides robust variability estimation of C-index values.

- The optional plot displays model performance with error bars ± MAD.

## See also

[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md),
[`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)
