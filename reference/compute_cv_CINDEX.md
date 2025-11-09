# Summarize and Visualize C-index Results from Survival Model Cross-Validation

This function aggregates and visualizes the C-index (concordance index)
results obtained from cross-validation of multiple survival models. Each
model should contain a `Resample_matrix` element with per-fold C-index
values. The function computes the median and MAD (median absolute
deviation) of the C-index for each model, identifies the top-performing
model, and optionally generates a bar plot summarizing model
performance.

## Usage

``` r
compute_cv_CINDEX(models, file_name = NULL, plot_results = TRUE)
```

## Arguments

- models:

  A named list of survival model objects, where each element corresponds
  to one fitted model. Each model must contain a `Resample_matrix` data
  frame with columns:

  `c_index`

  :   C-index value per resample (numeric).

  `Resample`

  :   Fold or resample identifier (e.g., "Fold1", "Fold2").

- file_name:

  Optional character string used to name the output PDF file saved under
  `"Results/CINDEX_CV_methods_<file_name>.pdf"`. If `NULL`, the file is
  not named explicitly.

- plot_results:

  Logical; if `TRUE` (default), generates and saves a PDF bar plot
  showing median C-index ± MAD per model.

## Value

A list with the following elements:

- `CINDEX_summary`:

  A tibble summarizing median and MAD per model.

- `All_folds`:

  A tibble with raw C-index values from all folds and models.

- `Top_model`:

  Character string naming the model with the highest median C-index.

## Details

For each model:

- The **median C-index** represents the typical discrimination
  performance across folds.

- The **MAD (Median Absolute Deviation)** measures variability in
  C-index values (robust equivalent of standard deviation).

The plot helps visualize and compare models based on survival prediction
performance, with error bars representing ± MAD.

## See also

[`aggregate_results_survival()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results_survival.md),
[`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)

## Examples

``` r
if (FALSE) { # \dontrun{
models <- list(
  cox_ph_survival = list(
    Resample_matrix = tibble::tibble(
      c_index = c(0.72, 0.75, 0.70),
      Resample = c("Fold1", "Fold2", "Fold3")
    )
  ),
  rand_forest_aorsf = list(
    Resample_matrix = tibble::tibble(
      c_index = c(0.80, 0.82, 0.78),
      Resample = c("Fold1", "Fold2", "Fold3")
    )
  )
)
compute_cv_CINDEX(models, file_name = "example")
} # }
```
