# Compute Cross-Validation Accuracy for ML Models

This function extracts cross-validated accuracy values from a list of
trained machine learning models, summarizes their median and standard
deviation, and optionally plots a bar chart or selects base models for
stacking.

## Usage

``` r
compute_cv_accuracy(
  models,
  file_name = NULL,
  base_models = FALSE,
  return = TRUE
)
```

## Arguments

- models:

  A named list of trained ML models, each with a `resample` element
  containing cross-validated accuracy.

- file_name:

  (Optional) Character string specifying the filename prefix for the
  saved accuracy plot (PDF format).

- base_models:

  Logical. If `TRUE`, the function selects and returns base models using
  [`choose_base_models()`](https://verapancaldilab.github.io/pipeML/reference/choose_base_models.md)
  for stacking.

- return:

  Logical. If `TRUE`, the function saves a barplot of the model accuracy
  values in the `Results/` directory.

## Value

A list containing:

- `Accuracy`: A data frame with the median and standard deviation of
  accuracy for each model.

- `Top_model`: A character string naming the model with the highest
  median accuracy.

- `Base_models` (optional): A character vector of selected base models
  if `base_models = TRUE`.

## Details

This function assumes that each model in the list has a `$resample`
component containing a column named `Accuracy`. It calculates the median
and standard deviation of accuracy for each model and creates a barplot
(if `return = TRUE`) with error bars.

If `base_models = TRUE`, it calls a helper function
[`choose_base_models()`](https://verapancaldilab.github.io/pipeML/reference/choose_base_models.md)
to select models for use in stacking.
