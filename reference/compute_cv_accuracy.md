# Internal: Compute Cross-Validation Accuracy for ML Models

Internal function to extract cross-validated accuracy from a list of
trained machine learning models, summarize their median and variability,
optionally generate a barplot, and select base models for stacking.

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

  Named list of trained ML models. Each model must contain a `$resample`
  data frame with a column named `Accuracy`.

- file_name:

  Optional character. Prefix for saving the accuracy barplot as a PDF in
  the `Results/` directory.

- base_models:

  Logical. If `TRUE`, selects base models using
  [`choose_base_models()`](https://verapancaldilab.github.io/pipeML/reference/choose_base_models.md)
  for stacking.

- return:

  Logical. If `TRUE`, saves a barplot of model accuracy values in the
  `Results/` directory.

## Value

A list containing:

- `Accuracy`: Data frame summarizing the median and MAD of accuracy for
  each model.

- `Top_model`: Character string with the model name having the highest
  median accuracy.

- `Base_models` (optional): Character vector of selected base models if
  `base_models = TRUE`.

## Details

The function assumes that each model contains a `$resample` component
with an `Accuracy` column. Median and MAD (median absolute deviation) of
accuracy are computed for each model. If `return = TRUE`, a PDF barplot
with error bars is created. When `base_models = TRUE`,
[`choose_base_models()`](https://verapancaldilab.github.io/pipeML/reference/choose_base_models.md)
is called to select models for stacking.
