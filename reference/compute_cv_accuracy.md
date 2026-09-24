# Internal: Compute Cross-Validation Accuracy for ML Models

Internal function to extract cross-validated accuracy from a list of
trained machine learning models, summarize their median and variability,
and optionally generate a barplot.

## Usage

``` r
compute_cv_accuracy(models, file_name = NULL, return = TRUE)
```

## Arguments

- models:

  Named list of trained ML models. Each model must contain a `$resample`
  data frame with a column named `Accuracy`.

- file_name:

  Optional character. Prefix for saving the accuracy barplot as a PDF in
  the `Results/` directory.

- return:

  Logical. If `TRUE`, saves a barplot of model accuracy values in the
  `Results/` directory.

## Value

A list containing:

- `Accuracy`: Data frame summarizing the median and MAD of accuracy for
  each model.

- `Top_model`: Character string with the model name having the highest
  median accuracy.

## Details

The function assumes that each model contains a `$resample` component
with an `Accuracy` column. Median and MAD (median absolute deviation) of
accuracy are computed for each model. If `return = TRUE`, a PDF barplot
with error bars is created.
