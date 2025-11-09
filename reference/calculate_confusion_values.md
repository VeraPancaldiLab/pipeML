# Calculate confusion values

This function computes the confusion matrix values (True Positives,
False Negatives, True Negatives, and False Positives) based on the given
metrics and the true target values.

## Usage

``` r
calculate_confusion_values(metrics, target)
```

## Arguments

- metrics:

  A data frame with metrics obtained using
  [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md),
  containing at least two columns: "Sensitivity" and "Specificity".

- target:

  A character vector containing the true values from the target
  variable. It should have the same length as the predictions.

## Value

A list containing four elements:

- `TP`: True Positives

- `FN`: False Negatives

- `TN`: True Negatives

- `FP`: False Positives
