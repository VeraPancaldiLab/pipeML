# Calculate precision values

This function calculates the precision of a model based on the provided
metrics and the true target values. Precision is defined as the ratio of
true positive predictions to all positive predictions (true positives +
false positives).

## Usage

``` r
calculate_precision(metrics, target)
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

A numeric vector representing the precision values. Precision is the
fraction of true positive predictions among all positive predictions.
