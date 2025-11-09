# Calculate recall values

This function calculates the recall (also known as sensitivity) of a
model based on the provided metrics and the true target values. Recall
is defined as the ratio of true positive predictions to all actual
positive instances (true positives + false negatives).

## Usage

``` r
calculate_recall(metrics, target)
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

A numeric vector representing the recall values. Recall is the fraction
of actual positive instances that were correctly predicted.
