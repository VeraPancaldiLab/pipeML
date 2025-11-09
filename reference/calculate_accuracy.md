# Calculates accuracy values from prediction

This function calculates the accuracy of a model based on the provided
metrics and the true target values. The accuracy is computed as the
ratio of correct predictions (both true positives and true negatives) to
the total number of predictions.

## Usage

``` r
calculate_accuracy(metrics, target)
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

A numeric vector representing the accuracy values. The result is the
fraction of correct predictions out of all predictions.
