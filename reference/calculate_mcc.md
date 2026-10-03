# Compute Matthews Correlation Coefficient (MCC) Score

The Matthews correlation coefficient (MCC) is a metric used for binary
classification problems. It takes into account true and false positives
and negatives, and is considered a balanced metric.

## Usage

``` r
calculate_mcc(metrics, target)
```

## Arguments

- metrics:

  A data frame with metrics obtained using
  [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md),
  containing at least two columns: "Sensitivity" and "Specificity" (one
  row per probability threshold).

- target:

  A vector of true class labels (`"yes"` / `"no"`).

## Value

A numeric vector with the MCC score (between -1 and 1) at each
probability threshold.
