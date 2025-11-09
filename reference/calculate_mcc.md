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

  A vector of predicted class labels or probabilities.

- target:

  A vector of true class labels.

## Value

The MCC score, a numeric value between -1 and 1.
