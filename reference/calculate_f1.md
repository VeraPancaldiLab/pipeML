# Compute F1 Score

The F1 score is the harmonic mean of precision and recall, and is used
to evaluate the balance between the two metrics. It is particularly
useful when the class distribution is imbalanced.

## Usage

``` r
calculate_f1(metrics, target)
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

A numeric vector with the F1 score (between 0 and 1) at each probability
threshold.
