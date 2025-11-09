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

  A vector of predicted class labels or probabilities.

- target:

  A vector of true class labels.

## Value

The F1 score, a numeric value between 0 and 1.
