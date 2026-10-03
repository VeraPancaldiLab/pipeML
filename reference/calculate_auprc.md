# Calculate AUC from Precision-Recall Curve

This function calculates the Area Under the Precision-Recall Curve
(AUPRC). It uses the trapezoidal rule to compute the AUPRC from the
recall and precision values. The curve is started at recall 0, with the
precision of its first point.

## Usage

``` r
calculate_auprc(recall, precision)
```

## Arguments

- recall:

  A numeric vector of recall values (sensitivity) from the
  precision-recall curve.

- precision:

  A numeric vector of precision values from the precision-recall curve.

## Value

The AUPRC score, a numeric value between 0 and 1.
