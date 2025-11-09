# Calculate AUC from ROC Curve

This function calculates the Area Under the Receiver Operating
Characteristic (ROC) curve. It uses the trapezoidal rule to compute the
AUC from the false positive rate (FPR) and sensitivity (true positive
rate).

## Usage

``` r
calculate_auroc(fpr, sensitivity)
```

## Arguments

- fpr:

  A numeric vector of false positive rates from the ROC curve.

- sensitivity:

  A numeric vector of sensitivities (true positive rates) from the ROC
  curve.

## Value

The AUC score, a numeric value between 0 and 1.
