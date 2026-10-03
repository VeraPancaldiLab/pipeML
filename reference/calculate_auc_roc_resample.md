# Internal: Calculate AUROC from Resample Predictions

Computes the Area Under the ROC Curve (AUROC) for a single
cross-validation resample. This function assumes binary classification
with the positive class labeled `"yes"`. Samples with the same predicted
probability enter the curve together, so the result does not depend on
the order of the samples.

## Usage

``` r
calculate_auc_roc_resample(obs, pred)
```

## Arguments

- obs:

  Vector of observed class labels (`"yes"` / `"no"`).

- pred:

  Numeric vector of predicted probabilities for the positive class
  `"yes"`.

## Value

Numeric value of AUROC. `NA` (with a warning) if the resample has
samples of a single class.
