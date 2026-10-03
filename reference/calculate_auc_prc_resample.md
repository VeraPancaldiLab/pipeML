# Internal: Calculate AUPRC from Resample Predictions

Computes the Area Under the Precision-Recall Curve (AUPRC) for a single
cross-validation resample. Assumes binary classification with positive
class `"yes"`. Samples with the same predicted probability enter the
curve together, so the result does not depend on the order of the
samples.

## Usage

``` r
calculate_auc_prc_resample(obs, pred)
```

## Arguments

- obs:

  Vector of observed class labels (`"yes"` / `"no"`).

- pred:

  Numeric vector of predicted probabilities for the positive class
  `"yes"`.

## Value

Numeric value of AUPRC. `NA` (with a warning) if the resample has
samples of a single class.
