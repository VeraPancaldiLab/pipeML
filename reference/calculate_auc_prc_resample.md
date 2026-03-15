# Internal: Calculate AUPRC from Resample Predictions

Computes the Area Under the Precision-Recall Curve (AUPRC) for a single
cross-validation resample. Assumes binary classification with positive
class `"yes"`.

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

Numeric value of AUPRC.
