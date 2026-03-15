# Internal: Calculate AUROC from Resample Predictions

Computes the Area Under the ROC Curve (AUROC) for a single
cross-validation resample. This function assumes binary classification
with the positive class labeled `"yes"`.

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

Numeric value of AUROC.
