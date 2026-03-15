# Internal: Calculate Accuracy and Kappa from Resample Predictions

Computes accuracy and Cohen's kappa statistic for a single
cross-validation resample. Ensures factors are aligned and calculates
expected agreement for kappa.

## Usage

``` r
calculate_accuracy_kappa_resample(obs, pred)
```

## Arguments

- obs:

  Vector of observed class labels (factor or character).

- pred:

  Vector of predicted class labels (factor or character), must match
  levels of `obs`.

## Value

Named numeric vector with `Accuracy_resample` and `Kappa_resample`.
