# Bootstrap-based AUROC and AUPRC Estimation (Internal)

Computes the bootstrap distribution of AUROC (Area Under the ROC Curve)
and AUPRC (Area Under the Precision-Recall Curve) for a set of
predictions. This function resamples the data with replacement and
computes the metrics for each bootstrap iteration, returning the mean,
95% confidence interval, and all bootstrap values.

## Usage

``` r
bootstrap_auc(predict, target, method, B = 1000, seed = 123, n_grid = 101)
```

## Arguments

- predict:

  Numeric vector or matrix of predicted values (scores).

- target:

  Numeric or factor vector of observed binary outcomes.

- method:

  Character string specifying the ML model or method name.

- B:

  Integer. Number of bootstrap iterations (default = 1000).

- seed:

  Integer. Random seed for reproducibility (default = 123).

- n_grid:

  Integer. Number of evenly spaced points in \[0, 1\] at which the
  pointwise curve bands are evaluated (default = 101).

## Value

A list with four elements:

- AUROC:

  List with `mean`, `lower`, `upper` 95% CI, and all `values` from
  bootstrap.

- AUPRC:

  List with `mean`, `lower`, `upper` 95% CI, and all `values` from
  bootstrap.

- ROC_band:

  Data frame with columns `fpr`, `lower`, `upper`: pointwise 95%
  bootstrap band for sensitivity at each false positive rate.

- PRC_band:

  Data frame with columns `recall`, `lower`, `upper`: pointwise 95%
  bootstrap band for precision at each recall level.
