# Bootstrap-based AUROC and AUPRC Estimation (Internal)

Computes the bootstrap distribution of AUROC (Area Under the ROC Curve)
and AUPRC (Area Under the Precision-Recall Curve) for a set of
predictions. This function resamples the data with replacement and
computes the metrics for each bootstrap iteration, returning the mean,
95% confidence interval, and all bootstrap values.

## Usage

``` r
bootstrap_auc(predict, target, method, B = 1000, seed = 123)
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

## Value

A list with two elements:

- AUROC:

  List with `mean`, `lower`, `upper` 95% CI, and all `values` from
  bootstrap.

- AUPRC:

  List with `mean`, `lower`, `upper` 95% CI, and all `values` from
  bootstrap.
