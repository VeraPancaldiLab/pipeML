# Bootstrap-based AUROC and AUPRC Estimation (Internal)

Computes the bootstrap distribution of AUROC (Area Under the ROC Curve)
and AUPRC (Area Under the Precision-Recall Curve) for a set of
predictions. This function resamples the data with replacement and
computes the metrics for each bootstrap iteration, returning the mean,
95% confidence interval, and all bootstrap values. The bootstrap is
stratified: the positive and negative samples are resampled separately,
keeping their numbers, so every resample contains both classes (as pROC
does for AUC intervals). The intervals and bands are the 2.5% and 97.5%
percentiles of the resamples.

## Usage

``` r
bootstrap_auc(predict, target, method, B = 1000, seed = 123, n_grid = 101)
```

## Arguments

- predict:

  Data frame of predicted class probabilities, with a column `yes`
  (probability of the positive class), one row per sample.

- target:

  Vector of observed outcomes (`"yes"` / `"no"`), in the order of the
  rows of `predict`.

- method:

  Character string specifying the ML model or method name.

- B:

  Integer. Number of bootstrap iterations (default = 1000).

- seed:

  Integer. Random seed of the bootstrap resamples (default = 123). It is
  only used inside the function: the random number state of the session
  is restored afterwards.

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
