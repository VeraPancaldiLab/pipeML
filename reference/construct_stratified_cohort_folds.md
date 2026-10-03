# Construct Stratified Cohort Folds for Cross-Validation

Generates stratified k-fold cross-validation indices for multiple
datasets while preserving class proportions within each dataset.
Supports multiple repeats.

## Usage

``` r
construct_stratified_cohort_folds(
  train_data,
  batch_id,
  target_id,
  k_folds,
  n_rep
)
```

## Arguments

- train_data:

  Data frame containing the training data.

- batch_id:

  Column name indicating cohort or batch membership for each sample.

- target_id:

  Column name of the target variable used for stratification.

- k_folds:

  Number of folds for cross-validation.

- n_rep:

  Number of repeated cross-validation runs.

## Value

A named list with the row indices of the training samples of each fold,
named `Fold<i>.Rep<j>` (`k_folds` x `n_rep` elements).

## Details

Each cohort is split into `k_folds` folds stratified by the target, and
the folds of all cohorts are merged, so every fold preserves the class
distribution within each cohort. Stops with an error if a cohort has
fewer samples than `k_folds`.
