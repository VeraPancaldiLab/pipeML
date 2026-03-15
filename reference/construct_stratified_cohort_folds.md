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

A named list of fold indices suitable for use in training and
evaluation.

## Details

Each fold preserves the class distribution within each cohort. Useful
for Leave-One-Dataset-Out (LODO) strategies or repeated stratified
k-fold cross-validation.
