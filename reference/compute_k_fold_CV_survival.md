# Nested Cross-Validation for Survival Models with Optional Custom Fold Construction

Performs *nested cross-validation* to evaluate and tune multiple
survival models using the **tidymodels** ecosystem. Supports both
standard event-stratified cross-validation and *Leave-One-Domain-Out
(LODO)* setups, enabling cohort-balanced model evaluation.
Hyperparameter grids are automatically constructed for each model type.

## Usage

``` r
compute_k_fold_CV_survival(
  df_features,
  df_outcome,
  outcome_col,
  event_col,
  k_folds,
  n_rep,
  ncores,
  LODO = FALSE,
  batch_id = NULL,
  file_name = NULL,
  fold_construction_fun = NULL,
  fold_construction_args_fixed = NULL,
  fold_construction_args_tunable = NULL
)
```

## Arguments

- df_features:

  A data frame of predictor variables (features).

- df_outcome:

  A data frame containing survival outcomes — typically including
  survival time and event indicator columns.

- outcome_col:

  Character string giving the name of the survival time column.

- event_col:

  Character string giving the name of the event indicator column
  (`0 = censored`, `1 = event`).

- k_folds:

  Integer. Number of folds for K-fold cross-validation (default = 5).

- n_rep:

  Integer. Number of repeated CV iterations (default = 1).

- ncores:

  Integer. Number of CPU cores to use for parallelization.

- LODO:

  Logical; if `TRUE`, performs Leave-One-Domain-Out cross-validation
  using `batch_id` to stratify samples by cohort.

- batch_id:

  Optional character string naming the column representing cohort or
  batch identifiers. Required if `LODO = TRUE`.

- file_name:

  Optional string specifying the suffix for the generated C-index
  summary PDF saved in the `"Results/"` directory.

- fold_construction_fun:

  Optional custom function for constructing data folds. Used to
  interface with external preprocessing workflows (e.g., CellTFusion).

- fold_construction_args_fixed:

  Optional list of fixed arguments passed to `fold_construction_fun()`.

- fold_construction_args_tunable:

  Optional list of tunable arguments passed to `fold_construction_fun()`
  during hyperparameter tuning.

## Value

A named list with the following elements:

- `Model`:

  The best-performing survival model retrained on the full dataset.

- `ML_Models`:

  All evaluated survival models with aggregated C-index results.

- `C_index_median`:

  Median C-index of the top-performing model.

- `Custom_output`:

  Optional list of custom outputs from fold construction.

## Details

Depending on the inputs, the function can:

1.  Build folds internally or accept custom folds from a user-defined
    function.

2.  Train survival models with or without hyperparameter tuning.

3.  Compute and aggregate the Concordance Index (C-index) across folds.

4.  Identify and retrain the top-performing model using optimal
    parameters.

Internally, the function:

- Merges predictors and outcomes into a single dataset.

- Creates stratified folds using **rsample**, either by event rate or by
  cohort × event combinations (if `LODO = TRUE`).

- Evaluates a predefined set of survival models: Cox PH, penalized Cox
  (glmnet), AFT (flexsurv), decision trees, bagged trees, and random
  forests.

- Aggregates the median and MAD of the C-index across resamples.

- Retrains the best-performing model with its optimal hyperparameters.

When a custom fold construction function is provided via
`fold_construction_fun`, the function handles folds in parallel, saves
intermediate results under `"Results/"`, and returns additional outputs
for advanced integration.

## See also

[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md),
[`compute_cv_CINDEX()`](https://verapancaldilab.github.io/pipeML/reference/compute_cv_CINDEX.md),
[`wrapper_train_best_hyperparams_survival()`](https://verapancaldilab.github.io/pipeML/reference/wrapper_train_best_hyperparams_survival.md)
