# Nested Cross-Validation for Survival Models (Internal)

Performs nested cross-validation for survival models using the
**tidymodels** ecosystem. Supports both standard event-stratified K-fold
CV and Leave-One-Domain-Out (LODO) setups. Hyperparameter grids are
automatically generated for each model type.

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
  return = FALSE,
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

  Data frame of predictor variables (features).

- df_outcome:

  Data frame of survival outcomes (time and event columns).

- outcome_col:

  Character. Name of the survival time column.

- event_col:

  Character. Name of the event indicator column (`0 = censored`,
  `1 = event`).

- k_folds:

  Integer. Number of folds for K-fold CV (default = 5).

- n_rep:

  Integer. Number of repeated CV iterations (default = 1).

- ncores:

  Integer. Number of CPU cores for parallelization.

- return:

  Logical. Whether to return the generated plots.

- LODO:

  Logical. If `TRUE`, performs Leave-One-Domain-Out CV using `batch_id`.

- batch_id:

  Character. Name of column representing cohort/batch. Required if
  `LODO = TRUE`.

- file_name:

  Optional string. Suffix for generated C-index summary PDF saved in
  `"Results/"`.

- fold_construction_fun:

  Optional custom function to construct data folds.

- fold_construction_args_fixed:

  Optional list of fixed arguments passed to `fold_construction_fun`.

- fold_construction_args_tunable:

  Optional list of tunable arguments passed to `fold_construction_fun`
  during hyperparameter tuning.

## Value

A named list containing:

- `Model`:

  The best-performing survival model retrained on the full dataset.

- `ML_Models`:

  All evaluated survival models with aggregated C-index results.

- `C_index_median`:

  Median C-index of the top-performing model.

- `Custom_output`:

  Optional outputs from the custom fold construction function.

## Details

The function can:

1.  Build folds internally or accept a custom fold construction
    function.

2.  Train multiple survival models with optional hyperparameter tuning.

3.  Compute and aggregate the Concordance Index (C-index) across folds.

4.  Retrain the top-performing model using its optimal hyperparameters.

The function internally:

- Merges predictors and outcomes.

- Creates stratified folds using **rsample**, either by event or by
  cohort × event (LODO).

- Evaluates predefined survival models: Cox PH, penalized Cox (glmnet),
  AFT (flexsurv), decision trees, bagged trees, and random forests.

- Aggregates the median and MAD of C-index across resamples.

- Retrains the top-performing model on the full dataset.

If `fold_construction_fun` is provided, the function handles folds in
parallel and returns additional outputs for advanced integration.

## See also

[`compute_ml_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_ml_survival.md),
[`get_default_hyperparams()`](https://verapancaldilab.github.io/pipeML/reference/get_default_hyperparams.md),
[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)
