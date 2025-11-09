# Nested Cross-Validation for Survival Models with Optional Custom Fold Construction

This function performs *nested cross-validation* to evaluate and tune
multiple survival models using the tidymodels ecosystem. It supports
both standard event-stratified cross-validation and
*Leave-One-Domain-Out (LODO)* setups, allowing cohort-balanced
performance evaluation. Hyperparameter grids are automatically
constructed for each model type.

## Usage

``` r
compute_k_fold_CV_survival(
  df_features,
  df_outcome,
  outcome_col,
  event_col,
  ml_options = list(nb_folds = 5, nb_repeats = 1, ncores = parallel::detectCores() - 1,
    LODO = FALSE, batch_id = NULL),
  file_name = NULL
)
```

## Arguments

- df_features:

  A data frame containing predictor variables (features).

- df_outcome:

  A data frame containing survival outcomes — typically including
  survival time and event indicator columns.

- outcome_col:

  Character string naming the survival time column.

- event_col:

  Character string naming the event indicator column (`0 = censored`,
  `1 = event`).

- ml_options:

  A named list of cross-validation and parallelization options:

  `nb_folds`

  :   Number of folds for K-fold cross-validation (default = 5).

  `nb_repeats`

  :   Number of repeated CV iterations (default = 1).

  `ncores`

  :   Number of CPU cores to use for parallelization.

  `LODO`

  :   Logical; if `TRUE`, performs Leave-One-Domain-Out stratification
      using `batch_id`.

  `batch_id`

  :   Name of the batch or cohort column (required if `LODO = TRUE`).

- file_name:

  Optional string specifying the suffix of the generated C-index summary
  PDF saved in `"Results/"`.

## Value

A named list containing:

- `Model`:

  The best-performing model retrained on full data.

- `ML_Models`:

  All evaluated survival models with aggregated C-index results.

- `C_index_median`:

  Median C-index of the top model.

- `Custom_output`:

  Optional custom fold output (if applicable).

## Details

Depending on the inputs, the function can:

1.  Build folds internally or accept custom folds from an external
    function.

2.  Train survival models with or without hyperparameter tuning.

3.  Compute and aggregate C-index (Concordance Index) across folds.

4.  Identify and retrain the top-performing model using optimal
    parameters.

Internally, the function:

- Combines predictors and outcomes into a unified dataset.

- Creates stratified folds using rsample, either by event rate or by
  cohort × event combinations (if LODO = TRUE).

- Evaluates a predefined set of survival models: Cox PH, penalized Cox
  (glmnet), AFT models, decision/bagged trees, and random forests.

- Aggregates median and MAD of C-index across resamples.

- Retrains the top model with its optimal hyperparameters.

The function can also interface with a user-defined
`fold_construction_fun()` to support custom preprocessing pipelines
(e.g., CellTFusion), handling folds in parallel and storing intermediate
results to disk.

## See also

[`aggregate_results_survival()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results_survival.md),
[`compute_cv_CINDEX()`](https://verapancaldilab.github.io/pipeML/reference/compute_cv_CINDEX.md),
`wrapper_train_best_hyperparams_survival()`
