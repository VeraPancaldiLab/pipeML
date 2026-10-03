# Repeated K-fold Cross-Validation for Survival Models (Internal)

Performs repeated k-fold cross-validation with hyperparameter tuning for
survival models using the **tidymodels** ecosystem. Supports both
standard event-stratified K-fold CV and Leave-One-Dataset-Out (LODO)
setups. Hyperparameter grids are automatically generated for each model
type.

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
  fold_construction_args_tunable = NULL,
  seed = 123,
  preprocess = TRUE
)
```

## Arguments

- df_features:

  Data frame of predictor variables (features). With `LODO = TRUE`, it
  also contains the cohort column named by `batch_id`.

- df_outcome:

  Data frame of survival outcomes (time and event columns).

- outcome_col:

  Character. Name of the survival time column. Must be `"time"`: the
  custom fold path, the final fit and the prediction functions use the
  columns `time` and `event` directly, so other names are not supported.

- event_col:

  Character. Name of the event indicator column (`0 = censored`,
  `1 = event`). Must be `"event"` (see `outcome_col`).

- k_folds:

  Integer. Number of folds for K-fold CV.

- n_rep:

  Integer. Number of repeated CV iterations.

- ncores:

  Integer. Number of CPU cores. Without `fold_construction_fun`, the
  folds run in parallel; with it and tunable arguments, the parameter
  combinations of each fold run in parallel. `NULL` or 1: sequential.

- return:

  Logical. Whether to save the plot of the cross-validation C-index in
  `"Results/"`.

- LODO:

  Logical. If `TRUE`, the folds are stratified by cohort (`batch_id`)
  and event.

- batch_id:

  Character. Name of the column of `df_features` with the cohort/batch.
  Required if `LODO = TRUE`. The column is removed after the folds are
  built, so the cohort is not a predictor.

- file_name:

  Optional string. Suffix for generated C-index summary PDF saved in
  `"Results/"`.

- fold_construction_fun:

  Optional custom function that builds the features of each fold
  (fold-aware features). It receives `data` (the features plus the
  `time` and `event` columns), `folds` and `bestune`: with
  `bestune = NULL` it saves each fold in `Results/fold_<fold name>.rds`;
  with `bestune` it returns the features built on all training samples
  (plus `time` and `event`), its custom output and the selected
  parameters. See the `fold_construction_fun` argument of
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  for the full contract.

- fold_construction_args_fixed:

  Optional list of fixed arguments passed to `fold_construction_fun`.

- fold_construction_args_tunable:

  Optional list of tunable arguments passed to `fold_construction_fun`
  during hyperparameter tuning.

- seed:

  Integer. Random seed set before the folds are drawn, so fold
  assignment and model fitting are reproducible. In the parallel
  custom-fold branch, each worker iteration is seeded from `seed`, the
  fold and the parameter configuration, so results do not depend on
  worker scheduling. `NULL` leaves the random number generator
  untouched.

- preprocess:

  Logical. If `TRUE` (default), near-zero variance and highly correlated
  features are removed with
  [`preprocess_features()`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md):
  once on all training samples before the cross-validation, or, with
  `fold_construction_fun`, on the training part of each fold and on the
  final training set.

## Value

A named list containing:

- `Model`:

  The selected model: a list with the cross-validation results
  (`Results_folds`, `Prediction_folds`, `Resample_matrix`), the selected
  hyperparameters (`bestTune`), the workflow fitted on all training
  samples (`Model_object`) and its training data (`trainingData`:
  features plus `time` and `event`).

- `ML_Models`:

  All survival models with their cross-validation results (`NULL` for a
  model excluded because no configuration could be fitted and
  evaluated).

- `C_index_median`:

  Median cross-validation C-index of the selected model.

- `Custom_output`:

  Only with `fold_construction_fun`: the custom output of the selected
  model, plus the selected feature parameters in `Parameters`.

## Details

The function can:

1.  Build folds internally or accept a custom fold construction
    function.

2.  Train multiple survival models with optional hyperparameter tuning.

3.  Compute and aggregate the Concordance Index (C-index) across folds.

4.  Train the selected model (best median C-index) on all training
    samples with its tuned hyperparameters. With
    `fold_construction_fun`, every model is trained on the features
    rebuilt on all training samples, and the best one is kept; without
    it, only the selected model is retrained.

The function internally:

- Merges predictors and outcomes.

- Creates stratified folds using **rsample**, either by event or by
  cohort x event (LODO).

- Preprocesses the features (`preprocess = TRUE`): once before the
  cross-validation, or, with `fold_construction_fun`, on each fold's
  training table and on the final training set.

- Builds the hyperparameter grids with
  [`get_default_hyperparams()`](https://verapancaldilab.github.io/pipeML/reference/get_default_hyperparams.md).
  With `fold_construction_fun`, the `mtry` of the forests is sized from
  the features built by the function (the fold with the fewest
  features), not from the input columns.

- Evaluates predefined survival models: Cox PH, penalized Cox (glmnet),
  AFT (flexsurv), decision trees, bagged trees, and random forests.

- Aggregates the median and MAD of C-index across resamples with
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md),
  which excludes every configuration with a failed fit or an `NA`
  C-index in a resample (a model without any configuration left is
  excluded).

- Trains the selected model on all training samples (see above). The
  function stops with a clear message if no model is left.

Without `fold_construction_fun`, the folds run in parallel when
`ncores > 1`. With `fold_construction_fun`, the folds are processed one
after the other; if it has tunable arguments, the parameter combinations
of each fold run in parallel. Additional outputs of the custom function
are returned in `Custom_output`.

## See also

[`compute_ml_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_ml_survival.md),
[`get_default_hyperparams()`](https://verapancaldilab.github.io/pipeML/reference/get_default_hyperparams.md),
[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)
