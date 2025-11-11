# Train the Best Survival Model Using Optimized Hyperparameters

Fits a survival model on the full training data using the optimal
hyperparameters obtained from nested cross-validation. This wrapper
ensures consistent retraining for different survival model types (Cox,
penalized Cox, AFT, tree-based, or ensemble models), and supports
preprocessing pipelines such as CellTFusion through a user-provided fold
construction function.

## Usage

``` r
wrapper_train_best_hyperparams_survival(
  train_data,
  optimized,
  ml_method,
  fold_construction_fun,
  fold_construction_args_fixed,
  outcome_col = "time",
  event_col = "event"
)
```

## Arguments

- train_data:

  A data frame containing the original training data used for
  cross-validation.

- optimized:

  A list output from `aggregate_results_survival()` or
  [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md),
  containing the best-tuned parameters (`bestTune`) and model
  performance summaries.

- ml_method:

  Character string specifying the survival model to train. Must be one
  of:

  - `"cox_ph_survival"` — Cox proportional hazards model.

  - `"proportional_hazards_glmnet"` — Penalized Cox (elastic net).

  - `"survreg_flexsurv"` — Parametric AFT model.

  - `"rand_forest_partykit"` — Random survival forest via `partykit`.

  - `"rand_forest_aorsf"` — Oblique random survival forest.

  - `"decision_tree_partykit"` — Single survival tree.

  - `"bag_tree_rpart"` — Bagged CART-based survival trees.

  - `"boost_tree_mboost"` — Gradient boosting for censored data.

- fold_construction_fun:

  A custom function used to construct folds and preprocessed data (e.g.,
  `prepare_CellTFusion_folds()`). Must accept arguments `data` and
  optionally `bestune`.

- fold_construction_args_fixed:

  A named list of fixed arguments to pass to `fold_construction_fun()`
  (e.g., paths, deconvolution matrices, etc.).

- outcome_col:

  Character string naming the survival time column (default = `"time"`).

- event_col:

  Character string naming the event indicator column (default =
  `"event"`).

## Value

A named list containing:

- `Model`:

  A list containing the fitted parsnip model object, resampling results,
  and tuning information.

- `training_set`:

  The final preprocessed training dataset used for fitting.

- `custom_output`:

  Additional data returned by the custom fold construction function
  (e.g., CellTFusion outputs or parameter tables).

## Details

This function performs the following steps:

1.  Extracts the optimal hyperparameters from the `optimized` object.

2.  Reconstructs the training dataset using the provided
    `fold_construction_fun()`, including any custom preprocessing or
    feature generation.

3.  Applies the optimal hyperparameters to the model specification.

4.  Fits the final model using the full training data.

If the selected model type has no tunable hyperparameters, the function
automatically detects this and proceeds with the default model
configuration.

## See also

[`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md),
`aggregate_results_survival()`,
[`compute_ml_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_ml_survival.md)
