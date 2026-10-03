# Train the Best Survival Model Using Optimized Hyperparameters

Fits a survival model on the full training data using the optimal
hyperparameters obtained from cross-validation. Used with a custom fold
construction function. This wrapper ensures consistent retraining for
different survival model types (Cox, penalized Cox, AFT, tree-based, or
ensemble models), and supports preprocessing pipelines such as
CellTFusion through a user-provided fold construction function.

## Usage

``` r
wrapper_train_best_hyperparams_survival(
  train_data,
  optimized,
  ml_method,
  fold_construction_fun,
  fold_construction_args_fixed,
  outcome_col = "time",
  event_col = "event",
  preprocess = TRUE
)
```

## Arguments

- train_data:

  A data frame containing the original training data used for
  cross-validation (features plus `time` and `event`), passed to
  `fold_construction_fun` as `data`.

- optimized:

  The cross-validation results of one model (one element of the output
  of
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)):
  `Results_folds`, `Prediction_folds`, `Resample_matrix` and `bestTune`.
  With tunable arguments, `bestTune` is a data frame with the selected
  feature parameters and model hyperparameters; without them, a list
  with the model hyperparameters and `fold_construction_args_fixed`.

- ml_method:

  Character string specifying the survival model to train. Must be one
  of:

  - `"cox_ph_survival"` - Cox proportional hazards model.

  - `"proportional_hazards_glmnet"` - Penalized Cox (elastic net).

  - `"survreg_flexsurv"` - Parametric AFT model.

  - `"rand_forest_partykit"` - Random survival forest via `partykit`.

  - `"rand_forest_aorsf"` - Oblique random survival forest.

  - `"decision_tree_partykit"` - Single survival tree.

  - `"bag_tree_rpart"` - Bagged CART-based survival trees.

  - `"boost_tree_mboost"` - Gradient boosting for censored data.

  `"rand_forest_partykit"` and `"boost_tree_mboost"` are supported here
  but are not run by default (they are commented out of `model_list` in
  [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)).

- fold_construction_fun:

  A custom function used to construct folds and preprocessed data (e.g.,
  `prepare_CellTFusion_folds()`). Must accept the arguments `data` and
  `bestune`. Called with `bestune`, it returns a list with the features
  built on all training samples plus `time` and `event`, its custom
  output, and the selected parameters.

- fold_construction_args_fixed:

  A named list of fixed arguments to pass to `fold_construction_fun()`
  (e.g., paths, deconvolution matrices, etc.).

- outcome_col:

  Character string naming the survival time column (default = `"time"`).

- event_col:

  Character string naming the event indicator column (default =
  `"event"`).

- preprocess:

  Logical. If `TRUE` (default), the training features are preprocessed
  with
  [`preprocess_features()`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md).

## Value

`NULL` if the final fit fails. Otherwise a named list containing:

- `Model`:

  A list containing the model name (`model`), the fitted workflow
  (`fitted`), the cross-validation results (`Results_folds`,
  `Prediction_folds`, `Resample_matrix`) and the selected model
  hyperparameters (`bestTune`; the feature parameters are in
  `custom_output$Parameters`).

- `training_set`:

  The final training set used for fitting (preprocessed if
  `preprocess = TRUE`).

- `custom_output`:

  Additional data returned by the custom fold construction function
  (e.g., CellTFusion outputs), plus its third element in `Parameters`.

## Details

This function performs the following steps:

1.  Extracts the optimal hyperparameters from the `optimized` object.

2.  Reconstructs the training dataset using the provided
    `fold_construction_fun()`, including any custom preprocessing or
    feature generation.

3.  Preprocesses the features
    ([`preprocess_features()`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md)),
    if `preprocess = TRUE`.

4.  Applies the optimal model hyperparameters to the model specification
    (the feature parameters are only used by `fold_construction_fun`).

5.  Fits the final model using the full training data. If the fit fails,
    a warning is given and `NULL` is returned, so the model is excluded.

If the selected model type has no tunable hyperparameters, the function
automatically detects this and proceeds with the default model
configuration.
