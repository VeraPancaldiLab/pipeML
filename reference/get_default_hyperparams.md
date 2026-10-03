# Get default hyperparameter grids for supported survival models

This function generates default hyperparameter grids for various
survival models compatible with the tidymodels framework. It uses the
default parameter ranges from the dials package and produces `levels`
evenly spaced values within those ranges for each tunable hyperparameter
(on the log scale for `penalty`). Only arguments that parsnip passes to
the engine are tuned. `min_n` is capped so that a node of the training
part of a fold can still be split, and `mtry` at the number of features;
values made equal by these caps are kept once.

## Usage

``` r
get_default_hyperparams(model_name, train_x = NULL, levels = 5, v = 5)
```

## Arguments

- model_name:

  Character string specifying the model name. Supported options include:

  - `"cox_ph_survival"` - Classic Cox proportional hazards model

  - `"proportional_hazards_glmnet"` - Penalized Cox (LASSO / Elastic
    Net)

  - `"survreg_flexsurv"` - Parametric accelerated failure time (AFT)

  - `"decision_tree_partykit"` - Single survival tree

  - `"bag_tree_rpart"` - Bagged CART survival trees

  - `"rand_forest_partykit"` - Random survival forest (ctree-based)

  - `"rand_forest_aorsf"` - Oblique random survival forest

  - `"boost_tree_mboost"` - Gradient boosting for survival

- train_x:

  Optional data frame or matrix of training predictors (features only).
  Its number of columns gives the range of `mtry` (number of features
  sampled at each split of the forests): without `train_x`, the forest
  grids have no `mtry` and the engine's default is used. Its number of
  rows, with `v`, caps `min_n` of the tree and the forests; without it,
  `min_n` is not capped.

- levels:

  Integer specifying how many values to generate per hyperparameter.
  Defaults to `5`. Must be at least 2.

- v:

  Integer. Number of folds for K-fold cross-validation (default = 5),
  used to cap `min_n` to the size of the training part of a fold.

## Value

A named list with one numeric vector of values per tuned argument (the
grid is all their combinations), or `NULL` for models without tunable
hyperparameters:

- `"proportional_hazards_glmnet"`: `penalty`, `mixture`.

- `"decision_tree_partykit"`: `tree_depth`, `min_n`.

- `"bag_tree_rpart"`: `trees`, the number of bagged trees (25 to 100).
  It is passed to the engine as `nbagg`
  ([`ipred::bagging()`](https://rdrr.io/pkg/ipred/man/bagging.html)) by
  [`compute_ml_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_ml_survival.md)
  and
  [`wrapper_train_best_hyperparams_survival()`](https://verapancaldilab.github.io/pipeML/reference/wrapper_train_best_hyperparams_survival.md).

- `"rand_forest_aorsf"`, `"rand_forest_partykit"`: `trees`, `min_n`,
  `mtry`.

- `"boost_tree_mboost"`: `trees`, `min_n`, `tree_depth`.

## Details

The function supports models such as Cox proportional hazards (regular
and penalized), parametric survival regression, decision trees, bagging,
random forests, and gradient boosting. For models without tunable
parameters (Cox proportional hazards and the AFT model), the function
returns `NULL`.

The helper function `vs()` internally calls
[`dials::value_seq()`](https://dials.tidymodels.org/reference/value_validate.html)
to generate evenly spaced sequences of parameter values across their
default ranges. For data-dependent parameters (like `mtry`),
[`dials::finalize()`](https://dials.tidymodels.org/reference/finalize.html)
is used to compute appropriate limits based on `train_x`.
