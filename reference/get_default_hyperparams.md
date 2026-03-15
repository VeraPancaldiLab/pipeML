# Get default hyperparameter grids for supported survival models

This function generates default hyperparameter grids for various
survival models compatible with the tidymodels framework. It uses the
default parameter ranges from the dials package and produces a sequence
of evenly spaced values within those ranges for each tunable
hyperparameter.

## Usage

``` r
get_default_hyperparams(model_name, train_x = NULL, levels = 5, v = 5)
```

## Arguments

- model_name:

  Character string specifying the model name. Supported options include:

  - `"cox_ph_survival"` – Classic Cox proportional hazards model

  - `"proportional_hazards_glmnet"` – Penalized Cox (LASSO / Elastic
    Net)

  - `"survreg_flexsurv"` – Parametric accelerated failure time (AFT)

  - `"decision_tree_partykit"` – Single survival tree

  - `"bag_tree_rpart"` – Bagged CART survival trees

  - `"rand_forest_partykit"` – Random survival forest (ctree-based)

  - `"rand_forest_aorsf"` – Oblique random survival forest

  - `"boost_tree_mboost"` – Gradient boosting for survival

- train_x:

  Optional data frame or matrix of training predictors. This is required
  for parameters that depend on the number of features, such as `mtry`,
  which determines the number of variables randomly sampled at each
  split in tree-based models.

- levels:

  Integer specifying how many values to generate per hyperparameter.
  Defaults to `5`. Must be at least 2.

- v:

  Integer. Number of folds for K-fold cross-validation (default = 5).

## Value

A named list of hyperparameter grids. Each element is a numeric vector
of sampled values for that parameter. Returns `NULL` for models without
tunable hyperparameters.

## Details

The function supports models such as Cox proportional hazards (regular
and penalized), parametric survival regression, decision trees, bagging,
random forests, and gradient boosting. For models without tunable
parameters (e.g., classic Cox or AFT models), the function returns
`NULL`.

The helper function `vs()` internally calls
[`dials::value_seq()`](https://dials.tidymodels.org/reference/value_validate.html)
to generate evenly spaced sequences of parameter values across their
default ranges. For data-dependent parameters (like `mtry`),
[`dials::finalize()`](https://dials.tidymodels.org/reference/finalize.html)
is used to compute appropriate limits based on `train_x`.
