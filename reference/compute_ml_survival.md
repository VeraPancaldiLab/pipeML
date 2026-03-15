# Train and Evaluate a Survival Model (Internal)

This internal function trains and evaluates a survival model using the
**tidymodels** framework. It supports Cox, penalized Cox, AFT, random
forests, bagged trees, and gradient boosting survival models.

## Usage

``` r
compute_ml_survival(
  df_train,
  df_test = NULL,
  outcome_col,
  event_col,
  model,
  models_hyperparameters,
  return_model = F
)
```

## Arguments

- df_train:

  A data frame containing training data including survival time, event
  indicator, and predictors.

- df_test:

  Optional data frame for testing. Must contain the same columns as
  `df_train`.

- outcome_col:

  Character. Name of the survival time column.

- event_col:

  Character. Name of the event indicator column (1 = event occurred, 0 =
  censored).

- model:

  Character. The type of survival model to train. Options:

  - `"cox_ph_survival"` – Cox proportional hazards model

  - `"proportional_hazards_glmnet"` – Penalized Cox (LASSO/Elastic Net)

  - `"survreg_flexsurv"` – Parametric AFT model

  - `"rand_forest_partykit"` – Random survival forest (ctree engine)

  - `"rand_forest_aorsf"` – Oblique random survival forest

  - `"decision_tree_partykit"` – Single survival tree

  - `"bag_tree_rpart"` – Bagged survival trees

  - `"boost_tree_mboost"` – Gradient boosting for survival data

- models_hyperparameters:

  Optional list of hyperparameter values to apply. Example:
  `list(list(trees = 500, min_n = 10))`. Defaults to `NULL` (use engine
  defaults).

- return_model:

  Logical. If `TRUE`, returns the fitted model along with predictions.
  Default is `FALSE`.

## Value

If `df_test` is provided:

- A tibble with columns `predictions` (predicted risk scores) and
  `c_index` (Concordance Index).

- If `return_model = TRUE`, a list with elements:

  - `Model` – fitted tidymodels workflow

  - `Metrics` – tibble with `predictions` and `c_index`. If `df_test` is
    `NULL`, the function returns only the fitted model object.

## Details

The function standardizes predictions across different engines and
evaluates model performance using the Concordance Index (C-index). Risk
scores are aligned such that higher values indicate higher risk.

The function uses the **parsnip** interface to define and fit survival
models. Workflow:

1.  Create model specification
    ([`parsnip::model_spec`](https://parsnip.tidymodels.org/reference/model_spec.html))
    based on `model` type.

2.  Apply optional hyperparameters via
    [`parsnip::set_args()`](https://parsnip.tidymodels.org/reference/set_args.html).

3.  Construct survival formula: `Surv(time, event) ~ .`.

4.  Fit model using
    [`parsnip::fit()`](https://generics.r-lib.org/reference/fit.html) on
    `df_train`.

5.  Evaluate performance using
    [`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)
    (C-index).

## See also

[`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)
