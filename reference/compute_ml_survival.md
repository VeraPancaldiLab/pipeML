# Train and Evaluate a Survival Model

This function trains and evaluates a single survival model using the
**tidymodels** framework. It supports a range of survival model types,
including Cox, penalized Cox, AFT, random forests, bagged trees, and
gradient boosting models.

## Usage

``` r
compute_ml_survival(
  df_train,
  df_test,
  outcome_col,
  event_col = NULL,
  model,
  models_hyperparameters
)
```

## Arguments

- df_train:

  A data frame containing the training data, including survival time,
  event indicator, and predictor variables.

- df_test:

  A data frame containing the test data with the same structure and
  columns as `df_train`.

- outcome_col:

  Character string specifying the name of the survival time column.

- event_col:

  Character string specifying the name of the event indicator column
  (`1 = event occurred`, `0 = censored`).

- model:

  Character string specifying which survival model to train. Supported
  values include:

  - `"cox_ph_survival"` – Cox proportional hazards model (survival)

  - `"proportional_hazards_glmnet"` – Penalized Cox model (LASSO/Elastic
    Net)

  - `"survreg_flexsurv"` – Parametric AFT model (flexsurv)

  - `"rand_forest_partykit"` – Random survival forest (ctree engine)

  - `"rand_forest_aorsf"` – Oblique random survival forest

  - `"decision_tree_partykit"` – Single survival tree

  - `"bag_tree_rpart"` – Bagged survival trees

  - `"boost_tree_mboost"` – Gradient boosting for survival data

- models_hyperparameters:

  A list containing model hyperparameter values to apply. Typically
  created from a tuning grid or optimization step, e.g.:
  `list(list(trees = 500, min_n = 10))`. If `NULL`, default parameters
  are used.

## Value

A tibble with two columns:

- `model`:

  The model name as a character string.

- `c_index`:

  Numeric value of the computed C-index on the test data.

## Details

The function automatically applies user-specified hyperparameters,
standardizes predictions from different engines, and evaluates model
performance using the **Concordance Index (C-index)** as the primary
metric. It ensures a consistent risk-score direction across models
(higher values indicate higher risk).

The function uses the **parsnip** interface from tidymodels to define,
train, and evaluate models.

Workflow:

1.  Defines the model specification
    ([`parsnip::model_spec`](https://parsnip.tidymodels.org/reference/model_spec.html)).

2.  Optionally applies user-specified hyperparameters via
    [`parsnip::set_args()`](https://parsnip.tidymodels.org/reference/set_args.html).

3.  Constructs a survival formula of the form `Surv(time, event) ~ .`.

4.  Fits the model using
    [`parsnip::fit()`](https://generics.r-lib.org/reference/fit.html) on
    the training data.

5.  Evaluates model performance using
    [`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md),
    which computes the C-index.

## See also

[`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)
