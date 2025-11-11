# Predict and Evaluate Survival Model Performance

This function generates predictions from a fitted survival model and
evaluates its performance using the Concordance Index (C-index). It
automatically handles different prediction output types supported by
various survival modeling engines and standardizes predictions into a
comparable numeric format.

## Usage

``` r
predict_and_evaluate_survival(
  model_fit,
  data,
  outcome_col = "time",
  event_col = "event"
)
```

## Arguments

- model_fit:

  A fitted survival model object (typically created via a `parsnip` or
  `workflow` model specification).

- data:

  A data frame containing the predictors and survival outcome variables
  used for prediction and evaluation.

- outcome_col:

  Character string specifying the name of the survival time column in
  `data`. Default is `"time"`.

- event_col:

  Character string specifying the name of the event indicator column in
  `data`. Default is `"event"`.

## Value

A list with two elements:

- `preds`:

  A tibble containing standardized numeric predictions.

- `c_index`:

  Numeric value representing the computed C-index.

## Details

The function attempts to generate predictions using multiple possible
prediction types, depending on what the model supports:

- `"linear_pred"` — Linear predictor or log hazard (e.g., Cox models),
  where higher values indicate higher risk.

- `"time"` — Expected survival time (e.g., AFT models), where higher
  values imply longer survival (automatically reversed).

- `"survival"` — Survival probability at a specific evaluation time,
  where higher values indicate better survival (automatically reversed).

It standardizes the prediction output into a tibble with one numeric
column `.pred`, ensuring consistency across engines. The function then
computes the Concordance Index (C-index) using
[`yardstick::concordance_survival_vec()`](https://yardstick.tidymodels.org/reference/concordance_survival.html)
to assess how well the model ranks predicted survival relative to actual
outcomes.

## See also

[`yardstick::concordance_survival_vec()`](https://yardstick.tidymodels.org/reference/concordance_survival.html),
[`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fitted_model <- parsnip::fit(
  parsnip::proportional_hazards(engine = "survival"),
  Surv(time, event) ~ .,
  data = lung
)
results <- predict_and_evaluate_survival(fitted_model, lung)
results$c_index
} # }
```
