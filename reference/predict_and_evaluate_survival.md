# Predict and Evaluate Survival Model Performance (Internal)

Generates predictions from a fitted survival model and evaluates
performance using the Concordance Index (C-index). Handles multiple
prediction output types from different survival engines and standardizes
predictions into a comparable numeric format.

## Usage

``` r
predict_and_evaluate_survival(
  model_fit,
  data,
  outcome_col = NULL,
  event_col = NULL,
  ci = TRUE
)
```

## Arguments

- model_fit:

  A fitted survival model object (typically from `parsnip` or
  `workflow`).

- data:

  Data frame with the predictors and, to compute the C-index, the
  survival outcome columns.

- outcome_col:

  Character string specifying the survival time column of `data`.
  Default `NULL`.

- event_col:

  Character string specifying the event indicator column of `data`.
  Default `NULL`. If `outcome_col` or `event_col` is `NULL`, the C-index
  is not computed (`NA`): only predictions are returned.

- ci:

  Logical. If `TRUE` (default), the 95% confidence interval of the
  C-index is computed with
  [`compute_cindex_ci()`](https://verapancaldilab.github.io/pipeML/reference/compute_cindex_ci.md)
  (1000 bootstrap resamples). `FALSE` computes the C-index only, as
  cross-validation does
  ([`compute_ml_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_ml_survival.md)),
  since its confidence interval is not used there.

## Value

A list containing:

- `preds`:

  Tibble with the risk scores in `.pred` (higher = higher risk), one row
  per sample of `data`.

- `c_index`:

  The C-index (`NA` without outcome columns, `NaN` if `data` has no
  event).

- `c_index_lower`:

  Lower bound of the 95% CI of the C-index (`NA` without outcome columns
  or with `ci = FALSE`).

- `c_index_upper`:

  Upper bound of the 95% CI of the C-index (same).

## Details

The function attempts predictions using multiple types depending on
model support:

- `"linear_pred"` - Linear predictor. `parsnip` returns it so that
  higher = longer survival (for Cox models, the log hazard with the sign
  changed); reversed internally.

- `"time"` - Expected survival time (higher = longer survival, reversed
  internally).

- `"survival"` - Survival probability at the median observed time of
  `data` (higher = better survival, reversed internally). Only tried
  when `outcome_col` and `event_col` are given; none of the default
  models needs it.

The C-index is computed before the predictions are reversed, since it
expects higher = longer survival. The returned predictions are risk
scores: higher = higher risk.

Standardizes output into a tibble with a single numeric `.pred` column.
Computes the C-index if outcome/event columns are provided, with its 95%
confidence interval if `ci = TRUE`
([`compute_cindex_ci()`](https://verapancaldilab.github.io/pipeML/reference/compute_cindex_ci.md):
percentile interval over 1000 bootstrap resamples of the samples,
stratified by event, seed 123, without changing the random number state
of the session).
