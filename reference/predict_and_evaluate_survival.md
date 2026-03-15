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
  event_col = NULL
)
```

## Arguments

- model_fit:

  A fitted survival model object (typically from `parsnip` or
  `workflow`).

- data:

  Data frame containing predictors and survival outcome variables.

- outcome_col:

  Character string specifying the survival time column (default =
  `"time"`).

- event_col:

  Character string specifying the event indicator column (default =
  `"event"`).

## Value

A list containing:

- `preds`:

  Tibble with standardized numeric predictions `.pred`.

- `c_index`:

  Numeric value of the computed C-index.

- `c_index_lower`:

  Lower bound of 95% CI for the C-index.

- `c_index_upper`:

  Upper bound of 95% CI for the C-index.

## Details

The function attempts predictions using multiple types depending on
model support:

- `"linear_pred"` — Linear predictor/log hazard (higher = higher risk).

- `"time"` — Expected survival time (higher = longer survival, reversed
  internally).

- `"survival"` — Survival probability at a median evaluation time
  (higher = better survival, reversed internally).

Standardizes output into a tibble with a single numeric `.pred` column.
Computes C-index using
[`compute_cindex_ci()`](https://verapancaldilab.github.io/pipeML/reference/compute_cindex_ci.md)
if outcome/event columns are provided.
