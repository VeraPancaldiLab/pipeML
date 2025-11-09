# Plot and save survival performance of a model on test data

This function groups individuals into risk strata (e.g.,
Low/Medium/High) based on their predicted risk scores from a fitted
survival model. It then plots Kaplan–Meier survival curves for each risk
group, including a log-rank test for group separation and displays the
concordance index (C-index) of the model on the test data.

## Usage

``` r
plot_survival_performance(
  df_test,
  c_index = NULL,
  n_groups = 3,
  file_name = NULL
)
```

## Arguments

- df_test:

  A data frame containing at least the following columns:

  time

  :   Observed survival or follow-up time (numeric).

  event

  :   Event indicator (1 = event occurred, 0 = censored).

  .pred

  :   Predicted risk score or linear predictor from the model (higher
      values indicate higher risk).

- c_index:

  Numeric. The concordance index (C-index) computed on the test data.

- n_groups:

  Integer. Number of risk groups to stratify by (default = 3). Typically
  3 groups correspond to "Low", "Medium", and "High" risk strata.

- file_name:

  Character (optional). If provided, the Kaplan–Meier plot will be saved
  to `"Results/Survival_KM_<file_name>.pdf"`.

## Value

Invisibly returns the `ggsurvplot` object for further customization, and
saves a PDF of the plot in the "Results/" directory if `file_name` is
provided.

## Details

Risk groups are defined by quantile-based cut points on the predicted
risk scores. The log-rank test is used to assess whether survival curves
differ significantly between risk strata. The C-index is displayed for
interpretability.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_survival_performance(
  df_test = test_data,
  c_index = 0.74,
  n_groups = 3,
  file_name = "cox_model_test"
)
} # }
```
