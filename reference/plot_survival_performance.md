# Plot Kaplan-Meier Curves by Predicted Risk Group

Stratifies individuals into risk groups based on predicted risk scores
from a fitted survival model, plots Kaplan-Meier survival curves per
risk group, performs a log-rank test, and displays the concordance index
(C-index) with confidence interval. Optionally saves the plot as a PDF
in "Results/".

## Usage

``` r
plot_survival_performance(df_test, prediction, n_groups = 3, file_name = NULL)
```

## Arguments

- df_test:

  Data frame with the observed outcome of the test samples, in columns
  `time` and `event` (1 = event, 0 = censored), in the order of the
  predictions.

- prediction:

  The output of
  [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
  with `task_type = "survival"` (risk scores in `$preds`, C-index and
  its confidence interval in `$c_index`, `$c_index_lower` and
  `$c_index_upper`).

- n_groups:

  Integer. Number of risk groups for stratification (default = 3).

- file_name:

  Optional character. If provided, saves the Kaplan-Meier plot to
  "Results/Survival_KM\_\<file_name\>.pdf".

## Value

Invisibly returns the `ggsurvplot` object for further customization.

## Details

Risk groups are defined by quantiles of the predicted risk scores.
Kaplan-Meier curves visualize survival per risk group, and a log-rank
test assesses differences. The C-index and its 95% confidence interval
are displayed in the plot subtitle.

## Examples

``` r
if (FALSE) { # \dontrun{
# res_survival: output of compute_features.training.ML(task_type = "survival")
pred <- compute_prediction(model = res_survival$Model,
                           test_data = X_test,
                           task_type = "survival",
                           time_var = time_test,
                           event_var = event_test)

plot_survival_performance(df_test = data.frame(time = time_test, event = event_test),
                          prediction = pred,
                          n_groups = 3,
                          file_name = "Example")
} # }
```
