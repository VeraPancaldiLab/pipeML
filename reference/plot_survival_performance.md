# Plot and Save Survival Performance of a Model (Internal)

Stratifies individuals into risk groups based on predicted risk scores
from a fitted survival model, plots Kaplan–Meier survival curves per
risk group, performs a log-rank test, and displays the concordance index
(C-index) with confidence interval. Optionally saves the plot as a PDF
in "Results/".

## Usage

``` r
plot_survival_performance(df_test, prediction, n_groups = 3, file_name = NULL)
```

## Arguments

- df_test:

  Data frame containing observed survival, event indicator and predicted
  risk score.

- prediction:

  List containing prediction results

- n_groups:

  Integer. Number of risk groups for stratification (default = 3).

- file_name:

  Optional character. If provided, saves the Kaplan–Meier plot to
  "Results/Survival_KM\_\<file_name\>.pdf".

## Value

Invisibly returns the `ggsurvplot` object for further customization.

## Details

Risk groups are defined by quantiles of the predicted risk scores.
Kaplan–Meier curves visualize survival per risk group, and a log-rank
test assesses differences. The C-index and its 95% confidence interval
are displayed in the plot subtitle.
