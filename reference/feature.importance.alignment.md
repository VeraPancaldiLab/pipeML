# Align Feature Importance Based on Direction of Association

Adjusts variable importance values according to the direction of
association with the target variable. Positive coefficients retain the
original importance sign, while negative coefficients flip it.

## Usage

``` r
feature.importance.alignment(model)
```

## Arguments

- model:

  A trained machine learning model object containing
  `result$Variable_importance`, `result$Cell_groups`, and
  `result$Model`.

## Value

A data frame of feature importance values with aligned signs in the
`final_importance` column.

## Details

This function fits a univariate logistic regression of each feature
against the outcome, then multiplies the original importance by 1 or -1
depending on the sign of the regression coefficient.
