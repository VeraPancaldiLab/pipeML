# Compute variable importance using SHAP values

Computes the variable importance for a machine learning model using SHAP
(SHapley Additive exPlanations) values.

## Usage

``` r
compute_variable.importance(model, stacking = FALSE, n_cores = 2)
```

## Arguments

- model:

  A trained machine learning model.

- stacking:

  A logical value indicating whether the model was trained using
  stacking (default is FALSE).

- n_cores:

  An integer specifying the number of workers to use for parallel
  computation (default is 2).

## Value

A data frame with SHAP values representing the variable importance for
each feature.

## Details

If `stacking` is TRUE, the function computes the SHAP values for each
base model in the stacked ensemble model and averages them. If
`stacking` is FALSE, the function computes the SHAP values for the
provided single machine learning model. The computed SHAP values are
returned as a data frame with features as rows and samples as columns.
