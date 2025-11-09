# Choose Top Base Models for Stacking Based on Accuracy or AUC Scores

This function selects three base models for stacking based on either
Accuracy or AUC metrics. It chooses the top models from different
categories (e.g., tree-based methods, linear models, instance-based
methods) according to the specified metric.

## Usage

``` r
choose_base_models(models, metric = "Accuracy")
```

## Arguments

- models:

  A list of trained machine learning models. Each model must contain a
  `resample` data frame with performance metrics (Accuracy, AUROC,
  AUPRC) from cross-validation.

- metric:

  A character string specifying the metric to use for model selection.
  Can be either "Accuracy", "AUROC", or "AUPRC". Default is "Accuracy".

## Value

A character vector containing the names of the top models selected based
on the specified metric.

## Examples

``` r
if (FALSE) { # \dontrun{
base_models = choose_base_models(models = ml_models, metric = "AUROC")
} # }
```
