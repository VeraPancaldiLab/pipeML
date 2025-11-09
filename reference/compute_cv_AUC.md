# Compute Cross-Validated AUC Values for Machine Learning Models

This function computes cross-validated AUROC and AUPRC scores for a list
of trained machine learning models. It can also save performance
barplots and optionally select base models for stacking.

## Usage

``` r
compute_cv_AUC(
  models,
  file_name = NULL,
  base_models = FALSE,
  AUC_type = "AUROC",
  return = TRUE
)
```

## Arguments

- models:

  A named list of trained machine learning models. Each model should
  contain a `resample` data frame with AUROC and AUPRC values from
  cross-validation.

- file_name:

  (Optional) Character string. Used as the prefix for the plot filenames
  if `save_plot = TRUE`.

- base_models:

  Logical. If `TRUE`, selects a subset of models as base learners for
  stacking using the
  [`choose_base_models()`](https://verapancaldilab.github.io/pipeML/reference/choose_base_models.md)
  function.

- AUC_type:

  Character. Either `"AUROC"` or `"AUPRC"`; determines which metric is
  used to select the top-performing model.

- return:

  Logical. Whether to return the results and generated plots.

## Value

A list containing:

- `AUROC`:

  A data frame with median and standard deviation of AUROC values for
  each model.

- `AUPRC`:

  A data frame with median and standard deviation of AUPRC values for
  each model.

- `Top_model`:

  The name of the model with the highest median value for the selected
  metric (`AUC_type`).

- `Base_models`:

  (Optional) A character vector of selected base models for stacking,
  returned if `base_models = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- compute_cv_AUC(
  models = ml_models,
  file_name = "Model_Performance",
  base_models = TRUE,
  AUC_type = "AUROC",
  save_plot = TRUE
)
} # }
```
