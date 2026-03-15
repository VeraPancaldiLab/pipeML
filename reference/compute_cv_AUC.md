# Internal: Compute Cross-Validated AUROC and AUPRC for ML Models

Internal function to summarize cross-validated AUROC and AUPRC values
from a list of trained machine learning models. Computes median and MAD
for each model, optionally generates barplots, and can select base
models for stacking.

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

  Named list of trained ML models. Each model must contain a `$resample`
  data frame with `AUROC` and `AUPRC` columns.

- file_name:

  Optional character string. Prefix for saving AUROC/AUPRC plots in the
  `Results/` directory.

- base_models:

  Logical. If `TRUE`, selects a subset of models as base learners for
  stacking using
  [`choose_base_models()`](https://verapancaldilab.github.io/pipeML/reference/choose_base_models.md).

- AUC_type:

  Character. Either `"AUROC"` or `"AUPRC"`, used to select the
  top-performing model.

- return:

  Logical. If `TRUE`, saves barplots of AUROC and AUPRC values in the
  `Results/` directory.

## Value

A list containing:

- `AUROC`:

  Data frame with median and MAD of AUROC for each model.

- `AUPRC`:

  Data frame with median and MAD of AUPRC for each model.

- `Top_model`:

  Character string: the model with the highest median value for the
  selected metric (`AUC_type`).

- `Base_models`:

  (Optional) Character vector of selected base models if
  `base_models = TRUE`.
