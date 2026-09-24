# Internal: Compute Cross-Validated AUROC and AUPRC for ML Models

Internal function to summarize cross-validated AUROC and AUPRC values
from a list of trained machine learning models. Computes median and MAD
for each model and optionally generates barplots.

## Usage

``` r
compute_cv_AUC(models, file_name = NULL, AUC_type = "AUROC", return = TRUE)
```

## Arguments

- models:

  Named list of trained ML models. Each model must contain a `$resample`
  data frame with `AUROC` and `AUPRC` columns.

- file_name:

  Optional character string. Prefix for saving AUROC/AUPRC plots in the
  `Results/` directory.

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
