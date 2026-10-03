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
  data frame with `AUROC` and `AUPRC` columns (and, to save the plots,
  `$trainingData` with the outcome in `.outcome`).

- file_name:

  Optional character string. Prefix for saving AUROC/AUPRC plots in the
  `Results/` directory.

- AUC_type:

  Character. Either `"AUROC"` or `"AUPRC"`, used to select the
  top-performing model.

- return:

  Logical. If `TRUE`, saves barplots of AUROC and AUPRC values in the
  `Results/` directory (`AUROC_CV_methods_<file_name>.pdf` and
  `AUPRC_CV_methods_<file_name>.pdf`): models sorted from best to worst,
  selected model highlighted, error bars of one MAD and a dashed line
  for a random classifier
  ([`plot_cv_metric()`](https://verapancaldilab.github.io/pipeML/reference/plot_cv_metric.md),
  also used for the survival C-index).

## Value

A list containing:

- `AUROC`:

  Data frame with the median (`Median_AUROC`) and MAD (`MAD_AUROC`) of
  AUROC for each model, sorted from best to worst.

- `AUPRC`:

  Data frame with the median (`Median_AUPRC`) and MAD (`MAD_AUPRC`) of
  AUPRC for each model, sorted from best to worst.

- `Top_model`:

  Character string: the model with the highest median value for the
  selected metric (`AUC_type`).
