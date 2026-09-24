# Get performance curves

This function generates and saves the Receiver Operating Characteristic
(ROC) curve and Precision-Recall curve based on the provided metrics. It
also includes the AUC values for both curves in the plot legends.

## Usage

``` r
get_curves(
  data,
  spec = "Specificity",
  sens = "Sensitivity",
  reca = "Recall",
  prec = "Precision",
  color,
  auc_roc,
  auc_prc,
  LODO = FALSE,
  file.name,
  width = 6,
  height = 6,
  roc_band = NULL,
  prc_band = NULL
)
```

## Arguments

- data:

  A data frame containing the prediction metrics.

- spec:

  The name of the column containing the specificity values.

- sens:

  The name of the column containing the sensitivity values.

- reca:

  The name of the column containing the recall values.

- prec:

  The name of the column containing the precision values.

- color:

  The name of the column containing the cohort names. Each cohort will
  have a corresponding color in the plot. Multiple cohorts will result
  in different curves.

- auc_roc:

  A list with elements `estimate`, `lower` and `upper` giving the AUROC
  and its confidence interval, as returned in
  `compute_prediction()$AUC$AUROC`. When `LODO = TRUE`, each element is
  a vector with one value per cohort, named to match the values of the
  `color` column.

- auc_prc:

  Same structure as `auc_roc`, for the AUPRC
  (`compute_prediction()$AUC$AUPRC`).

- LODO:

  Logical. If TRUE, the function assumes the data contains stacked
  predictions from multiple cohorts and assigns AUROC/AUPRC per cohort
  (default = FALSE).

- file.name:

  A character string used as the file name prefix for saving the plots.

- width:

  A numeric value for the width of plot

- height:

  A numeric value for the height of plot

- roc_band:

  Optional data frame with columns `fpr`, `lower`, `upper` (as in
  `compute_prediction()$Curve_bands$ROC`). If supplied and
  `LODO = FALSE`, it is drawn as a shaded pointwise confidence band
  around the ROC curve.

- prc_band:

  Optional data frame with columns `recall`, `lower`, `upper` (as in
  `compute_prediction()$Curve_bands$PRC`), drawn around the
  precision-recall curve in the same way.

## Value

Saves two PDF plots: one for the ROC curve and one for the
Precision-Recall curve in the "Results/" directory.
