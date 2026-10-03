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
  file.name = NULL,
  width = 6,
  height = 6,
  roc_band = NULL,
  prc_band = NULL
)
```

## Arguments

- data:

  A data frame containing the prediction metrics at each threshold, as
  returned in `compute_prediction()$Metrics` (sorted by decreasing
  predicted probability within each curve).

- spec:

  The name of the column containing the specificity values.

- sens:

  The name of the column containing the sensitivity values.

- reca:

  The name of the column containing the recall values.

- prec:

  The name of the column containing the precision values.

- color:

  The name of the column that identifies each curve (e.g. `"model"` for
  the output of
  [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md),
  or the column with the cohort names). Each value will have a
  corresponding color in the plot. Several curves (several values in
  this column) are only supported with `LODO = TRUE`.

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
  (default = FALSE). `auc_roc` and `auc_prc` must then hold named
  vectors, with the names of the cohorts.

- file.name:

  Optional character string added to the names of the saved plots.

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

No return value. Saves two PDF plots in the "Results/" directory:
`ROC_curve_<file.name>.pdf` for the ROC curve and
`PRC_curve_<file.name>.pdf` for the Precision-Recall curve
(`ROC_curve.pdf` and `PRC_curve.pdf` if `file.name` is `NULL`).

## Details

The ROC curve is drawn from the point (0, 0), and the precision-recall
curve from recall 0 with the precision of its first point, which are the
starting points used to calculate AUROC and AUPRC.

## Examples

``` r
if (FALSE) { # \dontrun{
# pred: output of compute_prediction() (classification)
get_curves(data = pred$Metrics,
           color = "model",
           auc_roc = pred$AUC$AUROC,
           auc_prc = pred$AUC$AUPRC,
           roc_band = pred$Curve_bands$ROC,
           prc_band = pred$Curve_bands$PRC,
           file.name = "Example")
} # }
```
