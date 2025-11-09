# Get performance curves

This function generates and saves the Receiver Operating Characteristic
(ROC) curve and Precision-Recall curve based on the provided metrics. It
also includes the AUC values for both curves in the plot legends.

## Usage

``` r
get_curves(data, spec, sens, reca, prec, color, auc_roc, auc_prc, file.name)
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

  A numeric value representing the AUC for the ROC curve.

- auc_prc:

  A numeric value representing the AUC for the Precision-Recall curve.

- file.name:

  A character string used as the file name prefix for saving the plots.

## Value

Saves two PDF plots: one for the ROC curve and one for the
Precision-Recall curve in the "Results/" directory.
