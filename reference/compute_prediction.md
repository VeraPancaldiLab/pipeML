# Compute Prediction Metrics for a Trained Machine Learning Model

This function computes prediction metrics for a given machine learning
model, including the confusion matrix, AUROC, AUPRC, and other
performance metrics such as Accuracy, Sensitivity, Specificity,
Precision, Recall, F1 score, and MCC. The function also determines the
optimal classification threshold based on a chosen metric (e.g.,
Accuracy, F1, or AUROC) and generates a confusion matrix plot.

## Usage

``` r
compute_prediction(
  model,
  test_data,
  target_var,
  trait.positive,
  stack = FALSE,
  file.name = NULL,
  maximize = "Accuracy",
  return = F
)
```

## Arguments

- model:

  The trained machine learning model returned from
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  or
  [`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md).

- test_data:

  A matrix or data frame containing the testing dataset (features only).

- target_var:

  A character vector of true target values for the test data (the
  observed labels).

- trait.positive:

  Value in `target_var` to be considered as the positive class.

- stack:

  Logical. If stacking was used during model training, this parameter
  should be set to TRUE in order to use the meta-learner for prediction.
  Default is FALSE.

- file.name:

  A character string to specify the filename for saving the confusion
  matrix plot (optional). If `NULL`, the plot is not saved.

- maximize:

  A character string indicating which metric to maximize when selecting
  the best threshold for the confusion matrix. Options include
  "Accuracy", "Precision", "Recall", "Specificity", "Sensitivity", "F1",
  or "MCC". Default is "Accuracy".

- return:

  Logical. Whether to return the results and generated plots.

## Value

A list containing:

- Metrics:

  A data frame with various performance metrics (Accuracy, Sensitivity,
  Specificity, Precision, Recall, F1 score, MCC) for each threshold.

- AUC:

  A list containing the AUROC and AUPRC values.

- Predictions:

  A data frame with the predicted probabilities for each class (e.g.,
  `yes` or `no`).

## Details

This function first generates predictions for the test dataset using the
trained machine learning model. It then calculates performance metrics
for a range of threshold values and selects the threshold that maximizes
the chosen metric (e.g., Accuracy, F1 score, etc.). The function returns
the metrics for the best threshold, including AUROC and AUPRC, and
produces a confusion matrix plot that compares predicted versus actual
labels.

The confusion matrix plot is saved as a PDF with the name
`Confusion_Matrix_<file.name>.pdf` if a valid `file.name` is provided.

## See also

[`confusionMatrix`](https://rdrr.io/pkg/caret/man/confusionMatrix.html),
[`varImp`](https://rdrr.io/pkg/caret/man/varImp.html),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
