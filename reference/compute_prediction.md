# Compute Prediction Metrics for a Trained Machine Learning Model

Computes prediction metrics for a trained machine learning model,
including the confusion matrix, AUROC, AUPRC, Accuracy, Sensitivity,
Specificity, Precision, Recall, F1 score, and MCC. For classification
tasks, it also determines the optimal classification threshold and
generates ROC, PRC, and confusion matrix plots. For survival analysis
tasks, it predicts risk scores and optionally generates Kaplan–Meier
plots.

## Usage

``` r
compute_prediction(
  model,
  test_data,
  target_var = NULL,
  trait.positive = NULL,
  task_type = "classification",
  time_var = NULL,
  event_var = NULL,
  stack = FALSE,
  file.name = NULL,
  return = FALSE
)
```

## Arguments

- model:

  The trained machine learning model returned from
  [`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md)
  or
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md).

- test_data:

  A data frame or matrix of predictor variables for the test set.

- target_var:

  Vector of true labels for the test set (classification only).

- trait.positive:

  Value in `target_var` representing the positive class (classification
  only).

- task_type:

  Character. Either `"classification"` or `"survival"`.

- time_var:

  Column or vector of survival/follow-up times (required for survival
  tasks).

- event_var:

  Column or vector of event indicators (1 = event, 0 = censored;
  required for survival tasks).

- stack:

  Logical. If TRUE, uses meta-learner predictions for stacked models
  (classification only).

- file.name:

  Character. Filename prefix for saving plots (optional). If NULL, plots
  are not saved.

- return:

  Logical. Whether to return metrics, predictions, and plots. Default =
  FALSE.

## Value

A list containing:

- `Metrics`:

  Data frame of performance metrics (Accuracy, Sensitivity, Specificity,
  Precision, Recall, F1 score, MCC) for each threshold (classification
  only).

- `AUC`:

  List containing AUROC and AUPRC values with optional bootstrap
  confidence intervals (classification only).

- `Predictions`:

  Data frame of predicted probabilities for each class (classification)
  or risk scores (survival).

## Details

For **classification**, the function:

1.  Uses the trained model (or meta-learner if `stack = TRUE`) to
    predict probabilities for the test data.

2.  Computes performance metrics across thresholds and selects the
    optimal threshold based on a chosen metric.

3.  Calculates AUROC and AUPRC and optionally bootstrapped confidence
    intervals.

4.  Generates ROC, PRC, and confusion matrix plots if `return = TRUE`
    and `file.name` is provided.

For **survival analysis**, the function:

1.  Predicts risk scores using the trained survival model.

2.  Optionally generates Kaplan–Meier plots stratified by predicted risk
    groups.

## See also

[`confusionMatrix`](https://rdrr.io/pkg/caret/man/confusionMatrix.html),
[`varImp`](https://rdrr.io/pkg/caret/man/varImp.html),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
