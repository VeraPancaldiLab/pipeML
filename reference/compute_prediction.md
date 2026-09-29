# Compute Prediction Metrics for a Trained Machine Learning Model

Applies a trained model to a test set and evaluates it. For
classification, it computes AUROC and AUPRC with bootstrap confidence
intervals, and Accuracy, Sensitivity, Specificity, Precision, Recall, F1
score and MCC at each probability threshold. For survival, it predicts
risk scores and computes the C-index.

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
  file.name = NULL,
  return = FALSE
)
```

## Arguments

- model:

  The trained model returned as `$Model` by
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  (e.g. `res$Model`).

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

- file.name:

  Character. File name prefix of the plots saved in `Results/`.

- return:

  Logical. Whether to save the plots in `Results/` (ROC and
  precision-recall curves for classification, Kaplan-Meier curves by
  predicted risk group for survival). Default = FALSE.

## Value

For classification, a list containing:

- `Metrics`:

  Data frame of performance metrics (Accuracy, Sensitivity, Specificity,
  Precision, Recall, F1 score, MCC) for each threshold.

- `AUC`:

  List with `AUROC` and `AUPRC`, each a list with the `estimate` on the
  test set and the `lower` and `upper` bounds of its 95% bootstrap
  confidence interval.

- `Predictions`:

  Data frame of predicted probabilities for each class.

- `Curve_bands`:

  List with `ROC` (columns `fpr`, `lower`, `upper`) and `PRC` (columns
  `recall`, `lower`, `upper`): pointwise 95% bootstrap bands around the
  ROC and precision-recall curves.

For survival, a list containing:

- `preds`:

  Predicted risk scores (higher values mean higher risk).

- `c_index`, `c_index_lower`, `c_index_upper`:

  C-index on the test set and its 95% confidence interval.

## Details

Confidence intervals of AUROC and AUPRC come from 1000 bootstrap
resamples of the test samples.

Survival models can predict a risk score, a survival time or a survival
probability. The last two are reversed, so that higher predictions
always mean higher risk.

## See also

[`get_curves`](https://verapancaldilab.github.io/pipeML/reference/get_curves.md),
[`plot_survival_performance`](https://verapancaldilab.github.io/pipeML/reference/plot_survival_performance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(data_example_classification)
X <- data_example_classification[, setdiff(colnames(data_example_classification), "target")]
y <- data_example_classification$target
set.seed(123)
train_idx <- caret::createDataPartition(y, p = 0.7, list = FALSE)

res <- compute_features.training.ML(features_train = X[train_idx, ],
                                    target_var = y[train_idx],
                                    task_type = "classification",
                                    trait.positive = "1",
                                    k_folds = 5,
                                    n_rep = 2)

pred <- compute_prediction(model = res$Model,
                           test_data = X[-train_idx, ],
                           target_var = y[-train_idx],
                           task_type = "classification",
                           trait.positive = "1")
pred$AUC
} # }
```
