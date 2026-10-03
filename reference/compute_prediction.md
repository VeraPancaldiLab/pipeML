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
  return = FALSE,
  n_groups = 2
)
```

## Arguments

- model:

  The trained model returned as `$Model` by
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  (e.g. `res$Model`).

- test_data:

  A data frame of predictor variables for the test set (for
  classification, a matrix is also accepted). It must contain all the
  features of the model; other columns are ignored for classification.

- target_var:

  Vector of true labels for the test set, in the order of the rows of
  `test_data` (classification only).

- trait.positive:

  Value in `target_var` representing the positive class (classification
  only).

- task_type:

  Character. Either `"classification"` or `"survival"`.

- time_var:

  Numeric vector of survival/follow-up times of the test samples, in the
  order of the rows of `test_data` (required for survival tasks).

- event_var:

  Numeric vector of event indicators of the test samples (1 = event, 0 =
  censored; required for survival tasks).

- file.name:

  Character. File name prefix of the plots saved in `Results/`.

- return:

  Logical. Whether to save the plots in `Results/` (ROC and
  precision-recall curves for classification, Kaplan-Meier curves by
  predicted risk group for survival). Default = FALSE.

- n_groups:

  Integer. Number of risk groups of the Kaplan-Meier plot (survival
  only, used with `return = TRUE`). Default = 2. See
  [`plot_survival_performance()`](https://verapancaldilab.github.io/pipeML/reference/plot_survival_performance.md).

## Value

For classification, a list containing:

- `Metrics`:

  Data frame of performance metrics (Accuracy, Sensitivity, Specificity,
  Precision, Recall, F1 score, MCC) for each threshold: one row per test
  sample, sorted by decreasing predicted probability.

- `AUC`:

  List with `AUROC` and `AUPRC`, each a list with the `estimate` on the
  test set and the `lower` and `upper` bounds of its 95% bootstrap
  confidence interval.

- `Predictions`:

  Data frame of predicted probabilities for each class (columns `no` and
  `yes`), in the order of the rows of `test_data`.

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
resamples of the test samples, stratified by class (positives and
negatives resampled separately, so every resample has both), and the
confidence interval of the C-index from 1000 bootstrap resamples
stratified by event. Both use the fixed seed 123 only inside the
bootstrap: the state of the random number generator of the session is
not changed.

Survival models can predict a linear predictor, a survival time or a
survival probability. The modelling library returns all three so that
higher values mean longer survival (also the linear predictor of Cox
models); they are reversed, so that higher predictions always mean
higher risk.

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
