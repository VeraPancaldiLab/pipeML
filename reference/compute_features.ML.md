# Train and evaluate machine learning models for classification or survival analysis

This function trains and evaluates machine learning models using
cross-validation on training data and then evaluates performance on
independent test data. It supports both **classification** and
**survival analysis** tasks, including hyperparameter tuning and
cohort-based (Leave-One-Dataset-Out, LODO) validation. For survival
models, it computes the **C-index** and generates Kaplan-Meier plots
stratified by predicted risk.

## Usage

``` r
compute_features.ML(
  features_train,
  features_test,
  coldata,
  task_type = c("classification", "survival"),
  trait = NULL,
  trait.positive = NULL,
  time_var = NULL,
  event_var = NULL,
  metric = "AUROC",
  k_folds = 10,
  n_rep = 5,
  LODO = FALSE,
  batch_id = NULL,
  file_name = NULL,
  ncores = NULL,
  return = FALSE,
  fold_construction_fun = NULL,
  fold_construction_args_fixed = NULL,
  fold_construction_args_tunable = NULL,
  seed = 123
)
```

## Arguments

- features_train:

  A data frame or matrix of predictor variables used for training (rows
  = samples, columns = features).

- features_test:

  A data frame or matrix of predictor variables used for testing.

- coldata:

  A data frame containing outcome information. Row names must match
  those of `features_train` and `features_test`.

- task_type:

  Character. Type of task: `"classification"` or `"survival"`.

- trait:

  Character. Column name in `coldata` used as the target variable
  (required for classification tasks).

- trait.positive:

  Value in `trait` that represents the positive class (classification
  only). Ensures all performance metrics and interpretability analyses
  consistently treat the correct class as positive.

- time_var:

  Character. Column name in `coldata` containing survival/follow-up time
  (required for survival tasks).

- event_var:

  Character. Column name in `coldata` indicating event occurrence (1 =
  event occurred, 0 = censored; required for survival tasks).

- metric:

  Character. Performance metric used for model tuning and selection:

  - Classification: `"AUROC"` (default) or `"AUPRC"`.

  - Survival: evaluated using concordance index (C-index).

- k_folds:

  Integer. Number of folds for cross-validation. Default: 10.

- n_rep:

  Integer. Number of repetitions for cross-validation. Default: 5.

- LODO:

  Logical. If `TRUE`, the cross-validation folds are stratified by
  cohort and outcome (see
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)).

- batch_id:

  Character. Column name in `coldata` with the cohort/batch of each
  sample (required if `LODO = TRUE`). The cross-validation folds are
  then stratified by cohort and outcome.

- file_name:

  Character. Base name used to save plots/results under `Results/`. For
  survival tasks, Kaplan-Meier plots are saved as
  `"Results/Survival_KM_<file_name>.pdf"`.

- ncores:

  Integer. Number of CPU cores for parallelization (cross-validation
  folds are processed in parallel). Default: `NULL` (sequential).

- return:

  Logical. Whether to save the plots in `Results/`. Default: `FALSE`.

- fold_construction_fun:

  Function. Optional custom function to construct cross-validation
  folds. Must accept a `bestune` argument internally to inject optimized
  hyperparameters. Used for both classification and survival.
  `features_test` is used as given: it must contain the features built
  by this function (e.g. the test samples projected onto the structure
  learned on the training samples).

- fold_construction_args_fixed:

  List. Fixed arguments passed to `fold_construction_fun` for both CV
  and final training.

- fold_construction_args_tunable:

  List. Arguments passed to `fold_construction_fun` defining
  hyperparameters to explore during CV.

- seed:

  Integer. Random seed for reproducible cross-validation (fold
  assignment and model fitting, including parallel runs). Default:
  `123`. Use `NULL` to leave the random number generator untouched. See
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  for details.

## Value

A named list, or `NULL` if no model could be trained (classification):

- Model:

  The output of
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  on the training set (the selected model is `$Model$Model`).

- AUC:

  Classification: AUROC and AUPRC on the test set, with bootstrap
  confidence intervals (see
  [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)).

- Metrics:

  Classification: threshold-based performance metrics on the test set.

- Prediction:

  Predicted class probabilities (classification) or risk scores
  (survival) of the test samples.

- Curve_bands:

  Classification: pointwise 95% bootstrap bands around the ROC and
  precision-recall curves.

- C_index:

  Survival: C-index on the test set.

## Details

For **classification tasks**, the function performs repeated k-fold
cross-validation with hyperparameter tuning, followed by evaluation on
the test set. ROC and PR curves are generated.

For **survival tasks**, it performs model selection using the C-index,
refits the best model on the full training data and evaluates the
C-index on the test set. With `return = TRUE`, Kaplan-Meier curves of
the test samples split at the median predicted risk are saved, with the
C-index and log-rank test p-value.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Classification ---
data(data_example_classification)
X <- data_example_classification[, setdiff(colnames(data_example_classification), "target")]
set.seed(123)
train_idx <- caret::createDataPartition(data_example_classification$target, p = 0.7, list = FALSE)
res <- compute_features.ML(features_train = X[train_idx, ],
                           features_test = X[-train_idx, ],
                           coldata = data_example_classification,
                           task_type = "classification",
                           trait = "target",
                           trait.positive = "1",
                           k_folds = 5,
                           n_rep = 2,
                           ncores = 2)
res$AUC

# --- Survival ---
data(data_example_survival)
X <- data_example_survival[, setdiff(colnames(data_example_survival), c("time", "status"))]
set.seed(123)
train_idx <- caret::createDataPartition(data_example_survival$status, p = 0.7, list = FALSE)
res_survival <- compute_features.ML(features_train = X[train_idx, ],
                                    features_test = X[-train_idx, ],
                                    coldata = data_example_survival,
                                    task_type = "survival",
                                    time_var = "time",
                                    event_var = "status",
                                    k_folds = 5,
                                    n_rep = 2,
                                    ncores = 2)
res_survival$C_index
} # }
```
