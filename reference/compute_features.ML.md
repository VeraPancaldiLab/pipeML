# Train and evaluate machine learning models for classification or survival analysis

This function trains and evaluates machine learning models using
cross-validation on training data and then tests performance on
independent test data. It supports both **classification** and
**survival analysis** tasks, including hyperparameter tuning, model
stacking, and cohort-based (Leave-One-Dataset-Out, LODO) validation. For
survival models, it computes the **C-index** and generates Kaplan–Meier
plots stratified by predicted risk.

## Usage

``` r
compute_features.ML(
  features_train,
  features_test,
  clinical,
  task_type = c("classification", "survival"),
  trait = NULL,
  trait.positive = NULL,
  time_var = NULL,
  event_var = NULL,
  metric = "Accuracy",
  stack = FALSE,
  k_folds = 10,
  n_rep = 5,
  LODO = FALSE,
  batch_id = NULL,
  file_name = NULL,
  ncores = NULL,
  maximize = "Accuracy",
  return = FALSE,
  fold_construction_fun = NULL,
  fold_construction_args_fixed = NULL,
  fold_construction_args_tunable = NULL
)
```

## Arguments

- features_train:

  A data frame or matrix of predictor variables used for training (rows
  as samples, columns as features).

- features_test:

  A data frame or matrix of predictor variables used for testing.

- clinical:

  A data frame containing clinical or outcome information. Row names
  must match those of `features_train` and `features_test`.

- task_type:

  Character. Type of task: either `"classification"` or `"survival"`.

- trait:

  Character. Name of the column in `clinical` used as the target
  variable (required for classification tasks).

- trait.positive:

  Character or numeric. Value in `trait` considered the positive class.

- time_var:

  Character. Name of the column in `clinical` containing survival or
  follow-up time (required for survival tasks).

- event_var:

  Character. Name of the column in `clinical` indicating event
  occurrence (1 = event occurred, 0 = censored; required for survival
  tasks).

- metric:

  Character. Performance metric for model tuning and selection.
  Supported options for classification: `"Accuracy"`, `"AUROC"`,
  `"AUPRC"`. For survival models, performance is evaluated using the
  concordance index (C-index).

- stack:

  Logical. Whether to perform model stacking (default = `FALSE`).

- k_folds:

  Integer. Number of folds for cross-validation (default = 10).

- n_rep:

  Integer. Number of repetitions for cross-validation (default = 5).

- LODO:

  Logical. If `TRUE`, performs Leave-One-Dataset-Out (LODO)
  cross-validation based on cohort identifiers.

- batch_id:

  Column name indicating where the cohort or batch membership for each
  sample is. Required if `LODO = TRUE`.

- file_name:

  Character. Base name used to save plots and results under the
  `Results/` directory. For survival tasks, this will be used to create
  a Kaplan–Meier plot named `"Results/Survival_KM_<file_name>.pdf"`.

- ncores:

  Integer. Number of CPU cores to use for parallelization. If not
  specified, defaults to `parallel::detectCores() - 1`.

- maximize:

  Character. Metric to maximize when selecting the optimal
  classification threshold (options: `"Accuracy"`, `"Precision"`,
  `"Recall"`, `"Specificity"`, `"Sensitivity"`, `"F1"`, or `"MCC"`).
  Default = `"Accuracy"`.

- return:

  Logical. Whether to return and save plots (default = `FALSE`).

- fold_construction_fun:

  Function. Optional custom function to construct cross-validation
  folds. Must accept a `bestune` argument internally for optimized
  hyperparameter injection.

- fold_construction_args_fixed:

  List. Fixed arguments passed to `fold_construction_fun`, used in both
  cross-validation and final model training.

- fold_construction_args_tunable:

  List. Tunable arguments passed to `fold_construction_fun`, defining
  hyperparameters to explore during cross-validation.

## Value

A named list containing:

- Model:

  The trained model or workflow (classification) or refitted best model
  (survival).

- Metrics:

  Performance metrics computed on the test data.

- AUC:

  For classification tasks, a list containing AUROC and AUPRC values.

- Prediction:

  Predicted class probabilities (classification) or risk scores
  (survival).

- CV_Results:

  Cross-validation results, including median and MAD of C-index
  (survival).

- Test_CINDEX:

  Concordance index for the test data (survival only).

- KM_Plot:

  Kaplan–Meier plot object (if `return = TRUE`).

## Details

For **classification tasks**, this function performs cross-validation
tuning based on the chosen performance metric (e.g., Accuracy, AUROC, or
AUPRC), followed by test-set evaluation and ROC/PR curve plotting.

For **survival analysis tasks**, it performs model selection using the
C-index, refits the best model on the full training data, evaluates the
test-set C-index, and plots Kaplan–Meier survival curves across
quantile-based risk strata (Low/Medium/High risk). The C-index and
log-rank test p-value are displayed on the plot.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Classification Example ---
results_classif <- compute_features.ML(
  features_train = X_train,
  features_test  = X_test,
  clinical       = clin_df,
  task_type      = "classification",
  trait          = "Response",
  trait.positive = "Responder",
  k_folds        = 5,
  n_rep          = 1,
  file_name      = "classification_example",
  return         = TRUE
)

# --- Survival Example ---
results_surv <- compute_features.ML(
  features_train = X_train,
  features_test  = X_test,
  clinical       = clin_df,
  task_type      = "survival",
  time_var       = "time",
  event_var      = "status",
  k_folds        = 5,
  n_rep          = 1,
  file_name      = "cox_survival_example",
  return         = TRUE
)
} # }
```
