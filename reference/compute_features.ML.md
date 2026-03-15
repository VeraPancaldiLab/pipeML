# Train and evaluate machine learning models for classification or survival analysis

This function trains and evaluates machine learning models using
cross-validation on training data and then evaluates performance on
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
  coldata,
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
  return = FALSE,
  fold_construction_fun = NULL,
  fold_construction_args_fixed = NULL,
  fold_construction_args_tunable = NULL
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

  Character. Column name in `clinical` used as the target variable
  (required for classification tasks).

- trait.positive:

  Value in `trait` that represents the positive class (classification
  only). Ensures all performance metrics and interpretability analyses
  consistently treat the correct class as positive.

- time_var:

  Character. Column name in `clinical` containing survival/follow-up
  time (required for survival tasks).

- event_var:

  Character. Column name in `clinical` indicating event occurrence (1 =
  event occurred, 0 = censored; required for survival tasks).

- metric:

  Character. Performance metric used for model tuning and selection:

  - Classification: `"Accuracy"`, `"AUROC"`, `"AUPRC"`.

  - Survival: evaluated using concordance index (C-index).

- stack:

  Logical. Perform model stacking (ensemble meta-learning). Default:
  `FALSE`.

- k_folds:

  Integer. Number of folds for cross-validation. Default: 10.

- n_rep:

  Integer. Number of repetitions for cross-validation. Default: 5.

- LODO:

  Logical. If `TRUE`, performs Leave-One-Dataset-Out cross-validation
  based on cohorts.

- batch_id:

  Column name indicating cohort or batch membership for each sample
  (required if `LODO = TRUE`).

- file_name:

  Character. Base name used to save plots/results under `Results/`. For
  survival tasks, Kaplan–Meier plots are saved as
  `"Results/Survival_KM_<file_name>.pdf"`.

- ncores:

  Integer. Number of CPU cores for parallelization. Default:
  `parallel::detectCores() - 1`.

- return:

  Logical. Whether to return and save plots/results. Default: `FALSE`.

- fold_construction_fun:

  Function. Optional custom function to construct cross-validation
  folds. Must accept a `bestune` argument internally to inject optimized
  hyperparameters.

- fold_construction_args_fixed:

  List. Fixed arguments passed to `fold_construction_fun` for both CV
  and final training.

- fold_construction_args_tunable:

  List. Arguments passed to `fold_construction_fun` defining
  hyperparameters to explore during CV.

## Value

A named list containing:

- Model:

  Trained model or workflow (classification) or refitted best model
  (survival).

- Metrics:

  Performance metrics computed on the test data.

- AUC:

  For classification tasks, a list containing AUROC and AUPRC values.

- Prediction:

  Predicted class probabilities (classification) or risk scores
  (survival).

- CV_Results:

  Cross-validation results, including median and MAD of C-index for
  survival tasks.

- Test_CINDEX:

  Concordance index on test data (survival only).

- KM_Plot:

  Kaplan–Meier plot object (if `return = TRUE`).

## Details

For **classification tasks**, the function performs repeated k-fold
cross-validation with hyperparameter tuning, followed by evaluation on
the test set. ROC and PR curves are generated.

For **survival tasks**, it performs model selection using the C-index,
refits the best model on the full training data, evaluates test-set
C-index, and plots Kaplan–Meier curves across quantile-based risk strata
(Low/Medium/High). The C-index and log-rank test p-value are displayed.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Classification Example ---
results_classif <- compute_features.ML(
  features_train = X_train,
  features_test  = X_test,
  coldata        = clin_df,
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
  coldata        = clin_df,
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
