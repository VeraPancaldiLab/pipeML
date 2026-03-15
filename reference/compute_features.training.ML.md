# Train machine learning or survival models with optional stacking and custom cross-validation

This function trains one or more machine learning models using repeated
k-fold cross-validation, with optional model stacking, feature
selection, and support for both classification and survival tasks. It
allows flexible cross-validation schemes, including:

- Standard stratified k-fold cross-validation

- Leave-One-Dataset-Out (LODO) stratified folds by cohort

- User-defined custom fold construction via a `fold_construction_fun`

## Usage

``` r
compute_features.training.ML(
  features_train,
  task_type = c("classification", "survival"),
  target_var = NULL,
  trait.positive = NULL,
  time_var = NULL,
  event_var = NULL,
  metric = NULL,
  stack = FALSE,
  k_folds = 10,
  n_rep = 5,
  LODO = FALSE,
  batch_var = NULL,
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

  A data frame with samples in rows and features in columns.

- task_type:

  Character. Prediction task type: `"classification"` or `"survival"`.

- target_var:

  Vector. Target variable for classification tasks.

- trait.positive:

  Value in `target_var` representing the positive class.

- time_var:

  Character. Name of the survival time variable (required for survival
  tasks).

- event_var:

  Character. Name of the event indicator (1 = event occurred, 0 =
  censored) for survival tasks.

- metric:

  Character. Performance metric for model selection and tuning.
  Supported values:

  - `"Accuracy"` — classification accuracy

  - `"AUROC"` — area under the ROC curve

  - `"AUPRC"` — area under the precision-recall curve

  - `"C-index"` — concordance index (for survival tasks)

- stack:

  Logical. Perform model stacking (ensemble meta-learning). Default:
  `FALSE`.

- k_folds:

  Integer. Number of folds for cross-validation. Default: 10.

- n_rep:

  Integer. Number of repetitions for repeated CV. Default: 5.

- LODO:

  Logical. If `TRUE`, constructs folds stratified by cohort (LODO
  scheme).

- batch_var:

  Character. Batch membership for each sample. Required if
  `LODO = TRUE`.

- file_name:

  Character. File name prefix used to save performance plots in
  `"Results/"`.

- ncores:

  Integer. Number of CPU cores for parallelization. Default:
  `parallel::detectCores() - 1`.

- return:

  Logical. Whether to return the trained models and plots. Default:
  `FALSE`.

- fold_construction_fun:

  Function. Optional user-defined function for fold construction. Must
  accept a `bestune` argument:

  - `bestune = NULL` — explore parameter grid across folds (parallelized
    via `foreach`).

  - `bestune provided` — rebuild features on the full dataset using
    optimized parameters.

  The function should save individual folds as `"Results/fold_*.rds"`
  with:

  - `train_data` — training data

  - `test_data` — testing data

  - `obs_test` — observed outcomes

  - `params` — parameters used (if applicable)

- fold_construction_args_fixed:

  List of arguments passed to `fold_construction_fun` that remain fixed
  across CV and final training.

- fold_construction_args_tunable:

  List of arguments passed to `fold_construction_fun` for hyperparameter
  tuning.

## Value

A list containing:

- Trained model(s) or meta-learner (if `stack = TRUE`)

- Features used for training

- Cross-validation performance results and plots

- Best hyperparameter configuration (if applicable)

## Details

The function supports both classification and survival analysis
pipelines via `task_type = "classification"` or
`task_type = "survival"`.

The function provides:

- Automatic feature preprocessing (e.g., correlation filtering,
  low-variance removal)

- Parallelized cross-validation across folds and repetitions

- Integration with custom model pipelines (e.g., CellTFusion,
  pathway-based deconvolution)

- Unified handling of both survival and classification models

When a custom fold constructor is provided, default k-fold logic is
bypassed, and results are computed using the pre-generated folds.
