# Train machine learning or survival models with custom cross-validation

This function trains one or more machine learning models using repeated
k-fold cross-validation, with optional feature selection, and support
for both classification and survival tasks. It allows flexible
cross-validation schemes, including:

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
  k_folds = 10,
  n_rep = 5,
  LODO = FALSE,
  batch_var = NULL,
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

  A data frame with samples in rows and features in columns.

- task_type:

  Character. Prediction task type: `"classification"` or `"survival"`.

- target_var:

  Vector. Target variable for classification tasks.

- trait.positive:

  Value in `target_var` representing the positive class.

- time_var:

  Numeric vector. Survival/follow-up time of each sample (required for
  survival tasks).

- event_var:

  Numeric vector. Event indicator of each sample (1 = event occurred, 0
  = censored; required for survival tasks).

- metric:

  Character. Performance metric for model selection and tuning
  (classification):

  - `"AUROC"` - area under the ROC curve (used when `NULL`)

  - `"AUPRC"` - area under the precision-recall curve

  Survival models are always selected by the concordance index
  (C-index); `metric` is ignored.

- k_folds:

  Integer. Number of folds for cross-validation. Default: 10.

- n_rep:

  Integer. Number of repetitions for repeated CV. Default: 5.

- LODO:

  Logical. If `TRUE`, the cross-validation folds are stratified by
  cohort and outcome, for Leave-One-Dataset-Out analyses (train on some
  cohorts, test on a left-out one). The cohort is not used as a
  predictor.

- batch_var:

  Vector. Cohort/batch of each sample. Required if `LODO = TRUE`.

- file_name:

  Character. File name prefix used to save performance plots in
  `"Results/"`.

- ncores:

  Integer. Number of CPU cores for parallelization (cross-validation
  folds are processed in parallel). Default: `NULL` (sequential).

- return:

  Logical. Whether to save the cross-validation performance plots in
  `"Results/"`. Default: `FALSE`.

- fold_construction_fun:

  Function. Optional user-defined function for fold construction. It is
  called with `data` (the training features plus the outcome: a `target`
  column coded `"no"`/`"yes"` for classification, `time` and `event`
  columns for survival), `folds` and `bestune`. It must remove the
  outcome columns before building features. Must accept a `bestune`
  argument:

  - `bestune = NULL` - explore parameter grid across folds (parallelized
    via `foreach`).

  - `bestune provided` - rebuild features on the full dataset using
    optimized parameters, and return a list with the features plus the
    outcome columns, any custom output, and `bestune`.

  The function should save individual folds as `"Results/fold_*.rds"`
  with:

  - `train_data` - training features plus the outcome columns

  - `test_data` - test features (plus `time` and `event` for survival)

  - `obs_test` - observed outcomes (classification)

  - `params` - parameters used (if applicable)

- fold_construction_args_fixed:

  List of arguments passed to `fold_construction_fun` that remain fixed
  across CV and final training.

- fold_construction_args_tunable:

  List of arguments passed to `fold_construction_fun` for hyperparameter
  tuning.

- seed:

  Integer. Random seed for reproducible cross-validation: it fixes the
  fold assignment and the randomness in model fitting (e.g. random
  forest, bagging, boosting), including when running in parallel with
  `ncores`. Default: `123`. Use `NULL` to leave the random number
  generator untouched. Randomness inside a user-supplied
  `fold_construction_fun` that runs its own parallel workers is not
  covered: seed those workers inside that function.

## Value

A named list, or `NULL` if no model could be trained:

- Model:

  The selected model, trained on all training samples. Classification: a
  caret `train` object (tuned hyperparameters in `$bestTune`,
  performance per resample in `$resample`). Survival: a list with the
  model name (`$model`), the fitted workflow (`$Model_object`), the
  tuned hyperparameters (`$bestTune`), the C-index per resample
  (`$Resample_matrix`) and the training data (`$trainingData`).

- ML_Models:

  All trained models.

- AUROC_median, AUPRC_median:

  Classification: median and MAD of the cross-validation AUROC and AUPRC
  of each model.

- C_index_median:

  Survival: median cross-validation C-index of the selected model.

- Custom_output:

  Only with `fold_construction_fun`: the custom output returned by that
  function on all training samples, and the selected feature parameters
  in `$Parameters`.

## Details

The function supports both classification and survival analysis
pipelines via `task_type = "classification"` or
`task_type = "survival"`.

Classification trains and tunes 11 algorithms with caret (`treebag`,
`rf`, `C5.0`, `glmnet`, lasso, ridge, `knn`, `rpart`, `svmRadial`,
`svmLinear`, `xgbTree`). Survival trains and tunes 6 models with parsnip
and censored (Cox PH, elastic-net Cox, parametric AFT, conditional
inference tree, bagged CART, oblique random survival forest). The best
model is selected by the cross-validation `metric` (classification) or
C-index (survival), and trained on all training samples with its tuned
hyperparameters.

When `fold_construction_fun` is provided, the features of each fold are
built by that function, and near-constant, highly correlated (\|r\| \>
0.9) and, for classification, class-constant features are removed from
the features built on each training part. The final model is trained on
the features the function builds on all training samples. See
[`vignette("a5_custom_folds", package = "pipeML")`](https://verapancaldilab.github.io/pipeML/articles/a5_custom_folds.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Classification ---
data(data_example_classification)
X <- data_example_classification[, setdiff(colnames(data_example_classification), "target")]
res <- compute_features.training.ML(features_train = X,
                                    target_var = data_example_classification$target,
                                    task_type = "classification",
                                    trait.positive = "1",
                                    metric = "AUROC",
                                    k_folds = 5,
                                    n_rep = 2,
                                    ncores = 2)
res$Model
res$AUROC_median

# --- Survival ---
data(data_example_survival)
X <- data_example_survival[, setdiff(colnames(data_example_survival), c("time", "status"))]
res_survival <- compute_features.training.ML(features_train = X,
                                             task_type = "survival",
                                             time_var = data_example_survival$time,
                                             event_var = data_example_survival$status,
                                             k_folds = 5,
                                             n_rep = 2,
                                             ncores = 2)
res_survival$Model$model
res_survival$C_index_median
} # }
```
