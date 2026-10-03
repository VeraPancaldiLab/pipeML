# Train machine learning or survival models with custom cross-validation

This function trains several machine learning models using repeated
k-fold cross-validation, with hyperparameter tuning, and supports both
classification and survival tasks. It allows flexible cross-validation
schemes, including:

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
  seed = 123,
  preprocess = TRUE
)
```

## Arguments

- features_train:

  A data frame with samples in rows and features in columns.

- task_type:

  Character. Prediction task type: `"classification"` or `"survival"`.

- target_var:

  Vector. Target variable for classification tasks: one value per
  sample, in the same order as the rows of `features_train`, without
  missing values.

- trait.positive:

  Value in `target_var` representing the positive class. All other
  values form the negative class.

- time_var:

  Numeric vector. Survival/follow-up time of each sample, in the same
  order as the rows of `features_train` (required for survival tasks).

- event_var:

  Vector. Event indicator of each sample (1 = event occurred, 0 =
  censored), in the same order as the rows of `features_train` (required
  for survival tasks).

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

  Vector. Cohort/batch of each sample, in the same order as the rows of
  `features_train`. Required if `LODO = TRUE`. For classification, each
  cohort needs at least `k_folds` samples.

- file_name:

  Character. File name prefix used to save performance plots in
  `"Results/"`.

- ncores:

  Integer. Number of CPU cores for parallelization (cross-validation
  folds are processed in parallel). Default: `NULL` (sequential). For
  classification with `fold_construction_fun`, the models are trained
  sequentially and `ncores` is not used.

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

  - `bestune = NULL` - build the features of every fold (for every
    combination of the tunable arguments, if any) and save them. The
    function can run its own parallel workers for this.

  - `bestune provided` - rebuild features on the full dataset using
    optimized parameters, and return a list with the features plus the
    outcome columns, any custom output, and `bestune` (with tunable
    arguments: a data frame with the selected values of these
    arguments).

  The function should save each fold as
  `"Results/fold_<fold name>.rds"`, a list with:

  - `train_data` - training features plus the outcome columns

  - `test_data` - test features (plus `time` and `event` for survival)

  - `obs_test` - observed outcomes of the test samples, in the order of
    the rows of `test_data` (classification)

  - `rowIndex` - row indices of the test samples in `data`

  - `fold_name` - name of the fold

  - `params` - data frame with the values of the tunable arguments used
    (only with tunable arguments)

  With tunable arguments, the file of a fold contains a list with one
  such element per combination.

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

- preprocess:

  Logical. If `TRUE` (default), near-constant and highly correlated
  (\|r\| \> 0.9) features are removed before training (see Details). The
  features must then be numeric. Use `FALSE` to train on the features as
  given.

## Value

A named list:

- Model:

  The selected model, trained on all training samples. Classification: a
  caret `train` object (tuned hyperparameters in `$bestTune`,
  performance per resample in `$resample`). Survival: a list with the
  fitted workflow (`$Model_object`), the tuned hyperparameters
  (`$bestTune`), the C-index per resample (`$Resample_matrix`, with the
  name of the model in its column `model`) and the training data
  (`$trainingData`).

- ML_Models:

  All trained models. Classification: models that predict the same value
  for all training samples are excluded.

- AUROC_median, AUPRC_median:

  Classification: median and MAD of the cross-validation AUROC and AUPRC
  of each model.

- C_index_median:

  Survival: median cross-validation C-index of the selected model.

- Custom_output:

  Only with `fold_construction_fun`: the custom output returned by that
  function on all training samples and, with tunable arguments, the
  selected feature parameters in `$Parameters`.

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

With `preprocess = TRUE`, near-constant features and, for each pair of
features with an absolute correlation above 0.9, one of the two are
removed. The outcome is not used. Without `fold_construction_fun`, this
is done once on all training samples before the cross-validation, so the
folds and the final model use the same features. The test samples are
not involved:
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
keeps the features of the final model.

When `fold_construction_fun` is provided, the features of each fold are
built by that function, and the preprocessing is done inside each fold:
on the features built on the training part, keeping the same features in
the held-out part. The final model is trained on the features the
function builds on all training samples, preprocessed in the same way.
See
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
unique(res_survival$Model$Resample_matrix$model) # selected model
res_survival$C_index_median
} # }
```
