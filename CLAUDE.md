# pipeML — Developer Context for Claude

## Package Overview

**pipeML** is an R package providing a modular, leakage-free machine
learning framework with fold-aware feature recomputation. Its core
innovation is custom cross-validation fold construction that allows
features depending on dataset context (enrichment scores, correlations,
network features) to be independently recomputed within each CV fold,
preventing information leakage.

- **Version:** 0.0.1 (active development)
- **License:** GPL \>= 3
- **R requirement:** \>= 4.3
- **GitHub:** <https://github.com/VeraPancaldiLab/pipeML>
- **Website:** <https://verapancaldilab.github.io/pipeML>
- **Authors:** Marcelo Hurtado (aut, cre), Vera Pancaldi (aut)
- **Target domain:** Biomedical/omics high-dimensional machine learning

------------------------------------------------------------------------

## Repository Structure

    R/
      pipeML-package.R       # Package metadata & namespace declarations
      data.R                 # Documentation for bundled example datasets
      machine_learning.R     # All implementation (~5,500 lines, core file)
    vignettes/
      pipeML.Rmd             # "Get started" landing page (intro, install, workflow chooser, citation)
      a1_classification.Rmd  # Articles, one per topic (same layout as multideconv)
      a2_survival.Rmd
      a3_shap.Rmd
      a4_lodo.Rmd
      a5_custom_folds.Rmd
      figures/               # Static figures shown by the articles (code chunks are eval = FALSE)
    data/                    # Bundled example datasets (.rda)
    man/                     # Auto-generated Roxygen docs (never edit manually)
    docs/                    # pkgdown website output (gitignored)
    inst/
      create_hex_logo.R      # Hex sticker generation
    .github/workflows/
      R-CMD-check.yaml       # CI: R package checks across OS/R versions
      pkgdown.yaml           # CI: Build and deploy pkgdown site
    DESCRIPTION              # Package metadata & dependencies
    NAMESPACE                # Exported functions (managed by Roxygen)
    _pkgdown.yml             # pkgdown website configuration
    pipeML.Rproj             # RStudio project config

------------------------------------------------------------------------

## Exported Functions (6 total)

All live in `R/machine_learning.R`.

| Function | Purpose |
|----|----|
| [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md) | Train models on training data with repeated k-fold CV |
| [`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md) | Combined train + predict workflow (training + testing) |
| [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md) | Generate predictions on test data using trained model (returns AUROC/AUPRC with bootstrap CIs and `Curve_bands`) |
| [`get_curves()`](https://verapancaldilab.github.io/pipeML/reference/get_curves.md) | ROC and Precision-Recall curves with pointwise bootstrap confidence bands |
| [`compute_shap_values()`](https://verapancaldilab.github.io/pipeML/reference/compute_shap_values.md) | SHAP values of the final model (trained on all training samples) for every training sample; takes only the trained model |
| [`plot_survival_performance()`](https://verapancaldilab.github.io/pipeML/reference/plot_survival_performance.md) | Kaplan-Meier curves stratified by predicted risk groups |

------------------------------------------------------------------------

## Supported ML Algorithms

**Classification (11, via caret):** `treebag`, `rf`, `C5.0`, `glmnet`,
`lasso` and `ridge` (both `glmnet` with fixed `alpha`), `knn`, `rpart`,
`svmRadial`, `svmLinear`, `xgbTree`

**Survival (6 active, via tidymodels/parsnip/censored; `model_list` in
[`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)):**
`cox_ph_survival` (Cox PH), `proportional_hazards_glmnet` (elastic-net
Cox), `survreg_flexsurv` (parametric AFT), `decision_tree_partykit`
(conditional inference tree), `bag_tree_rpart` (bagged CART),
`rand_forest_aorsf` (oblique random survival forest).
`rand_forest_partykit` and `boost_tree_mboost` are implemented but
commented out of `model_list`.

------------------------------------------------------------------------

## Key Dependencies

**Imports (must be installed):** `caret`, `doParallel`, `foreach`,
`dplyr`, `tidyr`, `tibble`, `purrr (>= 1.0.2)`, `ggplot2`, `reshape2`,
`survival`, `survminer`, `fastshap`, `dials`, `parsnip`, `rsample`,
`workflows`, `tune`, `yardstick`, `grDevices`, `parallel`, `stats`

**Suggests (optional, needed for specific algorithms):**
`testthat (>= 3.0.0)`, `knitr`, `rmarkdown`, `C50`, `randomForest`,
`glmnet`, `xgboost`, `kernlab`, `recipes`, `tidyverse`, `tidymodels`,
`censored`, `flexsurv`, `coin`, `aorsf`, `WGCNA`, `cowplot`, `matlib`,
`shapviz`

**Remotes (GitHub):** - `VeraPancaldiLab/multideconv` — custom
deconvolution package (optional) - `bgreenwell/fastshap@v0.3.0` —
`fastshap` was archived from CRAN on 2026-05-27, so CI (pak) can only
install it from GitHub. Pinned to the v0.3.0 release tag; its
`explain(object, X, pred_wrapper, newdata, nsim, adjust)` API and
default estimator are unchanged from CRAN 0.1.1.

------------------------------------------------------------------------

## Architecture & Key Design Patterns

### Leakage-Free CV

The key innovation:
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
and
[`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md)
accept custom fold construction functions. These functions receive a
`bestune` argument, allowing feature pipelines to recompute features
independently per fold — then be called again on the full training set
with best hyperparameters.

Contract of a fold construction function: - It receives `data` =
training features **plus the outcome**: a `target` column
(`"no"`/`"yes"`) for classification, `time` and `event` columns for
survival. It must remove them before building features and add them back
to what it returns. - Fold mode (`bestune = NULL`): saves
`Results/fold_<fold name>.rds` with `train_data` (features + outcome
columns), `test_data` (features; plus `time`/`event` for survival),
`obs_test` (classification), `rowIndex`, `fold_name`, and `params` when
there are tunable arguments. - Final mode (`bestune` given): returns
`list(features + outcome columns, custom output, bestune or the selected parameters)`. -
Every `Results/fold_*.rds` file is read as a resample, so both CV
functions delete leftover ones (from an interrupted run) before calling
the fold function, and delete the new ones after reading them. - Before
2026-09-29 survival fold functions received the features only and had to
get time/event through `fold_construction_args_fixed`. Nothing checks
that the returned features exclude the outcome:
`check_survival_fold_features()` (which stopped on a missing
`time`/`event` or a feature identical to them) was removed on 2026-10-03
at the user’s request; the fold function must follow the contract.

### Cross-Validation Strategies

- Default: Repeated stratified k-fold
- Multi-cohort: Leave-One-Dataset-Out (LODO) via
  [`construct_stratified_cohort_folds()`](https://verapancaldilab.github.io/pipeML/reference/construct_stratified_cohort_folds.md)

### Parallelization

- `doParallel` + `foreach %dopar%` for cross-fold parallelization
- XGBoost uses internal threading — external parallel is disabled to
  avoid contention
- Classification with `fold_construction_fun` trains the models
  sequentially: `ncores` is not used there (parallelism belongs to the
  fold function). The survival custom-fold branch with tunable arguments
  runs the parameter combinations of each fold in parallel.
- Survival custom-fold branch with tunable arguments: since 2026-10-03
  it runs sequentially (`%do%`) when `ncores` is `NULL` or 1, and
  creates one cluster for all folds when `ncores > 1`. Before, it always
  called `makeCluster(ncores)` inside the fold loop, so the default
  `ncores = NULL` made a 0-worker cluster and stopped with “subscript
  out of bounds”.
- Clusters are released also on error:
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html) after
  `makeCluster()` in
  [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
  and in the survival custom-fold branch, `tryCatch(finally = )` in the
  standard survival CV.

### Hyperparameter Tuning

- Metric-based optimization: AUROC or AUPRC (classification), C-index
  (survival)
- Grid search within CV folds → best params applied to full training
  data
- Standard classification path: caret tunes on Accuracy and fits
  `finalModel` with those values;
  [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
  re-tunes by `metric`
  ([`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md))
  and retrains every model on all training samples with that `bestTune`
  ([`caret_train()`](https://verapancaldilab.github.io/pipeML/reference/caret_train.md)
  with `trainControl(method = "none")`, as the custom-fold path does),
  re-attaching `$results`, `$pred`, `$resample` and `$bestTune`. Before
  2026-10-01 only `$bestTune` was replaced, so the returned model was
  still the Accuracy-tuned one.
- Custom-fold classification path:
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)
  selects `bestTune` on Accuracy; each branch must copy `res$bestTune`
  from
  [`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md)
  into the model before the final `train(method = "none")`. The branch
  without tunable arguments did not do this before 2026-10-01.
- treebag has no hyperparameters: its entry in the `hyperparams` lists
  is `"parameter"` (the column caret fills with `"none"`), not `NULL`.
  With `NULL`,
  [`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md)
  takes a branch that gives every resample the same AUROC (MAD 0).
- Custom-fold path: the `rf` `mtry` grid is sized from `grid_data`, the
  fold training table with the fewest features after preprocessing, read
  from the fold files before the loop. It must be the same in every fold
  (the sanity check requires each setting in all resamples), so it
  cannot be sized per fold; the raw `train_data` has an unrelated number
  of columns.
- Grids: the custom-fold classification path uses
  [`get_tune_grid()`](https://verapancaldilab.github.io/pipeML/reference/get_tune_grid.md);
  the standard path uses caret’s default grids (except lasso/ridge,
  which share the `lambda` values), so the two paths do not tune over
  the same values. Since 2026-10-03 `lambda` is log-spaced
  (`10^seq(-3, 0, length = 20)`, before `seq(0.001, 1, length = 20)`,
  where 18 of 20 values were \>= 0.05) and the custom `svmRadial`
  `sigma` is estimated from `grid_data` (quantiles of 1/\|x - x’\|^2
  over all sample pairs on scaled features, deterministic so the grid is
  the same in every fold; before fixed at 0.01/0.05/0.1, too large with
  many features).
- Survival grids
  ([`get_default_hyperparams()`](https://verapancaldilab.github.io/pipeML/reference/get_default_hyperparams.md)):
  only arguments that parsnip passes to the engine are tuned. Since
  2026-10-03 the partykit tree has no `cost_complexity` (dropped by
  parsnip: 5 identical trees per configuration before) and the
  (inactive) mboost grid has only `trees`, `min_n`, `tree_depth` (was
  5^7 = 78,125 configurations, with 3 dropped arguments and
  `loss_reduction` outside mincriterion’s 0-1 range). `bag_tree_rpart`
  tunes `trees` (number of bags, 25 to 100): `bag_tree()` has no such
  argument and `set_args(trees = )` is silently ignored
  ([`ipred::bagging()`](https://rdrr.io/pkg/ipred/man/bagging.html)
  built 25 trees for every value before 2026-10-03), so
  [`compute_ml_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_ml_survival.md)
  and
  [`wrapper_train_best_hyperparams_survival()`](https://verapancaldilab.github.io/pipeML/reference/wrapper_train_best_hyperparams_survival.md)
  pass it with `set_engine("rpart", nbagg = !!n_bags)`.
- [`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md)
  sorts `Prediction_folds` by resample, hyperparameters and decreasing
  `yes` before `bind_cols(metrics)`, because `metrics` comes out in that
  order. `$pred` is therefore sorted that way, not in caret’s order.

### Feature Preprocessing (internal `preprocess_features()`, `preprocess` argument, default `TRUE`)

- Takes the features only (no outcome column) and returns the kept
  features; every caller sets the outcome columns aside (`target`, or
  `time` and `event`) and adds them back with
  [`cbind()`](https://rdrr.io/r/base/cbind.html). Before 2026-10-03 it
  took the outcome too (`target_col`, `time_var`, `event_var`).
- Near-zero variance removal
- Collinearity filtering (correlation threshold)
- Where it runs (classification and survival), all switched off by
  `preprocess = FALSE`:
  - Standard path: once on all training samples before the CV
    ([`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
    after the `dataset` column is dropped;
    [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)
    before the hyperparameter grids, so `mtry` is sized from the kept
    features). Added on 2026-10-02 as the user’s decision: the filter is
    unsupervised, but the held-out fold’s feature values take part in
    it. Before, the standard path had no filter.
  - Custom folds: on each fold’s training table (per parameter
    combination), the held-out table is cut to the kept columns; and
    once on the final training table of all samples (classification: in
    [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md),
    after re-tuning by `metric`; survival: in
    [`wrapper_train_best_hyperparams_survival()`](https://verapancaldilab.github.io/pipeML/reference/wrapper_train_best_hyperparams_survival.md)).
    Both custom paths (with and without tunable arguments) build,
    preprocess and train the final table once per model. Before
    2026-10-03 the classification path with tunable arguments also
    called `wrapper_train_best_hyperparams_classification()` with the
    Accuracy `bestTune` (fold function on all samples, preprocessing and
    training), whose results were all discarded; that function is now
    commented out in `machine_learning.R` (kept for reference, no help
    page).
  - [`wrapper_train_best_hyperparams_survival()`](https://verapancaldilab.github.io/pipeML/reference/wrapper_train_best_hyperparams_survival.md)
    returns `NULL` with a warning if the final fit fails, so that model
    is excluded
    ([`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)
    stops with a clear message only if no model is left, in both paths).
    Its `Model$bestTune` holds the model hyperparameters only, as in
    classification; the feature parameters are in
    `custom_output$Parameters`.
  - Custom folds: `grid_data` is taken from the fold tables after
    preprocessing, since the models are trained on those. Classification
    sizes the `rf` grid from it; survival (since 2026-10-03) re-sizes
    the `mtry` of the forests (`rand_forest_aorsf`,
    `rand_forest_partykit`) from it. Before, the survival custom path
    sized `mtry` from the input columns; parsnip (`min_cols()`) then
    reset every `mtry` above the number of predictors to that number
    with a warning, so several grid values were the same model.
  - The test set is never involved:
    [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
    keeps the features of the final model.
- The outcome is not used: the step removing features with near-zero
  variance within a class (classification) was deleted on 2026-10-02,
  because it also removed features that separate the classes (e.g. a
  marker absent in one class)
- Stops if a feature is not numeric or if no feature is left

### Reproducibility (`seed` argument, default `123`)

- [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md),
  [`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md)
  and
  [`compute_shap_values()`](https://verapancaldilab.github.io/pipeML/reference/compute_shap_values.md)
  take `seed`; `NULL` leaves the RNG untouched.
- The training functions pass it to
  [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
  /
  [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md),
  which call `set.seed(seed)` right before drawing folds. caret draws
  per-resample seeds for its parallel workers from the main session’s
  RNG, so this also covers `ncores > 1`.
- `%dopar%` loops don’t share the main RNG, so each iteration seeds
  itself: `seed + fold_i` in the standard survival CV (folds run in
  parallel with `ncores > 1`, sequentially otherwise) and
  `seed + 1000L * fold_i + parameter_i` in the survival custom-fold
  branch. Results are therefore independent of worker scheduling and
  `ncores`.
  [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)
  also calls `set.seed(seed)` right before the final training: without
  `ncores` the per-fold seeds are set in the main session, with `ncores`
  in the workers, so before 2026-10-03 the final survival model
  (forests, bagged trees) depended on `ncores` while the CV results did
  not.
- [`compute_shap_values()`](https://verapancaldilab.github.io/pipeML/reference/compute_shap_values.md)
  passes `seed` to
  [`fastshap::explain()`](https://bgreenwell.github.io/fastshap/reference/explain.html)
  and runs sequentially: fastshap’s `parallel = TRUE` does not seed its
  workers, so its results change between runs.
- Randomness inside a user-supplied `fold_construction_fun` that runs
  its own parallel workers is not covered.

### SHAP

- `compute_shap_values(model_trained, task_type, seed)` explains the
  **final model** (trained on all training samples with `bestTune`) on
  all its training samples, which are also the background data.
  Classification: the caret `train` object and its `$trainingData`;
  survival: `$Model_object` and `$trainingData`.
- `$trainingData` always holds the final model’s features, so each
  feature has one fixed definition (custom-fold features are computed
  once on all training samples). For survival,
  [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md)
  stores the custom `training_set` on the custom path and `df_all`
  (minus `strata`) on the standard path.
- The output data frame (samples x features) carries the model’s average
  prediction in `attr(, "baseline")`; baseline + row sums = the sample’s
  prediction.
- Fold models are not saved during CV: SHAP doesn’t use them.

### Performance Curves (`compute_prediction()` → `get_curves()`)

- AUROC/AUPRC are stored as `list(estimate, lower, upper)`: `estimate`
  is the value on the full test set, and the CI comes from 1000
  bootstrap resamples
  ([`bootstrap_auc()`](https://verapancaldilab.github.io/pipeML/reference/bootstrap_auc.md),
  seed 123). Since 2026-10-03 both bootstraps are stratified:
  [`bootstrap_auc()`](https://verapancaldilab.github.io/pipeML/reference/bootstrap_auc.md)
  resamples positives and negatives separately (before, a resample with
  one class gave NaN and
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) stopped
  [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md),
  e.g. 20 test samples with 2 positives), and
  [`compute_cindex_ci()`](https://verapancaldilab.github.io/pipeML/reference/compute_cindex_ci.md)
  resamples events and censored samples separately (before, resamples
  without events were silently dropped by `na.rm = TRUE`).
- The shaded bands are pointwise 95% bootstrap bands.
  [`bootstrap_auc()`](https://verapancaldilab.github.io/pipeML/reference/bootstrap_auc.md)
  evaluates each resample’s curve on a 101-point grid with straight
  lines between consecutive points
  ([`roc_at_grid()`](https://verapancaldilab.github.io/pipeML/reference/roc_at_grid.md),
  [`prc_at_grid()`](https://verapancaldilab.github.io/pipeML/reference/prc_at_grid.md)),
  as
  [`get_curves()`](https://verapancaldilab.github.io/pipeML/reference/get_curves.md)
  draws the curves (`geom_line()`) and as AUROC/AUPRC are integrated
  (trapezoids), returned as `Curve_bands`. Before 2026-10-03 the ROC was
  read as a staircase (best sensitivity at FPR \<= x: band up to ~0.16
  below the drawn curve on diagonal segments from tied probabilities)
  and the PR curve at the next point (gap up to ~0.13 even without
  ties). `get_curves(roc_band, prc_band)` draws them under the curve;
  the LODO branch draws no band.
- Tied probabilities: the curves are built from cumulative TP/FP counts
  over samples sorted by decreasing probability
  ([`calculate_auc_roc_resample()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auc_roc_resample.md),
  [`calculate_auc_prc_resample()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auc_prc_resample.md),
  [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md)).
  Samples with the same probability get the counts at the end of their
  tie group (`stats::ave(tp, yes, FUN = max)`); without this,
  AUROC/AUPRC depend on the row order (a constant predictor scored 1 or
  0 with rows sorted by class).
  [`calculate_auroc()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auroc.md)
  starts the curve at (0, 0) and
  [`calculate_auprc()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auprc.md)
  at recall 0 with the precision of the first point, so a perfect model
  has AUPRC 1 and a constant one has AUPRC = prevalence.

### Outcome Alignment (`compute_features.ML()`)

- Labels are taken from `coldata` by row name
  (`coldata[rownames(features_train), , drop = FALSE]`), so they follow
  the feature rows’ order; missing samples are an error. Never subset
  with `%in%` and then attach columns by position — that silently
  scrambled labels whenever feature rows weren’t in `coldata` order.

------------------------------------------------------------------------

## Bundled Example Datasets

| Dataset | Content |
|----|----|
| `data_example_classification` | Breast Cancer Wisconsin (from mlbench) |
| `data_example_survival` | Lung cancer survival (from survival package) |
| `counts_example` | Gene expression matrix — Gide et al. 2019 melanoma cohort (4.4 MB) |
| `coldata_example` | Metadata with anti-PD-1 therapy response labels for Gide cohort |

Access via `data(dataset_name)` after
[`library(pipeML)`](https://verapancaldilab.github.io/pipeML).

------------------------------------------------------------------------

## Development Workflow

### Code Style

- 2-space indentation (configured in pipeML.Rproj)
- Roxygen2 for all exported function documentation
- Run `devtools::document()` after any signature or documentation
  changes (regenerates `man/` and `NAMESPACE`)
- Never edit `man/` or `NAMESPACE` directly

### Building & Checking

``` r

devtools::load_all()       # Load package in dev mode
devtools::document()       # Regenerate docs & NAMESPACE
devtools::check()          # Full R CMD check
devtools::build_vignettes() # Build vignettes
pkgdown::build_site()      # Rebuild docs website
```

### CI/CD (GitHub Actions)

- **R-CMD-check:** Runs on push/PR to main — tests macOS, Windows,
  Ubuntu across R devel/release/oldrel-1
- **pkgdown:** Auto-deploys website to `gh-pages` branch on push to main
  or release

### Testing

- `testthat` (v3) suite in `tests/testthat/`: helpers
  ([`preprocess_features()`](https://verapancaldilab.github.io/pipeML/reference/preprocess_features.md),
  hyperparameter grids, example data,
  [`ensure_caret()`](https://verapancaldilab.github.io/pipeML/reference/ensure_caret.md)),
  survival helpers (failed fits, empty `bestTune`, predictions without
  outcome,
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)
  exclusion of failed configurations) and end-to-end workflows
  (training, prediction, SHAP additivity/reproducibility; standard and
  custom folds).
- End-to-end tests are `skip_on_cran()`; the survival training test
  (several minutes) is also `skip_on_ci()`, so run `devtools::test()`
  locally before releasing.
- Tests run in a temporary working directory (`local_temp_wd()` in
  `helper-data.R`), since pipeML writes to `Results/`.
- The vignette code is `eval = FALSE`: run its chunks when changing
  user-facing behaviour. `a3_shap`’s survival section uses
  `res_survival` from `a2_survival`, so run the articles in order in one
  session.

------------------------------------------------------------------------

## Common Tasks

### Adding a new ML algorithm

1.  Add to
    [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
    /
    [`compute_custom_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_custom_k_fold_CV.md)
    in `machine_learning.R`
2.  Add a hyperparameter grid entry in
    [`get_tune_grid()`](https://verapancaldilab.github.io/pipeML/reference/get_tune_grid.md)
3.  Add its final `train(method = "none")` block in the custom path of
    [`compute_k_fold_CV()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV.md)
4.  Document in vignette and function `@param method` Roxygen docs

### Adding a new metric

1.  Add calculation helper (follow
    [`calculate_auroc()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auroc.md)
    pattern)
2.  Wire into
    [`calculate_cv_metrics()`](https://verapancaldilab.github.io/pipeML/reference/calculate_cv_metrics.md)
    and
    [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)

### Modifying the vignettes

The vignettes follow multideconv’s layout: `pipeML.Rmd` is a short
landing page and each topic is an article `aN_<topic>.Rmd`, listed in
both the navbar `articles` menu and the `articles:` index of
`_pkgdown.yml`. A new article must be added in both places. Run
`devtools::build_vignettes()` to test locally. The pkgdown CI will
publish on merge to main.

------------------------------------------------------------------------

## Notes & Gotchas

- The core implementation is a single large file (`machine_learning.R`,
  ~5,500 lines). Internal helpers are not exported — check NAMESPACE
  before assuming a function is public.

- caret must be *attached*, not just loaded: its `"knn"` model code
  calls `knn3()` without a namespace prefix, so a knn fit with
  `trainControl(method = "none")` fails with
  `could not find function "knn3"` otherwise. Classification entry
  points call the internal
  [`ensure_caret()`](https://verapancaldilab.github.io/pipeML/reference/ensure_caret.md).

- Bootstraps with a fixed seed
  ([`bootstrap_auc()`](https://verapancaldilab.github.io/pipeML/reference/bootstrap_auc.md),
  [`compute_cindex_ci()`](https://verapancaldilab.github.io/pipeML/reference/compute_cindex_ci.md),
  seed 123) save and restore the session’s random state
  ([`withr::local_preserve_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  then [`set.seed()`](https://rdrr.io/r/base/Random.html)), so
  [`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
  no longer resets the user’s RNG and survival CV fits no longer restart
  from the same state after each evaluation (before 2026-10-03 they
  called [`set.seed()`](https://rdrr.io/r/base/Random.html) directly).
  [`withr::local_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  is not used: in withr 3.0.2 it does not restore the state. Survival CV
  evaluates folds with `predict_and_evaluate_survival(ci = FALSE)`: the
  bootstrap CI (about 1.5 s per call, ~180 configurations per fold) was
  computed and discarded for every fit.

- Don’t add [`set.seed()`](https://rdrr.io/r/base/Random.html) calls
  inside helpers: they override the user’s `seed` mid-run
  ([`get_tune_grid()`](https://verapancaldilab.github.io/pipeML/reference/get_tune_grid.md)
  used to do this). Seed only at the entry points described under
  Reproducibility.

- Roxygen markdown is enabled (`Roxygen: list(markdown = TRUE)`): write
  `95%`, not `95\%` (it becomes `\\%` in the Rd, which comments out the
  rest of the line), and avoid `[a, b]`-style ranges in plain text
  (parsed as links). Run `devtools::check_man()` after doc edits.

- CRAN doesn’t allow non-CRAN required dependencies: `multideconv` and
  `fastshap` (both in `Remotes`) block a CRAN submission.

- Survival models need the `censored` package (in Suggests) to register
  parsnip’s “censored regression” engines. Every survival entry point
  calls the internal
  [`ensure_censored()`](https://verapancaldilab.github.io/pipeML/reference/ensure_censored.md),
  which loads its namespace; users don’t need
  [`library(censored)`](https://github.com/tidymodels/censored).
  Survival formulas must use `survival::Surv(...)`, not bare `Surv(...)`
  — the bare form only works when some other package happened to attach
  `survival`.

- Direction of survival predictions: parsnip/censored return every
  prediction type as “higher = longer survival”, including `linear_pred`
  of Cox models (`censored:::predict_linear_pred._coxph()` negates the
  log hazard by default, `increasing = TRUE`).
  [`yardstick::concordance_survival_vec()`](https://yardstick.tidymodels.org/reference/concordance_survival.html)
  expects that direction, so
  [`predict_and_evaluate_survival()`](https://verapancaldilab.github.io/pipeML/reference/predict_and_evaluate_survival.md)
  computes the C-index first and only then negates all predictions into
  risk scores (higher = higher risk), which is what `preds`, the
  Kaplan-Meier risk groups and the survival SHAP values use. Before
  2026-10-02 nothing was negated (the `pred_type` set inside the
  `tryCatch` handlers never reached the outer variable), so “High risk”
  labelled the longest survivors. The `pred_type` assignments are still
  there but unused.

- `multideconv` is a remote (GitHub) dependency — not on CRAN.
  Installation requires
  `remotes::install_github("VeraPancaldiLab/multideconv")`.

- SHAP computation via
  [`fastshap::explain()`](https://bgreenwell.github.io/fastshap/reference/explain.html)
  can be memory-intensive on large datasets.

- xgboost \>= 3 breaks caret’s own `"xgbTree"` model (the booster is an
  ALTREP object: `modelFit$xNames <- ...` fails with “ALTLIST classes
  must provide a Set_elt method”). All
  [`caret::train()`](https://rdrr.io/pkg/caret/man/train.html) calls
  that may get `"xgbTree"` go through the internal
  [`caret_train()`](https://verapancaldilab.github.io/pipeML/reference/caret_train.md),
  which substitutes
  [`xgbtree_model()`](https://verapancaldilab.github.io/pipeML/reference/xgbtree_model.md)
  (booster stored in `modelFit$booster`) and resets `$method` to
  `"xgbTree"`. Use
  [`caret_train()`](https://verapancaldilab.github.io/pipeML/reference/caret_train.md)
  for any new call site.

- XGBoost parallel contention: when using `doParallel`, XGBoost nthread
  is set to 1 internally to prevent nested parallelism crashes.

- The `docs/` directory is gitignored — pkgdown output is built and
  deployed by CI only.

- `compute_shap_values(model_trained, ...)` expects `res$Model` from
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  (for classification, the caret `train` object). It validates this and
  stops with a clear error if given the wrong object (e.g. `res` instead
  of `res$Model`).

- Survival MAD: since 2026-10-03
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md)
  (`c_index_mad`) and
  [`compute_cv_CINDEX()`](https://verapancaldilab.github.io/pipeML/reference/compute_cv_CINDEX.md)
  (`MAD_CINDEX`) use the scaled MAD
  ([`stats::mad()`](https://rdrr.io/r/stats/mad.html) default), as
  classification does; before they used `constant = 1`, so the
  same-looking error bars were ~1.5x narrower. The C-index plot uses the
  same helper as the classification plots
  ([`plot_cv_metric()`](https://verapancaldilab.github.io/pipeML/reference/plot_cv_metric.md)).

- Survival CV,
  [`aggregate_results()`](https://verapancaldilab.github.io/pipeML/reference/aggregate_results.md):
  a configuration is excluded if its fit failed or its C-index is `NA`
  (e.g. a held-out fold without events) in at least one resample, with a
  message, so all configurations are summarized over the same resamples.
  Before 2026-10-03 an `NA` C-index was silently dropped
  (`na.rm = TRUE`). Known, left as is: the C-index median/MAD are
  computed over per-sample rows (the fold’s C-index repeated on each
  test sample), which differs from per-resample values with unequal
  folds and an even number of resamples.

- Survival LODO:
  [`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
  /
  [`compute_features.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.ML.md)
  now pass `LODO` and `batch_id = "dataset"` to
  [`compute_k_fold_CV_survival()`](https://verapancaldilab.github.io/pipeML/reference/compute_k_fold_CV_survival.md),
  which builds cohort x event stratified folds and then drops `dataset`
  and `strata`, so the cohort is not a predictor. Before 2026-09-28,
  survival LODO never built cohort folds: text cohort labels crashed and
  numeric labels were silently used as a predictor.

## Known Issues / TODO

### Open bugs

- Cosmetic: on precision-recall plots, the curve’s final vertical drop
  at recall = 1 (to precision = prevalence) falls below the confidence
  band, which ends at the first threshold reaching full recall.

### Performance

- **[`compute_shap_values()`](https://verapancaldilab.github.io/pipeML/reference/compute_shap_values.md)’s
  cost depends on the model’s prediction speed**: it runs one
  [`fastshap::explain()`](https://bgreenwell.github.io/fastshap/reference/explain.html)
  call on the final model (`nsim = 100`, sequential). Measured on 73
  samples x 118 features: `xgbTree` ≈ 17 s, `glmnet` ≈ 15 s, `svmRadial`
  ≈ 75 s. `nsim` is hardcoded; exposing it would let callers trade
  precision for speed with slow-predict models.
