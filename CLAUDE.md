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
For survival,
[`check_survival_fold_features()`](https://verapancaldilab.github.io/pipeML/reference/check_survival_fold_features.md)
stops if the returned features contain a copy of `time` or `event`
(outcome not removed). Before 2026-09-29 survival fold functions
received the features only and had to get time/event through
`fold_construction_args_fixed`.

### Cross-Validation Strategies

- Default: Repeated stratified k-fold
- Multi-cohort: Leave-One-Dataset-Out (LODO) via
  [`construct_stratified_cohort_folds()`](https://verapancaldilab.github.io/pipeML/reference/construct_stratified_cohort_folds.md)

### Parallelization

- `doParallel` + `foreach %dopar%` for cross-fold parallelization
- XGBoost uses internal threading — external parallel is disabled to
  avoid contention

### Hyperparameter Tuning

- Metric-based optimization: AUROC or AUPRC (classification), C-index
  (survival)
- Grid search within CV folds → best params applied to full training
  data

### Feature Preprocessing (internal `preprocess_features()`)

- Near-zero variance removal
- Collinearity filtering (correlation threshold)
- Removal of features constant within any target class (classification
  only)

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
  seed 123).
- The shaded bands are pointwise 95% bootstrap bands.
  [`bootstrap_auc()`](https://verapancaldilab.github.io/pipeML/reference/bootstrap_auc.md)
  evaluates each resample’s curve on a 101-point grid
  ([`roc_at_grid()`](https://verapancaldilab.github.io/pipeML/reference/roc_at_grid.md):
  best sensitivity at FPR ≤ x;
  [`prc_at_grid()`](https://verapancaldilab.github.io/pipeML/reference/prc_at_grid.md):
  precision at the first threshold reaching recall x), returned as
  `Curve_bands`. `get_curves(roc_band, prc_band)` draws them under the
  curve; the LODO branch draws no band.

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
3.  Handle retraining in
    [`wrapper_train_best_hyperparams_classification()`](https://verapancaldilab.github.io/pipeML/reference/wrapper_train_best_hyperparams_classification.md)
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
