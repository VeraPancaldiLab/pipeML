# pipeML — Developer Context for Claude

## Package Overview

**pipeML** is an R package providing a modular, leakage-free machine learning framework with fold-aware feature recomputation. Its core innovation is custom cross-validation fold construction that allows features depending on dataset context (enrichment scores, correlations, network features) to be independently recomputed within each CV fold, preventing information leakage.

- **Version:** 0.0.1 (active development)
- **License:** GPL >= 3
- **R requirement:** >= 4.3
- **GitHub:** https://github.com/VeraPancaldiLab/pipeML
- **Website:** https://verapancaldilab.github.io/pipeML
- **Authors:** Marcelo Hurtado (aut, cre), Vera Pancaldi (aut)
- **Target domain:** Biomedical/omics high-dimensional machine learning

---

## Repository Structure

```
R/
  pipeML-package.R       # Package metadata & namespace declarations
  data.R                 # Documentation for bundled example datasets
  machine_learning.R     # All implementation (~5,500 lines, core file)
vignettes/
  pipeML.Rmd             # Main tutorial vignette
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
```

---

## Exported Functions (7 total)

All live in `R/machine_learning.R`.

| Function | Purpose |
|---|---|
| `compute_features.training.ML()` | Train models on training data with repeated k-fold CV |
| `compute_features.ML()` | Combined train + predict workflow (training + testing) |
| `compute_prediction()` | Generate predictions on test data using trained model (returns AUROC/AUPRC with bootstrap CIs and `Curve_bands`) |
| `get_curves()` | ROC and Precision-Recall curves with pointwise bootstrap confidence bands |
| `compute_shap_values()` | Per-fold SHAP values (each sample explained only by the fold models that held it out); takes only the trained model |
| `plot_shap_stability()` | Visualize SHAP importance stability across resamples (input: `compute_shap_values(..., return_resamples = TRUE)$shap_resamples`) |
| `plot_survival_performance()` | Kaplan-Meier curves stratified by predicted risk groups |

---

## Supported ML Algorithms

**Classification (11, via caret):**
`treebag`, `rf`, `C5.0`, `glmnet`, `lasso` and `ridge` (both `glmnet` with fixed `alpha`), `knn`, `rpart`, `svmRadial`, `svmLinear`, `xgbTree`

**Survival (6 active, via tidymodels/parsnip/censored; `model_list` in `compute_k_fold_CV_survival()`):**
`cox_ph_survival` (Cox PH), `proportional_hazards_glmnet` (elastic-net Cox), `survreg_flexsurv` (parametric AFT), `decision_tree_partykit` (conditional inference tree), `bag_tree_rpart` (bagged CART), `rand_forest_aorsf` (oblique random survival forest). `rand_forest_partykit` and `boost_tree_mboost` are implemented but commented out of `model_list`.

---

## Key Dependencies

**Imports (must be installed):**
`caret`, `doParallel`, `foreach`, `dplyr`, `tidyr`, `tibble`, `purrr (>= 1.0.2)`, `ggplot2`, `reshape2`, `survival`, `survminer`, `fastshap`, `dials`, `parsnip`, `rsample`, `workflows`, `tune`, `yardstick`, `grDevices`, `parallel`, `stats`

**Suggests (optional, needed for specific algorithms):**
`testthat (>= 3.0.0)`, `knitr`, `rmarkdown`, `C50`, `randomForest`, `glmnet`, `xgboost`, `kernlab`, `recipes`, `tidyverse`, `tidymodels`, `censored`, `flexsurv`, `coin`, `aorsf`, `WGCNA`, `cowplot`, `matlib`, `shapviz`

**Remotes (GitHub):**
- `VeraPancaldiLab/multideconv` — custom deconvolution package (optional)
- `bgreenwell/fastshap@v0.3.0` — `fastshap` was archived from CRAN on 2026-05-27, so CI (pak) can only install it from GitHub. Pinned to the v0.3.0 release tag; its `explain(object, X, pred_wrapper, newdata, nsim, adjust)` API and default estimator are unchanged from CRAN 0.1.1.

---

## Architecture & Key Design Patterns

### Leakage-Free CV
The key innovation: `compute_features.training.ML()` and `compute_features.ML()` accept custom fold construction functions. These functions receive a `bestune` argument, allowing feature pipelines to recompute features independently per fold — then be called again on the full training set with best hyperparameters.

### Cross-Validation Strategies
- Default: Repeated stratified k-fold
- Multi-cohort: Leave-One-Dataset-Out (LODO) via `construct_stratified_cohort_folds()`

### Parallelization
- `doParallel` + `foreach %dopar%` for cross-fold parallelization
- XGBoost uses internal threading — external parallel is disabled to avoid contention

### Hyperparameter Tuning
- Metric-based optimization: AUROC, AUPRC, Accuracy, C-index
- Grid search within CV folds → best params applied to full training data

### Feature Preprocessing (internal `preprocess_features()`)
- Near-zero variance removal
- Collinearity filtering (correlation threshold)
- Removal of features constant within any target class (classification only)

### Reproducibility (`seed` argument, default `123`)
- `compute_features.training.ML()`, `compute_features.ML()` and `compute_shap_values()` take `seed`; `NULL` leaves the RNG untouched.
- The training functions pass it to `compute_k_fold_CV()` / `compute_k_fold_CV_survival()`, which call `set.seed(seed)` right before drawing folds. caret draws per-resample seeds for its parallel workers from the main session's RNG, so this also covers `ncores > 1`.
- `%dopar%` loops don't share the main RNG, so each iteration seeds itself: `seed + match(resample, resamples)` in `compute_shap_values()`, and `seed + 1000L * fold_i + parameter_i` in the survival custom-fold branch. Results are therefore independent of worker scheduling and `n_cores`.
- Randomness inside a user-supplied `fold_construction_fun` that runs its own parallel workers is not covered.

### Fold Models & SHAP
- During CV the model of each fold is saved to `fold_models_dir` (default `Results/fold_models/<task_type>`), and only the files of the selected model and `bestTune` are kept.
- `compute_shap_values(model_trained, task_type, ...)` takes everything (training data, outcome, folds, tuned hyperparameters) from the model: caret's `$trainingData` / `$pred` for classification, `$trainingData` / `$Resample_matrix` for survival. Each sample is explained only by the fold models that held it out, then summarized across repeats by the median.
- Standard-CV folds are refitted if their saved model is missing; custom-fold models (`fold_construction_fun`) cannot be rebuilt, so a missing fold model is an error.

### Performance Curves (`compute_prediction()` → `get_curves()`)
- AUROC/AUPRC are stored as `list(estimate, lower, upper)`: `estimate` is the value on the full test set, and the CI comes from 1000 bootstrap resamples (`bootstrap_auc()`, seed 123).
- The shaded bands are pointwise 95% bootstrap bands. `bootstrap_auc()` evaluates each resample's curve on a 101-point grid (`roc_at_grid()`: best sensitivity at FPR ≤ x; `prc_at_grid()`: precision at the first threshold reaching recall x), returned as `Curve_bands`. `get_curves(roc_band, prc_band)` draws them under the curve; the LODO branch draws no band.

### Outcome Alignment (`compute_features.ML()`)
- Labels are taken from `coldata` by row name (`coldata[rownames(features_train), , drop = FALSE]`), so they follow the feature rows' order; missing samples are an error. Never subset with `%in%` and then attach columns by position — that silently scrambled labels whenever feature rows weren't in `coldata` order.

---

## Bundled Example Datasets

| Dataset | Content |
|---|---|
| `data_example_classification` | Breast Cancer Wisconsin (from mlbench) |
| `data_example_survival` | Lung cancer survival (from survival package) |
| `counts_example` | Gene expression matrix — Gide et al. 2019 melanoma cohort (4.4 MB) |
| `coldata_example` | Metadata with anti-PD-1 therapy response labels for Gide cohort |

Access via `data(dataset_name)` after `library(pipeML)`.

---

## Development Workflow

### Code Style
- 2-space indentation (configured in pipeML.Rproj)
- Roxygen2 for all exported function documentation
- Run `devtools::document()` after any signature or documentation changes (regenerates `man/` and `NAMESPACE`)
- Never edit `man/` or `NAMESPACE` directly

### Building & Checking
```r
devtools::load_all()       # Load package in dev mode
devtools::document()       # Regenerate docs & NAMESPACE
devtools::check()          # Full R CMD check
devtools::build_vignettes() # Build vignettes
pkgdown::build_site()      # Rebuild docs website
```

### CI/CD (GitHub Actions)
- **R-CMD-check:** Runs on push/PR to main — tests macOS, Windows, Ubuntu across R devel/release/oldrel-1
- **pkgdown:** Auto-deploys website to `gh-pages` branch on push to main or release

### Testing
- `testthat` (v3) is configured in DESCRIPTION but **no test suite exists yet**
- Current testing is done via manual vignette execution
- New functionality should eventually be covered by tests in `tests/testthat/`

---

## Common Tasks

### Adding a new ML algorithm
1. Add to `compute_k_fold_CV()` / `compute_custom_k_fold_CV()` in `machine_learning.R`
2. Add a hyperparameter grid entry in `get_tune_grid()`
3. Handle retraining in `wrapper_train_best_hyperparams_classification()`
4. Document in vignette and function `@param method` Roxygen docs

### Adding a new metric
1. Add calculation helper (follow `calculate_auroc()` pattern)
2. Wire into `calculate_cv_metrics()` and `compute_prediction()`

### Modifying the vignette
Edit `vignettes/pipeML.Rmd`. Run `devtools::build_vignettes()` to test locally. The pkgdown CI will publish on merge to main.

---

## Notes & Gotchas

- The core implementation is a single large file (`machine_learning.R`, ~5,500 lines). Internal helpers are not exported — check NAMESPACE before assuming a function is public.
- caret must be *attached*, not just loaded: its `"knn"` model code calls `knn3()` without a namespace prefix, so a knn fit with `trainControl(method = "none")` fails with `could not find function "knn3"` otherwise. Classification entry points call the internal `ensure_caret()`.
- Don't add `set.seed()` calls inside helpers: they override the user's `seed` mid-run (`get_tune_grid()` used to do this). Seed only at the entry points described under Reproducibility.
- Roxygen markdown is enabled (`Roxygen: list(markdown = TRUE)`): write `95%`, not `95\%` (it becomes `\\%` in the Rd, which comments out the rest of the line), and avoid `[a, b]`-style ranges in plain text (parsed as links). Run `devtools::check_man()` after doc edits.
- CRAN doesn't allow non-CRAN required dependencies: `multideconv` and `fastshap` (both in `Remotes`) block a CRAN submission.
- Survival models need the `censored` package (in Suggests) to register parsnip's "censored regression" engines. Every survival entry point calls the internal `ensure_censored()`, which loads its namespace; users don't need `library(censored)`. Survival formulas must use `survival::Surv(...)`, not bare `Surv(...)` — the bare form only works when some other package happened to attach `survival`.
- `multideconv` is a remote (GitHub) dependency — not on CRAN. Installation requires `remotes::install_github("VeraPancaldiLab/multideconv")`.
- SHAP computation via `fastshap::explain()` can be memory-intensive on large datasets.
- XGBoost parallel contention: when using `doParallel`, XGBoost nthread is set to 1 internally to prevent nested parallelism crashes.
- The `docs/` directory is gitignored — pkgdown output is built and deployed by CI only.
- `compute_shap_values(model_trained, ...)` expects `res$Model` from `compute_features.training.ML()` (for classification, the caret `train` object). It validates this and stops with a clear error if given the wrong object (e.g. `res` instead of `res$Model`).

## Known Issues / TODO

### Open bugs
- **LODO is ignored for survival and leaks the cohort label.** `compute_features.training.ML()` and `compute_features.ML()` add a `dataset` column when `LODO = TRUE`, but never pass `LODO`/`batch_id` to `compute_k_fold_CV_survival()`. No cohort-stratified folds are built, and `dataset` is fit as an ordinary predictor (`Surv(time, event) ~ .`) — silently.
- **`ncores` is ignored in standard survival CV.** `compute_k_fold_CV_survival()` only creates a cluster in the custom-fold branch; the standard branch loops sequentially over folds and hyperparameter grids (up to 125 combinations per model), so a 2-fold run can take 10+ minutes.
- **`compute_ml_survival()` swallows fitting errors**: it converts them into a warning and returns `NULL`, so failures go unnoticed (this is how the bare-`Surv()` bug stayed hidden). It also accepts a `fold_models_dir` argument it never uses.
- **`data_example_survival` codes events as 1/2** (from `survival::lung`), while the docs say 0/1. It works because `Surv()` accepts both.
- **No test suite.** `tests/testthat/` doesn't exist; verification is manual.
- Cosmetic: on precision-recall plots, the curve's final vertical drop at recall = 1 (to precision = prevalence) falls below the confidence band, which ends at the first threshold reaching full recall.

### Performance
- **`compute_shap_values()`'s cost is extremely method-dependent, and this is invisible
  to the caller until it's too late.** Measured empirically (melanoma LODO dataset,
  ~250-300 x ~15 NMF-factor features, `nsim = 100`, one `fastshap::explain()` call per
  CV resample): `glmnet` ≈ 16s/resample, `KNN` ≈ 25s/resample, but `svmRadial` ≈
  **930s/resample** — a ~58x slowdown, because SVM/KNN-family predict methods are
  computationally heavier per call (KNN recomputes distances to the full training set;
  SVM's kernel evaluation is per-support-vector) and `fastshap::explain()` calls the
  prediction function repeatedly (proportional to `nsim`) for every resample. With the
  default `k_folds x n_rep` producing up to 100 resamples, an SVM-selected model can
  turn what's normally a ~30min job into a ~26-hour one, with zero warning beforehand.
  Concrete improvements worth making:
  - **Expose `nsim`** as a `compute_shap_values()` parameter instead of the hardcoded
    100 (`machine_learning.R` line ~2964) — callers with a slow-predict method could
    trade precision for speed deliberately, instead of being stuck with a fixed cost
    multiplier they can't control.
  - **Expose which/how many resamples to explain**, rather than always looping over
    every one of `unique(model_trained$pred$Resample)` — explaining a representative
    subset (e.g. 10-20 of 100) would give an approximate-but-fast SHAP estimate on
    request.
  - **Print a per-resample time estimate after the first resample completes** (before
    committing to the rest of the loop) — currently there's no feedback at all until
    the whole thing finishes or the caller gives up waiting; even a single
    `cat(sprintf("First resample took %.1fs; estimated total: %.1fmin for %d
    resamples\n", ...))` after resample 1 would let users abort early with an informed
    decision instead of guessing.
  - Consider flagging known-slow methods (`svmRadial`, `svmLinear`, `knn`, and other
    instance-/kernel-based predictors) with a `message()` up front suggesting a
    reduced `nsim` or resample subset for those specifically.
