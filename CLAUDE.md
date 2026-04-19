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
  machine_learning.R     # All implementation (~5,700 lines, core file)
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

## Exported Functions (8 total)

All live in `R/machine_learning.R`.

| Function | Purpose |
|---|---|
| `feature.selection.boruta()` | Repeated Boruta feature selection with parallel support |
| `compute_features.training.ML()` | Train models on training data with repeated k-fold CV |
| `compute_features.ML()` | Combined train + predict workflow (training + testing) |
| `compute_prediction()` | Generate predictions on test data using trained model |
| `get_curves()` | ROC and Precision-Recall curves with confidence intervals |
| `compute_shap_values()` | SHAP feature importance values across resamples |
| `plot_shap_stability()` | Visualize SHAP importance stability across resamples |
| `plot_survival_performance()` | Kaplan-Meier curves stratified by predicted risk groups |

---

## Supported ML Algorithms

**Classification (11, via caret):**
`treebag`, `rf`, `C5.0`, `glmnet` (elastic/lasso/ridge), `knn`, `rpart`, `svmRadial`, `svmLinear`, `xgbTree`

**Survival (7, via tidymodels/parsnip/censored):**
Cox PH, Elastic Net Cox, AFT parametric, Conditional Inference Trees, Bagged CART, Random Survival Forests, Gradient Boosting (censored)

---

## Key Dependencies

**Imports (must be installed):**
`caret`, `Boruta`, `doParallel`, `foreach`, `dplyr`, `tidyr`, `tibble`, `purrr (>= 1.0.2)`, `ggplot2`, `reshape2`, `survival`, `survminer`, `fastshap`, `dials`, `parsnip`, `rsample`, `workflows`, `tune`, `yardstick`, `grDevices`, `parallel`, `stats`

**Suggests (optional, needed for specific algorithms):**
`testthat (>= 3.0.0)`, `knitr`, `rmarkdown`, `C50`, `randomForest`, `glmnet`, `xgboost`, `kernlab`, `recipes`, `tidymodels`, `censored`, `flexsurv`, `coin`, `aorsf`, `WGCNA`, `cowplot`, `matlib`

**Remotes (GitHub):**
`VeraPancaldiLab/multideconv` — custom deconvolution package (optional)

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
- Boruta iterations support parallel runs

### Hyperparameter Tuning
- Metric-based optimization: AUROC, AUPRC, Accuracy, C-index
- Grid search within CV folds → best params applied to full training data

### Model Stacking
- Ensemble meta-learning support
- Weighted feature importance from base models + meta-learner

### Feature Preprocessing (internal `preprocess_features()`)
- Near-zero variance removal
- Collinearity filtering (correlation threshold)
- Boruta selection with tentative feature handling

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

- The core implementation is a single large file (`machine_learning.R`, ~5,700 lines). Internal helpers are not exported — check NAMESPACE before assuming a function is public.
- `multideconv` is a remote (GitHub) dependency — not on CRAN. Installation requires `remotes::install_github("VeraPancaldiLab/multideconv")`.
- SHAP computation via `fastshap::explain()` can be memory-intensive on large datasets.
- XGBoost parallel contention: when using `doParallel`, XGBoost nthread is set to 1 internally to prevent nested parallelism crashes.
- The `docs/` directory is gitignored — pkgdown output is built and deployed by CI only.
