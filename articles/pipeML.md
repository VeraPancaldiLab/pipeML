# pipeML

``` r
library(pipeML)
#> Warning: replacing previous import 'dplyr::explain' by 'fastshap::explain' when
#> loading 'pipeML'
library(caret)
#> Loading required package: ggplot2
#> Loading required package: lattice
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(survival)
#> 
#> Attaching package: 'survival'
#> The following object is masked from 'package:caret':
#> 
#>     cluster
library(censored)
#> Loading required package: parsnip
```

This tutorial demonstrates how to use the `pipeML` package for training
and testing machine learning models. It introduces the main functions of
the pipeline and guides you through additional functions for
visualization.

## **Classification tasks**

Load example data

``` r
data = pipeML::data_example_classification
X <- data %>% dplyr::select(-target)
y <- data$target
```

For example purposes we make a Train/Test split

``` r
set.seed(123)

train_idx <- caret::createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

y_train <- y[train_idx]
y_test  <- y[-train_idx]
```

### **Feature selection**

Apply repeated feature selection using the Boruta algorithm:

``` r
y = as.factor(y)
res_boruta <- feature.selection.boruta(
  data = data,
  iterations = 20,
  fix = FALSE,
  doParallel = FALSE,
  threshold = 0.8,
  file_name = "Test",
  return = FALSE
)
```

Inspect the results:

``` r
head(res_boruta$Matrix_Importance)

cat("Confirmed features:\n", res_boruta$Confirmed)

cat("Tentative features:\n", res_boruta$Tentative)
```

To further assess tentative features, set `fix = TRUE` to rerun the
selection and confirm or reject them:

``` r
res_boruta <- feature.selection.boruta(
  data = data,
  iterations = 20,
  fix = TRUE,
  doParallel = FALSE,
  threshold = 0.8,
  file_name = "Test",
  return = FALSE
)
```

For faster execution, enable parallelization (ensure parameters
`doParallel` and `workers` are set):

``` r
res_boruta <- feature.selection.boruta(
  data = data,
  iterations = 20,
  fix = FALSE,
  doParallel = FALSE,
  workers = 2,
  threshold = 0.8,
  file_name = "Test",
  return = FALSE
)
```

### **Train machine learning models for classification tasks**

Train and tune models using repeated stratified k-fold cross-validation:

``` r
res <- compute_features.training.ML(features_train = X_train, 
                                    target_var = y_train,
                                    task_type = "classification",
                                    trait.positive = "1",
                                    metric = "AUROC",
                                    k_folds = 2,
                                    n_rep = 1,
                                    file_name = "Example_classification",
                                    return = F)
```

Access the best-trained model:

``` r
res$Model
```

View all trained and tuned machine learning models:

``` r
names(res$ML_Models)
```

``` r
knitr::include_graphics("figures/AUROC_classification.png")
```

![Figure 1. Models training
performance.](figures/AUROC_classification.png)

Figure 1. Models training performance.

### **Model prediction**

After training, users can predict on new data by using the
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
function. You can specify which metric to maximize when determining the
optimal classification threshold. Supported values for maximize include:
“Accuracy”, “Precision”, “Recall”, “Specificity”, “Sensitivity”, “F1”,
and “MCC”.

``` r
pred = compute_prediction(model = res$Model, 
                          test_data = X_test, 
                          target_var = y_test, 
                          task_type = "classification",
                          trait.positive = "1", 
                          file.name = "Example_classification", 
                          return = T)
```

Check predicitions

``` r
head(pred$Predictions)
```

Verify threshold-based predictions metrics

``` r
head(pred$Metrics)
```

If `return = TRUE`,
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
will saved the RO and PR curves in the `Results/` directory.

``` r
knitr::include_graphics("figures/ROC.png")
```

![Figure 2. ROC curve.](figures/ROC.png)

Figure 2. ROC curve.

``` r
knitr::include_graphics("figures/PR.png")
```

![Figure 3. PR curve.](figures/PR.png)

Figure 3. PR curve.

### **Compute SHAP Values**

SHAP values help interpret machine learning predictions by quantifying
feature contributions.

``` r
df = cbind(X_train, target = as.numeric(y_train)) ## compute_shap_values expects predictors and target in a single data.frame

shap_classification <- compute_shap_values(
  model_trained = res$Model, 
  data_train = df, 
  task_type = "classification", 
  target_col = "target", 
  trait.positive = "2", # Notice here that because of as.numeric() our target variable changed to 1 and 2 so we will consider 2 as our new 1.
  n_cores = 2,  
  file.name = "Example_classification"
)
```

Visualize feature importance and interactions using shapviz:

``` r
sv <- shapviz::shapviz(
  shap = as.matrix(shap_classification$shap_values),
  X = df[, setdiff(colnames(df), "target")]
)

# Global feature importance (bar plot)
shapviz::sv_importance(sv, kind = "bar") +
  ggplot2::ggtitle("Global Feature Importance")

# Beeswarm plot
shapviz::sv_importance(sv, kind = "beeswarm") +
  ggplot2::ggtitle("SHAP Beeswarm Summary")

# Feature dependence plots for top 6 features
top_features <- names(sort(colMeans(abs(shap_classification$shap_values)), decreasing = TRUE))[1:6]

for (f in top_features) {
  print(
    shapviz::sv_dependence(sv, v = f) +
    ggplot2::ggtitle(paste0("Dependence: ", f))
  )
}
```

To apply model stacking, set `stack = TRUE`:

``` r
res <- compute_features.training.ML(features_train = deconvolution, 
                                    target_var = traitData$Best.Confirmed.Overall.Response,
                                    task_type = "classification",
                                    trait.positive = "CR",
                                    metric = "AUROC",
                                    stack = TRUE,
                                    k_folds = 2,
                                    n_rep = 2,
                                    LODO = FALSE,
                                    file_name = "Test",
                                    ncores = 2,
                                    return = FALSE)
```

Inspect the base models used in stacking:

``` r
res$Model$Base_models
```

Access the meta-learner:

``` r
res$Model$Meta_learner
```

## **Survival tasks**

In addition to classification and regression tasks, `pipeML` supports
survival analysis, allowing users to train machine learning models that
predict time-to-event outcomes.

In survival analysis, the response variable is defined by two
components:

- Time: follow-up or survival time
- Event: indicator of whether the event occurred (1) or the observation
  was censored (0)

`pipeML` automates the training, evaluation, and interpretation of
survival models using cross-validation and SHAP-based explanations.

Load example dataset for survival

``` r
data = pipeML::data_example_survival
X <- data %>% dplyr::select(-time, -status)
time <- data$time
event <- data$status
```

Similar to the previous example, we will split our data in Train/Test

``` r
set.seed(123)
train_idx <- caret::createDataPartition(event, p = 0.7, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

time_train <- time[train_idx]
time_test  <- time[-train_idx]

event_train <- event[train_idx]
event_test  <- event[-train_idx]
```

### **Train machine learning models for survival tasks**

The function
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
can also train survival models by setting task_type = “survival”.

``` r
res_survival <- compute_features.training.ML(
  features_train = X_train, 
  task_type = "survival",
  time_var = time_train,
  event_var = event_train,
  k_folds = 2,
  n_rep = 1,
  file_name = "Example_survival",
  ncores = 1,
  return = TRUE
)
#> Creating stratified v-fold CV (stratified by event rate only)
#> Running fold_i 1 
#> Running  cox_ph_survival 
#> Running  proportional_hazards_glmnet 
#> Running  survreg_flexsurv
#> Warning in value[[3L]](cond): Model fitting failed for survreg_flexsurv :
#> Please install the flexsurv package to use this engine.
#> Running  decision_tree_partykit
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Running  bag_tree_rpart 
#> Running  rand_forest_aorsf
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Running fold_i 2 
#> Running  cox_ph_survival 
#> Running  proportional_hazards_glmnet 
#> Running  survreg_flexsurv
#> Warning in value[[3L]](cond): Model fitting failed for survreg_flexsurv :
#> Please install the flexsurv package to use this engine.
#> Running  decision_tree_partykit
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Warning in value[[3L]](cond): Model fitting failed for decision_tree_partykit :
#> there is no package called 'coin'
#> Running  bag_tree_rpart 
#> Running  rand_forest_aorsf
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Warning in value[[3L]](cond): Model fitting failed for rand_forest_aorsf :
#> Please install the aorsf package to use this engine.
#> Best ML model found:  proportional_hazards_glmnet 
#> Returning model trained
```

Access to the best model

``` r
names(res_survival$ML_Models)
#> [1] "cox_ph_survival"             "proportional_hazards_glmnet"
#> [3] "survreg_flexsurv"            "decision_tree_partykit"     
#> [5] "bag_tree_rpart"
res_survival$Model$Model_object
#> ══ Workflow [trained] ══════════════════════════════════════════════════════════
#> Preprocessor: Formula
#> Model: proportional_hazards()
#> 
#> ── Preprocessor ────────────────────────────────────────────────────────────────
#> Surv(time, event) ~ .
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> 
#> Call:  glmnet::glmnet(x = data_obj$x, y = data_obj$y, family = "cox",      weights = weights, alpha = alpha, lambda = lambda) 
#> 
#>     Df %Dev  Lambda
#> 1    8 0.00 205.200
#> 2    8 0.03 187.000
#> 3    8 0.03 170.400
#> 4    8 0.04 155.200
#> 5    8 0.04 141.400
#> 6    8 0.04 128.900
#> 7    8 0.05 117.400
#> 8    8 0.05 107.000
#> 9    8 0.06  97.490
#> 10   8 0.06  88.830
#> 11   8 0.07  80.940
#> 12   8 0.08  73.750
#> 13   8 0.08  67.200
#> 14   8 0.09  61.230
#> 15   8 0.10  55.790
#> 16   8 0.11  50.830
#> 17   8 0.12  46.320
#> 18   8 0.13  42.200
#> 19   8 0.15  38.450
#> 20   8 0.16  35.040
#> 21   8 0.17  31.930
#> 22   8 0.19  29.090
#> 23   8 0.21  26.500
#> 24   8 0.22  24.150
#> 25   8 0.24  22.000
#> 26   8 0.27  20.050
#> 27   8 0.29  18.270
#> 28   8 0.32  16.650
#> 29   8 0.34  15.170
#> 30   8 0.37  13.820
#> 31   8 0.40  12.590
#> 32   8 0.44  11.470
#> 33   8 0.47  10.450
#> 34   8 0.51   9.525
#> 35   8 0.55   8.679
#> 36   8 0.59   7.908
#> 37   8 0.64   7.206
#> 38   8 0.69   6.565
#> 39   8 0.74   5.982
#> 40   8 0.79   5.451
#> 41   8 0.85   4.966
#> 42   8 0.91   4.525
#> 43   8 0.97   4.123
#> 44   8 1.03   3.757
#> 45   8 1.09   3.423
#> 46   8 1.16   3.119
#> 
#> ...
#> and 54 more lines.
```

Check training metrics

``` r
head(res_survival$Model$Prediction_folds)
#>   predictions                       model Resample rowIndex penalty mixture
#> 1  0.13035757 proportional_hazards_glmnet    Fold1        2   1e-10       0
#> 2  1.22509805 proportional_hazards_glmnet    Fold1        9   1e-10       0
#> 3  0.92662674 proportional_hazards_glmnet    Fold1       12   1e-10       0
#> 4 -0.07545289 proportional_hazards_glmnet    Fold1       14   1e-10       0
#> 5  0.82133556 proportional_hazards_glmnet    Fold1       16   1e-10       0
#> 6  0.99870685 proportional_hazards_glmnet    Fold1       17   1e-10       0
#>     c_index
#> 1 0.6141016
#> 2 0.6141016
#> 3 0.6141016
#> 4 0.6141016
#> 5 0.6141016
#> 6 0.6141016
```

``` r
knitr::include_graphics("figures/cindex_survival.png")
```

![Figure 4. Models training performance.](figures/cindex_survival.png)

Figure 4. Models training performance.

### **Compute SHAP Values**

We will used the same function `compute_shap_values` computed before but
this time with the `task_type = "survival"`.

``` r
shap_survival <- compute_shap_values(
  model_trained = res_survival$Model,
  data_train = df,
  task_type = "survival",
  time_col = "time",
  event_col = "status",
  n_cores = 2,
  file.name = "Example_survival"
)
```

### **Predict Survival risk**

After training, predictions can be generated using
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md).

``` r
pred <- compute_prediction(
  model = res_survival$Model,
  test_data = X_test,
  task_type = "survival",
  time_var = time_test,
  event_var = event_test,
  file.name = "Example_survival",
  return = TRUE
)
```

Unlike classification models, survival models may return different types
of predictions depending on the model used. Some models predict a risk
score, while others predict expected survival time or survival
probability.

The
[`compute_prediction()`](https://verapancaldilab.github.io/pipeML/reference/compute_prediction.md)
function automatically handles these differences. Internally, it
attempts several prediction types and converts them into a standardized
risk score so that results can be compared across models:

- Risk score (**linear_pred**) – typically produced by Cox models.
  Higher values indicate higher predicted risk.
- Predicted survival time (**time**) – produced by some parametric
  survival models. Higher values indicate longer survival, so the values
  are internally reversed to represent risk.
- Survival probability (**survival**) – probability of surviving at a
  given time point. Higher probabilities correspond to lower risk, so
  these values are also internally reversed.

After this standardization step, predictions are always interpreted in
the same way:

**Higher prediction values correspond to higher predicted risk**

This allows `pipeML` to compute performance metrics such as the
concordance index (C-index) and to stratify patients into risk groups
for Kaplan–Meier visualization.

``` r
knitr::include_graphics("figures/KM.png")
```

![Figure 5. KM plot.](figures/KM.png)

Figure 5. KM plot.

### **Train and predict using the machine learning models**

If a separate testing dataset is already available, pipeML allows you to
train models and generate predictions in a single step using the
compute_features.ML() function. This function performs model training,
hyperparameter tuning, and evaluation on the training data, and then
applies the selected model to the test dataset.

Load data

``` r
data = pipeML::data_example_classification
X <- data %>% dplyr::select(-target)
y <- data$target

set.seed(123)

train_idx <- caret::createDataPartition(y, p = 0.8, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
```

Classification tasks

``` r
res <- compute_features.ML(features_train = X_train, 
                           features_test = X_test, 
                           coldata = data,
                           task_type = "classification",
                           trait = "target",
                           trait.positive = "1",
                           metric = "AUROC",
                           k_folds = 2,
                           n_rep = 1,
                           ncores = 2,
                           file_name = "Test",
                           return = FALSE)
```

Check training metrics

``` r
head(res$Model$Model$pred)
```

Check prediction metrics

``` r
head(res$Metrics)
```

Load data

``` r
data = pipeML::data_example_survival
X <- data %>% dplyr::select(-time, -status)
time <- data$time
event <- data$status

set.seed(123)
train_idx <- caret::createDataPartition(event, p = 0.7, list = FALSE)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

time_train <- time[train_idx]
time_test  <- time[-train_idx]

event_train <- event[train_idx]
event_test  <- event[-train_idx]
```

Survival tasks

``` r
res <- compute_features.ML(features_train = X_train, 
                           features_test = X_test, 
                           coldata = data,
                           task_type = "survival",
                           time_var = "time",
                           event_var = "status",
                           k_folds = 2,
                           n_rep = 1,
                           ncores = 2,
                           file_name = "Test",
                           return = FALSE)
```

### **Leaving one dataset out (LODO) analysis**

If your data includes features from multiple cohorts, `pipeML` provides
a flexible approach to perform Leave-One-Dataset-Out (LODO) analysis.
This is achieved by applying k-fold stratified sampling across batches,
ensuring that each fold preserves the batch structure while maintaining
class balance. For this set `LODO = TRUE` and provided column name
containing the batch information in the `batch_var` variable

Below, we demonstrate how to perform a LODO analysis using simulated
datasets:

Simulate datasets from different batches

``` r
set.seed(123)

# Simulate traitData with 3 cohorts
traitData <- data.frame(
  Sample = paste0("Sample", 1:90),
  Response = sample(c("R", "NR"), 90, replace = TRUE),
  Cohort = rep(paste0("Cohort", 1:3), each = 30),
  stringsAsFactors = FALSE
)
rownames(traitData) <- traitData$Sample

# Simulate counts matrix 
counts_all <- matrix(rnorm(1000 * 90), nrow = 1000, ncol = 90)
rownames(counts_all) <- paste0("Gene", 1:1000)
colnames(counts_all) <- traitData$Sample

# Simulate some example features 
features_all <- matrix(runif(90 * 15), nrow = 90, ncol = 15)
rownames(features_all) <- traitData$Sample
colnames(features_all) <- paste0("Feature", 1:15)
```

Perform LODO analysis

``` r
prediction = list()
i = 1
for (cohort in unique(traitData$Cohort)) {
    
  # Test cohort
  traitData_test = traitData %>%
    filter(Cohort == cohort)
  counts_test = counts_all[,colnames(counts_all)%in%rownames(traitData_test)]
  features_test = features_all[rownames(features_all)%in%rownames(traitData_test),]

  # Train cohort
  traitData_train = traitData %>%
    filter(Cohort != cohort)
  counts_train = counts_all[,colnames(counts_all)%in%rownames(traitData_train)]
  features_train = features_all[rownames(features_all)%in%rownames(traitData_train),]

  #### ML Training
  res = compute_features.training.ML(features_train = features_train, 
                                     target_var = traitData_train$Response, 
                                     task_type = "classification",
                                     trait.positive = "R", 
                                     metric = "AUROC", 
                                     k_folds = 2, 
                                     n_rep = 1, 
                                     LODO = TRUE,
                                     batch_var = traitData_train$Cohort, 
                                     ncores = 2, 
                                     return = F)

  #### ML predicting
  
  #### Testing
  pred = compute_prediction(model = res$Model, 
                            test_data = features_test, 
                            target_var = traitData_test$Response, 
                            task_type = "classification",
                            trait.positive = "R", 
                            return = F)
  
  #### Save results
  prediction[[i]] = pred
  names(prediction)[i] = cohort
  i = i + 1
}
```

For plotting we will make use of our `get_curves` function adapted in a
case for multiple cohorts. This function is implemented inside
`compute_prediction` to make the RO and PR curves.

Extract prediction metrics and join

``` r
roc_data <- lapply(names(prediction), function(cohort) {
  df <- prediction[[cohort]]$Metrics
  df$cohort <- cohort
  df
}) %>% dplyr::bind_rows()

auc_roc <- list(
  mean  = sapply(prediction, function(x) x$AUC$AUROC$mean),
  lower = sapply(prediction, function(x) x$AUC$AUROC$lower),
  upper = sapply(prediction, function(x) x$AUC$AUROC$upper)
)

auc_prc <- list(
  mean  = sapply(prediction, function(x) x$AUC$AUPRC$mean),
  lower = sapply(prediction, function(x) x$AUC$AUPRC$lower),
  upper = sapply(prediction, function(x) x$AUC$AUPRC$upper)
)
```

Plot RO and PR curves

``` r
get_curves(
  data = roc_data,
  color = "cohort",
  auc_roc = auc_roc,
  auc_prc = auc_prc,
  LODO = TRUE,
  file.name = "LODO_cohort_example",
  width = 9,
  height = 9
)
```

``` r
knitr::include_graphics("figures/RO_LODO.png")
```

![Figure 6. LODO RO curves.](figures/RO_LODO.png)

Figure 6. LODO RO curves.

``` r
knitr::include_graphics("figures/PR_LODO.png")
```

![Figure 7. LODO PR curves.](figures/PR_LODO.png)

Figure 7. LODO PR curves.

### **Leakage-aware custom cross-validation**

A central design principle of `pipeML` is to prevent information leakage
during model training and evaluation. In many machine learning
workflows, feature engineering steps are applied to the full dataset
before cross-validation, which can inadvertently introduce information
from the test folds into the training process. This leads to
overoptimistic performance estimates.

To address this, `pipeML` provides built-in support for custom fold
construction through the fold_construction_fun argument in
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md).
This mechanism allows feature engineering and preprocessing steps to be
recomputed independently within each cross-validation fold, ensuring
that test samples never influence the training process.

This capability is a core component of the `pipeML` pipeline, enabling
leakage-aware model development for datasets where features depend on
the full sample structure.

##### Why custom fold construction is important?

In many biological and high-dimensional datasets, features are not
independent variables but are derived from the data itself. Examples
include:

- correlation-based clustering
- dimensionality reduction (e.g., PCA)
- gene set enrichment or pathway scoring
- transcription factor activity inference
- aggregation of features across samples

If these transformations are applied to the entire dataset before
cross-validation, the test samples influence how the feature space is
constructed. As a result, the model indirectly “sees” information from
the test data during training.

By recomputing these steps within each fold, `pipeML` ensures that: -
training data are used to define the feature space - test samples are
projected onto the learned space without influencing it - model
evaluation reflects true out-of-sample performance

This approach closely mimics how the model would behave when applied to
completely unseen data, producing more realistic performance estimates.

Step 1 – Define a base feature function (example: WGCNA)

**Structure of the Base Feature Function**

The function is designed around two operational modes: - Training Mode
(modules is NULL): The function learns the feature structure from the
input dataset. - Projection Mode (modules provided): The function
applies a previously learned structure to new data.

Why the Function Returns Two Objects - features → always returned; the
transformed representation of the current dataset (training or test). -
modules (or structure) → returned only when learning from training data;
needed to project test data in future steps.

This ensures a leakage-aware workflow: - Training data defines the
feature space. - Test data is projected without altering the learned
structure.

Structure of base function

- data: features as rows, samples as columns
- structure: precomputed structure (e.g., clusters, components); NULL
  for training
- … : additional arguments specific to the algorithm used

``` r
compute_features_modular <- function(data, structure = NULL, ...) {
  

  # TRAINING MODE
  if (is.null(structure)) {
    
    # -------------------- REPLACE THIS BLOCK --------------------
    structure <- learn_structure(data, ...)   # user-defined function
    # -------------------- REPLACE THIS BLOCK --------------------
    
  }
  
  # PROJECT MODE
  
  # -------------------- REPLACE THIS BLOCK --------------------
  features <- project_data(data, structure, ...)   # user-defined function
  # -------------------- REPLACE THIS BLOCK --------------------

  
  return(list(features = features, structure = structure))
}
```

Here we illustrate a correlation-based feature computation using
Weighted Gene Co-expression Network Analysis (WGCNA). This is just an
example: in practice, you can use any feature computation that depends
on multiple samples, such as clustering, PCA, among others.

``` r
compute_features_modular <- function(counts, power, modules = NULL) {

  ## Just preprocessing (IGNORE)
  rownames(counts) <- gsub("-", ".", rownames(counts)) 
  datExpr <- t(counts) 
  cor <- WGCNA::cor

  # TRAINING MODE
  if (is.null(modules)) {
    
    # -------------------- REPLACE THIS BLOCK --------------------
    net <- WGCNA::blockwiseModules(datExpr, power = power)
    modules <- net$colors
    names(modules) <- colnames(datExpr)
    # -------------------- REPLACE THIS BLOCK --------------------
    
  }
  
  # PROJECT MODE 
  
  # -------------------- REPLACE THIS BLOCK --------------------
  module_features <- sapply(unique(modules), function(mod) {
    genes <- names(modules[modules == mod])
    pc <- prcomp(datExpr[, genes, drop = FALSE])
    pc$x[, 1]  
  })
  # -------------------- REPLACE THIS BLOCK --------------------

  ## Just formatting (IGNORE)
  colnames(module_features) <- paste0("Module_", seq_len(ncol(module_features)))
  
  
  return(list(features = as.matrix(module_features), structure = modules))
}
```

Step 2 – Make the function suitable for `pipeML`

We then need to extend and give the correct format to this function to
make it suitable for running across folds inside `pipeML`

This template provides a modular framework to prepare cross-validation
folds for `pipeML` in a leakage-aware way. It separates training vs
projection:

- Training mode: computes features and learns the data structure from
  the training folds
- Projection mode: applies the learned structure to held-out folds
  without influencing it.

Users can easily adapt this template by replacing their previous
`compute_features_modular` function

**Note**

In `pipeML` data corresponds to samples as rows and features as columns.
If your `compute_features_modular()` needs features as rows, make sure
to [`t()`](https://rdrr.io/r/base/t.html) inside this function.

The parameters `data`, `folds`, `besttune` are taking care inside pipeML
automatically once `fold_construction_fun` is set. Be sure to not change
the parameteres names or delete them!

- … : Additional parameters passed to your feature function

Make sure `compute_features_modular()` returns a matrix with the
features. If not, make sure to extract them before adding the target
column

``` r
prepare_custom_folds <- function(data, folds = NULL, bestune = NULL, ...) {
  
  if (!is.null(bestune)) {
    
    obs_train <- data$target
    data$target <- NULL
    
    
    # -------------------- REPLACE THIS BLOCK --------------------
    result <- compute_features_modular(data, ...)  # or t(data) if necessary
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
    train_features_final = result$features ### if necessary
    train_features_final$target <- obs_train
    
    custom_output <- result
    
    return(list(train_features_final, custom_output, bestune))
    
  } else {
    
    processed_folds <- list()
    
    for (i in seq_along(folds)) {
      
      train_idx <- folds[[i]]
      test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
      
      train_data <- data[train_idx, , drop = FALSE]
      obs_train <- train_data$target
      train_data$target <- NULL
      
      
      # -------------------- REPLACE THIS BLOCK --------------------
      train_result <- compute_features_modular(train_data, ...) # or t(train_data) if necessary
      # -------------------- REPLACE THIS BLOCK --------------------
      
      
      train_features = train_result$features ## if necessary
      train_features$target <- obs_train
      
      test_data <- data[test_idx, , drop = FALSE]
      obs_test <- test_data$target
      test_data$target <- NULL
      
      
      # -------------------- REPLACE THIS BLOCK --------------------
      test_features <- compute_features_modular(
        test_data, # or t(test_data) if necessary
        structure = train_result$structure,  
        ...
      )
      # -------------------- REPLACE THIS BLOCK --------------------
      
      
      test_features = test_features$features ## if necessary

      processed_folds[[i]] <- list(
        train_data = train_features,
        test_data  = test_features,
        obs_test   = obs_test,
        rowIndex   = test_idx,
        fold_name  = names(folds)[i]
      )
    }
    
    for (i in seq_along(processed_folds)) {
      filename <- file.path("Results", paste0("fold_", names(folds)[i], ".rds"))
      saveRDS(processed_folds[[i]], file = filename)
    }
    
    return(processed_folds)
  }
}
```

Here we illustrate how the function will look applying our
`compute_features_modular()` function

Notice that each time I call my function I am setting my additional
argument `power`

``` r
prepare_WGCNA_folds <- function(data, folds = NULL, bestune = NULL, power) {
  
  if (!is.null(bestune)) {
    
    obs_train <- data$target
    data$target <- NULL
        
    
    # -------------------- REPLACE THIS BLOCK --------------------
    wgcna_result <- compute_features_modular(t(data), power = power)
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
    train_cell_data_final <- as.data.frame(wgcna_result$features)
    train_cell_data_final$target <- obs_train
    
    custom_output <- wgcna_result
    
    return(list(train_cell_data_final, custom_output, bestune))
    
  } else {
    
    processed_folds <- list()
    
    for (i in seq_along(folds)) {
      
      train_idx <- folds[[i]]
      test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
      
      train_data <- data[train_idx, , drop = FALSE]
      obs_train <- train_data$target
      train_data$target <- NULL
      
      
      # -------------------- REPLACE THIS BLOCK --------------------
      train_result <- compute_features_modular(t(train_data), power = power)
      # -------------------- REPLACE THIS BLOCK --------------------
      
      
      train_features <- as.data.frame(train_result$features)
      train_features$target <- obs_train
      
      test_data <- data[test_idx, , drop = FALSE]
      obs_test <- test_data$target
      test_data$target <- NULL
      
      
      # -------------------- REPLACE THIS BLOCK --------------------
      test_features <- compute_features_modular(t(test_data), modules = train_result$modules)
      # -------------------- REPLACE THIS BLOCK --------------------
      
      
      test_features <- as.data.frame(test_features$features)
      
      processed_folds[[i]] <- list(
        train_data = train_features,
        test_data  = test_features,
        obs_test   = obs_test,
        rowIndex   = test_idx,
        fold_name  = names(folds)[i]
      )
    }
    
    for (i in seq_along(processed_folds)) {
      filename <- file.path("Results", paste0("fold_", names(folds)[i], ".rds"))
      saveRDS(processed_folds[[i]], file = filename)
    }
    
    return(processed_folds)
  }
}
```

Load data example (Gene expression counts)

``` r
counts = pipeML::counts_example
coldata = pipeML::coldata_example

set.seed(123)

train_idx <- caret::createDataPartition(coldata$Response, p = 0.8, list = FALSE)

counts_train <- counts[train_idx, ]
counts_test  <- counts[-train_idx, ]

coldata_train <- coldata[train_idx]
coldata_test  <- coldata[-train_idx]
```

Step 3 – Run Custom Cross-Validation

Once your custom fold function is ready, pass it to
[`compute_features.training.ML()`](https://verapancaldilab.github.io/pipeML/reference/compute_features.training.ML.md)
via fold_construction_fun.

Note the argument `fold_construction_args_fixed` corresponds to the
additional parameters passed in `prepare_custom_folds()` function (in
this example `prepare_WGCNA_folds()`). They should be set with the value
to use when running your function. If your function does not have any
parameters, you can ignore this argument and it will be set up to NULL.

``` r
res_custom   = compute_features.training.ML(features_train = t(counts_train),
                                            target_var     = coldata_train$Response,
                                            task_type      = "classification",
                                            trait.positive = "R",
                                            metric         = "AUROC",
                                            k_folds        = 2,
                                            n_rep          = 1,
                                            return         = FALSE,
                                            fold_construction_fun = prepare_WGCNA_folds,
                                            fold_construction_args_fixed = list(power=6))
```

Notice that `res_custom$Custom_output` contains the output features of
your base function in case these are needed (e.g. for prediction see
below)

``` r
res_custom$Custom_output
```

##### Tunable hyperparameters within custom fold functions

In some scenarios, the feature construction step may include parameters
whose values can influence model performance. For example, when
computing WGCNA, parameters such as soft-thresholding power, minimum
module size, module merging threshold, and module splitting sensitivity
may affect the resulting features and therefore the downstream model
performance.

To address this, `pipeML` supports hyperparameter tuning within your
custom fold functions. This allows users to identify which parameter
values lead to the best predictive performance.

Briefly, the user provides a grid of candidate parameter values, and
`pipeML` will:

- construct custom cross-validation folds for each parameter combination
- train and evaluate machine learning models for each configuration
- compare the resulting performance across folds and repetitions
- return the parameter values that maximize the selected performance
  metric

Before that, user needs to modify `compute_features_modular` and
`prepare_custom_folds()` to account for all these parameters.

For the `compute_features_modular` we are only going to add the tunable
parameters in our function call

``` r
compute_features_modular <- function(
    counts, 
    power,                  
    modules = NULL,    
    ## tunable parameters
    minModuleSize,         
    mergeCutHeight,      
    deepSplit,              
) {
  
  ## Just preprocessing (IGNORE)
  rownames(counts) <- gsub("-", ".", rownames(counts)) 
  datExpr <- t(counts) 
  cor <- WGCNA::cor
  
  # TRAINING MODE
  
  if (is.null(modules)) {
    
    
    # -------------------- REPLACE THIS BLOCK --------------------
    net <- WGCNA::blockwiseModules(
      datExpr,
      power = power,
      minModuleSize = minModuleSize,
      mergeCutHeight = mergeCutHeight,
      deepSplit = deepSplit
    )
    modules <- net$colors
    names(modules) <- colnames(datExpr)
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
  }
  
  # PROJECT MODE
  
  
  # -------------------- REPLACE THIS BLOCK --------------------
  module_features <- sapply(unique(modules), function(mod) {
    genes <- names(modules[modules == mod])
    pc <- prcomp(datExpr[, genes, drop = FALSE])
    pc$x[, 1]   
  })
  # -------------------- REPLACE THIS BLOCK --------------------
  
  
  ## Just formatting (IGNORE)
  colnames(module_features) <- paste0("Module_", seq_len(ncol(module_features)))
  
  return(list(features = as.matrix(module_features), structure = modules))
}
```

Then we will use this version of `prepare_custom_folds()` that have some
modifications to account for parameters combinations

Notice we have an additional parameter `ncores` that will control
paralellization for evaluating all your runs (don’t remove it!)

``` r
prepare_custom_folds_tuning <- function(data, folds = NULL, bestune = NULL, ncores = NULL, ...){
  
  if (!is.null(bestune)) {
    
    obs_train <- data$target
    data$target <- NULL
    
    # -------------------- REPLACE THIS BLOCK --------------------
    required_cols <- c()# list your params separated by a comma
    # -------------------- REPLACE THIS BLOCK --------------------
    
    best_params <- if (is.data.frame(bestune)) {
      if (all(required_cols %in% names(bestune))) {
        dplyr::select(bestune, dplyr::all_of(required_cols))
      }else{stop("Not all tunable params found. Verify your function.")}
    } else if (is.list(bestune)) {
      if (all(required_cols %in% names(bestune))) {
        tibble::as_tibble(bestune[required_cols])
      }else{stop("Not all tunable params found. Verify your function.")}
    } else {
      stop("`bestune` must be a data.frame or list.")
    }
    
    
    # -------------------- REPLACE THIS BLOCK --------------------
    res_final <- compute_features_modular(
      counts = data,  # or t(data) if necessary
      # here you will put your tunable parameters (just change param1, param2, ... for the names of your parameters)
      param1 = best_params$param1,
      param2 = best_params$param2,
      param3 = best_params$param3
      ...
    )
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
    train_features_final = res_final$features  ## if necessary
    train_features_final$target <- obs_train
    
    custom_output <- res_final
    
    return(list(train_features_final, custom_output, bestune))
    
  } else {

    # -------------------- REPLACE THIS BLOCK --------------------    
    custom_grid <- expand.grid(
      param1 = param1,
      param2 = param2,
      param3 = param3,
      ...,
      stringsAsFactors = FALSE
    )
    # -------------------- REPLACE THIS BLOCK --------------------
    
    if (is.null(ncores)) ncores <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    
    processed_folds <- foreach::foreach(i = seq_along(folds), .packages = c("dplyr"),
                                        .export = c("compute_features_modular")) %dopar% {  
                                          
      train_idx <- folds[[i]]
      test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
      
      train_data <- data[train_idx, , drop = FALSE]
      obs_train <- train_data$target
      train_data$target <- NULL
      
      fold_results <- lapply(seq_len(nrow(param_grid)), function(j) {
        params <- param_grid[j, ]
        
        
        # -------------------- REPLACE THIS BLOCK --------------------
        res_train <- compute_features_modular(
          counts = train_data,  # or t(train_data) if necessary
          param1 = params$param1,
          param2 = params$param2,
          param3 = params$param3,
          ...
        )
        # -------------------- REPLACE THIS BLOCK --------------------
        
        
        train_features <- as.data.frame(res_train$features) #if necessary
        train_features$target <- obs_train
        
        test_data <- data[test_idx, , drop = FALSE]
        obs_test <- test_data$target
        test_data$target <- NULL
        
        
        # -------------------- REPLACE THIS BLOCK --------------------
        test_features <- compute_features_modular(
          test_data, # or t(test_data) if necessary
          structure = res_train$structure,  
          ...
        )
        # -------------------- REPLACE THIS BLOCK --------------------
        
        
        test_features = test_features$features ## if necessary

        list(
          train_data = train_features,
          test_data  = test_features,
          obs_test   = obs_test,
          rowIndex   = test_idx,
          fold_name  = names(folds)[i],
          params     = params
        )
      })
      
      filename <- file.path("Results", paste0("fold_", names(folds)[i], ".rds"))
      saveRDS(processed_folds[[i]], file = filename)
      
      fold_results
      
    }
    
    parallel::stopCluster(cl)
    unregister_dopar()
    
  }
}
```

In our case it would be:

``` r
prepare_WGCNA_folds_modular <- function(
    data,
    folds = NULL,
    bestune = NULL,
    power,
    ncores = NULL,
    ### tunable parameters
    minModuleSize,
    mergeCutHeight,
    deepSplit
) {
  
  if (!is.null(bestune)) {
    
    obs_train <- data$target
    data$target <- NULL
    
    
    # -------------------- REPLACE THIS BLOCK --------------------
    required_cols <- c("minModuleSize", "mergeCutHeight", "deepSplit")
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
    best_params <- if (is.data.frame(bestune)) {
      if (all(required_cols %in% names(bestune))) {
        dplyr::select(bestune, dplyr::all_of(required_cols))
      }else{stop("Not all tunable params found. Verify your function.")}
    } else if (is.list(bestune)) {
      if (all(required_cols %in% names(bestune))) {
        tibble::as_tibble(bestune[required_cols])
      }else{stop("Not all tunable params found. Verify your function.")}
    } else {
      stop("`bestune` must be a data.frame or list.")
    }
    
    
    # -------------------- REPLACE THIS BLOCK --------------------
    res_final <- compute_features_modular(
      counts = t(data),
      power = power,
      ## tunable parameters
      minModuleSize = best_params$minModuleSize,
      mergeCutHeight = best_params$mergeCutHeight,
      deepSplit = best_params$deepSplit
    )
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
    train_cell_data_final <- as.data.frame(res_final$features)
    train_cell_data_final$target <- obs_train
    
    custom_output <- res_final
    
    return(list(train_cell_data_final, custom_output, best_params))
    
  } else {
    
    
    # -------------------- REPLACE THIS BLOCK --------------------
    custom_grid <- expand.grid(
      power = power,
      minModuleSize = minModuleSize,
      mergeCutHeight = mergeCutHeight,
      deepSplit = deepSplit,
      stringsAsFactors = FALSE
    )
    # -------------------- REPLACE THIS BLOCK --------------------
    
    
    if (is.null(ncores)) ncores <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    
    processed_folds <- foreach::foreach(i = seq_along(folds), .packages = c("dplyr"),
                                          .export = c("compute_features_modular")) %dopar% {
      
      train_idx <- folds[[i]]
      test_idx <- setdiff(seq_len(nrow(data)), train_idx)
      
      train_data <- data[train_idx, , drop = FALSE]
      obs_train <- train_data$target
      train_data$target <- NULL
      
      fold_results <- lapply(seq_len(nrow(custom_grid)), function(j) {

        params <- custom_grid[j, ]
        
        
        # -------------------- REPLACE THIS BLOCK --------------------
        wgcna_train <- compute_features_modular(
          counts = t(train_data),
          power = power,
          ## tunable parameters
          minModuleSize = params$minModuleSize,
          mergeCutHeight = params$mergeCutHeight,
          deepSplit = params$deepSplit
        )
        # -------------------- REPLACE THIS BLOCK --------------------
        
        
        train_features <- as.data.frame(wgcna_train$features) # if necessary
        train_features$target <- obs_train
        
        test_data <- data[test_idx, , drop = FALSE]
        obs_test <- test_data$target
        test_data$target <- NULL
      
        
        # -------------------- REPLACE THIS BLOCK --------------------
        wgcna_test <- compute_features_modular(
          counts = t(test_data),
          modules = wgcna_train$modules
        )
        # -------------------- REPLACE THIS BLOCK --------------------
        
        
        test_features <- as.data.frame(wgcna_test$features) ) # if necessary
        
        list(
          train_data = train_features,
          test_data = test_features,
          obs_test = obs_test,
          rowIndex = test_idx,
          fold_name = names(folds)[i],
          params = params
        )
      })
      
      filename <- file.path("Results", paste0("fold_", names(folds)[i], ".rds"))
      saveRDS(fold_results, file = filename)
      
      fold_results
    }
    
    parallel::stopCluster(cl)
    unregister_dopar() 
    
  }
}
```

To enable this functionality, the user must specify the argument
`fold_construction_args_tunable`, which contains the set of parameter
values to be evaluated

Note the argument `fold_construction_args_fixed` corresponds to the
additional parameters passed in `prepare_custom_folds()` function that
are not tunable (in our case `power` and `ncores`). This means these
parameters are going to have the same value across all runs (**fix**).
Be sure of setting `ncores` to a value that your computer can handle to
avoid crashing.

The tunable parameters define the search space explored during feature
computation. The total number of configurations corresponds to all
combinations of the provided parameter values. In this example:

- minModuleSize = c(20, 50, 100) → 3 values
- mergeCutHeight = c(0.15, 0.25, 0.35) → 3 values
- deepSplit = c(1, 2, 3) → 3 values

This results in **3 × 3 × 3 = 27** feature parameter combinations.

For each of these 27 configurations, the machine learning models are
trained and tuned. For example, if logistic regression with elastic net
(`glmnet`) is used, the model internally evaluates different values of
the hyperparameters `alpha` and `lambda`. If the model tests **10 alpha
values and 20 lambda values**, this results in **200 model
configurations** for each feature combination.

Therefore, the total number of evaluated models becomes:

**27 (feature configurations) × 200 (model hyperparameter combinations)
= 5400** model fits, which are further multiplied by the number of folds
and repetitions used in cross-validation.

``` r
res_params <- compute_features.training.ML(features_train = t(counts_train), 
                                          target_var     = coldata_train$Response,
                                          task_type = "classification",
                                          trait.positive = "R",
                                          metric = "AUROC",
                                          k_folds = 2,
                                          n_rep = 1,
                                          return = FALSE,
                                          fold_construction_fun = prepare_WGCNA_folds_modular,
                                          fold_construction_args_fixed = list(power=6, ncores = 2),
                                          fold_construction_args_tunable = list(
                                            minModuleSize = c(20, 50, 100),       
                                            mergeCutHeight = c(0.15, 0.25, 0.35),
                                            deepSplit = c(1, 2, 3)           
                                          ))
```

Step 4 – Prediction on test data

To apply the model to new data, compute the same type of features using
the learned modules from the training set.

``` r
test = compute_features_modular(counts_test, modules = res_custom$Custom_output$modules)
test_features = test$features
```

Prediction

``` r
pred <- compute_prediction(res_custom$Model,
                           test_features,
                           coldata_test$Response,
                           trait.positive = "R",
                           return = FALSE)
```

Why do we need `bestune` argument even if I don’t have hyperparams in my
function?

For this to work, your custom function must accept a `bestune` argument,
which is used internally to inject the optimized parameter values found
during the tuning step.

When `bestune` is NULL, the function assumes that tuning has not yet
been performed. A grid of candidate parameter values (defined in
`fold_construction_args_tunable`) is generated. For each fold, the
function iterates through all combinations of parameter values and
recomputes the features. This exploration step is parallelized across
folds using foreach and doParallel, allowing multiple folds to be
processed simultaneously. Parallelization reduces runtime considerably
when the parameter grid or number of folds is large.

When `bestune` is not NULL, it means the tuning process has already been
completed. The optimized parameter values are extracted from the bestune
object. Features are then recomputed once on the full training dataset
using these tuned parameters. This ensures the final model is trained
with the best parameter setting identified during cross-validation.

In summary, the `bestune` argument acts as a control switch:

- NULL → run parameter search with parallelized fold-level evaluation.
- non-NULL → lock in the tuned parameter values and rebuild the features
  for final training.

This design allows a single custom fold-construction function to handle
both hyperparameter tuning (exploration, parallelized) and final model
preparation (exploitation, single optimized run).

**NOTE:** `pipeML` is built on top of existing frameworks and makes
extensive use of the R packages `caret`, `parnsip`, `tidymodels`,
`parsnip` and `censored`. If you use `pipeML` in your work, please cite
our package along with these foundational packages.

## References
