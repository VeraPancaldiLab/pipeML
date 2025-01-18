# pipeML
A robust R machine learning pipeline for classification tasks and survival analysis.

### Install from source
Clone the repository:
```
git clone https://github.com/mhurtado13/pipeML
```

## Environment

pipeML uses [renv](https://rstudio.github.io/renv/index.html) for creating a reproducible r-environment that can be easily share between users, this, to allow reproducibility and avoid libraries conflicts. Setting it up will install the neccessary packages, along with their specific versions in an isolated environment. 

For this, open the project `pipeML.Rproj` and in the R console run:

```r
# To activate the R environment (if you are using it for the first time)
renv::activate()
# To download and install all the require libraries and packages (if you are using it for the first time)
renv::restore() 
```

Once all packages have been installed, you can start creating your own scripts but be sure to still be inside the .Rproj!

Note that this is a once-step only when running multideconv for the first time. For the following times, you will only need to open the `pipeML.Rproj` and you are ready to go!

If you are doing other analyses that require the installation of extra libraries not present in the environment, you can install them as usual but after that, make sure to run `renv::snapshot()` to update your environment.

Make sure to run `renv::deactivate()` when finishing running multideconv, to avoid conflicts whenever you start a different project.

For more information about how `renv` works, visit https://rstudio.github.io/renv/articles/renv.html.

## ML Pipeline

pipeML is a robust pipeline done in R to train and test machine learning models for classification tasks. It is develop for a fast and easy use, its complexity allows a robust implementation during each step of the pipeline (Figure 1). 

<p align="center">
 <img src="man/MLpipeline.png?raw=true" />
</p>

<p align="center"><i>
   Figure 1. Machine learning pipeline.
</i></p>

## Key features

* Stratified data split
* Iterative Boruta algorithm for feature selection
* Repeated k-fold cross validation (kCV)
* Hyperparameter tuning based on AUC or Accuracy
* Stratified k-fold construction
* Model stacking
* 13 Machine Learning methods implemented
  * Bagged CART
  * Random Forest (RF)
  * C50
  * Logistic regression (LG)
  * CART
  * Naive Bayes (NB)
  * Regularized Lasso
  * Ridge regression
  * Linear Discriminant Analysis (LDA)
  * Regularized Logistic Regression (Elastic net)
  * K-nearest neighbors (KNN)
  * Support vector machine with radial kernel (SVMr)
  * Support vector machine with linear kernel (SVMl)

## General usage

The main function to compute the ML pipeline is through the `compute.ML` function. Refer to `examples/mlpipeline.Rmd` to see a detailed explanation of the parameters inside the function. 

```r
res_ml = compute.ML(raw.counts, normalized = T, clinical, trait = "Response",trait.positive = "CR", partition = 0.8, metric = "AUC", stack = T, feature.selection = F,seed = 1234, doParallel = T,  workers = 2, file_name = "Test", return = T)
```

The previous function is a one-partition run meaning it will partitionate your data into train and test once. For reproducibility, a seed needs to be specified. For a robust analysis we recommend testing performance from different partitions of train/test. For doing this users can use the `compute.bootstrap.ML` function:

```r
compute.bootstrap.ML(raw.counts, normalized = T, clinical, trait = "Response", trait.positive = "YES", partition = 0.8, metric = "Accuracy", iterations = 20, feature.selection = F, stack = T, workers = 4, file.name = "Test", return = F)
```

ML will be saved in your working directory inside a `Results` folder. In order to pool the performance across all models we provide the function `get_pooled_roc_curves` that allows you to plot a boxplot with the AUC-ROC and AUC-PRC scores across all your ML models.

```r
get_pooled_roc_curves("Test", "Results/ML_models") #Get boxplot with AUC scores distribution across iterations
```

pipeML provides also a function to perform the Leaving-One-Dataset-Out (LODO) approach.

```r
res_ml = compute.LODO.ML(raw.counts, normalized = T, clinical, trait = "Response",trait.positive = "R", trait.out = "Cohort", out = "Dupont", metric = "Accuracy", stack = T, feature.selection = F, doParallel = T, workers = 4, file_name = "Test", return = F)
```

## Additional notes
pipeML is fully implemented and adapted as part of [CellTFusion](https://github.com/VeraPancaldiLab/CellTFusion) R package, an integration tool for cell deconvolution couple with transcription factor activities to deconvolve cell states of the tumor microenvironment. If you would like to use it in this context we invite you to visit our github repository.

## Citing pipeML

If you use pipeML in a scientific publication, we would appreciate citation to the :

```
XXXXX
```
