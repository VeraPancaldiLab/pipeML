# Plot Pooled AUROC and AUPRC Performance Curves

This function reads multiple `.rds` files containing machine learning
results, pools the AUROC and AUPRC metrics, and generates boxplots
summarizing performance across iterations. Median values are annotated
on the plots.

## Usage

``` r
get_pooled_roc_curves(file.name, folder_path)
```

## Arguments

- file.name:

  Character. Name to use when saving the plots (used as a prefix in the
  output file names).

- folder_path:

  Character. Path to the directory containing the `.rds` files with ML
  model results.

## Value

Saves two PDF files in the `Results/` directory:

- Boxplot of AUROC values with median annotation

- Boxplot of AUPRC values with median annotation

No value is returned to the R environment.

## Details

Each `.rds` file is expected to contain a list with a `result$AUC`
element that includes both `AUROC` and `AUPRC` values.
