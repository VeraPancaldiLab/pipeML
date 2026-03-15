# Internal: Plot Pooled AUROC and AUPRC Performance Curves

Internal function to read multiple `.rds` files containing machine
learning results, pool the AUROC and AUPRC metrics, and generate
boxplots summarizing performance across iterations. Median values are
annotated on the plots.

## Usage

``` r
get_pooled_roc_curves(file.name, folder_path)
```

## Arguments

- file.name:

  Character. Name used as a prefix when saving output plots.

- folder_path:

  Character. Path to the directory containing the `.rds` files with ML
  model results.

## Details

Each `.rds` file is expected to contain a list with a `result$AUC`
element, including both `AUROC` and `AUPRC` values. The function saves
two PDF files in the `Results/` directory:

- Boxplot of AUROC values with median annotation

- Boxplot of AUPRC values with median annotation

No value is returned to the R environment.
