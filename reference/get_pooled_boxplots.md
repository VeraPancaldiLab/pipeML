# Internal: Plot Pooled AUROC and AUPRC Boxplots Across Multiple Folders

Internal function to aggregate AUROC and AUPRC metrics from multiple
folders (e.g., different cohorts or models), and generate comparative
boxplots showing model performance across groups.

## Usage

``` r
get_pooled_boxplots(folder_paths, file_name, width = 12, height = 8)
```

## Arguments

- folder_paths:

  Character vector. Paths to folders containing `.rds` files with ML
  model results.

- file_name:

  Character. Prefix used when saving the resulting PDF plots.

- width:

  Numeric. Width of the saved plots in inches. Default is 12.

- height:

  Numeric. Height of the saved plots in inches. Default is 8.

## Details

Each `.rds` file should contain a list with a `result$AUC` element
including numeric values for both `AUROC` and `AUPRC`. Folder names are
used as grouping labels in the plots. Red dashed horizontal lines are
drawn at a reference value (0.7) for visual interpretation. Two PDF
files are saved in the `Results/` directory:

- `Boxplots_AUROC_performance_<file_name>.pdf`

- `Boxplots_AUPRC_performance_<file_name>.pdf`

No object is returned to the R environment.
