# Plot Pooled AUROC and AUPRC Boxplots Across Multiple Folders

This function aggregates AUROC and AUPRC metrics from multiple folders
(typically corresponding to different cohorts or models), and generates
boxplots comparing model performance across groups.

## Usage

``` r
get_pooled_boxplots(folder_paths, file_name, width = 12, height = 8)
```

## Arguments

- folder_paths:

  Character vector. Paths to folders containing `.rds` files with ML
  model performance results.

- file_name:

  Character. Prefix for the saved PDF files containing the plots.

- width:

  Numeric. Width of the saved plots in inches. Default is 12.

- height:

  Numeric. Height of the saved plots in inches. Default is 8.

## Value

Saves two PDF files in the `Results/` directory:

- `Boxplots_AUROC_performance_<file_name>.pdf`

- `Boxplots_AUPRC_performance_<file_name>.pdf`

No object is returned to the R environment.

## Details

Each `.rds` file is expected to contain a list object with a
`result$AUC` element, which includes numeric values for both `AUROC` and
`AUPRC`. Folder names are used as grouping labels in the plots.

Red dashed lines are drawn at a fixed reference value (e.g., 0.7) for
visual interpretation.
