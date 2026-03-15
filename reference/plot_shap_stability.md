# Plot SHAP Feature Importance Stability Across Resamples

This function visualizes the stability of SHAP feature importance values
across cross-validation resamples. It computes the mean absolute SHAP
value and standard deviation per feature, then generates a horizontal
bar plot with error bars representing variability.

## Usage

``` r
plot_shap_stability(shap_df, file.name)
```

## Arguments

- shap_df:

  A data frame of SHAP values with samples in rows and features in
  columns.

- file.name:

  Character string specifying the filename prefix for the saved PDF
  plot.

## Value

A horizontal bar plot saved as a PDF in the `Results/` directory. The
plot shows the mean absolute SHAP value per feature and error bars
indicating the standard deviation across resamples.

## Details

The function reshapes the SHAP values into long format, calculates the
mean absolute value and standard deviation per feature, and then creates
a ggplot2 horizontal bar chart. The plot is saved as a PDF with
filename: `"Results/SHAP_stability_resample_<file.name>.pdf"`.
