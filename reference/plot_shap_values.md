# Plot SHAP values

Plots the variable importance based on SHAP (SHapley Additive
exPlanations) values and saves the plot to the `Results/` directory.

## Usage

``` r
plot_shap_values(shap_df, ml_model, file_name, width = 10, height = 10)
```

## Arguments

- shap_df:

  A data frame or matrix containing SHAP values where rows represent
  samples and columns represent features.

- ml_model:

  A character string representing the name of the machine learning
  model.

- file_name:

  A character string representing the name of the file to save the plot
  in the `Results/` directory.

- width:

  A numeric value specifying the width of the plot in inches (default is
  10).

- height:

  A numeric value specifying the height of the plot in inches (default
  is 10).

## Value

A plot saved as a PDF file in the `Results/` directory showing the
variable importance of the machine learning model.

## Details

This function generates a bar plot of the SHAP values, where the
features are sorted by their mean SHAP value. The plot distinguishes
between features that increase the predicted outcome (colored in red)
and those that decrease the predicted outcome (colored in blue). The
plot is saved as a PDF file in the `Results/` directory, with the
filename specified by the user.
