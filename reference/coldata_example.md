# coldata_example

Metadata for the Gide et al. (2019) cohort, containing the response to
anti-PD-1 therapy. Each row corresponds to a patient/sample and columns
describe clinical or response information. This dataset is used as the
target variable for supervised learning tasks.

## Usage

``` r
coldata_example
```

## Format

A data frame with samples as rows and a column `Response` indicating
treatment outcome.

## Source

Gide T.N., Quek C., Menzies A.M., Tasker A.T., Shang P., Holst J.,
Madore J., Lim S.Y., Velickovic R., Wongchenko M., et al. Distinct
Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and
Anti-PD-1/Anti-CTLA-4 Combined Therapy. Cancer Cell. 2019;35:238–255.
doi: 10.1016/j.ccell.2019.01.003.

## Examples

``` r
data(coldata_example)
head(coldata_example)
#>                      Response Cohort
#> Sample_10_PD1_PRE          NR   Gide
#> Sample_10_ipiPD1_PRE       NR   Gide
#> Sample_11_PD1_PRE          NR   Gide
#> Sample_12_PD1_PRE          NR   Gide
#> Sample_12_ipiPD1_PRE       NR   Gide
#> Sample_13_PD1_PRE          NR   Gide
```
