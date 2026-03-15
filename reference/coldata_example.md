# coldata_example

Metadata for the Hugo et al. (2016) cohort, containing the response to
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

Hugo W., Zaretsky J.M., Sun L., Song C., Moreno B.H., Hu-Lieskovan S.,
Berent-Maoz B., Pang J., Chmielowski B., Cherry G., et al. Genomic and
Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic
Melanoma. Cell. 2016;165:35–44. doi: 10.1016/j.cell.2016.02.065

## Examples

``` r
data(coldata_example)
head(coldata_example)
#>                      Response Cohort
#> Sample_Pt1.baseline        NR   Hugo
#> Sample_Pt10.baseline       NR   Hugo
#> Sample_Pt12.baseline       NR   Hugo
#> Sample_Pt13.baseline        R   Hugo
#> Sample_Pt14.baseline       NR   Hugo
#> Sample_Pt15.baseline        R   Hugo
```
