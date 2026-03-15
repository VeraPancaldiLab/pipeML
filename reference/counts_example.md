# counts_example

Gene expression matrix from the Hugo et al. (2016) metastatic melanoma
cohort. Rows correspond to genes (HUGO gene symbols) and columns to
patient samples. This dataset is used to compute features for training
machine learning models.

## Usage

``` r
counts_example
```

## Format

A numeric matrix with genes as rows and samples as columns.

## Source

Hugo W., Zaretsky J.M., Sun L., Song C., Moreno B.H., Hu-Lieskovan S.,
Berent-Maoz B., Pang J., Chmielowski B., Cherry G., et al. Genomic and
Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic
Melanoma. Cell. 2016;165:35–44. doi: 10.1016/j.cell.2016.02.065

## Examples

``` r
data(counts_example)
head(counts_example)
#>         Sample_Pt1.baseline Sample_Pt10.baseline Sample_Pt12.baseline
#> A1BG             5.97269265           4.98093927            5.3008559
#> A1CF             0.00000000           0.00000000            0.0000000
#> A2M              8.29746657           9.51652738            9.5387321
#> A2ML1            0.02856915           0.05658353            0.0000000
#> A3GALT2          0.00000000           0.16349873            0.0976108
#> A4GALT           3.69821848           5.39025496            1.8599695
#>         Sample_Pt13.baseline Sample_Pt14.baseline Sample_Pt15.baseline
#> A1BG              5.67440415           5.16269343             6.426265
#> A1CF              0.00000000           0.01435529             0.000000
#> A2M               9.77865212           6.53480866             8.465526
#> A2ML1             0.00000000           0.00000000             0.000000
#> A3GALT2           0.07038933           0.11103131             0.000000
#> A4GALT            3.66220550           2.80115866             1.097611
#>         Sample_Pt19.baseline Sample_Pt2.baseline Sample_Pt20.baseline
#> A1BG                5.176323          9.07860450           5.65363331
#> A1CF                0.000000          4.09508049           0.01435529
#> A2M                 5.422906          9.00419237           7.88472001
#> A2ML1               0.000000          0.15055968           0.07038933
#> A3GALT2             0.000000          0.04264434           0.00000000
#> A4GALT              2.204767          0.54596837           1.49057013
#>         Sample_Pt22.baseline Sample_Pt23.baseline Sample_Pt25.baseline
#> A1BG              4.38612116           4.98185265           6.40922120
#> A1CF              0.00000000           0.00000000           0.00000000
#> A2M               7.53294029           8.55754029           8.89050749
#> A2ML1             0.13750352           0.01435529           0.01435529
#> A3GALT2           0.02856915           0.00000000           0.00000000
#> A4GALT            3.09423607           2.15704371           1.55090066
#>         Sample_Pt27.baseline Sample_Pt28.baseline Sample_Pt29.baseline
#> A1BG              4.94953493           3.50716035           4.89917563
#> A1CF              0.00000000           0.28688115           0.04264434
#> A2M               9.13857900           6.39420543           9.73172632
#> A2ML1             0.17632277           6.94801674           0.12432814
#> A3GALT2           0.08406426           0.05658353           0.11103131
#> A4GALT            1.69599381           2.27202319           1.60880924
#>         Sample_Pt31.baseline Sample_Pt32.baseline Sample_Pt35.baseline
#> A1BG              5.73226920            6.0740344             6.046797
#> A1CF              0.00000000            0.0000000             0.000000
#> A2M               6.12163712            8.1092561             6.921008
#> A2ML1             0.01435529            0.3219281             0.000000
#> A3GALT2           0.11103131            0.1110313             0.000000
#> A4GALT            3.18428029            5.0373822             1.490570
#>         Sample_Pt37.baseline Sample_Pt38.baseline Sample_Pt4.baseline
#> A1BG              5.73741637            4.4093909          4.90689060
#> A1CF              0.00000000            0.0000000          0.00000000
#> A2M               6.16611286            9.9298945          6.05701697
#> A2ML1             0.00000000            0.0000000          0.02856915
#> A3GALT2           0.04264434            0.1634987          0.02856915
#> A4GALT            0.91073266            3.1358632          0.69599381
#>         Sample_Pt5.baseline Sample_Pt6.baseline Sample_Pt7.baseline
#> A1BG               5.977967          3.25701062          6.34677964
#> A1CF               0.000000          0.01435529          0.00000000
#> A2M                6.049631          6.80374367          6.61117238
#> A2ML1              0.000000          0.00000000          0.01435529
#> A3GALT2            0.000000          0.00000000          0.00000000
#> A4GALT             1.049631          1.07724300          1.15055968
#>         Sample_Pt8.baseline Sample_Pt9.baseline
#> A1BG              5.2562559           5.2452675
#> A1CF              0.0000000           0.0000000
#> A2M               7.5156998           9.0395771
#> A2ML1             0.2141248           0.6040713
#> A3GALT2           0.0000000           0.0000000
#> A4GALT            2.5631581           2.5753123
```
