# Clinical data

Metadata containing the variable to predict

## Usage

``` r
traitData
```

## Format

Matrix with samples as rows and target variable as column

## Source

Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501

## Examples

``` r
data(traitData)
head(traitData)
#>                 Best.Confirmed.Overall.Response binaryResponse Enrollment.IC
#> SAM7f0d9cc7f001                              PD          SD/PD           IC1
#> SAM4305ab968b90                              PD          SD/PD           IC0
#> SAMcf018fee2acd                              PD          SD/PD           IC0
#> SAMcc4675f394a1                              PD          SD/PD           IC0
#> SAM49f9b2e57aa5                              PD          SD/PD           IC0
#> SAM2e7aa8fa0ab3                              PD          SD/PD           IC1
#>                 IC.Level TC.Level Immune.phenotype FMOne.mutation.burden.per.MB
#> SAM7f0d9cc7f001      IC1      TC0         excluded                           NA
#> SAM4305ab968b90      IC0      TC0           desert                           12
#> SAMcf018fee2acd      IC0      TC0         excluded                           NA
#> SAMcc4675f394a1      IC0      TC0             <NA>                           NA
#> SAM49f9b2e57aa5      IC0      TC0         excluded                            5
#> SAM2e7aa8fa0ab3      IC1      TC0         excluded                           14
#>                 Sex  Race Intravesical.BCG.administered Baseline.ECOG.Score
#> SAM7f0d9cc7f001   M WHITE                             N                   1
#> SAM4305ab968b90   M WHITE                             N                   1
#> SAMcf018fee2acd   M WHITE                             N                   0
#> SAMcc4675f394a1   M WHITE                             N                   1
#> SAM49f9b2e57aa5   M WHITE                             Y                   1
#> SAM2e7aa8fa0ab3   M WHITE                             N                   1
#>                 Tobacco.Use.History Met.Disease.Status         Sample.age
#> SAM7f0d9cc7f001               NEVER              Liver  more than 2 years
#> SAM4305ab968b90            PREVIOUS              Liver  more than 2 years
#> SAMcf018fee2acd            PREVIOUS           Visceral  more than 2 years
#> SAMcc4675f394a1            PREVIOUS              Liver (less than) 1 year
#> SAM49f9b2e57aa5            PREVIOUS              Liver          1-2 years
#> SAM2e7aa8fa0ab3            PREVIOUS           Visceral (less than) 1 year
#>                  Tissue Received.platinum Sample.collected.pre.platinum
#> SAM7f0d9cc7f001 bladder                 Y                             Y
#> SAM4305ab968b90 bladder                 Y                             N
#> SAMcf018fee2acd    lung                 Y                             Y
#> SAMcc4675f394a1   liver                 Y                             Y
#> SAM49f9b2e57aa5  kidney                 Y                             Y
#> SAM2e7aa8fa0ab3  kidney                 N                          <NA>
#>                 Neoantigen.burden.per.MB sizeFactor ANONPT_ID        os censOS
#> SAM7f0d9cc7f001                       NA  0.6498134     10220  4.632444      1
#> SAM4305ab968b90                0.3921569  1.3451211     10375  1.412731      1
#> SAMcf018fee2acd                       NA  1.7799627     10175 16.229979      1
#> SAMcc4675f394a1                       NA  1.1994070     10142  3.121150      1
#> SAM49f9b2e57aa5                1.4509804  0.9977972     10037  7.457906      1
#> SAM2e7aa8fa0ab3                2.7058824  1.2795319     10172  7.720739      1
#>                  Lund       Lund2 TCGA.Subtype
#> SAM7f0d9cc7f001 MS2b1 Infiltrated           II
#> SAM4305ab968b90  MS1b        UroA            I
#> SAMcf018fee2acd  MS1a        UroA            I
#> SAMcc4675f394a1  MS1a        UroA            I
#> SAM49f9b2e57aa5  MS1b        UroA            I
#> SAM2e7aa8fa0ab3  MS1b        UroA           II
```
