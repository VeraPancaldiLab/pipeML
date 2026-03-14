#' data_example_classification
#'
#' Example dataset for classification tasks.
#' Derived from the Breast Cancer Wisconsin dataset (`mlbench::BreastCancer`), cleaned by removing rows with missing values.
#'
#' @format A data frame with samples as rows and features as columns (e.g., cell or tissue measurements).
#'
#' @source `mlbench::BreastCancer` dataset
#'
#' @examples
#' data(data_example_classification)
#' head(data_example_classification)
"data_example_classification"

#' data_example_survival
#'
#' Example dataset for survival analysis.
#' Uses the lung cancer dataset from the `survival` package.
#'
#' @format A data frame with samples as rows and variables as columns (e.g., survival time, status, covariates)
#'
#' @source `survival::lung`
#'
#' @examples
#' data(data_example_survival)
#' head(data_example_survival)
"data_example_survival"

#' counts_example
#'
#' Gene expression matrix from the Hugo et al. (2016) metastatic melanoma cohort.
#' Rows correspond to genes (HUGO gene symbols) and columns to patient samples.
#' This dataset is used to compute features for training machine learning models.
#'
#' @format A numeric matrix with genes as rows and samples as columns.
#'
#' @source Hugo W., Zaretsky J.M., Sun L., Song C., Moreno B.H., Hu-Lieskovan S., Berent-Maoz B., Pang J., Chmielowski B., Cherry G., et al.
#' Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell. 2016;165:35–44. doi: 10.1016/j.cell.2016.02.065
#'
#' @examples
#' data(counts_example)
#' head(counts_example)
"counts_example"

#' coldata_example
#'
#' Metadata for the Hugo et al. (2016) cohort, containing the response to anti-PD-1 therapy.
#' Each row corresponds to a patient/sample and columns describe clinical or response information.
#' This dataset is used as the target variable for supervised learning tasks.
#'
#' @format A data frame with samples as rows and a column `Response` indicating treatment outcome.
#'
#' @source Hugo W., Zaretsky J.M., Sun L., Song C., Moreno B.H., Hu-Lieskovan S., Berent-Maoz B., Pang J., Chmielowski B., Cherry G., et al.
#' Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell. 2016;165:35–44. doi: 10.1016/j.cell.2016.02.065
#'
#' @examples
#' data(coldata_example)
#' head(coldata_example)
"coldata_example"
