#' data_example_classification
#'
#' Example dataset for classification tasks.
#' Derived from the Breast Cancer Wisconsin dataset (`mlbench::BreastCancer`), cleaned by removing rows with missing values.
#'
#' @format A data frame with 683 samples as rows and 10 columns: 9 cytological features (`Cl.thickness`,
#'   `Cell.size`, `Cell.shape`, `Marg.adhesion`, `Epith.c.size`, `Bare.nuclei`, `Bl.cromatin`, `Normal.nucleoli`,
#'   `Mitoses`) and the outcome `target` (1 = malignant, 0 = benign).
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
#' Uses the lung cancer dataset from the `survival` package (complete cases only).
#'
#' @format A data frame with 167 samples as rows and 10 columns: the survival time in days (`time`), the event
#'   indicator (`status`: 1 = death, 0 = censored) and 8 covariates.
#'
#' @source `survival::lung`, with `status` recoded from 1 = censored / 2 = dead to 0 = censored / 1 = death.
#'
#' @examples
#' data(data_example_survival)
#' head(data_example_survival)
"data_example_survival"

#' counts_example
#'
#' Gene expression matrix of melanoma samples, including the Gide et al. (2019) metastatic melanoma cohort.
#' Rows correspond to genes (HUGO gene symbols) and columns to patient samples.
#' This dataset is used to compute features for training machine learning models.
#'
#' @format A numeric matrix with 5000 genes as rows and 336 samples as columns, with normalized (non-integer)
#'   expression values. The 73 samples of the Gide cohort are the ones annotated in `coldata_example`
#'   (`counts_example[, rownames(coldata_example)]`).
#'
#' @source Gide T.N., Quek C., Menzies A.M., Tasker A.T., Shang P., Holst J., Madore J., Lim S.Y., Velickovic R., Wongchenko M., et al. Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy. Cancer Cell. 2019;35:238–255. doi: 10.1016/j.ccell.2019.01.003.
#'
#' @examples
#' data(counts_example)
#' head(counts_example)
"counts_example"

#' coldata_example
#'
#' Metadata for the Gide et al. (2019) cohort, containing the response to anti-PD-1 therapy.
#' Each row corresponds to a patient/sample and columns describe clinical or response information.
#' This dataset is used as the target variable for supervised learning tasks.
#'
#' @format A data frame with 73 samples as rows (row names: sample identifiers, matching column names of
#'   `counts_example`) and 2 columns: `Response`, the treatment outcome (`"R"` = responder, `"NR"` = non-responder),
#'   and `Cohort` (`"Gide"`).
#'
#' @source Gide T.N., Quek C., Menzies A.M., Tasker A.T., Shang P., Holst J., Madore J., Lim S.Y., Velickovic R., Wongchenko M., et al. Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy. Cancer Cell. 2019;35:238–255. doi: 10.1016/j.ccell.2019.01.003.
#'
#' @examples
#' data(coldata_example)
#' head(coldata_example)
"coldata_example"
