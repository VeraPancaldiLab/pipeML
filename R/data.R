#' Deconvolution matrix
#'
#' A matrix with the cell type deconvolution features to use as train data from the bulk RNAseq of Mariathasan et al. (2018)
#'
#' @format Matrix with samples as rows and deconvolution features as columns
#'
#' @examples
#' data(deconvolution)
#' head(deconvolution)
"deconvolution"

#' Raw counts
#'
#' Raw gene expression matrix from bulk RNAseq.
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(raw.counts)
#' head(raw.counts)
"raw.counts"

#' ML model trained
#'
#' Machine learning model trained
#'
#' @format Output from compute.features.ML()
#'
#' @examples
#' data(model_trained)
#' head(model_trained)
"model_trained"

#' Clinical data
#'
#' Metadata containing the variable to predict
#'
#' @format Matrix with samples as rows and target variable as column
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(traitData)
#' head(traitData)
"traitData"
