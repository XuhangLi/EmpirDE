#' Example Statistic Table
#'
#' A DataFrame object containing an example test statistic matrix. This dataset includes 10 genes across 300 conditions.
#' @format A data frame with 10 rows (genes) and 300 columns (conditions).
#' @return The corresponding DataFrame object
#' @usage data("example_stat_table")
"example_stat_table"

#' An WPS dataset: count table
#'
#' A gene-by-sample count table that includes 10 genes across 301 conditions based on a real WPS data
#' @format A gene-by-sample dataframe of read counts
#' @return The corresponding dataframe object
#' @usage data("countTable")
"countTable"

#' An WPS dataset: metadata table
#'
#' A meta-data table for the 301 samples in the example count table. WPS meta-data table should contain the following columns:
#' \describe{
#'  \item{\code{sampleID}}{column that corresponds to the column names in \code{countTable}.}
#'  \item{\code{covTreatment}}{column that is the covariate for the treatment (condition) and must have a 'control' condition.}
#'  \item{\code{covBatch}}{column that is the covariate for batches, which by default is the replicate batch.}
#'  \item{\code{libID}}{column that is ID of sequencing libraries. WPS DE analysis is first conducted at individual sequencing library level thus this ID is used to match samples pooled in the same library.}
#'  \item{\code{plate2}}{column that is a logical value indicating if the sample is cultured in the second 96-well plate. This column is only needed for WPS experiments that involves plate 2 confounding effects.}
#' }
#' @format A dataframe object for metadata information.
#' @return The corresponding dataframe object
#' @usage data("metaDataTable")
"metaDataTable"

#' An example formated DESeqDataSet object for testing control-independent DE analysis
#'
#' To run control-independent DE analysis individually, one should format an input DESeqDataSet with covariates included in the example:
#' \describe{
#'  \item{\code{RNAi}}{covariate describing treatment (condition) of samples. A 'control' condition is not required but this covariate should be the last term in the design formula and is being tested for by default.}
#'  \item{\code{batchLabel}}{(optional) covariate describing experimental batches and can be skipped.}
#'  \item{\code{plate2}}{(optional) covariate indicating if the sample is cultured in the second 96-well plate. Only for WPS experiments that involves plate 2 confounding effects.}
#' }
#' @format An example DESeqDataSet object
#' @return The corresponding DESeqDataSet object
#' @usage data("example_dds")
"example_dds"
