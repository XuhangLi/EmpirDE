#' Example Statistic Table
#'
#' A DataFrame object containing an example test statistic matrix.
#' This dataset includes 10 genes across 300 conditions.
#'
#' @format A data frame with 10 rows (genes) and 300 columns (conditions).
#' @usage data("example_stat_table")
# "example_stat_table"

#' Example WPS dataset
#'
#' An example WPS dataset containing two objects: (1) a gene-by-sample countTable that includes 10 genes across 301 conditions based on a real WPS data; (2) a meta-data table for the 301 samples.
#' The values in count table are raw read counts from WPS experiments.
#' Meta-data table must contain the following columns:
#' \code{sampleID} column that corresponds to the column names in \code{countTable}.
#' \code{covTreatment} column that is the covariate for the treatment (condition) and must have a 'control' condition.
#' \code{covBatch} column that is the covariate for batches, which by default is the replicate batch.
#' \code{libID} column that is ID of sequencing libraries. WPS DE analysis is first conducted at individual sequencing library level thus this ID is used to match samples pooled in the same library.
#' \code{plate2} column that is a logical value indicating if the sample is cultured in the second 96-well plate. This column is only needed for WPS experiments that involves plate 2 confounding effects.
#'
#' @format A gene-by-sample data frame for read counts and a data frame documenting the metadata information.
#' @usage data("WPS_example_data")
# "WPS_example_data"
