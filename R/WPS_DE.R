#' @title WPS Differential Expression (DE) analysis
#'
#' @description
#' Perform the two-pronged WPS DE analysis with an input dataset consistent with WPS data structure.
#' @param statTbl A gene-by-condition matrix of the DE test statistic (Wald statistic) from DEseq2
#' @param display Whether or not display the fitting of genes with heavy tails (default is FALSE)
#'
#' @return A list of the fitting results:
#' \describe{
#'   \item{\code{p_mat}}{The matrix of empirical p-values (gene-by-conditions)}
#'   \item{\code{nulls}}{A data frame for the fitted mean and standard deviations for each gene}
#'   \item{\code{not_fit}}{Genes that were not fitted because the number of non-NA test statistics is fewer than 100 (theoratical null is used)}
#'   \item{\code{qualityMetrics}}{A list of fitting quality metrics, including the genes whose test statistic distribution is heavy tailed (these tails were trimmed prior to fitting) and number of trimmed conditions for the fitting of each gene}
#' }
#'
#'
#' @export WPS_DE
#'
#' @author Xuhang Li
#' @examples
#' data(WPS_example_data)
#' result <- WPS_DE(countTable, metaDataTable)

# all non-base functions to be called by ::
# enter the project folder
# create function files in R
# run devtools::document() to document the new changes
# then run build - check
# when done, run Git - commit - push

WPS_DE <- function(countTable, metaDataTable) {

  # make sure the counttable and meta data table are aligned
  metaDataTable = metaDataTable[match(metaDataTable$sampleID, colnames(countTable)),]
  if (!identical(metaDataTable$sampleID, colnames(countTable))){
    stop('The meta-data table is not aligned with column names of the count table!')
  }

  # step1: perform the control dependent and independent DE analysis

  cat("\033[31mStarting control-(in)dependent DE analysis per library ...\033[39m\n")
  libs = unique(metaDataTable$libID)
  ctr_dep_DE_res = list()
  for (i in 1:length(libs)){
    # subset to the libary being analyzed
    subsetInd = metaDataTable$libID == libs[i]
    input_subset = countTable[,subsetInd]
    batchLabel_subset = metaDataTable$covBatch[subsetInd]
    batchLabel_subset = as.factor(batchLabel_subset)
    RNAi_subset = factor(metaDataTable$covTreatment[subsetInd], levels = c('control', setdiff(metaDataTable$covTreatment[subsetInd], 'control')))
    # create dds object
    coldata =  data.frame(batchLabel = batchLabel_subset, RNAi = RNAi_subset)
    rownames(coldata) = colnames(input_subset)
    if (any(colnames(metaDataTable) == 'plate2')){
      coldata$plate2 =  metaDataTable$plate2[subsetInd]
      }
    cat("Performing control-dependent DE analysis for library", libs[i],"\n")
    invisible(suppressMessages(dds <- DESeq2::DESeqDataSetFromMatrix(countData = input_subset,
                                                colData = coldata,
                                                design= ~ batchLabel + RNAi)))
    # pre-filtering
    keep <- rowSums(counts(dds)>=10) >= 1
    dds <- dds[keep,]
    # pre-calculating size factors
    dds = estimateSizeFactors(dds)

    # perform control-dependent DE analysis
    invisible(suppressMessages(DE_res_list <- control_dependent_DE(dds)))
    names(DE_res_list) = paste(libs[i], names(DE_res_list),sep = '_')
    ctr_dep_DE_res = c(ctr_dep_DE_res, DE_res_list)

    # perform control-independent DE analysis
    cat("Performing control-independent DE analysis for library", libs[i],"\n")
    cat("Fit main populations ...\n")
    invisible(suppressMessages(zMat <- fit_main_population(dds)))
    cat("Runing DEseq2 again ...\n")
    DE_res_list <- control_independent_DE(dds, zMat)






  }


}
