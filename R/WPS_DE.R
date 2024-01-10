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

WPS_DE <- function(countTable, metaDataTable, params = NULL) {
  # all libs to analyze
  libs = unique(metaDataTable$libID)

  # set parameters
  if (is.null(params)){
    pcutoffs = 0.005
    freqCutoff = length(libs)*0.25
    independentFilteringCutoff = 30
  }else{
    pcutoffs = params$p_out_cutoffs
    freqCutoff = params$core_ctr_outlier_cutoff
    independentFilteringCutoff = params$independentFilteringCutoff
  }

  # make sure the counttable and meta data table are aligned
  metaDataTable = metaDataTable[match(metaDataTable$sampleID, colnames(countTable)),]
  if (!identical(metaDataTable$sampleID, colnames(countTable))){
    stop('The meta-data table is not aligned with column names of the count table!')
  }

  # step1: perform the control dependent and independent DE analysis

  cat("\033[31mStarting control-(in)dependent DE analysis per library ...\033[39m\n")
  ctr_dep_DE_res = list()
  ctr_indep_DE_res = list()
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
    libs_w_ctr = c()
    if (any(dds$RNAi == 'control')){
      libs_w_ctr = c(libs_w_ctr, libs[i])
      invisible(suppressMessages(DE_res_list <- control_dependent_DE(dds)))
      names(DE_res_list) = paste(libs[i], names(DE_res_list),sep = '_')
      ctr_dep_DE_res = c(ctr_dep_DE_res, DE_res_list)
    }

    # perform control-independent DE analysis
    cat("Performing control-independent DE analysis for library", libs[i],"\n")
    cat("Fit main populations ...\n")
    invisible(suppressMessages(zMat <- fit_main_population(dds)))
    cat("Runing DEseq2 again ...\n")
    DE_res_list <- control_independent_DE(dds, zMat)
    names(DE_res_list) = paste(libs[i], names(DE_res_list),sep = '_')
    ctr_indep_DE_res = c(ctr_indep_DE_res, DE_res_list)
  }


  # step2: combine the control dependent and independent DE analysis
  cat("\033[31mCombining the two per library DE analysis...\033[39m\n")

  # find the control outlier genes
  res_list = find_control_outliers(ctr_dep_DE_res, libs_w_ctr, pcutoffs, freqCutoff) # cutoff could be calibrated by the vector like samples
  gene2rmList = res_list[[1]]

  # use control-independent DE analysis for control-outlier genes
  clean_DE_res = combine_DE_result(ctr_dep_DE_res, ctr_indep_DE_res, libs, pcutoffs, gene2rmList)


  # step3: perform empirical null modeling at the entire dataset level
  cat("\033[31mCorrecting test statistics based on emprical null...\033[39m\n")
  for (zz in 1:length(pcutoffs)){
    allgenes = c()
    for (i in 1:length(clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]])){
      allgenes = union(allgenes, rownames(clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[i]]))
    }
    RNAiName = names(clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]])

    # aggregate the raw tables
    pvalueTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(pvalueTbl) = allgenes
    colnames(pvalueTbl) = RNAiName
    pvalueTbl = as.data.frame(pvalueTbl)

    maxExpTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(maxExpTbl) = allgenes
    colnames(maxExpTbl) = RNAiName
    maxExpTbl = as.data.frame(maxExpTbl)

    statTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(statTbl) = allgenes
    colnames(statTbl) = RNAiName
    statTbl = as.data.frame(statTbl)

    for (i in 1:length(RNAiName)){
      tbl =  clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]]
      pvalueTbl[rownames(tbl),RNAiName[i]] = tbl$pvalue
      statTbl[rownames(tbl),RNAiName[i]] = tbl$stat
      maxExpTbl[rownames(tbl),RNAiName[i]] = apply(tbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
    }

    # there are some data cleaning mask (ie plate2 bad genes) in pvalue field. we transfer it to all the tables
    maxExpTbl[is.na(pvalueTbl)] = NA
    statTbl[is.na(pvalueTbl)] = NA
    # we remove genes without DE analysis output
    ind = matrixStats::rowAlls(is.na(pvalueTbl))
    statTbl = statTbl[!ind,]
    pvalueTbl = pvalueTbl[!ind,]
    maxExpTbl = maxExpTbl[!ind,]

    # fit the empirical null
    fit = fit_empirical_null(statTbl,display = F)

    # perform independent filtering and bi-directional p-adjustment
    p_adj_WPS = adjust_pvalues(emp_pmat = fit$p_mat,
                               maxExpTbl = maxExpTbl,
                               lowExpCutoff = independentFilteringCutoff)


    # update the DE tbls
    for (i in 1:length(RNAiName)){
      clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]]$empirical_pvalue =
                                      fit$p_mat[rownames(clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]]),
                                                  RNAiName[i]]
      clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]]$FDR =
                                      p_adj_WPS[rownames(clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]]),
                                                  RNAiName[i]]
      # rename the p-values produced by DEseq2
      colnames(clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]])[5] = 'pvalue_DESeq2'
      # remove the padj column
      clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]] =  clean_DE_res[[paste('cutoff_',as.character(pcutoffs[zz]),sep = '')]][[RNAiName[i]]][,-6]
    }
  }

  # final clean up
  if (length(pcutoffs) == 1){
    clean_DE_res = clean_DE_res[[1]]
  }

  return(clean_DE_res)

}
