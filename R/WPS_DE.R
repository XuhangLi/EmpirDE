#' @title WPS Differential Expression (DE) analysis
#'
#' @description
#' Perform the two-pronged WPS DE analysis with an input dataset consistent with WPS data structure.
#'
#' To run WPS DE analysis on custom datasets, the data should be collected over a few multiplexing libraries (or equivalent settings).
#' WPS DE first perform both control-dependent and -independent DE analysis on each individual library, followed by identification of
#' control-outlier genes in each library and fix such technical problem by combining with control-independent DE results. The second
#' stage of WPS DE combines conditions in all libraries and models the test statistic distribution based empirical null. This process
#' enventually produces the statistically rigorous FDR of each DE test.
#'
#' @param countTable A gene-by-sample matrix of read counts as the input dataset for DE analysis
#' @param metaDataTable A metadata table for the input dataset with the following columns:
#' \describe{
#'   \item{\code{sampleID}}{unique IDs for each sample that correspond to the column names in \code{countTable}.}
#'   \item{\code{covTreatment}}{covariate indicating experimental treatments to be tested for (e.g., RNAi conditions). Must have a 'control' condition if control-dependent DE analysis is deisred.}
#'   \item{\code{covBatch}}{covariate indicating experimental batches, which by default is the replicate batch. Other reasonable batch labels within each library may also be used.}
#'   \item{\code{libID}}{unique IDs for identifying the sequencing library of each sample. WPS DE analysis is first conducted at individual sequencing library level thus this ID is used to match samples pooled in the same library.}
#'   \item{\code{plate2}}{(optional) Logical values indicating if the sample is cultured in the second 96-well plate. This column is only needed for WPS experiments that involves plate 2 confounding effects.}
#' }
#' @param params Custom WPS DE paramerers (default is NULL).
#' \describe{
#'   \item{\code{pcutoffs}}{The p-outlier cutoff for selecting control-outlier genes. By default the value is 0.005 based on our benchmarking study. A single value or a numerical array can be supplied to titrate this parameter.}
#'   \item{\code{freqCutoff}}{The frequency cutoff for defining the core control-outlier genes. Default is 25 percent of the number of libraries.}
#'   \item{\code{independentFilteringCutoff}}{A numerical read count cutoff for conducting independent filtering prior to multiple testing adjustment. This cutoff is applied to the median normalized count of the control and treatment samples in comparison. The greater median should be higher than this cutoff to be included in the final result (otherwise masked to \code{NA}). Default is set to 30 based on our benchmarking study.}
#' }
#'
#' @return A list of the DE result tables for each condition. DE result table contains following columns:
#' \describe{
#'   \item{\code{baseMean}}{baseMean from DESeq2.}
#'   \item{\code{log2FoldChange_raw}}{log2FoldChange from DESeq2. Adding the suffix '_raw' is to emphasize this metric is directly from DESeq2 model without applying any shrinkage algorithm.}
#'   \item{\code{lfcSE}}{lfcSE from DESeq2.}
#'   \item{\code{stat}}{stat from DESeq2 (Wald testing statistic)}
#'   \item{\code{pvalue_DESeq2}}{original p-value produced by DESeq2}
#'   \item{\code{medianCount_RNAi}}{median normalized count of samples in the treatment condition. Used in the independent filtering.}
#'   \item{\code{medianCount_ctr}}{median normalized count of samples in the control condition. Used in the independent filtering.}
#'   \item{\code{DE_source}}{The DE result for this gene is based on which type of DE analysis, either control-dependent (vs. control) or independent (vs. control-independent null).}
#'   \item{\code{empirical_pvalue}}{Empirical p-values based on corrected test statistic. This is the final DE testing p-value of WPS DE framework.}
#'   \item{\code{FDR}}{False Discovery Rate (FDR) of the DE test. This is the final DE testing FDR of WPS DE framework.}
#' }
#'
#'
#' @export WPS_DE
#'
#' @author Xuhang Li
#' @examples
#' data("countTable")
#' data("metaDataTable")
#' result <- WPS_DE(countTable, metaDataTable)


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
    keep <- rowSums(DESeq2::counts(dds)>=10) >= 1
    dds <- dds[keep,]
    # pre-calculating size factors
    dds = DESeq2::estimateSizeFactors(dds)

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
