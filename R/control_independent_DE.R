#' @title Control-independent DE analysis
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



control_independent_DE <- function(dds_ori, adaZmat, zcutoff = 2.5){

  DE_res_list = list()

  # perform control-independent DE on a per-condition basis
  if (any(levels(dds_ori$RNAi) == 'control')){ # we allow input data even with missing controls
      myseq1 = 2:(length(levels(dds_ori$RNAi)))
    }else{
      myseq1 = 1:(length(levels(dds_ori$RNAi)))
  }
  for (j in myseq1){
    targetGene = levels(dds_ori$RNAi)[j]
    dds = dds_ori

    # modify the treatment coldata labels for control-independent DE
    if ('plate2' %in% names(SummarizedExperiment::colData(dds))){
      # Special treatment for plate2 sample - skipped for non-WPS data
      isPlate2 = unique(dds$plate2[dds$RNAi == targetGene])
      if (length(isPlate2) > 1){
        stop(paste('The plate2 annotation is not consistent within condition', targetGene))
      }
      if (isPlate2){
        dds_tmp = dds
        dds_tmp$RNAi_ori = dds_tmp$RNAi
        dds_tmp$RNAi = factor(as.character(dds_tmp$RNAi), levels = c('others',levels(dds_tmp$RNAi)))
        dds_tmp$RNAi[!(dds_tmp$RNAi %in% c('control',targetGene))] = 'others'
        dds_tmp$RNAi = factor(as.character(dds_tmp$RNAi), levels = c('others', 'control', targetGene))
        # remove extra plate2 in others (to avoid confounding)
        rmInd = dds_tmp$RNAi %in% 'others' & dds$plate2
        dds_tmp = dds_tmp[,!rmInd]
        dds <- dds_tmp
      }else{
        dds_tmp = dds
        dds_tmp$RNAi_ori = dds_tmp$RNAi
        dds_tmp$RNAi = factor(as.character(dds_tmp$RNAi), levels = c('others',levels(dds_tmp$RNAi)))
        dds_tmp$RNAi[!(dds_tmp$RNAi %in% c('control',targetGene))] = 'others'
        # remove all plate2 confounding samples (by default, control samples are plate2)
        rmInd = (dds_tmp$RNAi %in% 'control') | (dds_tmp$RNAi %in% 'others' & dds$plate2)
        dds_tmp = dds_tmp[,!rmInd]
        dds_tmp$RNAi = factor(as.character(dds_tmp$RNAi), levels = c('others',targetGene))
        dds <- dds_tmp
      }
    }else{ # for non-WPS data, we remove controls (as control-independent DE is for control-outlier genes)
      dds_tmp = dds
      dds_tmp$RNAi_ori = dds_tmp$RNAi
      dds_tmp$RNAi = factor(as.character(dds_tmp$RNAi), levels = c('others',levels(dds_tmp$RNAi)))
      dds_tmp$RNAi[!(dds_tmp$RNAi %in% c('control',targetGene))] = 'others'
      # remove control samples
      rmInd = dds_tmp$RNAi %in% 'control'
      dds_tmp$RNAi = factor(as.character(dds_tmp$RNAi), levels = c('others',targetGene))
      dds_tmp = dds_tmp[,!rmInd]
      dds <- dds_tmp
    }

    # clean up the outliers
    dds <- clean_outliers(dds, adaZmat, zcutoff)

    # removing outliers may reveal some actually lowly expressed genes, so refiltering low counts
    keep <- rowSums(counts(dds)>=10) >= 1
    dds <- dds[keep,]

    # run DEseq2
    invisible(suppressMessages(dds <- DESeq2::DESeq(dds,minReplicatesForReplace=7)))
    # filter results
    dds.filt = dds

    # Special treatment for plate2 sample - skipped for non-WPS data
    special_case = F
    if ('plate2' %in% names(SummarizedExperiment::colData(dds))){
      if (isPlate2 & any(dds.filt$RNAi %in% 'control')){# control exist and it is a plate2 confounded RNAi
        cat("\033[31mWPS DATA: detected fully confounded plate2 RNAi, filtering uncertain results to be conservative!\033[39m\n")
        special_case = T
        # define contrast - we look at two types
        mySample1 = paste('RNAi_',targetGene,'_vs_others',sep = '')
        mySample2 = paste('RNAi_control_vs_others',sep = '')
        mySample = mySample1

        res1 <- results(dds.filt,name = stringr::str_replace_all(mySample1,'-','.'),independentFiltering = F)
        colnames(res1)[2] = 'log2FoldChange_raw'

        res2 <- results(dds.filt,name = stringr::str_replace_all(mySample2,'-','.'),independentFiltering = F)
        colnames(res2)[2] = 'log2FoldChange_raw'

        # filter out the case where the FC direction of vector (when significant) is the same as target RNAi (confounded)
        # Definition of significance: we consider the significance as the RNAi of interest is not very different than vector, rather than checking the significance of vector itself
        # we mask the DE only if: (1) the DE is changing in the same direction as vector; (2) vector is significantly changed (p<0.01) and (3) the fold change of the DE is not very different than vector (2-fold in log scale)

        # Note for reproducibility
        # the DE analysis in metabolic WPS dataset uses shrunk FC in the following filter to be more conservative. We decided to skip the use of shrink FC to ensure computational speed for general application.
        # this general package is not designed to numericially reproduce metabolic WPS result!
        res1$pvalue[sign(res1$log2FoldChange_raw) == sign(res2$log2FoldChange_raw ) &
                      res2$pvalue < 0.01 &
                      abs(res1$log2FoldChange_raw) <= abs(res2$log2FoldChange_raw) * 2] = NA
        # after filtering, use res1
        res = res1

        # attach the counts information
        mVals1 = matrixStats::rowMedians(DESeq2::counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
        names(mVals1) = rownames(dds.filt)
        mVals2 = matrixStats::rowMedians(DESeq2::counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'others'])
        names(mVals2) = rownames(dds.filt)
        res = as.data.frame(res)
        res$medianCount_RNAi = mVals1[rownames(res)]
        res$medianCount_ctr = mVals2[rownames(res)]

        }
    }

    # otherwise it is either non-WPS setup or not confounded or no vector - we directly compare with other RNAi
    if (!special_case){
      # define the contrast
      mySample = paste('RNAi_',targetGene,'_vs_others',sep = '')
      res <- results(dds.filt,name = stringr::str_replace_all(mySample,'-','.'),independentFiltering = F)
      colnames(res)[2] = 'log2FoldChange_raw'

      # attach the counts information
      mVals1 = matrixStats::rowMedians(DESeq2::counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
      names(mVals1) = rownames(dds.filt)
      mVals2 = matrixStats::rowMedians(DESeq2::counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'others'])
      names(mVals2) = rownames(dds.filt)
      res = as.data.frame(res)
      res$medianCount_RNAi = mVals1[rownames(res)]
      res$medianCount_ctr = mVals2[rownames(res)]
    }

    DE_res_list[[mySample]] = res

  }
  return(DE_res_list)
}

