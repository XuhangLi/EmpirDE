#' @title Perform control dependent Differentially Expression (DE) analysis
#'
#' @description
#' This is a wrapper function to perform DE analysis in the cannonical way (control-dependent) using DEseq2.
#' This function is designed for internal use.
#' @param dds A DESeqDataSet object defined in \code{DESeq2} package
#'
#' @return A list of DE results (each element corresponds to the DE result of one condition, indexed by the condID).
#'
#'
#'
#' @author Xuhang Li


control_dependent_DE <- function(dds) {

  DE_res_list = list()
  dds <- DESeq2::DESeq(dds,minReplicatesForReplace=Inf)

  # reformat the results for custom filtering
  #list the DE for each individual RNAi
  for (j in 2:(length(levels(dds$RNAi)))){
    targetGene = levels(dds$RNAi)[j]
    mySample = paste('RNAi_',targetGene,'_vs_control',sep = '')

    # custom filtering setup
    dds.filt = dds
    res <- DESeq2::results(dds.filt,name = stringr::str_replace_all(mySample,'-','.'),independentFiltering = F) # turning off independent filtering and will impose it later on manually
    # in WPS DE, we used the simple FC estimates from GLM instead of shrunk FC that is often used in standard DE analysis
    # this is because WPS DE controls FDR through empirical null modeling, leaving the FC shrinkage unnecessary.
    # stat = res$stat
    # logFC_raw = res$log2FoldChange
    # res <- lfcShrink(dds.filt, coef = str_replace_all(mySample,'-','.'), type="apeglm",apeMethod="nbinomC", res=res)
    # res$stat = stat

    # to emphasize it is the unshrunk fc, we rename the column name
    colnames(res)[2] = 'log2FoldChange_raw'

    # calculate the median normalized counts for treatments and controls
    mVals1 = matrixStats::rowMedians(DESeq2::counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
    names(mVals1) = rownames(dds.filt)
    mVals2 = matrixStats::rowMedians(DESeq2::counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'control'])
    names(mVals2) = rownames(dds.filt)
    res = as.data.frame(res)
    res$medianCount_RNAi = mVals1[rownames(res)]
    res$medianCount_ctr = mVals2[rownames(res)]
    # save
    DE_res_list[[mySample]] = res
  }
  return(DE_res_list)
}
