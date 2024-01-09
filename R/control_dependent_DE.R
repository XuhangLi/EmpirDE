#' @title Perform control dependent Differentially Expression (DE) analysis
#'
#' @description
#' This is a wrapper function to perform DE analysis in the cannonical way (control-dependent) using DEseq2.
#' This function is designed for internal use.
#' @param dds A dds object from DEseq2
#'
#' @return A list of DE results (each element corresponds to the DE result of one condition, indexed by the condID):
#' \describe{
#'   \item{\code{DE_table}}{The data frame containing DE results from DEseq2 for one condition.}
#' }
#'
#'
#' @export control_dependent_DE
#'
#' @author Xuhang Li



control_dependent_DE <- function(dds) {


  library(dplyr)
  library(tibble)
  library(stringr)



  # filtering
  keep <- rowSums(counts(dds)>=10) >= 1
  dds <- dds[keep,]

  # library("BiocParallel")
  # register(MulticoreParam(2))
  dds <- DESeq(dds,minReplicatesForReplace=Inf)
  # resultsNames(dds) # lists the coefficients
  # plotDispEsts(dds)
  # assays(dds)[["cooks"]] #cooksCutoff in result function
  # par(mar=c(8,5,2,2))
  # boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

  #res <- results(dds, name = 'RNAi_nhr_20_vs_vector',alpha = 0.05)
  # or to shrink log fold changes association with condition:
  # res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
  # Note: We have sped up the apeglm method so it takes roughly about the same amount of time as normal, e.g. ~5 seconds for the pasilla dataset of ~10,000 genes and 7 samples. If fast shrinkage estimation of LFC is needed, but the posterior standard deviation is not needed, setting apeMethod="nbinomC" will produce a ~10x speedup, but the lfcSE column will be returned with NA. A variant of this fast method, apeMethod="nbinomC*" includes random starts.

  # write result
  #list the DE for each individual RNAi
  for (j in 2:(length(levels(RNAi_subset)))){
    targetGene = levels(RNAi_subset)[j]
    mySample = paste('RNAi_',targetGene,'_vs_x.vector',sep = '')
    # supply cooks filter metric (maximum cook for all samples related to tested covariate)
    dds.filt = dds

    # update base mean to the basemean between RNA and control
    # we use DEseq2's default setting to optimize power, but we define the base mean by only vectors plus the RNAi in query
    libWiseBaseMean = mcols(dds.filt)$baseMean
    mycounts = counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)]
    myconditions = dds.filt$RNAi[dds.filt$RNAi %in% c('x.vector',targetGene)]
    myWeights = rep(0,length(myconditions))
    myWeights[myconditions == 'x.vector'] = 0.5 / sum(myconditions == 'x.vector')
    myWeights[myconditions == targetGene] = 0.5 / sum(myconditions == targetGene)

    baseMean = apply(mycounts, 1, result <- function(x){result = weighted.mean(x,myWeights)})
    # baseVar = rowVars(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)])
    mcols(dds.filt)$baseMean = baseMean
    mcols(dds.filt)$baseVar = NA

    res <- results(dds.filt,name = str_replace_all(mySample,'-','.'),independentFiltering = F) # we may use "replace" instead of cook filter for met library
    stat = res$stat
    logFC_raw = res$log2FoldChange
    res <- lfcShrink(dds.filt, coef = str_replace_all(mySample,'-','.'), type="apeglm",apeMethod="nbinomC", res=res)
    res$stat = stat
    res$log2FoldChange_raw = logFC_raw

    # attach the counts information
    # filter by mean read counts of each component of contrast(will do in cutoff step)
    # batch effect is causing some trouble here, so we use median instead of mean; for met lib, we want to use mean
    # ==> change to median alway to ensure robustness
    mVals1 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
    names(mVals1) = rownames(dds.filt)
    mVals2 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'x.vector'])
    names(mVals2) = rownames(dds.filt)
    res = as.data.frame(res)
    res$medianCount_RNAi = mVals1[rownames(res)]
    res$medianCount_ctr = mVals2[rownames(res)]
    res$libWiseBaseMean = libWiseBaseMean # this will be used in the MERGE modeling! (otherwise the pseudocounts will not be comparable across conditions)
    # save
    write.csv(res,paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'.csv',sep = ''))
    # resOrdered <- res[order(res$padj),]
    # outputTbl = as.data.frame(subset(res, padj < 0.05))
    # if (nrow(outputTbl) > 0){
    #   outputTbl = cbind(data.frame(WBID = rownames(outputTbl)),outputTbl)
    #   rownames(outputTbl) = 1:nrow(outputTbl)
    #   outputTbl$RNAi = rep(targetGene,nrow(outputTbl))
    # }else{
    #   outputTbl = data.frame(WBID = 'NoHit')
    #   outputTbl$baseMean = NA
    #   outputTbl$log2FoldChange = NA
    #   outputTbl$lfcSE = NA
    #   outputTbl$pvalue = NA
    #   outputTbl$padj = NA
    #   outputTbl$medianCount_RNAi = NA
    #   outputTbl$medianCount_ctr = NA
    #   outputTbl$RNAi = targetGene
    #   outputTbl$libWiseBaseMean = NA
    #   rownames(outputTbl) = 1
    # }
    # outputTbl$batchID = batches[i]
    # mergeTbl = rbind(mergeTbl, outputTbl)
    print(mySample)
  }
  return()
}
