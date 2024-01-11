#' @title Clean up outlier genes in null population
#'
#' @description
#' This function uses the z-scores produced by robust main population fitting to identify the outlier conditions of each gene, and replace them with imputed inlier values.
#' This is a function for internal use only.
#' @param dds a DESeqDataSet object containing the input data
#' @param adaZmat z-score matrix corresponding to the data in \code{dds}, produced by \code{\link{fit_main_population}} function.
#' @param zcutoff z-score cutoff to define the outliers as compared with the inlier population (main/null population).
#' @param rand_seed random seed
#' @return an updated DESeqDataSet objective with all outliers in the null population replaced
#' \describe{
#'   \item{\code{dds}}{updated \code{dds}}
#' }
#'
#'
#' @author Xuhang Li




clean_outliers <- function(dds, adaZmat, zcutoff, rand_seed = 19951126){

  # first, find out the genes and samples to be replaced
  RNAi_others = unique(as.character(dds$RNAi_ori[dds$RNAi == 'others']))
  aveZmat = data.frame(row.names = rownames(dds))
  for (RNAiName in RNAi_others){
    aveZmat[RNAiName] = matrixStats::rowMedians(as.matrix(adaZmat[,dds$RNAi_ori %in% RNAiName]))
  }
  isOutlier = abs(aveZmat) > zcutoff

  # we skip the genes whose outlier content is more than 50%
  if (sum(rowSums(isOutlier) > 0.5 * ncol(isOutlier)) > 0.1 * nrow(isOutlier)){
    cat('Caution: outlier replacement was skipped for',sum(rowSums(isOutlier) > 0.5 * ncol(isOutlier)),'genes (>10%)\n')
  }
  isOutlier[rowSums(isOutlier) > 0.5 * ncol(isOutlier),] = FALSE
  if (sum(rowSums(isOutlier) <= 0.5 * ncol(isOutlier) & rowSums(isOutlier) > 0) > 0.9 * nrow(isOutlier)){
    cat('Caution: outlier replacement was performed for',sum(rowSums(isOutlier) <= 0.5 * ncol(isOutlier) & rowSums(isOutlier) > 0),'genes (>90%)\n')
  }
  # next, replace the outlier counts with imputed ones
  normCounts_imp = DESeq2::counts(dds,normalized=TRUE)
  normCounts_imp = normCounts_imp[,dds$RNAi == 'others']
  RNAiLabels = dds$RNAi_ori[dds$RNAi == 'others']
  newCounts <- DESeq2::counts(dds)
  # set outliers to NA
  for (i in 1:ncol(isOutlier)){
    normCounts_imp[isOutlier[,i],RNAiLabels %in% colnames(isOutlier)[i]] = NA
    newCounts[isOutlier[,i],dds$RNAi_ori %in% colnames(isOutlier)[i]] = NA
  }
  # impute the NAs
  set.seed(rand_seed)
  newCounts_tmp = newCounts[,colnames(normCounts_imp)]
  NAidx = which(is.na(newCounts_tmp))
  # bootstrap impute (sampling within batch)
  for (i in 1:nrow(normCounts_imp)){
    if (any(is.na(normCounts_imp[i,]))){
      tmp = normCounts_imp[i,]
      naInds = which(as.logical(is.na(tmp)))
      reps = as.character(SummarizedExperiment::colData(dds)[names(tmp)[naInds],'batchLabel'])
      tmp = tmp[!is.na(tmp)]
      reps2 = as.character(SummarizedExperiment::colData(dds)[names(tmp),'batchLabel'])
      for (rep in unique(reps)){
        normCounts_imp[i,naInds[reps == rep]] = as.numeric(sample(tmp[reps2==rep],sum(reps == rep),replace = T))
      }
    }
  }
  sizeF = DESeq2::sizeFactors(dds)[colnames(normCounts_imp)]
  sizeFmat = matrix(rep(sizeF,each=nrow(normCounts_imp)),nrow=nrow(normCounts_imp))
  replacementCounts = as.integer(normCounts_imp * sizeFmat)
  newCounts_tmp[NAidx] <- replacementCounts[NAidx]
  newCounts[, colnames(newCounts_tmp)] = newCounts_tmp

  # update dds
  DESeq2::counts(dds) <- newCounts
  return(dds)
}
