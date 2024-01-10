#' @title Clean up outlier genes in null population
#'
#' @description
#' To be completed
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
