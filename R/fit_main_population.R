#' @title Perform the robust fitting of the main population of each gene
#'
#' @description
#' To be updated!
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
#' @references This algorithm uses AdaTiSS published in Wang, Meng, Lihua Jiang, and Michael P. Snyder. "AdaTiSS: a novel data-Ada ptive robust method for identifying Ti ssue S pecificity S cores." Bioinformatics 37.23 (2021): 4469-4476.
#' @examples
#' data(WPS_example_data)
#' result <- WPS_DE(countTable, metaDataTable)

fit_main_population <- function(dds){

  # obtain the batch-free counts
  tmp = DESeq2::counts(dds, normalized=TRUE)
  normCounts = log2(tmp + 1)
  # remove batch effects in replicates
  normCounts_clean <- limma::removeBatchEffect(normCounts, batch = dds$batchLabel,design= model.matrix(~ dds$RNAi))

  # fit AdaTiss for every gene
  zMat = matrix(NA,nrow = nrow(normCounts_clean),ncol = ncol(normCounts_clean))
  rownames(zMat) = rownames(normCounts_clean)
  colnames(zMat) = colnames(normCounts_clean)
  case1 = c()
  case2 = c()
  for (i in 1:nrow(normCounts_clean)){
    y = normCounts_clean[i,]
    out = AdaReg(model.matrix(~1,data = as.data.frame(y)), y)
    pi0 = as.numeric(out$res.info['pi0.hat'])
    if (pi0 >= 0.7){# data fitted well on the main population (>70% inlier)
      zr1 = (y-out$beta.rob.fit)/sqrt(out$var.sig.gp.fit)
      zMat[i,] = zr1
      case1 = c(case1,i)
    }else{# poor fitting or too many zero
      # there may not be a single main population
      # we consider all samples together to be conservative

      # we found median of 10 is sufficient to filter out those lowly expressed ones, so
      # when expression is low (median <= 10), we consider the poor fitting as result of high variation, so we use mean for conservative
      # when expression is high (median > 10), we consider the poor fitting as result of high responses, so we use median for best power
      if (median(y) <= log2(10+1)){
        meanIn = mean(y)
        sdIn = sd(y)
        zr1 = (y-meanIn)/sdIn
        zMat[i,] = zr1
        case2 = c(case2,i)
      }else{
        meanIn = median(y)
        sdIn = mad(y)
        zr1 = (y-meanIn)/sdIn
        zMat[i,] = zr1
        case2 = c(case2,i)
      }

    }

    if (i %% 1000 == 0){
      print(paste('AdaTiss fitting for ',lib, ' ... ', 100*i/nrow(zMat),'%',sep = ''))
      print(paste('normal fit: ',length(case1)/i,sep = ''))
      print(paste('poor fit: ',length(case2)/i,sep = ''))
    }
  }

  return(zMat)
}


