#' @title Perform the robust fitting of the main population of each gene
#'
#' @description
#' This is a wrapper function that uses \code{AdaTiSS} algorithm (Wang et al. 2021) to fit the main population of the expression levels of a gene.
#' The mean and variance of this main population are used to calculate z-scores for each gene in each condition, which are further used in the control-independent
#' DE analysis to remove outliers of the main population.
#' @param dds A DESeqDataSet objective that contains the input data. When this function is used alone, the \code{dds} object should be constructed following the metadata format of WPS dataset (see examples).
#'
#' @return A matrix of fitted z-scores for each gene in each sample.
#'
#'
#' @export fit_main_population
#'
#' @author Xuhang Li
#' @references
#' This function uses AdaTiSS algorithm published by Wang et al.
#'
#' \cite{Wang, Meng, Lihua Jiang, and Michael P. Snyder. "AdaTiSS: a novel data-Ada ptive robust method for identifying Ti ssue S pecificity S cores." Bioinformatics 37.23 (2021): 4469-4476.}
#' @examples
#' data(example_dds)
#' adaZmat <- fit_main_population(example_dds)

fit_main_population <- function(dds){

  # obtain the batch-free counts
  tmp = DESeq2::counts(dds, normalized=TRUE)
  normCounts = log2(tmp + 1)
  # remove batch effects in replicates
  normCounts_clean <- limma::removeBatchEffect(normCounts, batch = dds$batchLabel,design= stats::model.matrix(~ dds$RNAi))

  # fit AdaTiss for every gene
  zMat = matrix(NA,nrow = nrow(normCounts_clean),ncol = ncol(normCounts_clean))
  rownames(zMat) = rownames(normCounts_clean)
  colnames(zMat) = colnames(normCounts_clean)
  case1 = c()
  case2 = c()
  for (i in 1:nrow(normCounts_clean)){
    y = normCounts_clean[i,]
    out = AdaReg(stats::model.matrix(~1,data = as.data.frame(y)), y)
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
      if (stats::median(y) <= log2(10+1)){
        meanIn = mean(y)
        sdIn = stats::sd(y)
        zr1 = (y-meanIn)/sdIn
        zMat[i,] = zr1
        case2 = c(case2,i)
      }else{
        meanIn = stats::median(y)
        sdIn = stats::mad(y)
        zr1 = (y-meanIn)/sdIn
        zMat[i,] = zr1
        case2 = c(case2,i)
      }

    }

    if (i %% 1000 == 0){
      print(paste('AdaTiss fitting', ' ... ', 100*i/nrow(zMat),'% (this may take a long time)',sep = ''))
      print(paste('normal fit: ',length(case1)/i,sep = ''))
      print(paste('poor fit: ',length(case2)/i,sep = ''))
    }
  }

  return(zMat)
}


