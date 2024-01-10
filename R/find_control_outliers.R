#' @title Find control-outlier genes for each library
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

find_control_outliers <- function(ctr_dep_DE_res, libs, pcutoffs, freqCutoff, FCtype = 'log2FoldChange_raw'){
  genelist_all_cutoffs = list()
  for (i in 1:length(pcutoffs)){
    genelist_all_cutoffs[[paste('cutoff_',as.character(pcutoffs[i]),sep = '')]] = list()
  }

  for (lib in libs){
    # grab all the DE res in this lib
    allres = names(ctr_dep_DE_res)[stringr::str_detect(names(ctr_dep_DE_res), paste('^',lib,'_',sep = ''))]
    allgenes = c()
    for (i in 1:length(allres)){
      allgenes = union(allgenes, rownames(ctr_dep_DE_res[[allres[i]]]))
    }
    tbl = data.frame(row.names = rownames(allgenes))
    for (thisRNAi in allres){
      tmp = ctr_dep_DE_res[[thisRNAi]]
      tbl[rownames(tmp),thisRNAi] = -log10(tmp$pvalue) * sign(tmp[,FCtype])
    }
    Qs = rowQuantiles(x = as.matrix(tbl),probs = c(0.75,0.5,0.25), na.rm = T)
    upperQ = Qs[,1]
    midQ = Qs[,2]
    lowerQ = Qs[,3]
    for (cutoff in pcutoffs){
      cutoff_lowerQ = cutoff * 10
      cutoff_midQ = cutoff
      badInd = which(
        (midQ > -log10(cutoff_midQ) &
           lowerQ > -log10(cutoff_lowerQ)) |
          (upperQ < log10(cutoff_lowerQ) &
             midQ < log10(cutoff_midQ)))
      badGenes = names(midQ)[badInd]

      badGenes = setdiff(badGenes,NA)
      genelist_all_cutoffs[[paste('cutoff_',as.character(cutoff),sep = '')]][[as.character(lib)]] = badGenes
    }
  }

  # next, define the core bad gene set
  for (z in 1:length(pcutoffs)){
    genelist2 = table(unlist(genelist_all_cutoffs[[z]]))
    coreset = names(genelist2)[genelist2 > freqCutoff]
    cat('At cutoff ', pcutoffs[z],', ',length(coreset),' core control-outlier genes identified',sep = '')
    coresetTbl = data.frame(gene = coreset)
    coresetTbl$freq = genelist2[coreset]

    for (i in 1:length(genelist_all_cutoffs[[z]])){
      genelist_all_cutoffs[[z]][[i]] = union(genelist_all_cutoffs[[z]][[i]], coreset)
    }
  }

  return(list(genelist_all_cutoffs,coresetTbl))
}
