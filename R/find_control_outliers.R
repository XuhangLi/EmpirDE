#' @title Find control-outlier genes for each library
#'
#' @description
#' This function identifies the control-outlier genes in each WPS library based on the results of control-dependent DE analysis.
#' This function is only designed for internal use.
#' @param ctr_dep_DE_res A list of control-dependent DE results
#' @param libs A character array of library IDs in the dataset
#' @param pcutoffs A single value or a numerical array of p-outlier cutoffs
#' @param freqCutoff A frequency cutoff to define the core set of control outlier gene. Default is 25% of total number of libraries.
#' @param FCtype A string defining the column name of logFC in the DE result table. Only for internal use.
#'
#'
#' @return A list of the control-outliers:
#' \describe{
#'   \item{\code{genelist_all_cutoffs}}{A list of control outlier genes in each library. When \code{pcutoffs} has more than one values, this list is first grouped by different cutoffs and then by each library.}
#'   \item{\code{coresetTbl}}{A data frame showing the name and frequency (i.e., appearing in # of libraries) of core control-outlier genes.}
#' }
#'
#'
#' @export
#'
#' @author Xuhang Li


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
