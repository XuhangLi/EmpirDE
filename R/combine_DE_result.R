#' @title Combine the control-dependent and -independent DE results to resolve control outlier genes
#'
#' @description
#' TBD
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

#lowExpCutoff = 30
combine_DE_result <- function(ctr_dep_DE_res, ctr_indep_DE_res, libs, cutoffs,gene2rmList, FCtype = 'log2FoldChange_raw'){

  clean_DE_res = list()
  for (zz in 1:length(libs)){
    lib = libs[zz]
    allres = names(ctr_indep_DE_res)[stringr::str_detect(names(ctr_indep_DE_res), paste('^',lib,'_',sep = ''))]

    cat("Processing ", lib,"\n")

    ctr_indep_rate = list()
    for (i in 1:length(cutoffs)){
      ctr_indep_rate[[paste('cutoff_',as.character(cutoffs[i]),sep = '')]] = NA
    }

    for (j in 1:length(allres)){
      # grab the two DE res for the current lib
      mySample2 = allres[j]
      mySample1 = stringr::str_replace(mySample2,'_vs_others$','_vs_control')
      mySample_clean = stringr::str_remove(mySample2,'_vs_others$')

      if(mySample1 %in% names(ctr_dep_DE_res)){
        outputTbl1 = ctr_dep_DE_res[[mySample1]]
        outputTbl1$DE_source = 'vs. control'
      }else{
        outputTbl1 = ctr_indep_DE_res[[mySample2]]
        outputTbl1$DE_source = 'vs. control-independent null'
      }
      outputTbl2 = ctr_indep_DE_res[[mySample2]]
      outputTbl2$DE_source = 'vs. control-independent null'

      # we merge the two DE result and readjust the FDR
      # first, overide with unique DE identified by ctr-indep res only (lower FDR and higher FC)
      # add unique rows in ctr-indep
      outputTbl1 = rbind(outputTbl1,outputTbl2[setdiff(rownames(outputTbl2),rownames(outputTbl1)),])
      commomG = intersect(rownames(outputTbl2),rownames(outputTbl1))
      # replace when p value is much smaller and fc is larger (more rep adds power and correct for some problematic control cases)
      repInd = outputTbl2[commomG,'pvalue'] < 0.01 * outputTbl1[commomG,'pvalue'] &
        abs(outputTbl2[commomG,FCtype]) > abs(outputTbl1[commomG,FCtype])
      repG = setdiff(commomG[repInd],NA)
      outputTbl1[repG,] = outputTbl2[repG,colnames(outputTbl1)]

      # replace the bad genes with ctr-indep result
      for (cutoff in cutoffs){
        gene2rm = gene2rmList[[paste('cutoff',cutoff,sep = '_')]][[lib]]
        ind1 = !(rownames(outputTbl1) %in% gene2rm)
        ind2 = rownames(outputTbl2) %in% gene2rm
        outputTbl_merge = rbind(outputTbl1[ind1, ],
                                outputTbl2[ind2, colnames(outputTbl1)])

        ctr_indep_rate[[paste('cutoff',cutoff,sep = '_')]] = c(ctr_indep_rate[[paste('cutoff',cutoff,sep = '_')]],
                                                               sum(outputTbl_merge$DE_source == 'vs. control-independent null')/nrow(outputTbl_merge))

        # mask FDR and wait for next step
        outputTbl_merge$padj = NA #p.adjust(outputTbl_merge$pvalue,method = 'BH')

        clean_DE_res[[paste('cutoff',cutoff,sep = '_')]][[mySample_clean]] = outputTbl_merge
      }
    }

    # print the rate of replacement
    for (i in 1:length(cutoffs)){
      ave = mean(ctr_indep_rate[[paste('cutoff_',as.character(cutoffs[i]),sep = '')]], na.rm = T) * 100
      cat(round(ave,1),'% genes used control-independent DE results\n',sep = '')
    }
  }
  return(clean_DE_res)
}


