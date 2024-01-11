#' @title Combine the control-dependent and -independent DE results
#'
#' @description
#' This function combines the DE results from both control-dependent and -independent DE analysis to resolve false positives related to control outlier genes.
#' This function is for internal-use only.
#' @param ctr_dep_DE_res A list of data frames that are control-dependent DE results for each condition.
#' @param ctr_indep_DE_res A list of data frames that are control-independent DE results for each condition.
#' @param libs A character array of library IDs.
#' @param cutoffs A single value or a numerical array of p-outlier cutoffs
#' @param gene2rmList A list of control outlier genes for each library, identified by \code{find_control_outliers} function.
#' @param FCtype A string defining the column name of logFC in the DE result table. Only for internal use.
#'
#' @return A list of clean DE results with each p-outlier cutoff.
#'
#'
#' @author Xuhang Li

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


