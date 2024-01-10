#' @title bi-directionally adjust the p-values for multiple testing
#'
#' @description
#' A wrapper function using \code{p.adjust} to adjust for multiple testing in occuring for both many genes and many testing conditions in large-scale experiments. This function is for internal use and must work with a maxExpTbl and lowExpCutoff for indepdendent filtering prior to adjustments.
#' @param emp_pmat A gene-by-condition matrix of p-values to be adjusted
#' @param maxExpTbl A gene-by-condition matrix of maximal median expression level between the control and case being contrasted in the statistical test. Used for independent filtering.
#' @param lowExpCutoff Normalized read count cutoff for independent filtering. Default is 30.
#'
#' @return A list of the fitting results:
#' \describe{
#'   \item{\code{dual_fdr_mat}}{A gene-by-condition matrix of adjusted p-values (empirical FDR). }
#' }
#'
#'
#' @export adjust_pvalues
#'
#' @author Xuhang Li


# all non-base functions to be called by ::
# enter the project folder
# create function files in R
# run devtools::document() to document the new changes
# then run build - check
# when done, run Git - commit - push


adjust_pvalues <- function(emp_pmat, maxExpTbl, lowExpCutoff = 30){ # performs independent filtering and two-dimensional p-adjustment
  row_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
  maxExpTbl_pass = maxExpTbl > lowExpCutoff
  maxExpTbl_pass[is.na(maxExpTbl_pass)] = F
  for (j in 1:nrow(row_fdr_mat)){
    row_fdr_mat[j,maxExpTbl_pass[j,]] = stats::p.adjust(emp_pmat[j, maxExpTbl_pass[j,]],method = 'BH')
  }
  col_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
  for (j in 1:ncol(col_fdr_mat)){
    col_fdr_mat[maxExpTbl_pass[,j],j] = stats::p.adjust(emp_pmat[maxExpTbl_pass[,j], j],method = 'BH')
  }
  dual_fdr_mat = pmax(col_fdr_mat, row_fdr_mat)
  rownames(dual_fdr_mat) = rownames(emp_pmat)
  colnames(dual_fdr_mat) = colnames(emp_pmat)
  return(dual_fdr_mat)
}
