#' @title WPS Differential Expression (DE) analysis
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

# all non-base functions to be called by ::
# enter the project folder
# create function files in R
# run devtools::document() to document the new changes
# then run build - check
# when done, run Git - commit - push

WPS_DE <- function(countTable, metaDataTable) {

  # step1: perform the control dependent DE analysis
  for (plateID in c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')){
    # load the target library set
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    batches = levels(batchLabel)
    library(DESeq2)
    library(limma)
    library(edgeR)
    #mergeTbl = data.frame()
    for (i in 1:length(batches)){
      if (!(batches[i] %in% c('met3_lib4','met2_lib5'))){
        # pick the control samples (in the batches belonging to target DE group)
        subsetInd = batchLabel == batches[i]
        # merge selected library with the extra controls
        input_subset = input[,subsetInd] # %>% rownames_to_column('gene') %>%
        # inner_join(input_extr_ctr %>% rownames_to_column('gene'), by = 'gene')
        # rownames(input_subset) = input_subset$gene
        # input_subset = input_subset[,-1]
        batchLabel_subset = colnames(input_subset)
        batchLabel_subset = str_extract(batchLabel_subset,'_rep._met[0-9]+_')
        batchLabel_subset = str_replace(batchLabel_subset,'_met._$','')
        batchLabel_subset = as.factor(batchLabel_subset)
        # batchLabel_subset = as.factor(c(as.character(batchLabel[subsetInd]), batchLabel_extr_ctr))
        RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
        # RNAi_subset[(1+length(RNAi_subset)):(length(batchLabel_extr_ctr)+length(RNAi_subset))] = 'x.vector'

        coldata =  data.frame(batchLabel = batchLabel_subset, RNAi = RNAi_subset)
        rownames(coldata) = colnames(input_subset)
        dds <- DESeqDataSetFromMatrix(countData = input_subset,
                                      colData = coldata,
                                      design= ~ batchLabel + RNAi)
         control_dependent_DE(dds)


}
