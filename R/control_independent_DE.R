library(dplyr)
library(tibble)
library(stringr)

cleanOutliers <- function(dds, adaZmat, zcutoff, imputeMethod){
  # some buggy places
  colnames(adaZmat) = str_replace_all(colnames(adaZmat),'_lin\\.7_','_lin-7_')
  colnames(adaZmat) = str_replace_all(colnames(adaZmat),'SPECIAL_acly\\.1_acly\\.2_','SPECIAL_acly-1_acly-2_')
  colnames(adaZmat) = str_replace_all(colnames(adaZmat),'SPECIAL_hxk\\.1_hxk\\.2_hxk\\.3_','SPECIAL_hxk-1_hxk-2_hxk-3_')
  # submatrix
  adaZmat =  adaZmat[rownames(dds), colnames(dds)]

  # first, find out the genes and samples to be replaced
  RNAi_others = unique(as.character(dds$RNAi_ori[dds$RNAi == 'others']))
  aveZmat = data.frame(row.names = rownames(dds))
  for (RNAiName in RNAi_others){
    aveZmat[RNAiName] = rowMedians(as.matrix(adaZmat[,dds$RNAi_ori %in% RNAiName]))
  }
  isOutlier = abs(aveZmat) > zcutoff

  # we skip the genes whose outlier content is more than 50%
  print(paste(str_extract(colnames(dds)[1],'met._lib.'),': outlier replacement was skipped for ',sum(rowSums(isOutlier) > 0.5 * ncol(isOutlier)),' genes.',sep = ''))
  isOutlier[rowSums(isOutlier) > 0.5 * ncol(isOutlier),] = FALSE
  print(paste(str_extract(colnames(dds)[1],'met._lib.'),': outlier replacement was performed for ',sum(rowSums(isOutlier) <= 0.5 * ncol(isOutlier) & rowSums(isOutlier) > 0),' genes.',sep = ''))
  # next, replace the outlier counts with imputed ones
  normCounts_imp = counts(dds,normalized=TRUE)
  normCounts_imp = normCounts_imp[,dds$RNAi == 'others']
  RNAiLabels = dds$RNAi_ori[dds$RNAi == 'others']
  newCounts <- counts(dds)
  # set outliers to NA
  for (i in 1:ncol(isOutlier)){
    normCounts_imp[isOutlier[,i],RNAiLabels %in% colnames(isOutlier)[i]] = NA
    newCounts[isOutlier[,i],dds$RNAi_ori %in% colnames(isOutlier)[i]] = NA
  }
  # impute the NAs
  if (imputeMethod == 'mean'){
    NAidx = which(is.na(newCounts))
    nullMean <- rowMeans(normCounts_imp,na.rm = T)
    replacementCounts = as.integer(outer(nullMean, sizeFactors(dds), "*"))
    newCounts[NAidx] <- replacementCounts[NAidx]
  }else if (imputeMethod == 'random'){
    set.seed(19951126)
    newCounts_tmp = newCounts[,colnames(normCounts_imp)]
    NAidx = which(is.na(newCounts_tmp))
    # bootstrap impute
    for (i in 1:nrow(normCounts_imp)){
      if (any(is.na(normCounts_imp[i,]))){
        tmp = normCounts_imp[i,]
        naInds = is.na(tmp)
        tmp = tmp[!is.na(tmp)]
        normCounts_imp[i,naInds] = as.numeric(sample(tmp,sum(is.na(normCounts_imp[i,])),replace = T))
      }
    }
    sizeF = sizeFactors(dds)[colnames(normCounts_imp)]
    sizeFmat = matrix(rep(sizeF,each=nrow(normCounts_imp)),nrow=nrow(normCounts_imp))
    replacementCounts = as.integer(normCounts_imp * sizeFmat)
    newCounts_tmp[NAidx] <- replacementCounts[NAidx]
    newCounts[, colnames(newCounts_tmp)] = newCounts_tmp
  }else if (imputeMethod == 'knn'){
    # knn may be intrisically not suitable as the outlier may coexpresses with outliers
    newCounts_tmp = newCounts[,colnames(normCounts_imp)]
    NAidx = which(is.na(newCounts_tmp))
    tmp = impute::impute.knn(normCounts_imp,rowmax = 0.8, rng.seed = 19951126, maxp = Inf)
    normCounts_imp = tmp$data
    sizeF = sizeFactors(dds)[colnames(normCounts_imp)]
    sizeFmat = matrix(rep(sizeF,each=nrow(normCounts_imp)),nrow=nrow(normCounts_imp))
    replacementCounts = as.integer(normCounts_imp * sizeFmat)
    newCounts_tmp[NAidx] <- replacementCounts[NAidx]
    newCounts[, colnames(newCounts_tmp)] = newCounts_tmp
  }else if (imputeMethod == 'random_strat'){
    set.seed(19951126)
    newCounts_tmp = newCounts[,colnames(normCounts_imp)]
    NAidx = which(is.na(newCounts_tmp))
    # bootstrap impute
    for (i in 1:nrow(normCounts_imp)){
      if (any(is.na(normCounts_imp[i,]))){
        tmp = normCounts_imp[i,]
        naInds = which(as.logical(is.na(tmp)))
        reps = str_extract(names(tmp)[naInds],'_rep._')
        tmp = tmp[!is.na(tmp)]
        reps2 = str_extract(names(tmp),'_rep._')
        for (rep in unique(reps)){
          normCounts_imp[i,naInds[reps == rep]] = as.numeric(sample(tmp[reps2==rep],sum(reps == rep),replace = T))
        }
      }
    }
    sizeF = sizeFactors(dds)[colnames(normCounts_imp)]
    sizeFmat = matrix(rep(sizeF,each=nrow(normCounts_imp)),nrow=nrow(normCounts_imp))
    replacementCounts = as.integer(normCounts_imp * sizeFmat)
    newCounts_tmp[NAidx] <- replacementCounts[NAidx]
    newCounts[, colnames(newCounts_tmp)] = newCounts_tmp
  }else if (imputeMethod == 'scImpute'){
    stop('unsuported method!')
    # not working! estimate too small counts; not suitable
  }else{
    stop('unsuported method!')
  }


  # update dds
  counts(dds) <- newCounts
  dds
}

# global para
zcutoff = 2.5
imputeMethods = c('random_strat')


# load the plate2 samples
plate2RNAi = read.csv('./../input_data/metaData/plate2_samples.csv')
plate2RNAi$RNAi_name_final = str_replace(plate2RNAi$RNAi_name_final,'-','_')
plate2RNAi$RNAi_name_final = paste('x.',str_replace(plate2RNAi$RNAi_name_final,' +','_'),sep = '')
plate2RNAi$library = str_replace(plate2RNAi$library,'-','_')

# plateID = 'met1'
# c('met10','met9','met8','met7','met6','met5','met4','met3','met2','met1')
for (imputeMethod in imputeMethods){
  for (plateID in c('met13','met11','met10','met9','met8','met7','met6','met5','met4','met3','met2','met1')){
    # load the target library set
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    batches = levels(batchLabel)
    library(DESeq2)
    library(limma)
    library(edgeR)
    #mergeTbl = data.frame()
    for (i in 1:length(batches)){
      adaZmat = read.csv(paste('0_adaZ_score/outputs/adaZ_',batches[i],'.csv',sep = ''),row.names = 1)
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

      # update if it is met2_lib5
      if (batches[i] == 'met2_lib5'){
        # add back in the extra replicate
        loadRData <- function(fileName){
          #loads an RData file, and returns it
          load(fileName)
          newdata = list(input, batchLabel, RNAi)
          return(newdata)
        }
        extra_data <- loadRData('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_extra_met2_lib5.Rdata')
        # merge
        rep3_input = extra_data[[1]]
        rep3_RNAi = extra_data[[3]]
        keepind = str_detect(colnames(rep3_input),'_rep3_extra_met2_lib5$')
        rep3_input = rep3_input[,keepind]
        rep3_RNAi = rep3_RNAi[keepind]

        keepind = !str_detect(colnames(input_subset),'_rep3_met2_lib5$')
        input_subset = input_subset[,keepind]
        coldata = coldata[keepind,]

        input_subset <- input_subset %>% rownames_to_column('gene') %>%
          full_join(rep3_input %>% rownames_to_column('gene'), by = 'gene')
        rownames(input_subset) = input_subset$gene
        input_subset = input_subset[,-1]
        input_subset[is.na(input_subset)] = 0
        coldata2 = data.frame(row.names = colnames(rep3_input))
        coldata2$batchLabel = '_rep3'
        coldata2$RNAi = as.character(rep3_RNAi)
        coldata = rbind(coldata, coldata2)

        #also update the adaZmat
        adaZmat2 = read.csv(paste('0_adaZ_score/outputs/adaZ_extra_met2_lib1.csv',sep = ''),row.names = 1)
        adaZmat = adaZmat %>% rownames_to_column('gene') %>%
          full_join(adaZmat2 %>% rownames_to_column('gene'), by = 'gene')
        rownames(adaZmat) = adaZmat$gene
        adaZmat = adaZmat[,-1]
        adaZmat[is.na(adaZmat)] = 0
      }

      plate2RNAi_tmp = plate2RNAi[plate2RNAi$library %in% batches[i],]
      #list the DE for each individual RNAi
      if (any(levels(RNAi_subset) == 'x.vector')){
        myseq1 = 2:(length(levels(RNAi_subset)))
      }else{
        myseq1 = 1:(length(levels(RNAi_subset)))
      }
      for (j in myseq1){
        targetGene = levels(RNAi_subset)[j]

        if (targetGene %in% plate2RNAi$RNAi_name_final){
          coldata_tmp = coldata
          coldata_tmp$RNAi_ori = coldata_tmp$RNAi
          coldata_tmp$RNAi = factor(coldata_tmp$RNAi, levels = c('others',levels(coldata_tmp$RNAi)))
          coldata_tmp$RNAi[!(coldata_tmp$RNAi %in% c('x.vector',targetGene))] = 'others'
          input_subset_tmp = input_subset
          #coldata_tmp$bacPlate = 'plate1'
          #coldata_tmp$bacPlate[coldata_tmp$RNAi_ori %in% plate2RNAi_tmp$RNAi_name_final] = 'plate2'
          #coldata_tmp$bacPlate[coldata_tmp$RNAi_ori %in% 'x.vector'] = 'plate2'
          coldata_tmp$RNAi = factor(as.character(coldata_tmp$RNAi), levels = c('others', 'x.vector', targetGene))
          # remove extra plate2 in others (to avoid confounding)
          rmInd = coldata_tmp$RNAi %in% 'others' & coldata_tmp$RNAi_ori %in% plate2RNAi_tmp$RNAi_name_final
          input_subset_tmp = input_subset[,!rmInd]
          coldata_tmp = coldata_tmp[!rmInd,]

          dds <- DESeqDataSetFromMatrix(countData = input_subset_tmp,
                                        colData = coldata_tmp,
                                        design= ~ batchLabel + RNAi)
        }else{
          coldata_tmp = coldata
          coldata_tmp$RNAi_ori = coldata_tmp$RNAi
          coldata_tmp$RNAi = factor(coldata_tmp$RNAi, levels = c('others',levels(coldata_tmp$RNAi)))
          coldata_tmp$RNAi[!(coldata_tmp$RNAi %in% c('x.vector',targetGene))] = 'others'
          # remove all plate2 confounding samples
          rmInd = (coldata_tmp$RNAi %in% 'x.vector') |
            (coldata_tmp$RNAi %in% 'others' & coldata_tmp$RNAi_ori %in% plate2RNAi_tmp$RNAi_name_final)
          input_subset_tmp = input_subset[,!rmInd]
          coldata_tmp = coldata_tmp[!rmInd,]
          dds <- DESeqDataSetFromMatrix(countData = input_subset_tmp,
                                        colData = coldata_tmp,
                                        design= ~ batchLabel + RNAi)
        }

        # clean up the outliers
        # filtering
        keep <- rowSums(counts(dds)>=10) >= 1
        dds <- dds[keep,]
        dds = estimateSizeFactors(dds)
        dds <- cleanOutliers(dds, adaZmat, zcutoff, imputeMethod)
        # refiltering low counts
        keep <- rowSums(counts(dds)>=10) >= 1
        dds <- dds[keep,]

        # library("BiocParallel")
        # register(MulticoreParam(2))
        dds <- DESeq(dds,minReplicatesForReplace=7)
        # supply cooks filter metric (maximum cook for all samples related to tested covariate)
        dds.filt = dds


        if (targetGene %in% plate2RNAi$RNAi_name_final & any(dds.filt$RNAi %in% 'x.vector')){# vector exist and it is a plate2 confounded RNAi
          # update base mean to the basemean between RNA and control
          # we use DEseq2's default setting to optimize power, but we define the base mean by only vectors plus the RNAi in query
          libWiseBaseMean = mcols(dds.filt)$baseMean
          mycounts = counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('others','x.vector',targetGene)]
          myconditions = dds.filt$RNAi[dds.filt$RNAi %in% c('others','x.vector',targetGene)]
          myWeights = rep(0,length(myconditions))
          myWeights[myconditions == 'others'] = 1/3 / sum(myconditions == 'others')
          myWeights[myconditions == targetGene] = 1/3 / sum(myconditions == targetGene)
          myWeights[myconditions == 'x.vector'] = 1/3 / sum(myconditions == 'x.vector')

          baseMean = apply(mycounts, 1, result <- function(x){result = weighted.mean(x,myWeights)})
          # baseVar = rowVars(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)])
          mcols(dds.filt)$baseMean = baseMean
          mcols(dds.filt)$baseVar = NA
          mySample1 = paste('RNAi_',targetGene,'_vs_others',sep = '')
          mySample2 = paste('RNAi_x.vector_vs_others',sep = '')
          mySample = mySample1
          res1 <- results(dds.filt,name = str_replace_all(mySample1,'-','.'),independentFiltering = F) # we may use "replace" instead of cook filter for met library
          stat_res1 = res1$stat
          logFC_raw = res1$log2FoldChange
          res1 <- lfcShrink(dds.filt, coef = str_replace_all(mySample1,'-','.'), type="apeglm",apeMethod="nbinomC", res=res1)
          res1$stat = stat_res1
          res1$log2FoldChange_raw = logFC_raw

          res2 <- results(dds.filt,name = str_replace_all(mySample2,'-','.'),independentFiltering = F) # we may use "replace" instead of cook filter for met library
          stat_res2 = res2$stat
          logFC_raw = res2$log2FoldChange
          res2 <- lfcShrink(dds.filt, coef = str_replace_all(mySample2,'-','.'), type="apeglm",apeMethod="nbinomC", res=res2)
          res2$stat = stat_res2
          res2$log2FoldChange_raw = logFC_raw
          # filter out the case where the FC direction of vector (when significant [this is important to max power when we uniformly replace the high frequency outlier genes in all library]) is the same as target RNAi (confounded)
          # 03232022: We revised the definition of significance: we consider the significance as the RNAi of interest is not very different than vector, rather than checking the significance of vector itself
          # we mask the DE only if: (1) the DE is changing in the same direction as vector; (2) vector is significantly changed (p<0.01) and (3) the fold change of the DE is not very different than vector (2-fold in log scale)
          res1$pvalue[sign(res1$log2FoldChange) == sign(res2$log2FoldChange) &
                        res2$pvalue < 0.01 &
                        abs(res1$log2FoldChange) <= abs(res2$log2FoldChange) * 2] = NA
          res = res1
          # attach the counts information
          # filter by mean read counts of each component of contrast(will do in cutoff step)
          # batch effect is causing some trouble here, so we use median instead of mean; for met lib, we want to use mean
          # ==> change to median alway to ensure robustness
          mVals1 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
          names(mVals1) = rownames(dds.filt)
          mVals2 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'others'])
          names(mVals2) = rownames(dds.filt)
          res = as.data.frame(res)
          res$medianCount_RNAi = mVals1[rownames(res)]
          res$medianCount_ctr = mVals2[rownames(res)]
          res$libWiseBaseMean = libWiseBaseMean # this will be used in the MERGE modeling! (otherwise the pseudocounts will not be comparable across conditions)
        }else{# either not confounded or no vector, we directly compare with other RNAi to maximize power
          # update base mean to the basemean between RNA, others and vector
          # we use DEseq2's default setting to optimize power, but we define the base mean by only vectors plus the RNAi in query
          libWiseBaseMean = mcols(dds.filt)$baseMean
          mycounts = counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('others',targetGene)]
          myconditions = dds.filt$RNAi[dds.filt$RNAi %in% c('others',targetGene)]
          myWeights = rep(0,length(myconditions))
          myWeights[myconditions == 'others'] = 0.5 / sum(myconditions == 'others')
          myWeights[myconditions == targetGene] = 0.5 / sum(myconditions == targetGene)

          baseMean = apply(mycounts, 1, result <- function(x){result = weighted.mean(x,myWeights)})
          # baseVar = rowVars(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi %in% c('x.vector',targetGene)])
          mcols(dds.filt)$baseMean = baseMean
          mcols(dds.filt)$baseVar = NA
          mySample = paste('RNAi_',targetGene,'_vs_others',sep = '')
          res <- results(dds.filt,name = str_replace_all(mySample,'-','.'),independentFiltering = F) # we may use "replace" instead of cook filter for met library
          stat = res$stat
          logFC_raw = res$log2FoldChange
          res <- lfcShrink(dds.filt, coef = str_replace_all(mySample,'-','.'), type="apeglm",apeMethod="nbinomC", res=res)
          res$stat = stat
          res$log2FoldChange_raw = logFC_raw

          # attach the counts information
          # filter by mean read counts of each component of contrast(will do in cutoff step)
          # batch effect is causing some trouble here, so we use median instead of mean; for met lib, we want to use mean
          # ==> change to median alway to ensure robustness
          mVals1 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == targetGene])
          names(mVals1) = rownames(dds.filt)
          mVals2 = rowMedians(counts(dds.filt, normalized=TRUE)[,dds.filt$RNAi == 'others'])
          names(mVals2) = rownames(dds.filt)
          res = as.data.frame(res)
          res$medianCount_RNAi = mVals1[rownames(res)]
          res$medianCount_ctr = mVals2[rownames(res)]
          res$libWiseBaseMean = libWiseBaseMean # this will be used in the MERGE modeling! (otherwise the pseudocounts will not be comparable across conditions)
        }

        # save
        write.csv(res,paste('output/raw_DE_output/DE_table_one2all_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'_imputMethod_',imputeMethod,'.csv',sep = ''))
        # resOrdered <- res[order(res$padj),]
        # outputTbl = as.data.frame(subset(res, padj < 0.05))
        # if (nrow(outputTbl) > 0){
        #   outputTbl = cbind(data.frame(WBID = rownames(outputTbl)),outputTbl)
        #   rownames(outputTbl) = 1:nrow(outputTbl)
        #   outputTbl$RNAi = rep(targetGene,nrow(outputTbl))
        # }else{
        #   outputTbl = data.frame(WBID = 'NoHit')
        #   outputTbl$baseMean = NA
        #   outputTbl$log2FoldChange = NA
        #   outputTbl$lfcSE = NA
        #   outputTbl$pvalue = NA
        #   outputTbl$padj = NA
        #   outputTbl$medianCount_RNAi = NA
        #   outputTbl$medianCount_ctr = NA
        #   outputTbl$RNAi = targetGene
        #   outputTbl$libWiseBaseMean = NA
        #   rownames(outputTbl) = 1
        # }
        # outputTbl$batchID = batches[i]
        # mergeTbl = rbind(mergeTbl, outputTbl)
        print(mySample)
      }
    }
    #mergeTbl = mergeTbl[order(mergeTbl["RNAi"], mergeTbl["log2FoldChange"]),]
    #write.csv(mergeTbl,file = paste('output/DE_one2all_master_table_FDR005_',plateID,'.csv',sep = ''))
  }
}
writeLines(capture.output(sessionInfo()), "sessionInfo_DEseq2_one2all.txt")
