#lowExpCutoff = 30
cutoffFC <- function(plateID,suffix,imputeM){
  if(imputeM == ''){
    mergeTbl = read.csv(paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR005_',plateID,'.csv',sep = ''),row.names = 1);

  }else{
    mergeTbl = read.csv(paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR005_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = 1);

  }

  # first, cut off low-expression DE
  # mergeTbl = mergeTbl[(mergeTbl$medianCount_RNAi >= lowExpCutoff | mergeTbl$medianCount_ctr >= lowExpCutoff) | mergeTbl$WBID == 'NoHit', ]

  # second, label out the gene type
  relCV = read.csv('output/median_deviation_of_CV_to_the_medianCV.csv',row.names = 1)
  #highVarGenes_broad = read.csv('output/conditionwise_cv_outlierGenes_2MADfilter.csv')
  iCELGenes = read.csv('./../input_data/otherTbls/WormPaths_Tables/wormPathTable.csv')
  metabolicGenes = read.csv('./../input_data/otherTbls/allMetGenes.txt')
  mergeTbl$isICEL = mergeTbl$WBID %in% iCELGenes$WormBase.ID
  mergeTbl$isMetabolic = mergeTbl$WBID %in% metabolicGenes$WBName
  mergeTbl$relCV = relCV$medianRelCV[match(mergeTbl$WBID, rownames(relCV))]

  # (TBD), label out the development confounded genes

  ###iCEL pathway annotation of each nhr DE genes
  library(stringr)
  library(tidyr)
  library(dbplyr)

  iCELGenes = read.csv('./../input_data/otherTbls/WormPaths_Tables/wormPathTable.csv')
  mergeTbl$LEVEL1 = iCELGenes$LEVEL.1[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
  mergeTbl$LEVEL1[is.na(mergeTbl$LEVEL1)] = 'Non-iCEL'
  mergeTbl$LEVEL2 = iCELGenes$LEVEL.2[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
  mergeTbl$LEVEL2[is.na(mergeTbl$LEVEL2)] = 'Non-iCEL'
  mergeTbl$LEVEL3 = iCELGenes$LEVEL.3[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
  mergeTbl$LEVEL3[is.na(mergeTbl$LEVEL3)] = 'Non-iCEL'
  mergeTbl$LEVEL4 = iCELGenes$LEVEL.4[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
  mergeTbl$LEVEL4[is.na(mergeTbl$LEVEL4)] = 'Non-iCEL'

  RNAi = unique(mergeTbl$RNAi)
  RNAi_ID = str_replace(RNAi, '^RNAi','')


  DEtbl = mergeTbl
  ETCgenes = iCELGenes$WormBase.ID[grep('ELECTRON TRANSPORT CHAIN', iCELGenes$LEVEL.2)]
  DE_ETC = DEtbl[is.element(DEtbl$WBID,ETCgenes), ]
  #hist(DE_ETC$log2FoldChange,breaks = 30)
  #abline(v = -log2(2))

  #hist(DEtbl$log2FoldChange,breaks = 30)
  #abline(v = -log2(2))

  a = DE_ETC[abs(DE_ETC$log2FoldChange) >= log2(2),]
  length(unique(a$RNAi))
  length(unique(DE_ETC$RNAi))

  #hist(table(DE_ETC$RNAi),breaks = 30)
  #hist(table(a$RNAi),breaks = 30)

  print(sum(abs(DEtbl$log2FoldChange) >= log2(2),na.rm = T) / nrow(DEtbl))
  print(sum(abs(DE_ETC$log2FoldChange) >= log2(2),na.rm = T) / nrow(DE_ETC))
  sum(abs(DE_ETC$log2FoldChange) >= log2(2),na.rm = T)
  sum(abs(DEtbl$log2FoldChange) >= log2(2),na.rm = T)

  DEtbl = mergeTbl
  DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(1.5) & DEtbl$padj <= 0.05) | DEtbl$WBID == 'NoHit', ]
  # add back no-hit
  noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
  noHitRNAi_batch = c()
  noHitRNAi_new = c()# if a noHit RNAi appears twice
  if (length(noHitRNAi) > 0){
    for (i in 1:length(noHitRNAi)){
      noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[i]]))
      noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[i])
    }
  }
  if (length(noHitRNAi_new) > 0){
    noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
    noHit$baseMean = NA
    noHit$log2FoldChange = NA
    noHit$log2FoldChange_raw = NA
    noHit$lfcSE= NA
    noHit$stat = NA
    noHit$pvalue = NA
    noHit$padj = NA
    noHit$medianCount_RNAi = NA
    noHit$medianCount_ctr = NA
    noHit$libWiseBaseMean = NA
    noHit$RNAi= noHitRNAi_new
    noHit$batchID = noHitRNAi_batch
    noHit$isICEL = NA
    noHit$isMetabolic = NA
    noHit$relCV = NA
    noHit$LEVEL1 = NA
    noHit$LEVEL2 = NA
    noHit$LEVEL3 = NA
    noHit$LEVEL4 = NA
    DEtbl = rbind(DEtbl, noHit)
  }
  WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
  DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
  DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
  if (imputeM == ''){
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR005_FC1.5_',plateID,'.csv',sep = ''),row.names = F)

  }else{
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR005_FC1.5_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)

  }
  hist(log10(table(DEtbl$RNAi)+1))

  DEtbl = mergeTbl
  DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(1.5) & DEtbl$padj <= 0.01) | DEtbl$WBID == 'NoHit', ]# add back no-hit
  # add back no-hit
  noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
  noHitRNAi_batch = c()
  noHitRNAi_new = c()# if a noHit RNAi appears twice
  if (length(noHitRNAi) > 0){
    for (i in 1:length(noHitRNAi)){
      noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[i]]))
      noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[i])
    }
  }
  if (length(noHitRNAi_new) > 0){
    noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
    noHit$baseMean = NA
    noHit$log2FoldChange = NA
    noHit$log2FoldChange_raw = NA
    noHit$lfcSE= NA
    noHit$stat = NA
    noHit$pvalue = NA
    noHit$padj = NA
    noHit$medianCount_RNAi = NA
    noHit$medianCount_ctr = NA
    noHit$libWiseBaseMean = NA
    noHit$RNAi= noHitRNAi_new
    noHit$batchID = noHitRNAi_batch
    noHit$isICEL = NA
    noHit$isMetabolic = NA
    noHit$relCV = NA
    noHit$LEVEL1 = NA
    noHit$LEVEL2 = NA
    noHit$LEVEL3 = NA
    noHit$LEVEL4 = NA
    DEtbl = rbind(DEtbl, noHit)
  }
  WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
  DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
  DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
  if (imputeM == ''){
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR001_FC1.5_',plateID,'.csv',sep = ''),row.names = F)

  }else{
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR001_FC1.5_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)

  }
  hist(log10(table(DEtbl$RNAi)+1))

  # IN THIS THERESHOLD, we also additionally cutoff the basemean (weighted) to ensure the robustness
  DEtbl = mergeTbl
  DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(2) & DEtbl$padj <= 0.05 & DEtbl$baseMean >= 10) | DEtbl$WBID == 'NoHit', ]
  # add back no-hit
  noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
  noHitRNAi_batch = c()
  noHitRNAi_new = c()# if a noHit RNAi appears twice
  if (length(noHitRNAi) > 0){
    for (i in 1:length(noHitRNAi)){
      noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[i]]))
      noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[i])
    }
  }
  if (length(noHitRNAi_new) > 0){
    noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
    noHit$baseMean = NA
    noHit$log2FoldChange = NA
    noHit$log2FoldChange_raw = NA
    noHit$lfcSE= NA
    noHit$stat = NA
    noHit$pvalue = NA
    noHit$padj = NA
    noHit$medianCount_RNAi = NA
    noHit$medianCount_ctr = NA
    noHit$libWiseBaseMean = NA
    noHit$RNAi= noHitRNAi_new
    noHit$batchID = noHitRNAi_batch
    noHit$isICEL = NA
    noHit$isMetabolic = NA
    noHit$relCV = NA
    noHit$LEVEL1 = NA
    noHit$LEVEL2 = NA
    noHit$LEVEL3 = NA
    noHit$LEVEL4 = NA
    DEtbl = rbind(DEtbl, noHit)
  }
  WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
  DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
  DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
  if (imputeM == ''){
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR005_FC2_',plateID,'.csv',sep = ''),row.names = F)

  }else{
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR005_FC2_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)

  }
  hist(log10(table(DEtbl$RNAi)+1))

  DEtbl = mergeTbl
  DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(2) & DEtbl$padj <= 0.01 & DEtbl$baseMean >= 10) | DEtbl$WBID == 'NoHit', ]
  # add back no-hit
  noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
  noHitRNAi_batch = c()
  noHitRNAi_new = c()# if a noHit RNAi appears twice
  if (length(noHitRNAi) > 0){
    for (i in 1:length(noHitRNAi)){
      noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[i]]))
      noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[i])
    }
  }
  if (length(noHitRNAi_new) > 0){
    noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
    noHit$baseMean = NA
    noHit$log2FoldChange = NA
    noHit$log2FoldChange_raw = NA
    noHit$lfcSE= NA
    noHit$stat = NA
    noHit$pvalue = NA
    noHit$padj = NA
    noHit$medianCount_RNAi = NA
    noHit$medianCount_ctr = NA
    noHit$libWiseBaseMean = NA
    noHit$RNAi= noHitRNAi_new
    noHit$batchID = noHitRNAi_batch
    noHit$isICEL = NA
    noHit$isMetabolic = NA
    noHit$relCV = NA
    noHit$LEVEL1 = NA
    noHit$LEVEL2 = NA
    noHit$LEVEL3 = NA
    noHit$LEVEL4 = NA
    DEtbl = rbind(DEtbl, noHit)
  }
  WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
  DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
  DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
  if (imputeM == ''){
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR001_FC2_',plateID,'.csv',sep = ''),row.names = F)

  }else{
    write.csv(DEtbl,file = paste('output/raw_DE_output_tables/DE_',suffix,'master_table_FDR001_FC2_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)

  }
  hist(log10(table(DEtbl$RNAi)+1))
  #
  # # cutoff the ANOVA type DE
  # mergeTbl = read.csv(paste('output/DE_master_table_ANNOVA_DEseq2_FDR005_',plateID,'.csv',sep = ''),row.names = 1);
  # mergeTbl = mergeTbl[(mergeTbl$maxMedianCount > lowExpCutoff & abs(mergeTbl$log2FoldChange) >= log2(2) & mergeTbl$padj <= 0.01)
  #                     | mergeTbl$WBID == 'NoHit', ]
  # write.csv(mergeTbl,file = paste('output/DE_master_table_ANNOVA_DEseq2_FDR005_FC2_',plateID,'.csv',sep = ''),row.names = F)
  #
}
mergeDE_v2 <- function(plateIDs,cutoffs,gene2rmList2, writeDetail,imputeM,FCtype){
  for (plateID in plateIDs){
    mergeTbl_list = list()
    for (cutoff in cutoffs){
      mergeTbl_list[[paste('cutoff',cutoff,sep = '_')]] = data.frame()
    }
    # load the target library set
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    batches = levels(batchLabel)
    for (i in 1:length(batches)){
      subsetInd = batchLabel == batches[i]
      RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
      if (any(levels(RNAi_subset) == 'x.vector')){
        myseq1 = 2:(length(levels(RNAi_subset)))
      }else{
        myseq1 = 1:(length(levels(RNAi_subset)))
      }
      for (j in myseq1){
        targetGene = levels(RNAi_subset)[j]
        mySample1 = paste('RNAi_',targetGene,'_vs_x.vector',sep = '')
        mySample2 = paste('RNAi_',targetGene,'_vs_others',sep = '')
        if(batches[i] != 'met3_lib4'){
          outputTbl1 = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_',mySample1,'_',batches[i],'.csv',sep = ''),row.names = 1)
          outputTbl1$DE_source = 'vs. vector'
        }else{
          outputTbl1 = read.csv(paste('output/raw_DE_output/DE_table_one2all_lfcShrink_RNAi_raw_',mySample2,'_',batches[i],'_imputMethod_',imputeM,'.csv',sep = ''),row.names = 1)
          outputTbl1 = outputTbl1[,c("baseMean","log2FoldChange","log2FoldChange_raw","lfcSE","stat","pvalue","padj","medianCount_RNAi","medianCount_ctr","libWiseBaseMean")]
          outputTbl1$DE_source = 'vs. other RNAi (same bacPlate)'
        }
        outputTbl2 = read.csv(paste('output/raw_DE_output/DE_table_one2all_lfcShrink_RNAi_raw_',mySample2,'_',batches[i],'_imputMethod_',imputeM,'.csv',sep = ''),row.names = 1)
        outputTbl2$DE_source = 'vs. other RNAi (same bacPlate)'
        # apply the existing mask
        #outputTbl1$pvalue[is.na(outputTbl1$padj)] = NA
        #outputTbl2$pvalue[is.na(outputTbl2$padj)] = NA

        # we merge the two DE result and readjust the FDR
        # first, overide with unique DE identified by one-all only (lower FDR and higher FC)
        # add unique rows in one2all
        outputTbl1 = rbind(outputTbl1,outputTbl2[setdiff(rownames(outputTbl2),rownames(outputTbl1)),])
        commomG = intersect(rownames(outputTbl2),rownames(outputTbl1))
        # replace when p value is much smaller and fc is larger (more rep adds power and correct for some problematic vector cases)
        repInd = outputTbl2[commomG,'pvalue'] < 0.01 * outputTbl1[commomG,'pvalue'] &
          abs(outputTbl2[commomG,FCtype]) > abs(outputTbl1[commomG,FCtype])
        repG = setdiff(commomG[repInd],NA)
        outputTbl1[repG,] = outputTbl2[repG,colnames(outputTbl1)]
        # replace when one2one is marked by NA
        # repG2 = intersect(rownames(outputTbl1)[is.na(outputTbl1$pvalue)],rownames(outputTbl2))
        # outputTbl1[repG2,] = outputTbl2[repG2,]
        #
        print(paste('percentage of DE overrided (or added): ',sum(outputTbl1$DE_source == 'vs. other RNAi (same bacPlate)')/nrow(outputTbl1)))

        # replace the bad genes with one2all result
        for (cutoff in cutoffs){
          gene2rm = gene2rmList2[[paste('cutoff',cutoff,sep = '_')]][[plateID]][[batches[i]]]
          ind1 = !(rownames(outputTbl1) %in% gene2rm)
          ind2 = rownames(outputTbl2) %in% gene2rm
          outputTbl_merge = rbind(outputTbl1[ind1, ],
                                  outputTbl2[ind2, colnames(outputTbl1)])
          nfilt = sum(ind1)
          nrescue = sum(ind2)
          print(paste('percentage of bad gene: ',nrescue/(nfilt+nrescue)))

          #redo FDR calculation
          outputTbl_merge$padj = p.adjust(outputTbl_merge$pvalue,method = 'BH')
          if (writeDetail){
            write.csv(outputTbl_merge,paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',cutoff,'_RNAi_', targetGene,'_',batches[i],'_imputMethod_',imputeM,'_',FCtype,'.csv',sep = ''))
          }

          #collect the FDR005 table
          resOrdered <- outputTbl_merge[order(outputTbl_merge$padj),]
          outputTbl_merge = as.data.frame(subset(outputTbl_merge, padj < 0.05))
          if (nrow(outputTbl_merge) > 0){
            outputTbl_merge = cbind(data.frame(WBID = rownames(outputTbl_merge)),outputTbl_merge)
            rownames(outputTbl_merge) = 1:nrow(outputTbl_merge)
            outputTbl_merge$RNAi = rep(targetGene,nrow(outputTbl_merge))
          }else{
            outputTbl_merge = data.frame(WBID = 'NoHit')
            outputTbl_merge$baseMean = NA
            outputTbl_merge$log2FoldChange = NA
            outputTbl_merge$log2FoldChange_raw = NA
            outputTbl_merge$lfcSE = NA
            outputTbl_merge$stat = NA
            outputTbl_merge$pvalue = NA
            outputTbl_merge$padj = NA
            outputTbl_merge$medianCount_RNAi = NA
            outputTbl_merge$medianCount_ctr = NA
            outputTbl_merge$RNAi = targetGene
            outputTbl_merge$libWiseBaseMean = NA
            outputTbl_merge$DE_source = 'both'
            rownames(outputTbl_merge) = 1
          }
          outputTbl_merge$batchID = batches[i]
          mergeTbl_list[[paste('cutoff',cutoff,sep = '_')]] = rbind(mergeTbl_list[[paste('cutoff',cutoff,sep = '_')]], outputTbl_merge)
        }
        print(paste(targetGene,batches[i]))
      }
    }

    # save the merged tables with diff fc and fdr cutoff
    for (cutoff in cutoffs){
      suffix = paste('merged_clean_cutoff_',cutoff,'_',sep = '')
      mergeTbl = mergeTbl_list[[paste('cutoff',cutoff,sep = '_')]]
      # second, label out the gene type
      relCV = read.csv('output/median_deviation_of_CV_to_the_medianCV.csv',row.names = 1)
      #highVarGenes_broad = read.csv('output/conditionwise_cv_outlierGenes_2MADfilter.csv')
      iCELGenes = read.csv('./../input_data/otherTbls/WormPaths_Tables/wormPathTable.csv')
      metabolicGenes = read.csv('./../input_data/otherTbls/allMetGenes.txt')
      mergeTbl$isICEL = mergeTbl$WBID %in% iCELGenes$WormBase.ID
      mergeTbl$isMetabolic = mergeTbl$WBID %in% metabolicGenes$WBName
      mergeTbl$relCV = relCV$medianRelCV[match(mergeTbl$WBID, rownames(relCV))]

      DEtbl = mergeTbl
      WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
      DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
      DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
      if (length(cutoffs)>1){
        write.csv(DEtbl,file = paste('output/cutoffTitration/DE_',suffix,'master_table_FDR005_',plateID,'_imputMethod_',imputeM,'_',FCtype,'.csv',sep = ''),row.names = F)

      }else{
        # write.csv(DEtbl,file = paste('output/clean_DE_output_tables/DE_',suffix,'master_table_FDR005_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      }


      # save other cutoffed tbles

      # DEtbl = mergeTbl
      # DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(1.5) & DEtbl$padj <= 0.05) | DEtbl$WBID == 'NoHit', ]
      # # add back no-hit
      # noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
      # noHitRNAi_batch = c()
      # noHitRNAi_new = c()# if a noHit RNAi appears twice
      # if (length(noHitRNAi) > 0){
      #   for (ii in 1:length(noHitRNAi)){
      #     noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[ii]]))
      #     noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[ii])
      #   }
      # }
      # if (length(noHitRNAi_new) > 0){
      #   noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
      #   noHit$baseMean = NA
      #   noHit$log2FoldChange = NA
      #   noHit$lfcSE= NA
      #   noHit$pvalue = NA
      #   noHit$padj = NA
      #   noHit$medianCount_RNAi = NA
      #   noHit$medianCount_ctr = NA
      #   noHit$libWiseBaseMean = NA
      #   noHit$RNAi= noHitRNAi_new
      #   noHit$batchID = noHitRNAi_batch
      #   noHit$DE_source = 'both'
      #   noHit$isICEL = NA
      #   noHit$isMetabolic = NA
      #   noHit$relCV = NA
      #   DEtbl = rbind(DEtbl, noHit)
      # }
      # WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
      # DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
      # DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
      # if (length(cutoffs)>1){
      #   write.csv(DEtbl,file = paste('output/cutoffTitration/DE_',suffix,'master_table_FDR005_FC1.5_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      #
      # }else{
      #   write.csv(DEtbl,file = paste('output/clean_DE_output_tables/DE_',suffix,'master_table_FDR005_FC1.5_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      # }
      #
      # DEtbl = mergeTbl
      # DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(1.5) & DEtbl$padj <= 0.01) | DEtbl$WBID == 'NoHit', ]
      # # add back no-hit
      # noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
      # noHitRNAi_batch = c()
      # noHitRNAi_new = c()# if a noHit RNAi appears twice
      # if (length(noHitRNAi) > 0){
      #   for (ii in 1:length(noHitRNAi)){
      #     noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[ii]]))
      #     noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[ii])
      #   }
      # }
      # if (length(noHitRNAi_new) > 0){
      #   noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
      #   noHit$baseMean = NA
      #   noHit$log2FoldChange = NA
      #   noHit$lfcSE= NA
      #   noHit$pvalue = NA
      #   noHit$padj = NA
      #   noHit$medianCount_RNAi = NA
      #   noHit$medianCount_ctr = NA
      #   noHit$libWiseBaseMean = NA
      #   noHit$RNAi= noHitRNAi_new
      #   noHit$batchID = noHitRNAi_batch
      #   noHit$DE_source = 'both'
      #   noHit$isICEL = NA
      #   noHit$isMetabolic = NA
      #   noHit$relCV = NA
      #   DEtbl = rbind(DEtbl, noHit)
      # }
      # WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
      # DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
      # DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
      # if (length(cutoffs)>1){
      #   write.csv(DEtbl,file = paste('output/cutoffTitration/DE_',suffix,'master_table_FDR001_FC1.5_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      #
      # }else{
      #   write.csv(DEtbl,file = paste('output/clean_DE_output_tables/DE_',suffix,'master_table_FDR001_FC1.5_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      # }
      # DEtbl = mergeTbl
      # DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(2) & DEtbl$padj <= 0.05 & DEtbl$baseMean >= 10) | DEtbl$WBID == 'NoHit', ]
      # # add back no-hit
      # noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
      # noHitRNAi_batch = c()
      # noHitRNAi_new = c()# if a noHit RNAi appears twice
      # if (length(noHitRNAi) > 0){
      #   for (ii in 1:length(noHitRNAi)){
      #     noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[ii]]))
      #     noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[ii])
      #   }
      # }
      # if (length(noHitRNAi_new) > 0){
      #   noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
      #   noHit$baseMean = NA
      #   noHit$log2FoldChange = NA
      #   noHit$lfcSE= NA
      #   noHit$pvalue = NA
      #   noHit$padj = NA
      #   noHit$medianCount_RNAi = NA
      #   noHit$medianCount_ctr = NA
      #   noHit$libWiseBaseMean = NA
      #   noHit$RNAi= noHitRNAi_new
      #   noHit$batchID = noHitRNAi_batch
      #   noHit$DE_source = 'both'
      #   noHit$isICEL = NA
      #   noHit$isMetabolic = NA
      #   noHit$relCV = NA
      #   DEtbl = rbind(DEtbl, noHit)
      # }
      # WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
      # DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
      # DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
      # if (length(cutoffs)>1){
      #   write.csv(DEtbl,file = paste('output/cutoffTitration/DE_',suffix,'master_table_FDR005_FC2_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      #
      # }else{
      #   write.csv(DEtbl,file = paste('output/clean_DE_output_tables/DE_',suffix,'master_table_FDR005_FC2_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      # }
      #
      # DEtbl = mergeTbl
      # DEtbl = DEtbl[(abs(DEtbl$log2FoldChange) >= log2(2) & DEtbl$padj <= 0.01 & DEtbl$baseMean >= 10) | DEtbl$WBID == 'NoHit', ]
      # # add back no-hit
      # noHitRNAi = setdiff(unique(mergeTbl$RNAi), unique(DEtbl$RNAi))
      # noHitRNAi_batch = c()
      # noHitRNAi_new = c()# if a noHit RNAi appears twice
      # if (length(noHitRNAi) > 0){
      #   for (ii in 1:length(noHitRNAi)){
      #     noHitRNAi_batch = c(noHitRNAi_batch, unique(mergeTbl$batchID[mergeTbl$RNAi == noHitRNAi[ii]]))
      #     noHitRNAi_new = c(noHitRNAi_new,noHitRNAi[ii])
      #   }
      # }
      # if (length(noHitRNAi_new) > 0){
      #   noHit = data.frame(WBID = rep('NoHit',length(noHitRNAi_new)))
      #   noHit$baseMean = NA
      #   noHit$log2FoldChange = NA
      #   noHit$lfcSE= NA
      #   noHit$pvalue = NA
      #   noHit$padj = NA
      #   noHit$medianCount_RNAi = NA
      #   noHit$medianCount_ctr = NA
      #   noHit$libWiseBaseMean = NA
      #   noHit$RNAi= noHitRNAi_new
      #   noHit$batchID = noHitRNAi_batch
      #   noHit$DE_source = 'both'
      #   noHit$isICEL = NA
      #   noHit$isMetabolic = NA
      #   noHit$relCV = NA
      #   DEtbl = rbind(DEtbl, noHit)
      # }
      # WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
      # DEtbl$Gene_name = WBID$Public.Name[match(DEtbl$WBID,WBID$WormBase.Gene.ID)]
      # DEtbl = DEtbl[,c('RNAi','Gene_name',setdiff(colnames(DEtbl),c('Gene_name','RNAi','lfcSE')))]
      # if (length(cutoffs)>1){
      #   write.csv(DEtbl,file = paste('output/cutoffTitration/DE_',suffix,'master_table_FDR001_FC2_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      #
      # }else{
      #   write.csv(DEtbl,file = paste('output/clean_DE_output_tables/DE_',suffix,'master_table_FDR001_FC2_',plateID,'_imputMethod_',imputeM,'.csv',sep = ''),row.names = F)
      # }
    }
  }
}
findgenes_pvalue_merge2 <- function(plateIDs, cutoffs,type,writeCore, freqCutoff,FCtype){
  # IN THE FINAL SETUP, WE DIDNT USE FREQUENCY CUTOFF FUNCTION, SO THERE IS NO CORE SET. WE SET THE CUTOFF TO INF TO DISABLIZE THE FUNCTION
  genelist_all_cutoffs = list()
  for (i in 1:length(cutoffs)){
    genelist_all_cutoffs[[paste('cutoff_',as.character(cutoffs[i]),sep = '')]] = list()
    for (j in 1:length(plateIDs)){
      genelist_all_cutoffs[[paste('cutoff_',as.character(cutoffs[i]),sep = '')]][[plateIDs[j]]] = list()
    }
  }
  for (plateID in plateIDs){
    library(limma)
    library(edgeR)
    library(DESeq2)
    library(stringr)
    # load the FC matrix from DE
    load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
    for (lib in levels(batchLabel)){
      if (lib != 'met3_lib4'){
        RNAis = droplevels(RNAi[batchLabel == lib])
        RNAis = setdiff(levels(RNAis),'x.vector')
        tbl = data.frame(row.names = rownames(input))
        for (thisRNAi in RNAis){
          tmp = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_RNAi_',thisRNAi,'_vs_x.vector_',lib,'.csv',sep = ''))
          tbl[tmp$X,thisRNAi] = -log10(tmp$pvalue) * sign(tmp[,FCtype])
        }
        Qs = rowQuantiles(x = as.matrix(tbl),probs = c(0.75,0.5,0.25), na.rm = T)
        upperQ = Qs[,1]
        midQ = Qs[,2]
        lowerQ = Qs[,3]
        for (cutoff in cutoffs){
          cutoff_lowerQ = cutoff * 10
          #cutoff_midQ = cutoff * 2
          #cutoff_upperQ = cutoff
          cutoff_midQ = cutoff

          if(type =='pvalue'){
            # badInd = which(
            #   (upperQ > -log10(cutoff_upperQ) &
            #   midQ > -log10(cutoff_midQ) &
            #   lowerQ > -log10(cutoff_lowerQ)) |
            #   (lowerQ < log10(cutoff_upperQ) &
            #    midQ < log10(cutoff_midQ) &
            #   upperQ < log10(cutoff_lowerQ)))
            badInd = which(
              (midQ > -log10(cutoff_midQ) &
                 lowerQ > -log10(cutoff_lowerQ)) |
                (upperQ < log10(cutoff_lowerQ) &
                   midQ < log10(cutoff_midQ)))
            badGenes = names(midQ)[badInd]
          }else if(type =='fdr'){
            stop('not done')
            badGenes = tbl$X[tbl$padj<cutoff]
          }
          badGenes = setdiff(badGenes,NA)
          genelist_all_cutoffs[[paste('cutoff_',as.character(cutoff),sep = '')]][[plateID]][[as.character(lib)]] = badGenes
        }
      }else{
        for (cutoff in cutoffs){
          genelist_all_cutoffs[[paste('cutoff_',as.character(cutoff),sep = '')]][[plateID]][[as.character(lib)]] = c()
        }
      }
    }
  }

  for (z in 1:length(cutoffs)){
    genelist2 = table(unlist(genelist_all_cutoffs[[z]]))
    coreset = names(genelist2)[genelist2 > freqCutoff]
    print(paste('at cutoff ', cutoffs[z],', the core set size is ',length(coreset),sep = ''))
    if (writeCore){
      coresetTbl = data.frame(gene = coreset)
      coresetTbl$freq = genelist2[coreset]
      write.csv(coresetTbl, file = paste('output/coreset_bad_genes_p_cutoff_',cutoffs[z],'_',FCtype,'.csv',sep = ''))
    }
    #  NOTE: we found that the pvalues for plate-wise coreset genes showed tendency for smaller p values within the plate; but the
    # whole-dataset coreset has no such tendency; therefore, it is true that if a gene is bad in multiple library in a plate, it is more
    # or less bad in all libraries of that plate; but even if a gene is bad in many libraries of the 60 libs, it doesnt mean this gene is
    # more or less bad in all 60 libraries.
    for (i in 1:length(genelist_all_cutoffs[[z]])){
      # genelist2 = table(unlist(genelist_all_cutoffs[[z]][[i]]))
      # coreset = names(genelist2)[genelist2 >= freqCutoff]
      # print(paste('at cutoff ', cutoffs[z],', the core set size of ',names(genelist_all_cutoffs[[z]])[i],' is ',length(coreset),sep = ''))
      # if (writeCore){
      #   # add in the pvalue table
      #   coresetTbl = data.frame(gene = coreset)
      #   coresetTbl$freq = genelist2[coreset]
      #   write.csv(coresetTbl, file = paste('output/coreset_bad_genes_p_cutoff_',cutoffs[z],'_',names(genelist_all_cutoffs[[z]])[i],'.csv',sep = ''))
      # }
      for (j in 1:length(genelist_all_cutoffs[[z]][[i]])){
        genelist_all_cutoffs[[z]][[i]][[j]] = union(genelist_all_cutoffs[[z]][[i]][[j]], coreset)
      }
    }
  }

  return(genelist_all_cutoffs)
}
####################### part 0: without cleaning results: independent filtering ################
# titrate cutoff
# titrating on all conditions merged is not good (optimal cutoff too small, so biased) and also complicated
# we do the titration on individual RNAi (so same as DEseq2), however, we use a uniform cutoff based on our evaluation
#mergeTbl = data.frame()
library(DESeq2)
library(edgeR)
seq1 = 0:100
N_hit_max = data.frame(row.names = seq1)
N_hit_mean = data.frame(row.names = seq1)
for (plateID in c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')){
  # load the target library set
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  batches = levels(batchLabel)
  for (i in 1:length(batches)){
    if (batches[i] != 'met3_lib4'){
      subsetInd = batchLabel == batches[i]
      RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
      for (j in 2:(length(levels(RNAi_subset)))){
        targetGene = levels(RNAi_subset)[j]
        mySample = paste('RNAi_',targetGene,'_vs_x.vector',sep = '')
        outputTbl = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'.csv',sep = ''))
        pVector = outputTbl$pvalue
        maxExp = apply(outputTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
        baseMean = outputTbl$baseMean
        N_hit_max_tmp = numeric()
        N_hit_mean_tmp = numeric()
        for (lowExpCutoff in seq1){
          tmp1 = p.adjust(pVector[maxExp >= lowExpCutoff],method = 'BH')
          N_hit_max_tmp = c(N_hit_max_tmp,sum(tmp1 < 0.05))
          tmp2 = p.adjust(pVector[baseMean >= lowExpCutoff],method = 'BH')
          N_hit_mean_tmp = c(N_hit_mean_tmp,sum(tmp2 < 0.05))
        }
        N_hit_max[targetGene] = N_hit_max_tmp
        N_hit_mean[targetGene] = N_hit_mean_tmp
        #mergeTbl = rbind(mergeTbl, outputTbl)
        print(paste(mySample,batches[i]))
      }
    }
  }
}
quantile(N_hit_max[1,])
hist(log2(1+as.numeric(N_hit_max[1,])),breaks = 100)
# check the median N_hit
median_RNAi_max = rowMedians(as.matrix(N_hit_max))
median_RNAi_mean = rowMedians(as.matrix(N_hit_mean))
mean_RNAi = rowMeans(as.matrix(N_hit_mean))
quant_RNAi = apply((as.matrix(N_hit_mean)), 1, function(x){quantile(x,0.75)})
plot(0:100, median_RNAi_max,col = 'blue')
points(0:100, median_RNAi_mean,col = 'red')
abline(v = 30)
#abline(h = 0.95)

# apply filters and write out the new tables
lowExpCutoff = 30
for (plateID in c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')){
  mergeTbl = data.frame()
  # load the target library set
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  batches = levels(batchLabel)
  for (i in 1:length(batches)){
    if (batches[i] != 'met3_lib4'){
      subsetInd = batchLabel == batches[i]
      RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
      for (j in 2:(length(levels(RNAi_subset)))){
        targetGene = levels(RNAi_subset)[j]
        mySample = paste('RNAi_',targetGene,'_vs_x.vector',sep = '')
        outputTbl = read.csv(paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'.csv',sep = ''),row.names = 1)
        pVector = outputTbl$pvalue
        maxExp = apply(outputTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)

        padj = p.adjust(pVector[maxExp >= lowExpCutoff],method = 'BH')
        outputTbl$padj[maxExp >= lowExpCutoff] = padj
        outputTbl$padj[maxExp < lowExpCutoff] = NA

        write.csv(outputTbl,paste('output/raw_DE_output/DE_table_lfcShrink_RNAi_alpha005_',mySample,'_',batches[i],'.csv',sep = ''))
        resOrdered <- outputTbl[order(outputTbl$padj),]
        outputTbl = as.data.frame(subset(outputTbl, padj < 0.05))
        if (nrow(outputTbl) > 0){
          outputTbl = cbind(data.frame(WBID = rownames(outputTbl)),outputTbl)
          rownames(outputTbl) = 1:nrow(outputTbl)
          outputTbl$RNAi = rep(targetGene,nrow(outputTbl))
        }else{
          outputTbl = data.frame(WBID = 'NoHit')
          outputTbl$baseMean = NA
          outputTbl$log2FoldChange = NA
          outputTbl$log2FoldChange_raw = NA
          outputTbl$lfcSE = NA
          outputTbl$stat = NA
          outputTbl$pvalue = NA
          outputTbl$padj = NA
          outputTbl$medianCount_RNAi = NA
          outputTbl$medianCount_ctr = NA
          outputTbl$RNAi = targetGene
          outputTbl$libWiseBaseMean = NA
          rownames(outputTbl) = 1
        }
        outputTbl$batchID = batches[i]
        mergeTbl = rbind(mergeTbl, outputTbl)
        print(paste(mySample,batches[i]))
      }
    }
  }
  mergeTbl = mergeTbl[order(mergeTbl["RNAi"], mergeTbl["log2FoldChange"]),]
  write.csv(mergeTbl,file = paste('output/raw_DE_output_tables/DE_master_table_FDR005_',plateID,'.csv',sep = ''))
}

# one2all
seq1 = 0:100
N_hit_max = data.frame(row.names = seq1)
N_hit_mean = data.frame(row.names = seq1)
for (plateID in c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')){
  # load the target library set
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  batches = levels(batchLabel)
  for (i in 1:length(batches)){
    if (batches[i] != 'met3_lib4'){
      subsetInd = batchLabel == batches[i]
      RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
      for (j in 2:(length(levels(RNAi_subset)))){
        targetGene = levels(RNAi_subset)[j]
        mySample = paste('RNAi_',targetGene,'_vs_others',sep = '')
        outputTbl = read.csv(paste('output/raw_DE_output/DE_table_one2all_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'_imputMethod_random_strat.csv',sep = ''))
        pVector = outputTbl$pvalue
        maxExp = apply(outputTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
        baseMean = outputTbl$baseMean
        N_hit_max_tmp = numeric()
        N_hit_mean_tmp = numeric()
        for (lowExpCutoff in seq1){
          tmp1 = p.adjust(pVector[maxExp >= lowExpCutoff],method = 'BH')
          N_hit_max_tmp = c(N_hit_max_tmp,sum(tmp1 < 0.05,na.rm = T))
          tmp2 = p.adjust(pVector[baseMean >= lowExpCutoff],method = 'BH')
          N_hit_mean_tmp = c(N_hit_mean_tmp,sum(tmp2 < 0.05,na.rm = T))
        }
        N_hit_max[targetGene] = N_hit_max_tmp
        N_hit_mean[targetGene] = N_hit_mean_tmp
        #mergeTbl = rbind(mergeTbl, outputTbl)
        print(paste(mySample,batches[i]))
      }
    }
  }
}
quantile(N_hit_max[1,])
hist(log2(1+as.numeric(N_hit_max[1,])),breaks = 100)
# check the median N_hit
median_RNAi_max = rowMedians(as.matrix(N_hit_max))
median_RNAi_mean = rowMedians(as.matrix(N_hit_mean))
mean_RNAi = rowMeans(as.matrix(N_hit_mean))
quant_RNAi = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.90)})
plot(0:100, quant_RNAi,col = 'blue')
#plot(0:100, quant_RNAi,col = 'red')
abline(v = 30)
#abline(h = 0.95)

# apply filters and write out the new tables
lowExpCutoff = 30
for (plateID in c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')){
  # load the target library set
  mergeTbl = data.frame()
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  batches = levels(batchLabel)
  for (i in 1:length(batches)){
    subsetInd = batchLabel == batches[i]
    RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
    if (any(levels(RNAi_subset) == 'x.vector')){
      myseq1 = 2:(length(levels(RNAi_subset)))
    }else{
      myseq1 = 1:(length(levels(RNAi_subset)))
    }
    for (j in myseq1){
      targetGene = levels(RNAi_subset)[j]
      mySample = paste('RNAi_',targetGene,'_vs_others',sep = '')
      outputTbl = read.csv(paste('output/raw_DE_output/DE_table_one2all_lfcShrink_RNAi_raw_',mySample,'_',batches[i],'_imputMethod_random_strat.csv',sep = ''),row.names = 1)
      pVector = outputTbl$pvalue
      maxExp = apply(outputTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)

      outputTbl$padj = NA
      padj = p.adjust(pVector[maxExp >= lowExpCutoff],method = 'BH')
      outputTbl$padj[maxExp >= lowExpCutoff] = padj
      outputTbl$padj[maxExp < lowExpCutoff] = NA

      write.csv(outputTbl,paste('output/raw_DE_output/DE_table_one2all_lfcShrink_RNAi_alpha005_',mySample,'_',batches[i],'_imputMethod_random_strat.csv',sep = ''))
      resOrdered <- outputTbl[order(outputTbl$padj),]
      outputTbl = as.data.frame(subset(outputTbl, padj < 0.05))
      if (nrow(outputTbl) > 0){
        outputTbl = cbind(data.frame(WBID = rownames(outputTbl)),outputTbl)
        rownames(outputTbl) = 1:nrow(outputTbl)
        outputTbl$RNAi = rep(targetGene,nrow(outputTbl))
      }else{
        outputTbl = data.frame(WBID = 'NoHit')
        outputTbl$baseMean = NA
        outputTbl$log2FoldChange = NA
        outputTbl$log2FoldChange_raw = NA
        outputTbl$lfcSE = NA
        outputTbl$stat = NA
        outputTbl$pvalue = NA
        outputTbl$padj = NA
        outputTbl$medianCount_RNAi = NA
        outputTbl$medianCount_ctr = NA
        outputTbl$RNAi = targetGene
        outputTbl$libWiseBaseMean = NA
        rownames(outputTbl) = 1
      }
      outputTbl$batchID = batches[i]
      mergeTbl = rbind(mergeTbl, outputTbl)
      print(paste(mySample,batches[i]))
    }
  }
  mergeTbl = mergeTbl[order(mergeTbl["RNAi"], mergeTbl["log2FoldChange"]),]
  write.csv(mergeTbl,file = paste('output/raw_DE_output_tables/DE_one2all_master_table_FDR005_',plateID,'_imputMethod_random_strat.csv',sep = ''))
}

######################### part 0: without cleaning results: cutoff DE by FC #####################
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
for (plateID in plateIDs){
  cutoffFC(plateID,'','')
}

plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
for (plateID in plateIDs){
  cutoffFC(plateID,'one2all_','random_strat')
}
#################### part I: cleaning up confounding effects ######################
# NOTE: this is an entire modeling pipeline of the whole dataset (12 plates)
# so we need to think about what to include what to exclude whenever running it with new data inputs
# now, we decide to do a full evaluation of cleaning strength (controled by p-value cutoff for vector outlier calling)
# this evaluation is done after the empirical modeling
# we do this because we found that the false positive calls in vectorlikes are always well controled at conservative FDR/FC cutoff (i.e., FC 1.5 FDR 0.1)
# However, the reproducibility of DE calls is sensitive to the cleaning strength. this is due to a trade-off between false DE calls introduced by
# extreme vector outliers and false DE calls introduced by insufficient empirical modeling of variance (that was masked by too much cleaning)

FCtypes = c("log2FoldChange","log2FoldChange_raw")
cutoff_all = c(1, 0.5, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001,0.00075, 0.0005, 1e-04, 1e-5, 1e-10, 1e-30, 0) # better revise this in the final run
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11', 'met13')
for (FCtype in FCtypes){
  for (cutoffs in cutoff_all){
    ######################### (1) save the result by different cleaning strength #####################
    # overide the original DE result with the one-all comparison result - replace all z-score outliers in one-all comparison
    plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
    # cutoffs = 0.005
    writeCore = T
    freqCutoff = length(plateIDs)*6*0.25
    gene2rmList2 = findgenes_pvalue_merge2(plateIDs,cutoffs,'pvalue',writeCore, freqCutoff,FCtype) # cutoff could be calibrated by the vector like samples
    a = gene2rmList2[[1]]
    N_bad = c()
    for (i in 1:11){
      for (j in 1:length(a[[i]])){
        b = length(a[[i]][[j]])
        names(b) = names(a[[i]])[j]
        N_bad = c(N_bad,b )

      }
    }
    N_bad
    hist(N_bad)
    writeDetail = T
    mergeDE_v2(plateIDs,cutoffs,gene2rmList2,writeDetail,'random_strat',FCtype)
  }
}
