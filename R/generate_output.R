


###
cutoff_all = c(1, 0.5, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001,0.00075, 0.0005, 1e-04, 1e-5, 1e-10, 1e-30, 0) # better revise this in the final run

strvector = c('1','0.5','0.2','0.1','0.075','0.05', '0.025','0.01', '0.0075','0.005','0.0025', '0.001','0.00075','5e-04','1e-04', '1e-05', '1e-10','1e-30', '0')
optIndFilter = str %in% c('1','0.005','0') # enabling this will take ~1day to finish; move this back to the loop if we want to run it!
FCtypes = c("log2FoldChange","log2FoldChange_raw")

for (FCtype in FCtypes){
  #str = '0.005'
  for (str in strvector){
    library(DESeq2)
    ############# make the pvalue and fc matrix ######
    plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
    allConditions = c()
    allgenes = c()
    for (plateID in plateIDs){
      # load the FC matrix from DE
      load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
      allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
      allgenes = union(allgenes, rownames(input))
    }
    RNAiName = c()
    batchID = c()
    for (i in 1:length(allConditions)){
      RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
      batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])

    }
    # remove vectors
    batchID = batchID[RNAiName != 'x.vector']
    RNAiName = RNAiName[RNAiName != 'x.vector']

    # aggregate the raw tables for cutoff titrating
    pvalueTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(pvalueTbl) = allgenes
    colnames(pvalueTbl) = paste(RNAiName,'_',batchID,sep = '')
    pvalueTbl = as.data.frame(pvalueTbl)

    maxExpTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(maxExpTbl) = allgenes
    colnames(maxExpTbl) = paste(RNAiName,'_',batchID,sep = '')
    maxExpTbl = as.data.frame(maxExpTbl)

    statTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(statTbl) = allgenes
    colnames(statTbl) = paste(RNAiName,'_',batchID,sep = '')
    statTbl = as.data.frame(statTbl)

    sourceTbl = matrix(NA, nrow = length(allgenes), ncol = length(RNAiName))
    rownames(sourceTbl) = allgenes
    colnames(sourceTbl) = paste(RNAiName,'_',batchID,sep = '')
    sourceTbl = as.data.frame(sourceTbl)
    for (i in 1:length(RNAiName)){
      tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',RNAiName[i],'_',batchID[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''),row.names = 1)
      pvalueTbl[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')] = tbl$pvalue
      statTbl[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')] = tbl$stat
      maxExpTbl[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')] = apply(tbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
      sourceTbl[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')] = tbl$DE_source
      if (i %% 100 == 0){
        print(i)
      }
    }
    # analyze the entanglement of one2all DE input
    sourceCount = rowSums(sourceTbl == 'vs. vector',na.rm = T) / rowSums(!is.na(sourceTbl))
    hist(sourceCount)

    maxExpTbl_bak = maxExpTbl
    pvalueTbl_bak = pvalueTbl
    statTbl_bak = statTbl

    # cutoff for minimal measurement
    maxExpTbl = maxExpTbl_bak
    pvalueTbl = pvalueTbl_bak
    statTbl = statTbl_bak

    # there are some data cleaning mask (ie plate2 bad genes) in pvalue field. we transfer it to all the tables
    maxExpTbl[is.na(pvalueTbl)] = NA
    statTbl[is.na(pvalueTbl)] = NA
    # we remove genes without DE analysis output
    ind = rowAlls(is.na(pvalueTbl))
    statTbl = statTbl[!ind,]
    pvalueTbl = pvalueTbl[!ind,]
    maxExpTbl = maxExpTbl[!ind,]

    # we require 100 observation to calculate the empirical null
    keep = (rowSums(!is.na(pvalueTbl)) >= 100)
    ther_null_genes = rownames(pvalueTbl)[!keep]

    ################# empirical null modeling based on Efron's implementation ############
    library(locfdr)
    emp_pmat = matrix(NA, nrow = nrow(statTbl), ncol = ncol(statTbl))
    #emp_padjmat = matrix(NA, nrow = nrow(statTbl), ncol = ncol(statTbl))

    compromisedInd = c()
    ntrim = c()
    emp_mean = c()
    emp_sd = c()
    emp_ind = c()
    for (i in 1:nrow(statTbl)){
      if (keep[i]){ # we consider empircial null
        data = as.numeric(statTbl[i,])
        nonNAid = which(!is.na(data))
        data = data[nonNAid]

        # the fitting of mixed density function will either fail or skewed if there are extreme and discrete outliers
        # always exclude extreme discrete outliers to avoid bugs
        qts = quantile(data, c(0.01,0.99))
        upper_loose = qts[2] + 3*mad(data)
        lower_loose = qts[1] - 3*mad(data)
        # trimming rules
        trimInd = data > lower_loose & data < upper_loose
        data_trim = data[trimInd]
        ntrim = c(ntrim, sum(!trimInd))

        # show the fitting for super heavy tails
        if(sum(data > lower_loose & data < upper_loose) > 0.995 *length(data)){
          showplot = 0
        }else{
          showplot = 1
          compromisedInd = c(compromisedInd, i)
        }

        brk = length(data_trim) %/% 8
        success = FALSE
        step = 0
        while(!success & step < 100){
          fit = try({lfdr_model_trim = locfdr(data_trim, bre = brk,plot = showplot, type = 0)},silent = T)
          if (class(fit) == 'try-error'){
            brk = brk + 1
            step = step + 1
          }else{
            f0_mean = lfdr_model_trim$fp0['mlest','delta']
            f0_sd = lfdr_model_trim$fp0['mlest','sigma']
            emp_mean = c(emp_mean, f0_mean)
            emp_sd = c(emp_sd, f0_sd)
            emp_ind = c(emp_ind, i)
            # empirical pvalue
            lowerT = pnorm(data, mean = f0_mean, sd = f0_sd, lower.tail = TRUE, log.p = FALSE)
            upperT = pnorm(data, mean = f0_mean, sd = f0_sd, lower.tail = FALSE, log.p = FALSE)
            p_data = rowMins(cbind(lowerT, upperT)) * 2

            # p_data = rep(NA, length(data))
            # p_data[!trimInd] = 0
            # p_data[trimInd] = ePval
            # padj_data = p.adjust(p_data,method = 'BH')
            success = TRUE
          }
        }
        if (step == 100){
          print(paste(i,'maximum break tuning step reaches'))
          # the fitting is failed. we do some compromise
          # we use median and mad as best guess of empirical null
          f0_mean = median(data_trim)
          f0_sd = mad(data_trim)
          emp_mean = c(emp_mean, f0_mean)
          emp_sd = c(emp_sd, f0_sd)
          emp_ind = c(emp_ind, i)
          # empirical pvalue
          lowerT = pnorm(data, mean = f0_mean, sd = f0_sd, lower.tail = TRUE, log.p = FALSE)
          upperT = pnorm(data, mean = f0_mean, sd = f0_sd, lower.tail = FALSE, log.p = FALSE)
          p_data = rowMins(cbind(lowerT, upperT)) * 2
          hist(data_trim)
          curve(length(data_trim) * dnorm(x, f0_mean, f0_sd), add=TRUE, col="red", lwd=2)
          # p_data = rep(NA, length(data))
          # p_data[!trimInd] = 0
          # p_data[trimInd] = ePval
          # padj_data = p.adjust(p_data,method = 'BH')
          success = TRUE


        }

        emp_pmat[i,nonNAid] = p_data
      }else{ # we consider theretical null
        f0_mean = 0
        f0_sd = 1
        data_trim = as.numeric(statTbl[i,])
        # empirical pvalue
        lowerT = pnorm(data_trim, mean = f0_mean, sd = f0_sd, lower.tail = TRUE, log.p = FALSE)
        upperT = pnorm(data_trim, mean = f0_mean, sd = f0_sd, lower.tail = FALSE, log.p = FALSE)
        emp_pmat[i,] = rowMins(cbind(lowerT, upperT)) * 2
      }
      #emp_padjmat[i,nonNAid] = padj_data
      if (i %% 100 == 0){
        print(i)
      }
    }
    hist(ntrim)
    length(compromisedInd)
    hist(emp_mean)
    hist(emp_sd)
    fittingInfo = data.frame(gene = rownames(statTbl)[emp_ind], mean = emp_mean, sd = emp_sd)
    write.csv(fittingInfo,paste('output/empirical_fitting_result_p',str,'_',FCtype,'.csv',sep = ''))
    # save
    rownames(emp_pmat) = rownames(statTbl)
    colnames(emp_pmat) = colnames(statTbl)
    write.csv(emp_pmat, paste('output/empirical_pvalue_matrix_p',str,'_',FCtype,'.csv',sep = ''))

    if (optIndFilter){
      # optimize the independent filtering cutoff for FDR correction
      cutoffs = 0:100
      N_hit_rowfdr = matrix(nrow = length(cutoffs),ncol = ncol(emp_pmat))
      N_hit_colfdr = matrix(nrow = length(cutoffs),ncol = ncol(emp_pmat))
      N_hit_dualfdr1 = matrix(nrow = length(cutoffs),ncol = ncol(emp_pmat))
      # N_hit_dualfdr2 = matrix(nrow = length(cutoffs),ncol = ncol(emp_pmat))
      # N_hit_matfdr = matrix(nrow = length(cutoffs),ncol = ncol(emp_pmat))

      for (i in 1:length(cutoffs)){
        lowExpCutoff = cutoffs[i]
        row_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
        maxExpTbl_pass = maxExpTbl > lowExpCutoff
        maxExpTbl_pass[is.na(maxExpTbl_pass)] = F
        for (j in 1:nrow(row_fdr_mat)){
          row_fdr_mat[j,maxExpTbl_pass[j,]] = p.adjust(emp_pmat[j, maxExpTbl_pass[j,]],method = 'BH')
        }
        col_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
        for (j in 1:ncol(col_fdr_mat)){
          col_fdr_mat[maxExpTbl_pass[,j],j] = p.adjust(emp_pmat[maxExpTbl_pass[,j], j],method = 'BH')
        }
        # mat_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
        # mat_fdr_mat[maxExpTbl_pass] = p.adjust(emp_pmat[maxExpTbl_pass],method = 'BH')

        N_hit_rowfdr[i,] = colSums(row_fdr_mat < 0.05,na.rm = T)
        N_hit_colfdr[i,] = colSums(col_fdr_mat < 0.05,na.rm = T)
        N_hit_dualfdr1[i,] = colSums(col_fdr_mat < 0.05 & row_fdr_mat < 0.05,na.rm = T) # maximum < 0.05
        # N_hit_dualfdr2[i,] = colSums((1-(1-col_fdr_mat)*(1-row_fdr_mat)) < 0.05,na.rm = T) # maximum < 0.05
        # N_hit_matfdr[i,] = colSums(mat_fdr_mat < 0.05,na.rm = T)
        #mergeTbl = rbind(mergeTbl, outputTbl)
        print(i)
      }
      # theoritically, the most meaningful combined FDR is the N_hit_dualfdr1.
      # the meaning of it is: for the DE called within each RNAi, the FDR is controlled at 5%; at the same time, for each gene
      # analyzed across ~1000 conditions, the DE called also is controlled under 5% FDR
      # we dont do a whole-matrix level FDR control because that will only ensure the overall FDR is within 5%, but clearly the
      # local FDR for conditions with less DE will be much higher, while with more DE being much lower.
      # in fact, what we want is to control FDR both for each gene and each condition, thereby, we use N_hit_dualfdr1.
      N_hit_max0 = N_hit_dualfdr1
      # N_hit in vectorlikes
      vectorlikeIDs = read.csv('output/vectorlikeIDs.csv')
      # also add in the real vector
      vectorlikeIDs = vectorlikeIDs$x
      realvectorIDs = c('x.REALVECTOR_vector_met11_lib1','x.REALVECTOR_vector_met11_lib2',
                        'x.REALVECTOR_vector_met11_lib3','x.REALVECTOR_vector_met11_lib5','x.REALVECTOR_vector_met11_lib6')

      N_hit_max = N_hit_max0[,colnames(emp_pmat) %in% vectorlikeIDs]
      mean_RNAi = apply((as.matrix(N_hit_max)), 1, function(x){mean(x,trim = .05)})
      quant_RNAi1 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.80)})
      quant_RNAi2 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.90)})
      quant_RNAi3 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.95)})
      #median_RNAi_max = median_RNAi_max/max(median_RNAi_max)
      #median_RNAi_mean = median_RNAi_mean/max(median_RNAi_mean)
      #mean_RNAi = mean_RNAi/max(mean_RNAi)
      #quant_RNAi = quant_RNAi/max(quant_RNAi)
      pdf(paste('figures/independent_filtering_optimization_cleaning_p',str,'_',FCtype,'.pdf',sep = ''),height = 10,width = 10)

      plot(0:100, quant_RNAi1/max(quant_RNAi1),col = 'red',ylim = c(0,1),lty = 2, type = 'l')
      lines(0:100, quant_RNAi2/max(quant_RNAi2),col = 'yellow',lty = 2)
      lines(0:100, quant_RNAi3/max(quant_RNAi3),col = 'blue',lty = 2)
      lines(0:100, mean_RNAi/max(mean_RNAi),col = 'black',lty = 2)

      N_hit_max = N_hit_max0
      mean_RNAi = apply((as.matrix(N_hit_max)), 1, function(x){mean(x,trim = .1)})
      quant_RNAi1 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.60)})
      quant_RNAi2 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.70)})
      quant_RNAi3 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.80)})
      #median_RNAi_max = median_RNAi_max/max(median_RNAi_max)
      #median_RNAi_mean = median_RNAi_mean/max(median_RNAi_mean)
      #mean_RNAi = mean_RNAi/max(mean_RNAi)
      #quant_RNAi = quant_RNAi/max(quant_RNAi)
      lines(0:100, quant_RNAi1/max(quant_RNAi1),col = 'red',lty = 1)
      lines(0:100, quant_RNAi2/max(quant_RNAi2),col = 'yellow',lty = 1)
      lines(0:100, quant_RNAi3/max(quant_RNAi3),col = 'blue',lty = 1)
      lines(0:100, mean_RNAi/max(mean_RNAi),col = 'black',lty = 1)
      abline(v = 30)
      abline(v = 25)

      N_hit_max = N_hit_max0[,colnames(emp_pmat) %in% realvectorIDs]
      mean_RNAi = apply((as.matrix(N_hit_max)), 1, function(x){mean(x,trim = .05)})
      quant_RNAi1 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.80)})
      quant_RNAi2 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.90)})
      quant_RNAi3 = apply((as.matrix(N_hit_max)), 1, function(x){quantile(x,0.95)})
      lines(0:100, quant_RNAi1/max(quant_RNAi1),col = 'red',lty = 3)
      lines(0:100, quant_RNAi2/max(quant_RNAi2),col = 'yellow',lty = 3)
      lines(0:100, quant_RNAi3/max(quant_RNAi3),col = 'blue',lty = 3)
      lines(0:100, mean_RNAi/max(mean_RNAi),col = 'black',lty = 3)
      dev.off()
    }
    # save p.adj
    lowExpCutoff = 30
    row_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
    maxExpTbl_pass = maxExpTbl > lowExpCutoff
    maxExpTbl_pass[is.na(maxExpTbl_pass)] = F
    for (j in 1:nrow(row_fdr_mat)){
      row_fdr_mat[j,maxExpTbl_pass[j,]] = p.adjust(emp_pmat[j, maxExpTbl_pass[j,]],method = 'BH')
    }
    col_fdr_mat = matrix(NA, nrow = nrow(emp_pmat), ncol = ncol(emp_pmat))
    for (j in 1:ncol(col_fdr_mat)){
      col_fdr_mat[maxExpTbl_pass[,j],j] = p.adjust(emp_pmat[maxExpTbl_pass[,j], j],method = 'BH')
    }
    dual_fdr_mat = pmax(col_fdr_mat, row_fdr_mat)
    rownames(row_fdr_mat) = rownames(statTbl)
    colnames(row_fdr_mat) = colnames(statTbl)
    rownames(col_fdr_mat) = rownames(statTbl)
    colnames(col_fdr_mat) = colnames(statTbl)
    rownames(dual_fdr_mat) = rownames(statTbl)
    colnames(dual_fdr_mat) = colnames(statTbl)
    write.csv(dual_fdr_mat, paste('output/empirical_dualFDR_matrix_p',str,'_',FCtype,'.csv',sep = ''))

    # update the DE tbls
    for (i in 1:length(RNAiName)){
      tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',RNAiName[i],'_',batchID[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''),row.names = 1)
      # commGenes = intersect(rownames(tbl),rownames(emp_pmat))
      if(length(setdiff(rownames(tbl),rownames(emp_pmat)))>0){
        stop('what happened?')
      }
      tbl$empirical_pvalue = emp_pmat[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')]
      tbl$dualFDR = dual_fdr_mat[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')]
      tbl$rowFDR = row_fdr_mat[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')]
      tbl$colFDR = col_fdr_mat[rownames(tbl),paste(RNAiName[i],'_',batchID[i],sep = '')]

      # tbl$emp_padj_row = 0 # no data assume no filter
      # commGenes = intersect(rownames(tbl),rownames(emp_padjmat))
      # tbl[commGenes,'emp_padj_row'] = emp_padjmat[commGenes,str_replace_all(paste(RNAiName[i],'_',batchID[i],sep = ''),'-','.')]
      # # calculate final FDR (2d)
      # tbl$FDR2d_thr_emp = 1 - (1-tbl$padj) * (1-tbl$emp_padj_row)
      # tbl$emp_padj_col = p.adjust(tbl$emp_pvalue_row, method = 'BH')
      # tbl$FDR2d_emp_emp = 1 - (1-tbl$emp_padj_col) * (1-tbl$emp_padj_row)
      write.csv(tbl, paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',RNAiName[i],'_',batchID[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
      print(i)
    }

    #################### re-cutoff optimal threshold: using the vectorlike samples to find best FC and FDR cutoff for DE calling  ######################
    library(stringr)
    plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
    allConditions = c()
    for (plateID in plateIDs){
      # load the FC matrix from DE
      load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
      allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
    }
    RNAiName = c()
    batchID = c()
    for (i in 1:length(allConditions)){
      RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
      batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])

    }

    # consider all vectorlike??
    vectorlikes = RNAiName[str_detect(RNAiName, '^x.VECTORLIKE_')
                           | str_detect(RNAiName, '^x.REALVECTOR_')
                           # | str_detect(allConditions, '^x.RCBVECTOR_')
                           # | str_detect(allConditions, '^x.NOSIGNAL_')
                           ]
    batch_vectorlikes = batchID[str_detect(RNAiName, '^x.VECTORLIKE_')
                                | str_detect(RNAiName, '^x.REALVECTOR_')]
    # we exclude the vectorlike that produce single stranded RNA
    # vectorlikes = setdiff(vectorlikes,c('x.VECTORLIKE_daao_1','x.VECTORLIKE_F53F4.10','x.VECTORLIKE_T26E3.7','x.VECTORLIKE_R151.2','x.VECTORLIKE_F13B12.4','x.VECTORLIKE_acs_3'))
    # # we exclude one vectorlike that has suspicous mapping ()
    # vectorlikes = setdiff(vectorlikes,c('x.VECTORLIKE_jhdm_1')) # it roughly maps to hrde-1 and down reg it for 2fc
    # good short - expression quantified, and expression not changed; insertion is short
    goodshort = c('x.SHORT_fdps_1','x.SHORT_lars_1','x.SHORT_mecr_1','x.SHORT_mtm_3','x.SHORT_alh_3','x.SHORT_bli_4','x.SHORT_C15H9.5')
    batch_goodshort = batchID[match(goodshort,RNAiName)]
    vectorlikes = c(vectorlikes, goodshort)
    batch_vectorlikes = c(batch_vectorlikes,batch_goodshort)

    # aggregate the raw tables for cutoff titrating
    mergeTbl = data.frame()
    for (i in 1:length(vectorlikes)){
      tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',vectorlikes[i],'_',batch_vectorlikes[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
      tbl$RNAi = vectorlikes[i]
      tbl$batchID = batch_vectorlikes[i]
      mergeTbl = rbind(mergeTbl, tbl)
    }
    maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
    # N_DE_list = list()
    fdr_seq = c(0.00001, 0.0001,0.001,seq(0.01,0.2,0.01), seq(0.25,1,0.05))
    fc_seq = seq(1,3,0.1)

    N_DE_quant = array(rep(NA, length(fdr_seq) * length(fc_seq)),
                       c(length(fdr_seq), length(fc_seq)))
    for (FDR in fdr_seq){
      for (FC in fc_seq){
        passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < FDR),]
        passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
        N_DE = table(passTbl$ID)
        zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
        tmp = rep(0, length(zeros))
        names(tmp) = zeros
        N_DE = c(N_DE, tmp)
        N_DE_quant[fdr_seq %in% FDR, fc_seq %in% FC] = quantile(N_DE,0.9) #
      }
    }

    library(pheatmap)
    #dev.off()
    pdf(paste('figures/2d_cutoff_titration_57vectorlikes_p',str,'_',FCtype,'.pdf',sep = ''),height = 10,width = 10)
    pheatmap(N_DE_quant,labels_row = fdr_seq, labels_col = fc_seq,cluster_rows = F,cluster_cols = F,display_numbers = T,
             breaks = seq(0,160,1.6)) #[,,eFDR_seq %in% eFDR]

    # evaluate FC cutoff at a fixed FDR
    # we end up with fdr and eFDR 0.05
    N_DE_fdr01 = matrix(NA, nrow = length(fc_seq), ncol = length(vectorlikes))
    colnames(N_DE_fdr01) = paste(vectorlikes, batch_vectorlikes)
    for (FC in fc_seq){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < 0.1),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_fdr01[fc_seq %in% FC,paste(vectorlikes, batch_vectorlikes)] = N_DE[paste(vectorlikes, batch_vectorlikes)]
    }
    plot(fc_seq, N_DE_fdr01[,1],type = 'o',ylim = c(0, 60))
    for (i in 2:ncol(N_DE_fdr01)){
      points(fc_seq, N_DE_fdr01[,i],type = 'o')
    }
    abline(h=5, col = 'red')
    dev.off()

    # do this evaluation with real vectors
    vectorlikes = RNAiName[str_detect(RNAiName, '^x.REALVECTOR_')
                           # | str_detect(allConditions, '^x.RCBVECTOR_')
                           # | str_detect(allConditions, '^x.NOSIGNAL_')
                           ]
    batch_vectorlikes = batchID[str_detect(RNAiName, '^x.REALVECTOR_')]

    # aggregate the raw tables for cutoff titrating
    mergeTbl = data.frame()
    for (i in 1:length(vectorlikes)){
      tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',vectorlikes[i],'_',batch_vectorlikes[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
      tbl$RNAi = vectorlikes[i]
      tbl$batchID = batch_vectorlikes[i]
      mergeTbl = rbind(mergeTbl, tbl)
    }
    maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
    # N_DE_list = list()
    fdr_seq = c(0.00001, 0.0001,0.001,seq(0.01,0.2,0.01), seq(0.25,1,0.05))
    fc_seq = seq(1,3,0.1)

    N_DE_quant = array(rep(NA, length(fdr_seq) * length(fc_seq)),
                       c(length(fdr_seq), length(fc_seq)))
    for (FDR in fdr_seq){
      for (FC in fc_seq){
        passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < FDR),]
        passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
        N_DE = table(passTbl$ID)
        zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
        tmp = rep(0, length(zeros))
        names(tmp) = zeros
        N_DE = c(N_DE, tmp)
        N_DE_quant[fdr_seq %in% FDR, fc_seq %in% FC] = mean(N_DE) #
      }
    }

    library(pheatmap)
    #dev.off()
    pdf(paste('figures/2d_cutoff_titration_realVector_p',str,'_',FCtype,'.pdf',sep = ''),height = 10,width = 10)
    pheatmap(N_DE_quant,labels_row = fdr_seq, labels_col = fc_seq,cluster_rows = F,cluster_cols = F,display_numbers = T,
             breaks = seq(0,160,1.6)) #[,,eFDR_seq %in% eFDR]

    # evaluate FC cutoff at a fixed FDR
    # we end up with fdr and eFDR 0.05
    N_DE_fdr01 = matrix(NA, nrow = length(fc_seq), ncol = length(vectorlikes))
    colnames(N_DE_fdr01) = paste(vectorlikes, batch_vectorlikes)
    for (FC in fc_seq){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FC) & mergeTbl$dualFDR < 0.1),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_fdr01[fc_seq %in% FC,paste(vectorlikes, batch_vectorlikes)] = N_DE[paste(vectorlikes, batch_vectorlikes)]
    }
    plot(fc_seq, N_DE_fdr01[,1],type = 'o',ylim = c(0, 60))
    for (i in 2:ncol(N_DE_fdr01)){
      points(fc_seq, N_DE_fdr01[,i],type = 'o')
    }
    abline(h=5, col = 'red')
    dev.off()
  }
}


#==> based on real vector, minimal FC valid is 1.3; 1.5 is sufficiently good
# write clean, final table
str = 0.005
# FDRcutoff = 0.05
# FDRcutoff = 0.2
FDRcutoff = 0.2
#eFDRcutoff = 0.05
FCcutoff = 1.5
FCtype = 'log2FoldChange_raw'

# write out the cutoffed table
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
allConditions = c()
allgenes = c()
for (plateID in plateIDs){
  # load the FC matrix from DE
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
  allgenes = union(allgenes, rownames(input))
}
RNAiName = c()
batchID = c()
for (i in 1:length(allConditions)){
  RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
  batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])

}
# remove vectors
batchID = batchID[RNAiName != 'x.vector']
RNAiName = RNAiName[RNAiName != 'x.vector']

batchID = batchID['x.vector' != RNAiName]
RNAiName = RNAiName['x.vector' != RNAiName]
mergeTbl = data.frame()
for (i in 1:length(RNAiName)){
  tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',RNAiName[i],'_',batchID[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
  tbl$RNAi = RNAiName[i]
  tbl$batchID = batchID[i]
  tbl = tbl[!is.na(tbl$dualFDR),]
  maxExp = apply(tbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
  tbl = tbl[(abs(tbl[,FCtype]) > log2(FCcutoff) & tbl$dualFDR < FDRcutoff), ]
  if(nrow(tbl)>0){
    mergeTbl = rbind(mergeTbl, tbl)
  }
}
colnames(mergeTbl)[1] = 'WBID'
# add annotations
iCELGenes = read.csv('./../input_data/otherTbls/WormPaths_Tables/wormPathTable.csv')
metabolicGenes = read.csv('./../input_data/otherTbls/allMetGenes.txt')
mergeTbl$isICEL = mergeTbl$WBID %in% iCELGenes$WormBase.ID
mergeTbl$isMetabolic = mergeTbl$WBID %in% metabolicGenes$WBName
relCV = read.csv('output/median_deviation_of_CV_to_the_medianCV.csv',row.names = 1)
mergeTbl$relCV = relCV$medianRelCV[match(mergeTbl$WBID, rownames(relCV))]
library(stringr)
library(tidyr)
library(dbplyr)

mergeTbl$LEVEL1 = iCELGenes$LEVEL.1[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
mergeTbl$LEVEL1[is.na(mergeTbl$LEVEL1)] = 'Non-iCEL'
mergeTbl$LEVEL2 = iCELGenes$LEVEL.2[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
mergeTbl$LEVEL2[is.na(mergeTbl$LEVEL2)] = 'Non-iCEL'
mergeTbl$LEVEL3 = iCELGenes$LEVEL.3[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
mergeTbl$LEVEL3[is.na(mergeTbl$LEVEL3)] = 'Non-iCEL'
mergeTbl$LEVEL4 = iCELGenes$LEVEL.4[match(mergeTbl$WBID, iCELGenes$WormBase.ID)]
mergeTbl$LEVEL4[is.na(mergeTbl$LEVEL4)] = 'Non-iCEL'

# add back no-hit
noHitRNAi = setdiff(paste(RNAiName, batchID), unique(paste(mergeTbl$RNAi, mergeTbl$batchID)))
noHitRNAi_batch = c()
noHitRNAi_new = c()# if a noHit RNAi appears twice
if (length(noHitRNAi) > 0){
  for (i in 1:length(noHitRNAi)){
    noHitRNAi_batch = c(noHitRNAi_batch, strsplit(noHitRNAi[i],' ')[[1]][[2]])
    noHitRNAi_new = c(noHitRNAi_new,strsplit(noHitRNAi[i],' ')[[1]][[1]])
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
  noHit$DE_source= NA
  noHit$RNAi= noHitRNAi_new
  noHit$batchID = noHitRNAi_batch
  noHit$isICEL = NA
  noHit$isMetabolic = NA
  noHit$relCV = NA
  noHit$LEVEL1 = NA
  noHit$LEVEL2 = NA
  noHit$LEVEL3 = NA
  noHit$LEVEL4 = NA
  noHit$empirical_pvalue = NA
  noHit$dualFDR = NA
  noHit$colFDR = NA
  noHit$rowFDR = NA
  mergeTbl = rbind(mergeTbl, noHit)
}
WBID = read.table('./../input_data/otherTbls/WBIDtbl.txt',header = T,sep = '\t')
mergeTbl$Gene_name = WBID$Public.Name[match(mergeTbl$WBID,WBID$WormBase.Gene.ID)]
mergeTbl = mergeTbl[,c('RNAi','Gene_name',setdiff(colnames(mergeTbl),c('Gene_name','RNAi','lfcSE')))]
write.csv(mergeTbl,file = paste('output/DE_merged_clean_pcutoff_',str,'_master_table_FDR2d',FDRcutoff,'_FC',FCcutoff,'_FCtype_',FCtype,'_ALL.csv',sep = ''),row.names = F) #'_eFDR',eFDRcutoff




# we found that the cleaniness of vectorlike/realvectors are robust to level of cleaning (p-cutoff) and to the choice of raw and shrink fold change
# but the reproducibility (power vs. FP in real treatments) is sensitive to level of cleaning.
# we titrate the optimal cleaning threshold here
library(ggrepel)
library(ggplot2)
# note: another cleaning parameter, the core cutoff, was not titrated. This parameter will interact with the p-value cutoff in the
# FP removing effect. we have the option to lower the core cutoff to make the smaller p-cutoff more effective
# maybe we can try a few cases without a full 2-D titration. (ps. a quick check found that 20% cutoff is not as good (slightly) as the 25%, so we dont change and stop further optimization)

# load all RNAi that has been tested multiple time, and check them all
RNAisheet = read.csv('./../1_QC_dataCleaning/outputs/RNAi_summary_final_dataset.csv')
RNAiIDsheet = read.csv('./../1_QC_dataCleaning/outputs/RNAi_batch_lookup.csv')
RNAiIDsheet$fileName = RNAiIDsheet$RNAi_ID
# some L2 labeling is messy we correct the problem we found
RNAisheet$independentExp[RNAisheet$RNAi_ID == 'x.elo_3'] = 1
RNAisheet$independentExp[RNAisheet$RNAi_ID == 'x.elo_3_L2'] = 2
RNAiIDsheet$RNAi_ID[RNAiIDsheet$batch_ID == 'met1_lib2' & RNAiIDsheet$RNAi_ID =='x.elo_3'] = 'x.elo_3_L2'
RNAisheet$independentExp[RNAisheet$RNAi_ID == 'x.elo_5'] = 1
RNAisheet$independentExp[RNAisheet$RNAi_ID == 'x.elo_5_L2'] = 2
RNAiIDsheet$RNAi_ID[RNAiIDsheet$batch_ID == 'met1_lib3' & RNAiIDsheet$RNAi_ID =='x.elo_5'] = 'x.elo_5_L2'
# we skip F45H10.2 because they are done by different RNAi and different stages

candidates = RNAisheet[RNAisheet$independentExp > 1,]
exclude = c('x.vha_5','x.dld_1') # these three were done with different RNAi sequence or was identified to have problematic RNAi
exclude_sample = c('x.acdh_1_met1_lib1', # exclude this because it only has two replicates and showed significant power-diff the three-rep data (2-rep gives less confident DE calls and more unstable FC assignment)
                   'x.nduf_6_met11_lib6' # excldde this because it is seeded in a different stage
)
# be careful on x.nduf_6, it is a posterier identified RNAi condition
candidates = candidates[!(candidates$RNAi_ID %in% exclude), ]


# these parameters are subjective to robustness analysis
FDRcutoff = 0.2 # not sensitive; larger makes pattern clearer
FCcutoff = 1.5 # not sensitive
#lowExpCutoff = 30
notRecap_cutoff = log2(1.1) # we consider these as strictly not reproduced; if this gets too large, it will be subjective to a power issue


pvector = c(1, 0.5, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001,0.00075, 0.0005, 1e-04, 1e-5, 1e-10, 1e-30, 0)
strvector = c('1','0.5','0.2','0.1','0.075','0.05', '0.025','0.01', '0.0075','0.005','0.0025', '0.001','0.00075','5e-04','1e-04', '1e-05', '1e-10','1e-30', '0')

library(lsa)
pair_A = c()
pair_B = c()
for (i in 1:nrow(candidates)){
  targetGene = candidates$RNAi_ID[i]
  samples = RNAiIDsheet[RNAiIDsheet$RNAi_ID == targetGene,]
  samples$fullID = paste(samples$RNAi_ID, samples$batch_ID,sep = '_')
  samples = samples[!(samples$fullID %in% exclude_sample),]
  for (j in 1:(nrow(samples)-1)){
    for (k in (j+1):nrow(samples)){
      pair_A = c(pair_A, paste(samples$RNAi_ID[j],'_',samples$batch_ID[j],sep = ''))
      pair_B = c(pair_B, paste(samples$RNAi_ID[k],'_',samples$batch_ID[k],sep = ''))
    }
  }
}

FCtypes = c("log2FoldChange","log2FoldChange_raw")
for (FCtype in FCtypes){
  pearson_logP_mat = c()
  pearson_logFC_mat = c()
  pearson_logFDR_mat = c()
  cosine_logFDR_mat = c()
  PCC_cleanFC_mat = c()
  Jacc_DE_call_mat= c()
  ave_N_DE_mat = c()
  Union_N_DE_mat = c()
  N_not_reproduce_mat = c()
  estFDR_mat = c()

  for (i in 1:length(strvector)){
    str = strvector[i]

    pearson_empPval = c()
    pearson_FC = c()
    pearson_empFDR = c()
    cosine_empFDR = c()
    PCC_FC_clean = c()
    Jacc_empFDR= c()

    N_DE_empFDR = c() # average N_DE (power)
    Union_DE_empFDR = c()
    N_not_reproduce = c()
    estFDR = c()

    for (i in 1:nrow(candidates)){
      targetGene = candidates$RNAi_ID[i]
      samples = RNAiIDsheet[RNAiIDsheet$RNAi_ID == targetGene,]
      samples$fullID = paste(samples$RNAi_ID, samples$batch_ID,sep = '_')
      samples = samples[!(samples$fullID %in% exclude_sample),]
      for (j in 1:(nrow(samples)-1)){
        for (k in (j+1):nrow(samples)){
          tbl_rep1 = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',samples$fileName[j],'_',samples$batch_ID[j],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
          tbl_rep2 = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',samples$fileName[k],'_',samples$batch_ID[k],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
          # avoid NaN in log calculation
          tbl_rep1$empirical_pvalue[tbl_rep1$empirical_pvalue==0] = 1e-300
          tbl_rep2$empirical_pvalue[tbl_rep2$empirical_pvalue==0] = 1e-300
          tbl_rep1$dualFDR[tbl_rep1$dualFDR==0] = 1e-300
          tbl_rep2$dualFDR[tbl_rep2$dualFDR==0] = 1e-300

          tbl_rep2$oriP = tbl_rep1$empirical_pvalue[match(tbl_rep2$X,tbl_rep1$X)]
          tbl_rep2$oriFC = tbl_rep1[match(tbl_rep2$X,tbl_rep1$X),FCtype]
          tbl_rep2$oriFDR = tbl_rep1$dualFDR[match(tbl_rep2$X,tbl_rep1$X)]

          pearson_empPval[length(pearson_empPval)+1] = cor(-log10(tbl_rep2$oriP) * sign(tbl_rep2$oriFC),
                                                           -log10(tbl_rep2$empirical_pvalue) * sign(tbl_rep2[,FCtype]),use = 'pairwise.complete.obs')
          pearson_empFDR[length(pearson_empFDR)+1] = cor(-log10(tbl_rep2$oriFDR) * sign(tbl_rep2$oriFC),
                                                         -log10(tbl_rep2$dualFDR) * sign(tbl_rep2[,FCtype]),use = 'pairwise.complete.obs')
          pearson_FC[length(pearson_FC)+1] = cor((tbl_rep2$oriFC),
                                                 (tbl_rep2[,FCtype]),use = 'pairwise.complete.obs')
          x = -log10(tbl_rep2$oriFDR) * sign(tbl_rep2$oriFC)
          y = -log10(tbl_rep2$dualFDR) * sign(tbl_rep2[,FCtype])
          paircomplete = !is.na(x) & !is.na(y)
          cosine_empFDR[length(cosine_empFDR)+1] = cosine(x[paircomplete], y[paircomplete])

          # DE calls
          DEcalls_rep1 = tbl_rep1$X[which(tbl_rep1$dualFDR < FDRcutoff & abs(tbl_rep1[,FCtype]) > log2(FCcutoff))]
          DEcalls_rep2 = tbl_rep2$X[which(tbl_rep2$dualFDR < FDRcutoff & abs(tbl_rep2[,FCtype]) > log2(FCcutoff))]
          N_DE_empFDR = c(N_DE_empFDR, mean(c(length(DEcalls_rep1), length(DEcalls_rep2))))

          # correlation of called DEG
          geneSet = intersect(union(DEcalls_rep1, DEcalls_rep2),
                              intersect(tbl_rep2$X, tbl_rep1$X))
          Union_DE_empFDR = c(Union_DE_empFDR, length(geneSet))

          Jacc_empFDR = c(Jacc_empFDR,
                          length(intersect(DEcalls_rep1, DEcalls_rep2))/ length(union(DEcalls_rep1,DEcalls_rep2)))
          # cosine distance of DE calls
          if (length(geneSet)>0){
            PCC_FC_clean[length(PCC_FC_clean)+1] = cor((tbl_rep2$oriFC[tbl_rep2$X %in% geneSet]),
                                                       (tbl_rep2[tbl_rep2$X %in% geneSet,FCtype]))
            # select the strictly not reproduced genes
            n = sum(tbl_rep2$X %in% geneSet & (
              abs(tbl_rep2[,FCtype])< notRecap_cutoff |
                abs(tbl_rep2$oriFC)< notRecap_cutoff |
                sign(tbl_rep2$oriFC) != sign(tbl_rep2[,FCtype])
            ),na.rm = T)
            N_not_reproduce = c(N_not_reproduce, n)
            estFDR = c(estFDR, n/length(geneSet))
          }else{
            PCC_FC_clean[length(PCC_FC_clean)+1] = 1
            N_not_reproduce = c(N_not_reproduce, 0)
            estFDR = c(estFDR, 0)
          }

        }
      }
    }
    pearson_logP_mat = rbind(pearson_logP_mat, pearson_empPval)
    pearson_logFC_mat = rbind(pearson_logFC_mat, pearson_FC)
    pearson_logFDR_mat = rbind(pearson_logFDR_mat, pearson_empFDR)
    cosine_logFDR_mat = rbind(cosine_logFDR_mat, cosine_empFDR)
    PCC_cleanFC_mat = rbind(PCC_cleanFC_mat, PCC_FC_clean)
    Jacc_DE_call_mat= rbind(Jacc_DE_call_mat, Jacc_empFDR)
    ave_N_DE_mat = rbind(ave_N_DE_mat, N_DE_empFDR)
    Union_N_DE_mat = rbind(Union_N_DE_mat, Union_DE_empFDR)
    N_not_reproduce_mat = rbind(N_not_reproduce_mat, N_not_reproduce)
    estFDR_mat = rbind(estFDR_mat, estFDR)

  }

  library(tune)
  library(ggplot2)
  library(ggrepel)
  pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_',FCtype,'.pdf',sep = ''))
  for (i in 1: ncol(estFDR_mat)){
    df = data.frame(x = Union_N_DE_mat[,i] - N_not_reproduce_mat[,i],y =  estFDR_mat[,i], strvector)
    print(ggplot(df, aes(log2(1+x), y,label = strvector)) +
            geom_point(color = "red")+geom_path(linetype = "dashed",)+
            xlab('# reproducible DEG calls (estimated power)')+ # this power was defiend by number of union DEG calls minus strictly unreproducible calls
            ylab('estimated total false positive call rate') + # this FDR was defined by [strictly unreproducible calls] / [number of union DEG calls]
            xlim(0,11)+ylim(0,1)+
            ggtitle(paste(pair_A[i],pair_B[i]))+
            geom_text_repel(max.overlaps =Inf, size = 5))
    df$x = df$x/max(df$x)
    df$x[is.na(df$x)] = 1
    print(ggplot(df, aes(x, y,label = strvector)) +
            geom_point(color = "red")+geom_path(linetype = "dashed",)+
            xlab('estimated recall of highest power (%)')+
            ylab('estimated total false positive call rate') +
            xlim(0,1)+ylim(0,1)+
            ggtitle(paste(pair_A[i],pair_B[i]))+
            geom_text_repel(max.overlaps =Inf, size = 5))
    df = data.frame(x = Union_N_DE_mat[,i] - N_not_reproduce_mat[,i],y =  N_not_reproduce_mat[,i], strvector)
    print(ggplot(df, aes(x, y,label = strvector)) +
            geom_point(color = "red")+geom_path(linetype = "dashed",)+
            xlab('# reproducible DEG calls (estimated TP)')+  # this power was defiend by number of union DEG calls minus strictly unreproducible calls
            ylab('# nonreproducible DEG calls (estimated FP)') + # this FP was defined by #[strictly unreproducible calls]
            ggtitle(paste(pair_A[i],pair_B[i]))+
            geom_text_repel(max.overlaps =Inf, size = 5)+ coord_obs_pred())
  }
  dev.off()
  library(stringr)

  pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_FP_',FCtype,'.pdf',sep = ''),height = 14,width = 14)
  df = data.frame(x = pvector, y =  N_not_reproduce_mat[,1], ann = str_remove(pair_A[1],'_met.+_lib.$'))
  for (i in 2: ncol(estFDR_mat)){
    tmp = data.frame(x = pvector, y =  N_not_reproduce_mat[,i], ann = str_remove(pair_A[i],'_met.+_lib.$'))
    df = rbind(df, tmp)
  }
  print(ggplot(df, aes(x, y,group = ann)) +
          geom_point(color = "red")+geom_path(linetype = "solid",)+
          xlab('p cutoff')+
          ylab('estimated # of false positive calls') + facet_grid(ann ~ ., scales="free_y") +
          scale_x_log10()
  )
  library(matrixStats)
  heatTbl = scale(N_not_reproduce_mat, center = F, scale = colMaxs(N_not_reproduce_mat,na.rm = T ))
  heatTbl[is.na(heatTbl)] = 0
  rowlab = str_remove(pair_A,'_met.+_lib.$')
  rowlab = paste(colMaxs(N_not_reproduce_mat), rowlab)
  pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab,cellwidth = 25, cellheight = 10)
  dev.off()


  pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_TP_',FCtype,'.pdf',sep = ''),height = 14,width = 14)
  df = data.frame(x = pvector, y =  Union_N_DE_mat[,1] - N_not_reproduce_mat[,1], ann = str_remove(pair_A[1],'_met.+_lib.$'))
  for (i in 2: ncol(estFDR_mat)){
    tmp = data.frame(x = pvector, y =  Union_N_DE_mat[,i] - N_not_reproduce_mat[,i], ann = str_remove(pair_A[i],'_met.+_lib.$'))
    df = rbind(df, tmp)
  }
  print(ggplot(df, aes(x, y,group = ann)) +
          geom_point(color = "red")+geom_path(linetype = "solid",)+
          xlab('p cutoff')+
          ylab('estimated # of true positive calls') + facet_grid(ann ~ ., scales="free_y") +
          scale_x_log10()
  )
  library(matrixStats)
  heatTbl = scale(Union_N_DE_mat - N_not_reproduce_mat, center = F, scale = colMaxs(Union_N_DE_mat - N_not_reproduce_mat,na.rm = T))
  heatTbl[is.na(heatTbl)] = 0
  rowlab = str_remove(pair_A,'_met.+_lib.$')
  rowlab = paste(colMaxs(Union_N_DE_mat - N_not_reproduce_mat), rowlab)
  pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab,cellwidth = 25, cellheight = 10)
  dev.off()

  # these three are not good metrc as they are not sensitive to signal or confounded by noises
  # pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_logP_PCC.pdf',sep = '')) # not a good metric
  # library(matrixStats)
  # heatTbl = pearson_logP_mat
  # rowlab = str_remove(pair_A,'_met.+_lib.$')
  # pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab)
  # dev.off()
  # pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_logFC_PCC.pdf',sep = '')) # not a good metric
  # library(matrixStats)
  # heatTbl = pearson_logFC_mat
  # rowlab = str_remove(pair_A,'_met.+_lib.$')
  # pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab)
  # dev.off()
  # pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_logFDR_PCC.pdf',sep = '')) # not a good metric
  # library(matrixStats)
  # heatTbl = pearson_logFDR_mat
  # rowlab = str_remove(pair_A,'_met.+_lib.$')
  # pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab)
  # dev.off()
  # pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_logFDR_Cosine.pdf',sep = '')) # not a good metric
  # library(matrixStats)
  # heatTbl = cosine_logFDR_mat
  # heatTbl = scale(heatTbl, center = F, scale = colMaxs(abs(heatTbl)))
  # rowlab = str_remove(pair_A,'_met.+_lib.$')
  # pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab)
  # dev.off()

  pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_clean_logFC_PCC_',FCtype,'.pdf',sep = '')) # not a good metric
  library(matrixStats)
  heatTbl = PCC_cleanFC_mat
  heatTbl = scale(heatTbl, center = F, scale = colMaxs(abs(heatTbl),na.rm = T))
  heatTbl[is.na(heatTbl)] = 0
  rowlab = str_remove(pair_A,'_met.+_lib.$')
  pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab)
  dev.off()

  pdf(paste('figures/0_DE_QA/benchmark_independent_repeats_titrate_cleaning_cutoff_DE_jaccard_',FCtype,'.pdf',sep = '')) # not a good metric
  library(matrixStats)
  heatTbl = Jacc_DE_call_mat
  heatTbl[is.na(heatTbl)] = 1
  heatTbl = scale(heatTbl, center = F, scale = colMaxs(abs(heatTbl),na.rm = T))
  heatTbl[is.na(heatTbl)] = 0
  rowlab = str_remove(pair_A,'_met.+_lib.$')
  rowlab = paste(colMaxs(Union_N_DE_mat),round(colMaxs(Jacc_DE_call_mat,na.rm = T),2), rowlab)
  pheatmap::pheatmap(t(heatTbl),cluster_cols = F, labels_col = strvector, labels_row = rowlab)
  dev.off()

}
dev.off()


######## visualize the titration on vectorlike-FP with fixed cutoff, too ######
cutoff_all = c(1, 0.5, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001,0.00075, 0.0005, 1e-04, 1e-5, 1e-10, 1e-30, 0) # better revise this in the final run
strvector = c('1','0.5','0.2','0.1','0.075','0.05', '0.025','0.01', '0.0075','0.005','0.0025', '0.001','0.00075','5e-04','1e-04', '1e-05', '1e-10','1e-30', '0')
FDRcutoffs = c(0.9,0.5,0.4,0.3,0.2,0.1,0.1,0.2,0.1,0.2,0.01) #
FCcutoffs = c(1.5,1.5,1.5,1.5,1.5,1.5,2,2,1,1,1) #

library(stringr)
plateIDs = c('met1','met2','met3','met4','met5','met6','met7','met8','met9','met10','met11','met13')
allConditions = c()
for (plateID in plateIDs){
  # load the FC matrix from DE
  load(paste('./../1_QC_dataCleaning/outputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  allConditions = c(allConditions, unique(paste(RNAi,batchLabel)))
}
RNAiName = c()
batchID = c()
for (i in 1:length(allConditions)){
  RNAiName = c(RNAiName, strsplit(allConditions[i],' ')[[1]][[1]])
  batchID = c(batchID, strsplit(allConditions[i],' ')[[1]][[2]])

}

# load a set of manually curated vectorlike conditions
vectorlikeTbl = read.csv('./../../MetabolicLibrary/input_data/metaData/manually_curated_vectorlike_conditions.csv')
vectorlikeTbl$vectorlike_RNAi = str_remove(vectorlikeTbl$NTP_condID,' met[0-9]+_lib.$')
vectorlikes = RNAiName[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                       | str_detect(RNAiName, '^x.REALVECTOR_')]
batch_vectorlikes = batchID[RNAiName %in% vectorlikeTbl$vectorlike_RNAi
                            | str_detect(RNAiName, '^x.REALVECTOR_')]
# some of the realvectors are technical replicates, we remove them
rmInd = str_detect(vectorlikes, '^x.REALVECTOR_') & batch_vectorlikes %in% c('met11_lib1','met11_lib2','met11_lib3','')
vectorlikes = vectorlikes[!rmInd]
batch_vectorlikes = batch_vectorlikes[!rmInd]


# do this evaluation with real vectors
realvectors = RNAiName[str_detect(RNAiName, '^x.REALVECTOR_')]
batch_realvectors = batchID[str_detect(RNAiName, '^x.REALVECTOR_')]
rmInd = str_detect(realvectors, '^x.REALVECTOR_') & batch_realvectors %in% c('met11_lib1','met11_lib2','met11_lib3','')
realvectors = realvectors[!rmInd]
batch_realvectors = batch_realvectors[!rmInd]


FCtypes = c("log2FoldChange","log2FoldChange_raw")

for (FCtype in FCtypes){
  #str = '0.005'
  N_DE_quant = list()
  N_DE_mean = list()
  for(zz in 1:length(FDRcutoffs)){
    N_DE_quant[[zz]] = vector() #
    N_DE_mean[[zz]] = vector() #
  }
  for (str in strvector){
    library(DESeq2)

    # aggregate the raw tables for cutoff titrating
    mergeTbl = data.frame()
    for (i in 1:length(vectorlikes)){
      tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',vectorlikes[i],'_',batch_vectorlikes[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
      tbl$RNAi = vectorlikes[i]
      tbl$batchID = batch_vectorlikes[i]
      mergeTbl = rbind(mergeTbl, tbl)
    }
    maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)

    for(zz in 1:length(FDRcutoffs)){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FCcutoffs[zz]) & mergeTbl$dualFDR < FDRcutoffs[zz]),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(vectorlikes, batch_vectorlikes), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_quant[[zz]] = c(N_DE_quant[[zz]], quantile(N_DE,0.9)) #
    }

    # aggregate the raw tables for cutoff titrating
    mergeTbl = data.frame()
    for (i in 1:length(realvectors)){
      tbl = read.csv(paste('output/clean_DE_output/DE_table_lfcShrink_merged_clean_cutoff_',str,'_RNAi_',realvectors[i],'_',batch_realvectors[i],'_imputMethod_random_strat_',FCtype,'.csv',sep = ''))
      tbl$RNAi = realvectors[i]
      tbl$batchID = batch_realvectors[i]
      mergeTbl = rbind(mergeTbl, tbl)
    }
    maxExp = apply(mergeTbl[,c('medianCount_RNAi','medianCount_ctr')],1,max)
    # N_DE_list = list()

    for(zz in 1:length(FDRcutoffs)){
      passTbl = mergeTbl[which(abs(mergeTbl[,FCtype]) > log2(FCcutoffs[zz]) & mergeTbl$dualFDR < FDRcutoffs[zz]),]
      passTbl$ID = paste(passTbl$RNAi, passTbl$batchID)
      N_DE = table(passTbl$ID)
      zeros = setdiff(paste(realvectors, batch_realvectors), names(N_DE))
      tmp = rep(0, length(zeros))
      names(tmp) = zeros
      N_DE = c(N_DE, tmp)
      N_DE_mean[[zz]] = c(N_DE_mean[[zz]], mean(N_DE)) #
    }

    print(str)
  }
  for(zz in 1:length(FDRcutoffs)){
    pdf(paste('figures/0_DE_QA/benchmark_vectorlikes_titrate_cleaning_cutoff_DE_',FCtype,'_FDR',FDRcutoffs[zz],'_FC',FCcutoffs[zz],'.pdf',sep = '')) # not a good metric
    df = data.frame(cutoff_all, N_DE_quant = N_DE_quant[[zz]])
    print(ggplot(df, aes(log10(cutoff_all), N_DE_quant,label = cutoff_all)) +
            geom_point(color = "red")+geom_line()+
            xlab('p_outlier cutoff')+
            ylab('90% quantile #DE in vectorlikes') +
            geom_text_repel(max.overlaps =Inf, size = 3))

    df = data.frame(cutoff_all, N_DE_mean = N_DE_mean[[zz]])
    print(ggplot(df, aes(log10(cutoff_all), N_DE_mean,label = cutoff_all)) +
            geom_point(color = "red")+geom_line()+
            xlab('p_outlier cutoff')+
            ylab('average #DE in vectorlikes') +
            geom_text_repel(max.overlaps =Inf, size = 3))
    dev.off()
  }

}



