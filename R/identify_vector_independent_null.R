Args <- commandArgs()
ind = as.numeric(Args[6])

source('AdaReg.R')
library(edgeR)
library(stringr)
library(limma)
library(edgeR)
library(DESeq2)
# this zscore is called AdaTiSS score https://www.biorxiv.org/content/10.1101/869404v1.full.pdf
calcAdaRegZ <- function(plateID, lib){
  # load the FC matrix from DE
  load(paste('./inputs/cleaned_merged_raw_data_',plateID,'.Rdata',sep = ''))
  if (lib != 'extra_met2_lib1'){
    # pick the control samples (in the batches belonging to target DE group)
    subsetInd = batchLabel == lib
    input_subset = input[,subsetInd] 
    batchLabel_subset = colnames(input_subset)
    batchLabel_subset = str_extract(batchLabel_subset,'_rep._met[0-9]+_')
    batchLabel_subset = str_replace(batchLabel_subset,'_met._$','')
    batchLabel_subset = as.factor(batchLabel_subset)
    RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
    coldata =  data.frame(batchLabel = batchLabel_subset, RNAi = RNAi_subset)
    rownames(coldata) = colnames(input_subset)
    
    dds <- DESeqDataSetFromMatrix(countData = input_subset,
                                  colData = coldata,
                                  design= ~ batchLabel + RNAi)
    # filtering 
    keep <- rowSums(counts(dds)>=10) >= 1
    dds <- dds[keep,]
    # calculating z-scores 
    dds = estimateSizeFactors(dds)
    tmp = counts(dds, normalized=TRUE)
    normCounts = log2(tmp + 1)
    # remove batch effects in replicates 
    normCounts_clean <- limma::removeBatchEffect(normCounts, batch = dds$batchLabel,design= model.matrix(~ dds$RNAi))
    
  }else{
    subsetInd = str_detect(colnames(input),'_rep3_extra_met2_lib5')
    input_subset = input[,subsetInd] 
    RNAi_subset = dropEmptyLevels(RNAi[subsetInd])
    coldata =  data.frame(RNAi = RNAi_subset)
    rownames(coldata) = colnames(input_subset)
    
    dds <- DESeqDataSetFromMatrix(countData = input_subset,
                                  colData = coldata,
                                  design= ~ RNAi)
    # filtering 
    keep <- rowSums(counts(dds)>=10) >= 1
    dds <- dds[keep,]
    # calculating z-scores 
    dds = estimateSizeFactors(dds)
    tmp = counts(dds, normalized=TRUE)
    normCounts_clean = log2(tmp + 1)
  }

  # fit every gene 
  zMat = matrix(NA,nrow = nrow(normCounts_clean),ncol = ncol(normCounts_clean))
  rownames(zMat) = rownames(normCounts_clean)
  colnames(zMat) = colnames(normCounts_clean)
  case1 = c()
  case2 = c()
  for (i in 1:nrow(normCounts_clean)){
    y = normCounts_clean[i,]
    out = AdaReg(model.matrix(~1,data = as.data.frame(y)), y)
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
      if (median(y) <= log2(10+1)){
        meanIn = mean(y)
        sdIn = sd(y)
        zr1 = (y-meanIn)/sdIn
        zMat[i,] = zr1
        case2 = c(case2,i)
      }else{
        meanIn = median(y)
        sdIn = mad(y)
        zr1 = (y-meanIn)/sdIn
        zMat[i,] = zr1
        case2 = c(case2,i)
      }

    }
    
    if (i %% 1000 == 0){
      print(paste('AdaTiss fitting for ',lib, ' ... ', 100*i/nrow(zMat),'%',sep = ''))
      print(paste('normal fit: ',length(case1)/i,sep = ''))
      print(paste('poor fit: ',length(case2)/i,sep = ''))
    }
  }
  
  write.csv(zMat,paste('outputs/adaZ_',lib,'.csv',sep = ''))
  return(zMat)
}



taskList = read.csv('taskList.csv',colClasses = 'character')

calcAdaRegZ(taskList$plateID[ind], taskList$lib[ind])





# inspection of the adaReg algorithm
# i = 6
# x1 = c()
# x2 = c()
# x3 = c()
# for (i in 1:6){
#   a = counts(dds,normalize = T, replace = F)[geneName,dds$batchLabel == levels(dds$batchLabel)[i]]
#   #y = scale(log2(1+a),center = median(log2(1+a)),scale = F)
#   y = log2(1+a)
#   
#   # try center the three reps seperately!!
#   
#   reps = names(a)
#   reps = str_extract(reps,'_rep._met[0-9]+_')
#   reps = str_replace(reps,'_met._$','')
#   reps = as.factor(reps)
# 
#   out = AdaReg(model.matrix(~1,data = as.data.frame(y)), y)#+ reps
#   
#   #out = AdaReg(model.matrix(~dropEmptyLevels(dds$RNAi[dds$batchLabel == levels(dds$batchLabel)[i]])), y,0)#+ reps
#   #out = AdaReg(model.matrix(~'x.vector' == dropEmptyLevels(dds$RNAi[dds$batchLabel == levels(dds$batchLabel)[i]]) ), y)#+ reps
#   #dev.new()
#   #hist(out$x.res, freq=FALSE, main="Histogram of residuals")
#   
#   zr1 = (y-out$beta.rob.fit)/sqrt(out$var.sig.gp.fit)
#   x1 = c(x1,mean(zr1[1:6]))
#   dev.new()
#   hist(zr1, freq=FALSE, main="Histogram of adareg")
#   curve(out$res.info["pi0.hat"] * dnorm(x, 0, 1), add=TRUE, col="red", lwd=2)
#   
#   zr2 = (y- median(y))/(1.486*mad(y))
#   x2 = c(x2,mean(zr2[1:6]))
#   dev.new()
#   hist(zr2, freq=FALSE, main="Histogram of rob z")
#   curve(1 * dnorm(x, 0, 1), add=TRUE, col="red", lwd=2)
#   
#   zr3 = (y - mean(y))/sd(y)
#   x3 = c(x3,mean(zr3[1:6]))
#   dev.new()
#   hist(zr3, freq=FALSE, main="Histogram of zscores")
#   curve(1 * dnorm(x, 0, 1), add=TRUE, col="red", lwd=2)
#   
# }
# x1
# x2
# x3
# plot(x1,x2)
# text(x1,x2, levels(batchLabel))
# plot(x1,x3)
# text(x1,x3, levels(batchLabel))
# plot(x2,x3)
# text(x2,x3, levels(batchLabel))
# 
# 
# 
# 
# 
# 
# 
# adareg.result = out
# adareg.coef = adareg.result$beta.rob.fit
# 
# adareg.res = adareg.result$x.res
# 
# res.fit = adareg.result$res.info
# hist(adareg.res, freq=FALSE, main="Histogram of residuals")
# 
# curve(res.fit["pi0.hat"] * dnorm(x, res.fit["mu0.hat"], res.fit["sd0.hat"]), min(adareg.res, na.rm=TRUE), max(adareg.res, na.rm=TRUE), add=TRUE, col="red", lwd=2)
# curve(res.fit["pi0.hat"] * dnorm(x, mean(y), sd(y)), min(adareg.res, na.rm=TRUE), max(adareg.res, na.rm=TRUE), add=TRUE, col="blue", lwd=2)
# curve(res.fit["pi0.hat"] * dnorm(x, res.fit["mu0.hat"], (1.486*mad(y))), min(adareg.res, na.rm=TRUE), max(adareg.res, na.rm=TRUE), add=TRUE, col="red", lwd=2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# out2 = AdaReg(model.matrix(~rep(1,51))[,1,drop = F], y)
# zr4 = (y-out2$res.info['mu0.hat'])/out2$res.info['sd0.hat']
# 
# hist(zr2, freq=FALSE, main="Histogram of zscores")
# curve(1 * dnorm(x, 0, 1), add=TRUE, col="red", lwd=2)
# 
# 
# 
