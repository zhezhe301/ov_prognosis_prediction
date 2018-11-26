#### Step1 Data Collection and Feature Pre-Selection

## TCGA (GDC) data collection
# patients info. The 23rd column was the treatment outcome represented by 1 and 2 
# where 1 means sensensitive and 2 means resistant 
# setwd('../../')
# please change file path on your local machine
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- core[complete.cases(core[,21]),]
core <- core[complete.cases(core[,22]),]
core <- as.data.frame(core)

# The file 'RNASeq_pvalue.csv' was the Fisher test for pre-selection molecular features 
# total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
# core <- total
# core <- core[complete.cases(core[,21]),]
# core <- core[complete.cases(core[,22]),]
# core <- as.data.frame(core)
# RNAseqall <- read.csv('tcga.RNAseq.csv',row.names = 1)
# RNAseq <- RNAseqall[which(rownames(RNAseqall)%in%rownames(core)),]
# 
# txoutcome <- core[rownames(RNAseq),23]
# table <- cbind(txoutcome, RNAseq)
# for (i in 1:length(colnames(table))) {
#   table[which(is.na(table[,i])==TRUE),i] <- median(table[,i],na.rm=T)
# }
# table[,2:length(colnames(table))] <- log2(table[,2:length(colnames(table))] +1)
# pmat <- matrix(nrow= length(colnames(table)), ncol=2)
# for (i in 2:length(colnames(table))) {
#   pmat[i,1] <- wilcox.test(table[,i]~table[,1])$p.value
# }
# rownames(pmat) <- colnames(table)[1:length(colnames(table))]
# write.csv(pmat, 'RNASeq_pvalue.csv')
list1 <- read.csv('RNASeq_pvalue.csv',row.names = 1)
list1<- as.matrix(list1)
list2 <- sort(list1[,1])
list3 <- list2[1:1282]  #<0.05, event :184
list3 <- list3[1:184]
# TCGA data
RNAseqall <- read.csv('tcga.RNAseq.csv',row.names = 1)
RNAseq <- RNAseqall[which(rownames(RNAseqall)%in%rownames(core)), intersect(names(list3),colnames(RNAseqall))]
txoutcome <- core[rownames(RNAseq),23]
table <- cbind(txoutcome, RNAseq)
for (i in 1:length(colnames(table))) {
  table[which(is.na(table[,i])==TRUE),i] <- median(table[,i],na.rm=T)
}
table[,2:length(colnames(table))] <- log2(table[,2:length(colnames(table))] +1) 

# save the files for next step in Python 
nrFolds <- 5
folds <- sample(rep_len(1:nrFolds, nrow(table)))
setwd('./2.tx.outcome.model/2.RNAseq')
for (k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- table[-fold,]
  data.test <- table[fold,]

  write.csv(data.train,paste('classifier.trainset.5fold',k,'.csv',sep=''))
  write.csv(data.test,paste('classifier.testset.5fold',k,'.csv',sep=''))
}

## PCT (or PCA method)to further reduce the dimensions of the data (file name of code: DNAseq.pca.5fold.py)
table1 <- table[,-1]
table1 <- as.matrix(table1)
pca <- princomp(table1,cor = TRUE)
#summary(pca)
## 1-45:80%, 1-56:85%, 1-71:90%, 1-93 :95%
pca.table <- table1 %*% pca$loadings[,1:54]
pca.table <- cbind(txoutcome, pca.table)

# save the files for next step in Python 
nrFolds <- 5
folds <- sample(rep_len(1:nrFolds, nrow(pca.table)))
for (k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- pca.table[-fold,]
  data.test <- pca.table[fold,]
  write.csv(data.train,paste('classifier.trainset.pca.5fold',k,'.csv',sep = ''))
  write.csv(data.test,paste('classifier.testset.pca.5fold',k,'.csv',sep=''))
}


#### Step2 Model Training and Predictions

### five-fold cross validation with and without PCA see in Python

### the C-index calculation
library(Hmisc) 

mat <- matrix(nrow= 100, ncol= 6)
for (k in 1:100) {
  f <- paste('dt_classifier.txoutcome.5fold', k, '.txt', sep="")
  treeclass <- read.table(f,header=F)
  colnames(treeclass) <- c('actual','pred')
  cindex.dt <- rcorr.cens(treeclass[,2], treeclass[,1])["C Index"]
  mat[k,1] <- cindex.dt
  
  f <- paste('rf_classifier.txoutcome.5fold', k, '.txt', sep="")
  rfclass <- read.table(f,header=F)
  colnames(rfclass) <- c('actual','pred')
  cindex.rf <- rcorr.cens(rfclass[,2], rfclass[,1])["C Index"] 
  mat[k,2] <- cindex.rf
  
  f <- paste('dnn_classifier.txoutcome.5fold', k, '.txt', sep="")
  dnnclass <- read.table(f,header=F)
  colnames(dnnclass) <- c('actual','pred')  
  cindex.dnnclass <- rcorr.cens(dnnclass[,2], dnnclass[,1])["C Index"] 
  mat[k,3] <- cindex.dnnclass
  
  f <- paste('dt_classifier.txoutcome.pca.5fold', k, '.txt', sep="")
  treeclass <- read.table(f,header=F)
  colnames(treeclass) <- c('actual','pred')
  cindex.dt <- rcorr.cens(treeclass[,2], treeclass[,1])["C Index"]
  mat[k,4] <- cindex.dt
  
  f <- paste('rf_classifier.txoutcome.pca.5fold', k, '.txt', sep="")
  rfclass <- read.table(f,header=F)
  colnames(rfclass) <- c('actual','pred')
  cindex.rf <- rcorr.cens(rfclass[,2], rfclass[,1])["C Index"] 
  mat[k,5] <- cindex.rf
  
  f <- paste('dnn_classifier.txoutcome.pca.5fold', k, '.txt', sep="")
  dnnclass <- read.table(f,header=F)
  colnames(dnnclass) <- c('actual','pred')  
  cindex.dnnclass <- rcorr.cens(dnnclass[,2], dnnclass[,1])["C Index"] 
  mat[k,6] <- cindex.dnnclass  
}
colnames(mat) <- c('dt.outcome.5fold','rf.outcome.5fold','dnnclass.outcome.5fold','dt.outcome.pca.5fold','rf.outcome.pca.5fold','dnnclass.outcome.pca.5fold')
write.csv(mat, 'mat.100python.txoutcome.csv')







