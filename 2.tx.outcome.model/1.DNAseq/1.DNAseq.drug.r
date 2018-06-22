#### Step1 Data Collection and Feature Pre-Selection

## TCGA (GDC) data collection
# patients info. The 23rd column was the treatment outcome represented by 1 and 2 
# where 1 means sensensitive and 2 means resistant 
setwd('/Users//zhangzhe//Project//ov_prognosis_drugprecision//result//1.')
total <- read.table('totalpatient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- core[complete.cases(core[,21]),]
core <- core[complete.cases(core[,22]),]
core <- as.data.frame(core)

# the file 'DNAseq.pvalue.csv' was the Fisher test for pre-selection molecular features 
########################################################################
# setwd('/data1/users/luym/project/tcga/1.DNASeq/old')
# outcome_total <- read.table('outcome_summary.txt', header=T,row.names=1)
# outcome_core <- outcome_total[complete.cases(outcome_total[,2]),]
# sen <- row.names(outcome_core)[which(outcome_core[,2]=='Sensitive')]
# res <- row.names(outcome_core)[which(outcome_core[,2] %in% c('Resistant','Persistent','Refractory'))]
# 
# 
# DNASeq_total <- read.table('DNASeq_summary.txt',header=T,sep='\t',quote=NULL, comment='')
# DNASeq_total[which(DNASeq_total[,1]%in%sen), 9] <-'Sensitive'
# DNASeq_total[which(DNASeq_total[,1]%in%res), 9] <-'Resistant'
# 
# DNASeq_core <- DNASeq_total[complete.cases(DNASeq_total[,9]),]
# DNASeq_core <- DNASeq_core[which(DNASeq_core[,2]=='ov'),]
# count <- DNASeq_core[which(DNASeq_core[,9]=='Sensitive'),]
# len_sen <- length(unique(count[,1]))
# count <- DNASeq_core[which(DNASeq_core[,9]=='Resistant'),]
# len_res <- length(unique(count[,1]))
# 
# 
# hugo <- as.data.frame(summary(DNASeq_core[,3],maxsum = 50000))
# res_mat <- matrix(NA,ncol=11,nrow=length(hugo[,1]))
# for (i in 1:length(hugo[,1])) {
#   gene_all <- DNASeq_core[which(DNASeq_core[,3]==row.names(hugo)[i]),]
#   gene_sen <- gene_all[which(gene_all[,9]=='Sensitive'),]
#   yes <- length(unique(gene_sen[,1]))
#   no <- len_sen - yes
#   p_value <- data.frame(yes, no)
#   res_mat[i,1] <- yes
#   res_mat[i,2] <- no
#   gene_res <- gene_all[which(gene_all[,9]=='Resistant'),]
#   yes <- length(unique(gene_res[,1]))
#   no <- len_res - yes
#   p_value[2,1] <- yes
#   p_value[2,2] <- no
#   res_mat[i,3] <- yes
#   res_mat[i,4] <- no
#   res_mat[i,9] <- chisq.test(p_value)$p.value
#   res_mat[i,10]<- fisher.test(p_value)$p.value
#   a <- p_value[3,1]+p_value[4,1]
#   b <- p_value[1,1]
#   c <- p_value[3,2]+p_value[4,2]
#   d <- p_value[1,2]
#   res_mat[i,11] <- a/(a+b)*(c+d)/c
# }
# hugo <- cbind(hugo, res_mat)
# colnames(hugo) <- c('Genecount', 'SenYes', 'SenNo','ResYes','ResNo', 'PerYes','PerNo','RefYes','RefNo','p_Value','Fisher_p','RR')
# write.csv(hugo,"DNAseq.pvalue.csv",quote=F)
setwd('~/Project/ov_prognosis_drugprecision/result/5.tx.pvalue/')
list1 <- read.csv('DNAseq.pvalue.csv',row.names = 1)
list1<- as.matrix(list1)
list2 <- sort(list1[,11])
list3 <- list2[1:43]  # p<0.05, event :266

# TCGA data
setwd('~/Project/ov_prognosis_drugprecision/result/')
DNAseqall <- read.table('ov_DNASeq_summ2.txt',header=T, sep='\t',quote='',row.names=1)
DNAseq <- DNAseqall[which(rownames(DNAseqall)%in%rownames(core)),intersect(names(list3),colnames(DNAseqall))]
txoutcome <- core[rownames(DNAseq),23]
table <- cbind(txoutcome, DNAseq)

# save the files for next step in Python 
nrFolds <- 5
folds <- sample(rep_len(1:nrFolds, nrow(table)))
for (k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- table[-fold,]
  data.test <- table[fold,]
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/')
  write.csv(data.train,paste('classifier.trainset.5fold',k,'.csv',sep=''))
  write.csv(data.test,paste('classifier.testset.5fold',k,'.csv',sep=''))
}


## PCT (or PCA method)to further reduce the dimensions of the data (file name of code: DNAseq.pca.5fold.py)
table1 <- table[,-1]
table1 <- as.matrix(table1)
pca <- princomp(table1,cor = TRUE)
#summary(pca)
## 1-22:80%, 1-26:85%, 1-29:90%, 1-33 :95%
pca.table <- table1 %*% pca$loadings[,1:40]
pca.table <- cbind(txoutcome, pca.table)

# save the files for next step in Python 
nrFolds <- 5
folds <- sample(rep_len(1:nrFolds, nrow(pca.table)))
for (k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- pca.table[-fold,]
  data.test <- pca.table[fold,]
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/')
  write.csv(data.train,paste('classifier.trainset.pca.5fold',k,'.csv',sep = ''))
  write.csv(data.test,paste('classifier.testset.pca.5fold',k,'.csv',sep=''))
}


#### Step2 Model Training and Predictions

### five-fold cross validation with and without PCA see in Python

### the C-index calculation
library(Hmisc) 
mat <- matrix(nrow= 100, ncol= 6)
for (k in 1:100) {
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/dt_classifier/')
  f <- paste('dt_classifier.txoutcome.5fold', k, '.txt', sep="")
  treeclass <- read.table(f,header=F)
  colnames(treeclass) <- c('actual','pred')
  cindex.dt <- rcorr.cens(treeclass[,2], treeclass[,1])["C Index"]
  mat[k,1] <- cindex.dt
  
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/rf_classifier/')
  f <- paste('rf_classifier.txoutcome.5fold', k, '.txt', sep="")
  rfclass <- read.table(f,header=F)
  colnames(rfclass) <- c('actual','pred')
  cindex.rf <- rcorr.cens(rfclass[,2], rfclass[,1])["C Index"] 
  mat[k,2] <- cindex.rf
  
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/dnn_classifier/')
  f <- paste('dnn_classifier.txoutcome.5fold', k, '.txt', sep="")
  dnnclass <- read.table(f,header=F)
  colnames(dnnclass) <- c('actual','pred')  
  cindex.dnnclass <- rcorr.cens(dnnclass[,2], dnnclass[,1])["C Index"] 
  mat[k,3] <- cindex.dnnclass
  
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/dt_classifier/')
  f <- paste('dt_classifier.txoutcome.pca.5fold', k, '.txt', sep="")
  treeclass <- read.table(f,header=F)
  colnames(treeclass) <- c('actual','pred')
  cindex.dt <- rcorr.cens(treeclass[,2], treeclass[,1])["C Index"]
  mat[k,4] <- cindex.dt
  
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/rf_classifier/')
  f <- paste('rf_classifier.txoutcome.pca.5fold', k, '.txt', sep="")
  rfclass <- read.table(f,header=F)
  colnames(rfclass) <- c('actual','pred')
  cindex.rf <- rcorr.cens(rfclass[,2], rfclass[,1])["C Index"] 
  mat[k,5] <- cindex.rf
  
  setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq/dnn_classifier/')
  f <- paste('dnn_classifier.txoutcome.pca.5fold', k, '.txt', sep="")
  dnnclass <- read.table(f,header=F)
  colnames(dnnclass) <- c('actual','pred')  
  cindex.dnnclass <- rcorr.cens(dnnclass[,2], dnnclass[,1])["C Index"] 
  mat[k,6] <- cindex.dnnclass  
}
colnames(mat) <- c('dt.outcome.5fold','rf.outcome.5fold','dnnclass.outcome.5fold','dt.outcome.pca.5fold','rf.outcome.pca.5fold','dnnclass.outcome.pca.5fold')
setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction/1.DNAseq')
write.csv(mat, 'mat.100python.txoutcome.csv')



