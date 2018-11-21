### Correlation plot between GDC platform and ICGC platform. 
## (A) gene expression data (RNA expression) 
setwd('/data1/users/luym/project/ov_prediction/data')
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
setwd('/data1/users/luym/project/ov_prediction/data')
RNAseqall <- read.csv('tcga.RNAseq.csv',row.names = 1)
RNAseq <- RNAseqall[which(rownames(RNAseqall)%in%rownames(core)),]
vital_status <- core[rownames(RNAseq),3]
days <- core[rownames(RNAseq),4]
table <- cbind(days, vital_status, RNAseq)
library("survival")
attach(table)
table$days <- as.numeric(table$days)
table[,2] <- as.numeric(table[,2])
list1 <- c()
for (i in colnames(table)[3:20531]) {
  z <- Surv(table$days,table$vital_status==2)
  y <- coxph(z~table[,i])
  y.s <- summary(y)
  list1 <- append(list1,y.s$logtest['pvalue'])
}
names(list1) <- colnames(table)[3:20531]
list2 <- sort(list1)[1:1685] ####p<0.05
list3 <- sort(list1)[1:184]
setwd('/data1/users/luym/project/ov_prediction/result/2.RNAseq')
write.table(list3,'coxRNAseq.txt',sep='\t')
list3 <-   read.table('coxRNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}
setwd('/data1/users/luym/project/ov_prediction/data')
testset <- read.csv('icgc.RNAseq.csv')
rownames(testset) <- testset[,2]
testset <- testset[,-c(1,2)]
surv.info <- read.table('icgc.patient.txt',header=T,sep='\t',row.names = 1)
surv.info <- as.data.frame(surv.info)
vital_status <- surv.info[rownames(testset),1]
days <- surv.info[rownames(testset),6]
testcore <- cbind(days,vital_status,testset)
intersectname_cox <- intersect(colnames(table),colnames(testset))
trainset.aftercox <- table[,intersectname_cox]
for (i in 1:length(colnames(trainset.aftercox))) {
  trainset.aftercox[which(is.na(trainset.aftercox[,i])==TRUE),i] <- median(trainset.aftercox[,i],na.rm=T)
}
testset.aftercox <- testset[,intersectname_cox]
for (i in 1:length(colnames(testset.aftercox))) {
  testset.aftercox[which(is.na(testset.aftercox[,i])==TRUE),i] <- median(testset.aftercox[,i],na.rm=T)
}
setwd('~/Project/ov_prognosis_drugprecision/result/9.correlation_two_database/')
setEPS()
postscript('RNAseq_correlation5.eps',width = 8,height = 8)
x <- apply(trainset.aftercox,2,median)
y <- apply(testset.aftercox,2,median)
y <- y[-which(x==0)]
x <- x[-which(x==0)]
x <- x[-which(y==0)]
y <- y[-which(y==0)]
x <- log(x,2)
y <- log(y,2)
cor.test(x,y)
plot(x,y,col='darkgreen',pch=16,cex=0.3,ylim = c(-3,15))
dev.off()  

## (B) DNA methylation data
setwd('/data1/users/luym/project/ov_prediction/data')
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
setwd('/data1/users/luym/project/ov_prediction/data')
DNAmethall <- read.table('tgca.DNAmethylation.txt',header=T, sep='\t',quote='',row.names=1)
DNAmeth <- DNAmethall[which(rownames(DNAmethall)%in%rownames(core)),]
vital_status <- core[rownames(DNAmeth),3]
days <- core[rownames(DNAmeth),4]
table <- cbind(days, vital_status, DNAmeth)
table <- table[,-(which(is.na(table[1,])==TRUE))]
library("survival")
attach(table)
table$days <- as.numeric(table$days)
table[,2] <- as.numeric(table[,2])
list1 <- c()
for (i in colnames(table)[3:14103]) {
  z <- Surv(table$days,table$vital_status==2)
  y <- coxph(z~table[,i])
  y.s <- summary(y)
  list1 <- append(list1,y.s$logtest['pvalue'])
}
names(list1) <- colnames(table)[3:14103]
list2 <- sort(list1)[1:362] ####p<0.05
list3 <- sort(list1)[1:344]
setwd('/data1/users/luym/project/ov_prediction/result/4.DNAmethy')
write.table(list3,'coxDNAmeth.txt',sep='\t')
setwd('/data1/users/luym/project/ov_prediction/result/4.DNAmethy')
list3 <-   read.table('coxDNAmeth.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}
setwd('/data1/users/luym/project/ov_prediction/data')
testset <- read.table('icgc.DNAmethylation.txt',sep='\t',header=T,row.names = 1)
surv.info <- read.table('icgc.patient.txt',header=T,sep='\t',row.names = 1)
surv.info <- as.data.frame(surv.info)
vital_status <- surv.info[rownames(testset),1]
days <- surv.info[rownames(testset),6]
testcore <- cbind(days,vital_status,testset)
for (i in 1:length(colnames(testcore))) {
  testcore[which(is.na(testcore[,i])==TRUE),i] <- median(testcore[,i],na.rm=T)
}
intersectname_cox <- intersect(colnames(table),colnames(testset))
trainset.aftercox <- table[,intersectname_cox]
for (i in 1:length(colnames(trainset.aftercox))) {
  trainset.aftercox[which(is.na(trainset.aftercox[,i])==TRUE),i] <- median(trainset.aftercox[,i],na.rm=T)
}

testset.aftercox <- testset[,intersectname_cox]
for (i in 1:length(colnames(testset.aftercox))) {
  testset.aftercox[which(is.na(testset.aftercox[,i])==TRUE),i] <- median(testset.aftercox[,i],na.rm=T)
}
setwd('~/Project/ov_prognosis_drugprecision/result/9.correlation_two_database/')
setEPS()
postscript('DNAmethy_correlation.eps',width = 8,height = 8)
x <- apply(trainset.aftercox,2,median)
y <- apply(testset.aftercox,2,median)
cor.test(x,y)
plot(x,y,col='darkgreen',pch=16,cex=0.2)
dev.off()

