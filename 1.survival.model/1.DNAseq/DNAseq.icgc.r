#### Step1 Data Collection and Feature Pre-Selection
## TCGA (also known as GDC) data collection
# setwd('../../')
# please change file path on your local machine
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)

DNAseqall <- read.table('tcga.DNAseq.txt',header=T, sep='\t',quote='',row.names=1)
DNAseq <- DNAseqall[which(rownames(DNAseqall)%in%rownames(core)),]
vital_status <- core[rownames(DNAseq),3]
days <- core[rownames(DNAseq),4] ## 1==alive, 2==dead #####
table <- cbind(days, vital_status, DNAseq)
table <- table[,-(which(is.na(table[1,])==TRUE))]
## Cox proportional hazards regression identifying the most significant features
library(survival)
attach(table)
table$days <- as.numeric(table$days)
table[,2] <- as.numeric(table[,2])
list1 <- c()
for (i in colnames(table)[3:12363]) {
  z <- Surv(table$days,table$vital_status==2)
  y <- coxph(z~table[,i])
  y.s <- summary(y)
  list1 <- append(list1,y.s$logtest['pvalue'])
}
names(list1) <- colnames(table)[3:12363]
list2 <- sort(list1)[1:892] ####p<0.05, events :266
list3 <- sort(list1)[(1+22):(266+22)]

write.table(list3,'coxDNAseq.txt',sep='\t')
list3 <- read.table('coxDNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=TRUE)
}
## ICGC data collection
testset <- read.csv('icgc.DNAseq.csv')
rownames(testset) <- testset[,2]
testset <- testset[,-c(1,2,3)]
surv.info <- read.table('icgc.patient.txt',header=T,sep='\t',row.names = 1)
surv.info <- as.data.frame(surv.info)
vital_status <- surv.info[rownames(testset),1]
days <- surv.info[rownames(testset),6]
testcore <- cbind(days,vital_status,testset)
testcore[,1] <- as.numeric(testcore[,1])
testcore[,2] <- as.numeric(testcore[,2])
# for (i in 1:length(colnames(testcore))) {
#   testcore[which(is.na(testcore[,i])==TRUE),i] <- median(testcore[,i],na.rm=T)
# }


## intersect between TCGA and ICGC
intersectname_cox <- intersect(rownames(list3),colnames(testset))
trainset.aftercox <- cbind(table[,c(1,2)],table[,intersectname_cox])
for (i in 1:length(colnames(trainset.aftercox))) {
  trainset.aftercox[which(is.na(trainset.aftercox[,i])==TRUE),i] <- median(trainset.aftercox[,i],na.rm=T)
}

testset.aftercox <- cbind(testcore[,c(1,2)],testcore[,intersectname_cox])
for (i in 1:length(colnames(testset.aftercox))) {
  testset.aftercox[which(is.na(testset.aftercox[,i])==TRUE),i] <- median(testset.aftercox[,i],na.rm=T)
}
setwd('./1.survival.model/1.DNAseq')
write.csv(trainset.aftercox,'trainset.icgc.csv')
write.csv(testset.aftercox,'testset.icgc.csv')

#### Step2 Model Training and Predictions

### ICGC independent validation without PCA 

## LASSO and Cox * 100 times
mat.100 <- matrix(nrow=100, ncol=2)
library(Hmisc) 
library(randomForestSRC)
library("glmnet")
library("survival")
for (k in 1:100) {
  repeat{
    x <- as.matrix(table3[,-c(1,2)])
    y <- as.matrix(table3[,c(1,2)])
    y[which(y[,2]==1),] <- 0
    y[which(y[,2]==2),] <- 1
    y[,1] <- table3[,1]
    colnames(y) <- c('time', 'status')
    fit <- glmnet(x,y,family='cox')
    plot(fit)
    
    fit_cv=cv.glmnet(x,y,family="cox",nfolds = 10)
    best_lambda <- fit_cv$lambda.min
    fit.coef <- coef(fit,s= fit_cv$lambda.min)
    listgene_lasso <- rownames(fit.coef)[which(fit.coef != 0)]
    if (length(listgene_lasso) >0) {
      break
    }
  }
 	max.dev.index     <- which.max(fit$dev.ratio)
 	optimal.lambda <- fit$lambda[max.dev.index]
 	optimal.beta  <- fit$beta[,max.dev.index]
 	nonzero.coef <- abs(optimal.beta)>0
 	selectedXname <- names(nonzero.coef)[which(nonzero.coef, arr.ind = TRUE)]
 	intersectname <- intersect(selectedXname,colnames(testcore))
 	selectedtestX <- testcore[,intersectname]
 	selectedtestX <- as.matrix(selectedtestX)
 	# selectedtestX <- scale(selectedtestX,scale=F)
 	selectedBeta <- optimal.beta[intersectname]
 	testcore$days <- as.numeric(testcore$days)
 	testcore$vital_status <- as.numeric(testcore$vital_status)
 	coxph.model <- coxph(Surv(testcore$days,testcore$vital_status) ~selectedtestX,init=scale(selectedBeta),iter=0)
	predict <- survfit(coxph.model,newdata= testcore[,intersectname])
	xxx <- summary(predict)
  predictval <- xxx$table[,'median']
  rcorr.cens(predictval, Surv(testcore$days, testcore$vital_status))["C Index"]
  cindex <- 1-rcorr.cens(predict(coxph.model,newdata=testcore[,intersectname],type="risk"), Surv(testcore$days, testcore$vital_status))["C Index"]
  if (cindex >=0.5) {
    mat.100[k,1] <- cindex
  } else {
    mat.100[k,1] <- 1 - cindex
  }
  

	## random forest * 100 times
	intersectname_cox <- intersect(rownames(list3),colnames(testset))
	trainset.aftercox <- cbind(table3[,c(1,2)],table3[,intersectname_cox])
	for (i in 1:length(colnames(trainset.aftercox))) {
	  trainset.aftercox[which(is.na(trainset.aftercox[,i])==TRUE),i] <- median(trainset.aftercox[,i],na.rm=T)
	}

	testset.aftercox <- cbind(testcore[,c(1,2)],testcore[,intersectname_cox])
	for (i in 1:length(colnames(testset.aftercox))) {
	  testset.aftercox[which(is.na(testset.aftercox[,i])==TRUE),i] <- median(testset.aftercox[,i],na.rm=T)
	}
	trainset.aftercox$vital_status <- trainset.aftercox$vital_status-1
	testset.aftercox$vital_status <- testset.aftercox$vital_status-1
	rf <- rfsrc( Surv(days, vital_status) ~ . , data = trainset.aftercox, ntree = 200,nsplit = 2)
	imp.rf <- sort(rf$importance) #variables importance
	rf.pred <- predict(rf,newdata = testset.aftercox)
	cindex2 <- 1-rcorr.cens(rf.pred$predicted , Surv(testset.aftercox$days, testset.aftercox$vital_status))["C Index"]
	### cindex:0.4663653
	if (cindex2 >0.5) {
	  mat.100[k,2] <- cindex2
	} else {
	  mat.100[k,2] <- 1 - cindex2
	}
	
}

colnames(mat.100) <- c('lasso','rf.surv')
write.csv(mat.100, 'mat.100.icgc.csv')

