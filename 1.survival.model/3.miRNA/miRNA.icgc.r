#### Step1 Data Collection and Feature Pre-Selection

## TCGA (GDC) data collection
# setwd('../../')
# please change file path on your local machine
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
miRNAseqall <- read.csv('tcga.miRNAseq.csv',row.names = 1)
miRNAseq <- miRNAseqall[which(rownames(miRNAseqall)%in%rownames(core)),]
vital_status <- core[rownames(miRNAseq),3]
days <- core[rownames(miRNAseq),4]
##### 1==alive, 2==dead #####
table <- cbind(days, vital_status, miRNAseq)

## Cox proportional hazards regression identifying the most significant features
library("survival")
attach(table)
table$days <- as.numeric(table$days)
table[,2] <- as.numeric(table[,2])
list1 <- c()
for (i in colnames(table)[3:1870]) {
  z <- Surv(table$days,table$vital_status==2)
  y <- coxph(z~table[,i])
  y.s <- summary(y)
  list1 <- append(list1,y.s$logtest['pvalue'])
}
names(list1) <- colnames(table)[3:1870]
list2 <- sort(list1)[1:163] ####p<0.05
## events :300
list3 <- sort(list1)[1:163]
write.table(list3,'coxmiRNAseq.txt',sep='\t')
## ICGC data collection
list3 <-   read.table('coxmiRNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}




## intersect between TCGA and ICGC
testset <- read.table('icgc.miRNAseq.txt',sep = '\t',quote='',header = TRUE,row.names = 1)
surv.info <- read.table('icgc.patient.txt',header=T,sep='\t',row.names = 1)
surv.info <- as.data.frame(surv.info)
vital_status <- surv.info[rownames(testset),1]
days <- surv.info[rownames(testset),6]
testcore <- cbind(days,vital_status,testset)
for (i in 1:length(colnames(testcore))) {
  testcore[which(is.na(testcore[,i])==TRUE),i] <- median(testcore[,i],na.rm=T)
}

table3[3:length(colnames(table3))] <- log(table3[3:length(colnames(table3))]+1)
testcore[3:length(colnames(testcore))] <- log(testcore[3:length(colnames(testcore))] +1)



###### write csv to save after cox ######
intersectname_cox <- intersect(rownames(list3),colnames(testset))
trainset.aftercox <- cbind(table[,c(1,2)],table[,intersectname_cox])
for (i in 1:length(colnames(trainset.aftercox))) {
  trainset.aftercox[which(is.na(trainset.aftercox[,i])==TRUE),i] <- median(trainset.aftercox[,i],na.rm=T)
}

testset.aftercox <- cbind(testcore[,c(1,2)],testcore[,intersectname_cox])
for (i in 1:length(colnames(testset.aftercox))) {
  testset.aftercox[which(is.na(testset.aftercox[,i])==TRUE),i] <- median(testset.aftercox[,i],na.rm=T)
}
## random forest * 100 times
setwd('./1.survival.model/3.miRNA') 
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
	x <- as.matrix(table3[,-c(1,2)])
	y <- as.matrix(table3[,c(1,2)])
	y[which(y[,2]==1),] <- 0
	y[which(y[,2]==2),] <- 1
	y[,1] <- table3[,1]
	colnames(y) <- c('time', 'status')
	fit <- glmnet(x,y,family='cox')
	plot(fit)
	plot(fit,xvar="lambda",label=TRUE)
	fit_cv=cv.glmnet(x,y,family="cox",nfolds = 10)
	best_lambda <- fit_cv$lambda.min
	fit.coef <- coef(fit,s= fit_cv$lambda.min)
	listgene_lasso <- rownames(fit.coef)[which(fit.coef != 0)]
	max.dev.index     <- which.max(fit$dev.ratio)
	optimal.lambda <- fit$lambda[max.dev.index]
	optimal.beta  <- fit$beta[,max.dev.index]
	nonzero.coef <- abs(optimal.beta)>0
	selectedXname <- names(nonzero.coef)[which(nonzero.coef, arr.ind = TRUE)]
	intersectname <- intersect(selectedXname,colnames(testcore))
	selectedtestX <- testcore[,intersectname]
	selectedtestX <- as.matrix(selectedtestX)
	selectedBeta <- optimal.beta[intersectname]
	coxph.model <- coxph(Surv(testcore$days,testcore$vital_status) ~selectedtestX,init=selectedBeta,iter=0) 
	predict <- survfit(coxph.model,newdata= testcore[,intersectname])
	xxx <- summary(predict)
	predictval <- xxx$table[,'median']
	rcorr.cens(predictval, Surv(testcore$days, testcore$vital_status))["C Index"]
	cindex <- 1-rcorr.cens(predict(coxph.model,newdata=testcore[,intersectname],type="risk"), Surv(testcore$days, testcore$vital_status))["C Index"]
	### cindex :0.5062483 
	mat.100[k,1] <- cindex

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
	cindex <- rcorr.cens(rf.pred$predicted , Surv(testset.aftercox$days, testset.aftercox$vital_status))["C Index"]
	### cindex:0.4663653 
	mat.100[k,2] <- cindex
}

colnames(mat.100) <- c('lasso','rf.surv')
write.csv(mat.100, 'mat.100.icgc.csv')





