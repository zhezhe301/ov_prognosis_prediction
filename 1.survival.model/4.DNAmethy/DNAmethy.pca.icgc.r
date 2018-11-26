#### Step1 Data Collection and Feature Pre-Selection

## TCGA (GDC) data collection
# setwd('../../')
# please change file path on your local machine
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
DNAmethall <- read.table('tcga.DNAmethylation.txt',header=T, sep='\t',quote='',row.names=1)
DNAmeth <- DNAmethall[which(rownames(DNAmethall)%in%rownames(core)),]
vital_status <- core[rownames(DNAmeth),3]
days <- core[rownames(DNAmeth),4]
##### 1==alive, 2==dead #####
table <- cbind(days, vital_status, DNAmeth)
table <- table[,-(which(is.na(table[1,])==TRUE))]

## Cox proportional hazards regression identifying the most significant features
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
## events :344
list3 <- sort(list1)[1:344]
write.table(list3,'coxDNAmeth.txt',sep='\t')
list3 <-   read.table('coxDNAmeth.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}

## ICGC data collection
testset <- read.table('icgc.DNAmethylation.txt',sep='\t',header=T,row.names = 1)
surv.info <- read.table('icgc.patient.txt',header=T,sep='\t',row.names = 1)
surv.info <- as.data.frame(surv.info)
vital_status <- surv.info[rownames(testset),1]
days <- surv.info[rownames(testset),6]
testcore <- cbind(days,vital_status,testset)
for (i in 1:length(colnames(testcore))) {
  testcore[which(is.na(testcore[,i])==TRUE),i] <- median(testcore[,i],na.rm=T)
}

## intersect between TCGA and ICGC
int.name <- intersect(colnames(testcore), rownames(list3))
table3 <- cbind(table3[,c(1,2)],table3[,int.name])
testcore <- cbind(testcore[,c(1,2)],testcore[,int.name])
## PCT (or PCA method)to further reduce the dimensions of the data
table.train.pca <- table3[,-c(1,2)]
table.train.pca <- as.matrix(table.train.pca)
pca <- princomp(table.train.pca,cor = TRUE)
#summary(pca)
## 1-89:80%, 1-111:85%, 1-140:90%, 1-183:95%
new.data.train <- table.train.pca %*% pca$loadings[,1:89]
new.data.train <- cbind(table3[,c(1,2)],new.data.train)

table.test.pca <- testcore[,-c(1,2)]
table.test.pca <- as.matrix(table.test.pca)
new.data.test <- table.test.pca %*% pca$loadings[,1:89]
new.data.test <- cbind(testcore[,c(1,2)],new.data.test)


###### write csv to save after cox ######
trainset.aftercox <- new.data.train
testset.aftercox <- new.data.test
setwd('./1.survival.model/4.DNAmethy') 
write.csv(trainset.aftercox,'trainset.pca.icgc.csv')
write.csv(testset.aftercox,'testset.pca.icgc.csv')



#### Step2 Model Training and Predictions

### ICGC independent validation with PCA 

## LASSO and Cox * 100 times
mat.100 <- matrix(nrow=100, ncol=2)
library(Hmisc) 
library(randomForestSRC)
library("glmnet")
for (k in 1:100) {
  x <- as.matrix(new.data.train[,-c(1,2)])
  y <- as.matrix(new.data.train[,c(1,2)])
  y[which(y[,2]==1),] <- 0
  y[which(y[,2]==2),] <- 1
  y[,1] <- new.data.train[,1]
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
  intersectname <- intersect(selectedXname,colnames(new.data.test))
  selectedtestX <- new.data.test[,intersectname]
  selectedtestX <- as.matrix(selectedtestX)
  selectedBeta <- optimal.beta[intersectname]
  coxph.model <- coxph(Surv(new.data.test$days,new.data.test$vital_status) ~selectedtestX,init=selectedBeta,iter=0) 
  predict <- survfit(coxph.model,newdata= new.data.test[,intersectname])
  xxx <- summary(predict)
  predictval <- xxx$table[,'median']
  cindex1 <- 1-rcorr.cens(predict(coxph.model,newdata=new.data.test[,intersectname],type="risk"), Surv(new.data.test$days, new.data.test$vital_status))["C Index"]
  mat.100[k,1] <- cindex1
  
  
  
  ## random forest * 100 times 
  
  trainset.aftercox <- new.data.train
  testset.aftercox <- new.data.test
  trainset.aftercox$vital_status <- trainset.aftercox$vital_status-1
  testset.aftercox$vital_status <- testset.aftercox$vital_status-1
  rf <- rfsrc( Surv(days, vital_status) ~ . , data = trainset.aftercox, ntree = 200,nsplit = 2)
  imp.rf <- sort(rf$importance) #variables importance
  rf.pred <- predict(rf,newdata = testset.aftercox)
  
  cindex2<- rcorr.cens(-rf.pred$predicted , Surv(testset.aftercox$days, testset.aftercox$vital_status))["C Index"]
  mat.100[k,2]<- cindex2
  
}
colnames(mat.100) <- c('lasso','rf.surv')
write.csv(mat.100, 'mat.100.pca.icgc.csv')











