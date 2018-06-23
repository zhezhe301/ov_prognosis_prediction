#### Step1 Data Collection and Feature Pre-Selection

## TCGA (GDC) data collection
setwd('/data1/users/luym/project/ov_prediction/data')
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)

setwd('/data1/users/luym/project/ov_prediction/data')
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
setwd('/data1/users/luym/project/ov_prediction/result/3.miRNAseq')
write.table(list3,'coxmiRNAseq.txt',sep='\t')

setwd('/data1/users/luym/project/ov_prediction/result/3.miRNAseq')
list3 <-   read.table('coxmiRNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}

table3[3:length(colnames(table3))] <- log(table3[3:length(colnames(table3))]+1)
#### Step2 Model Training and Predictions

### five-fold cross validation without PCA 

## LASSO and Cox * 100 times
library("glmnet")
library(Hmisc) 
library(randomForestSRC)
mat.100 <- matrix(nrow= 100, ncol=2)
for (i in 1:20) {
  nrFolds <- 5
  folds <- sample(rep_len(1:nrFolds, nrow(table3)))
  
  for (k in 1:nrFolds) {
    fold <- which(folds == k)
    data.train <- table3[-fold,]
    data.test <- table3[fold,]
    setwd('/data1/users/luym/project/ov_prediction/result/3.miRNA')
    write.csv(data.train,'trainset.5fold.csv')
    write.csv(data.test,'testset.5fold.csv')
    
    x <- as.matrix(data.train[,-c(1,2)])
    y <- as.matrix(data.train[,c(1,2)])
    y[which(y[,2]==1),] <- 0
    y[which(y[,2]==2),] <- 1
    y[,1] <- data.train[,1]
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
    selectedtestX <- data.test[,selectedXname]
    selectedtestX <- as.matrix(selectedtestX)
    selectedBeta <- optimal.beta[selectedXname]
    coxph.model <- coxph(Surv(data.test$days,data.test$vital_status) ~selectedtestX,init=selectedBeta,iter=0) 
    predict <- survfit(coxph.model,newdata= data.test[,selectedXname])
    xxx <- summary(predict)
    predictval <- xxx$table[,'median']
    rcorr.cens(predictval, Surv(data.test$days,data.test$vital_status))["C Index"]
    cindex <- 1-rcorr.cens(predict(coxph.model,newdata=data.test[,selectedXname],type="risk"), Surv(data.test$days,data.test$vital_status))["C Index"]
    mat.100[5*(i-1)+k,1] <- cindex
    ###
    
    
    trainset.aftercox <- data.train
    testset.aftercox <- data.test
    trainset.aftercox$vital_status <- trainset.aftercox$vital_status-1
    testset.aftercox$vital_status <- testset.aftercox$vital_status-1
    rf <- rfsrc( Surv(days, vital_status) ~ . , data = trainset.aftercox, ntree = 200,nsplit = 2)
    rf$xvar.names # picked variables
    imp.rf <- sort(rf$importance) #variables importance
    rf$predicted.oob
    rf.pred <- predict(rf,newdata = testset.aftercox)

    cindex <- rcorr.cens(-rf.pred$predicted , Surv(testset.aftercox$days, testset.aftercox$vital_status))["C Index"]
    mat.100[5*(i-1)+k,2] <- cindex 
  }
}
colnames(mat.100) <- c('lasso','rf.surv')
write.csv(mat.100, 'mat.100.5fold.csv')
## random forest * 100 times
library("glmnet")
library(Hmisc) 
library(randomForestSRC)
mat.100 <- matrix(nrow= 100, ncol=2)
for (i in 1:20) {
  nrFolds <- 5
  folds <- sample(rep_len(1:nrFolds, nrow(table3)))
  
  for (k in 1:nrFolds) {
    fold <- which(folds == k)
    data.train <- table3[-fold,]
    data.test <- table3[fold,]
    setwd('/data1/users/luym/project/ov_prediction/result/3.miRNAseq')
    write.csv(data.train,'trainset.pca.5fold.csv')
    write.csv(data.test,'testset.pca.5fold.csv')
    
    x <- as.matrix(data.train[,-c(1,2)])
    y <- as.matrix(data.train[,c(1,2)])
    y[which(y[,2]==1),] <- 0
    y[which(y[,2]==2),] <- 1
    y[,1] <- data.train[,1]
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
    selectedtestX <- data.test[,selectedXname]
    selectedtestX <- as.matrix(selectedtestX)
    selectedBeta <- optimal.beta[selectedXname]
    coxph.model <- coxph(Surv(data.test$days,data.test$vital_status) ~selectedtestX,init=selectedBeta,iter=0) 
    predict <- survfit(coxph.model,newdata= data.test[,selectedXname])
    xxx <- summary(predict)
    predictval <- xxx$table[,'median']
    rcorr.cens(predictval, Surv(data.test$days,data.test$vital_status))["C Index"]
    cindex <- 1-rcorr.cens(predict(coxph.model,newdata=data.test[,selectedXname],type="risk"), Surv(data.test$days,data.test$vital_status))["C Index"]
    mat.100[5*(i-1)+k,1] <- cindex
    ###
    
    
    trainset.aftercox <- data.train
    for (j in 1:length(colnames(trainset.aftercox))) {
      trainset.aftercox[which(is.na(trainset.aftercox[,j])==TRUE),j] <- min(trainset.aftercox[,j],na.rm=TRUE)
    }
    testset.aftercox <- data.test
    for (j in 1:length(colnames(testset.aftercox))) {
      testset.aftercox[which(is.na(testset.aftercox[,j])==TRUE),j] <- min(testset.aftercox[,j],na.rm=TRUE)
    }
    trainset.aftercox$vital_status <- trainset.aftercox$vital_status-1
    testset.aftercox$vital_status <- testset.aftercox$vital_status-1
    rf <- rfsrc( Surv(days, vital_status) ~ . , data = trainset.aftercox, ntree = 200,nsplit = 2)
    imp.rf <- sort(rf$importance) #variables importance
    rf.pred <- predict(rf,newdata = testset.aftercox)
    

    cindex <- rcorr.cens(-rf.pred$predicted , Surv(testset.aftercox$days, testset.aftercox$vital_status))["C Index"]
    mat.100[5*(i-1)+k,2] <- cindex 
  }
}
colnames(mat.100) <- c('lasso','rf.surv')
write.csv(mat.100, 'mat.100.5fold.csv')











