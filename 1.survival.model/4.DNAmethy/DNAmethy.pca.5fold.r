#### Step1 Data Collection and Feature Pre-Selection

## TCGA (GDC) data collection
setwd('/data1/users/luym/project/ov_prediction/data')
total <- read.table('totalpatient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)

setwd('/data1/users/luym/project/ov_prediction/data')
DNAmethall <- read.table('ovDNAmethylation_summary_3.txt',header=T, sep='\t',quote='',row.names=1)
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
setwd('/data1/users/luym/project/ov_prediction/result/4.DNAmethy')
write.table(list3,'coxDNAmeth.txt',sep='\t')

setwd('/data1/users/luym/project/ov_prediction/result/4.DNAmethy')
list3 <-   read.table('coxDNAmethy.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}
## PCT (or PCA method)to further reduce the dimensions of the data
table.train.pca <- table3[,-c(1,2)]
table.train.pca <- as.matrix(table.train.pca)
pca <- princomp(table.train.pca,cor = TRUE)
#summary(pca)
## 1-89:80%, 1-111:85%, 1-140:90%, 1-183:95%
new.data.train <- table.train.pca %*% pca$loadings[,1:89]
new.data.train <- cbind(table3[,c(1,2)],new.data.train)

#### Step2 Model Training and Predictions

### five-fold cross validation with PCA 

## LASSO and Cox * 100 times
library("glmnet")
library(Hmisc) 
library(randomForestSRC)
mat.100 <- matrix(nrow= 100, ncol=2)
for (i in 1:20) {
  nrFolds <- 5
  folds <- sample(rep_len(1:nrFolds, nrow(new.data.train)))
  
  for (k in 1:nrFolds) {
    fold <- which(folds == k)
    data.train <- new.data.train[-fold,]
    data.test <- new.data.train[fold,]
    setwd('/data1/users/luym/project/ov_prediction/result/4.DNAmethy')
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
    
    ## random forest * 100 times    
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
write.csv(mat.100, 'mat.100.pca.5fold.csv')







