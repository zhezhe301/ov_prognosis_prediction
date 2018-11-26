#### Step1 Data Collection and Feature Pre-Selection
## TCGA (GDC) data collection
# setwd('../../')
# please change file path on your local machine
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
DNAseqall <- read.table('tcga.DNAseq.txt',header=T, sep='\t',quote='',row.names=1)
DNAseq <- DNAseqall[which(rownames(DNAseqall)%in%rownames(core)),]
vital_status <- core[rownames(DNAseq),3]
days <- core[rownames(DNAseq),4]
table <- cbind(days, vital_status, DNAseq)
table <- table[,-(which(is.na(table[1,])==TRUE))]
## Cox proportional hazards regression identifying the most significant features
library("survival")
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

## PCT (or PCA method)to further reduce the dimensions of the data
table.train.pca <- table3[,-c(1,2)]
table.train.pca <- as.matrix(table.train.pca)
pca <- princomp(table.train.pca,cor = TRUE)
#summary(pca)
## 1-57:80%, 1-70:85%, 1-85:90%, 1-109 :95%
new.data.train <- table.train.pca %*% pca$loadings[,1:57]
new.data.train <- cbind(table3[,c(1,2)],new.data.train)

#### Step2 Model Training and Predictions

### five-fold cross validation with PCA 

## LASSO and Cox * 100 times
setwd('./1.survival.model/1.DNAseq')
library("glmnet")
library(Hmisc) 
library(randomForestSRC)
mat.100 <- matrix(nrow= 100, ncol=2)
for (i in 1:20) {
  nrFolds <- 5
  folds <- sample(rep_len(1:nrFolds, nrow(new.data.train)))
  
  for (k in 1:nrFolds) {
    repeat {
      fold <- which(folds == k)
      data.train <- new.data.train[-fold,]
      data.test <- new.data.train[fold,]

      write.csv(data.train,'trainset.pca.5fold.csv')
      write.csv(data.test,'testset.pca.5fold.csv')
      
      x <- as.matrix(data.train[,-c(1,2)])
      x <- scale(x )
      y <- as.matrix(data.train[,c(1,2)])
      y[which(y[,2]==1),] <- 0
      y[which(y[,2]==2),] <- 1
      y[,1] <- data.train[,1]
      colnames(y) <- c('time', 'status')
      fit <- glmnet(x,y,family='cox')
      # plot(fit)
      # plot(fit,xvar="lambda",label=TRUE)
      fit_cv=cv.glmnet(x,y,family="cox",nfolds = 10)
      best_lambda <- fit_cv$lambda.min
      fit.coef <- coef(fit,s= fit_cv$lambda.min)
      listgene_lasso <- rownames(fit.coef)[which(fit.coef != 0)]
      if (length(listgene_lasso) >0) {
        break
      }
    }
    
    
    optimal.beta  <- fit$beta[,which(fit$lambda==best_lambda)]
    nonzero.coef <- abs(optimal.beta)>0
    selectedXname <- names(nonzero.coef)[which(nonzero.coef, arr.ind = TRUE)]
    intersectname <- intersect(selectedXname,colnames(data.test))
    selectedBeta <- optimal.beta[intersectname]
    selectedtestX <- data.test[,intersectname]
    selectedtestX <- as.matrix(selectedtestX)
    selectedtestX <- scale(selectedtestX)
    coxph.model <- coxph(Surv(data.test$days,data.test$vital_status) ~selectedtestX,init=selectedBeta,iter=0)
    predict <- survfit(coxph.model,newdata= as.data.frame(data.test[,intersectname]))
    xxx <- summary(predict)
    predictval <- xxx$table[,'median']
    cindex <- rcorr.cens(predictval, Surv(data.test$days,data.test$vital_status))["C Index"]
    if (cindex >= 0.5) {
      mat.100[5*(i-1)+k,1] <- cindex
    } else {
      mat.100[5*(i-1)+k,1] <- 1-cindex
    }
    
    
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
     
    if (cindex >= 0.5) {
      mat.100[5*(i-1)+k,2] <- cindex
    } else {
      mat.100[5*(i-1)+k,2] <- 1- cindex
    } 
  }
}
colnames(mat.100) <- c('lasso','rf.surv')
write.csv(mat.100, 'mat.100.pca.5fold.csv')





