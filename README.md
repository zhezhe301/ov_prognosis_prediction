# OV_Prognosis_Prediction

Assessing the clinical utility of multi-omics data for predicting serous ovarian cancer prognosis.


## Prerequisites

Basically, we use R and Python programming language to write these scripts. Required packages are free to get, and specificly, TensorFlow should follow the instruction on https://www.tensorflow.org/tutorials/, and circos plot should follow the tutorials on http://circos.ca/documentation/tutorials/ to install them.

## Getting Started

These instructions will get you how to run on your local machine for development and testing purposes. There are four parts of our codes.

## 1. Part1 survival.model

This part provided prognosis prediction for survival model. There are four folders (*1.DNAseq*, *2.RNAseq*, *3.miRNA*, *4.DNAmethy*), and in each of them, there are four files written by R scripts. *XX.5fold.r* and *XX.pca.5fold.r* are used for creating files for figure3A, while *XX.icgc.r* and *XX.pca.icgc.r* are used for creating files for figure3B. We will take *DNAseq.5fold.r* as an example for step-by-step instruction on running the code.

Make sure you have downloaded the following files: *tcga.patient.txt*, *tcga.DNAseq.txt*, and you have installed following R packages: *survival*, *glmnet*, *Hmisc*, and *randomForestSRC*. Remember change file path on your local machine.

Step1. Data collection and feature pre-selection

```
# TCGA (GDC) data collection
setwd('/data1/users/luym/project/ov_prediction/data')
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
setwd('/data1/users/luym/project/ov_prediction/data')
DNAseqall <- read.table('tcga.DNAseq.txt',header=T, sep='\t',quote='',row.names=1)
DNAseq <- DNAseqall[which(rownames(DNAseqall)%in%rownames(core)),]
vital_status <- core[rownames(DNAseq),3]
days <- core[rownames(DNAseq),4]##### 1==alive, 2==dead #####
table <- cbind(days, vital_status, DNAseq)
table <- table[,-(which(is.na(table[1,])==TRUE))]
# Cox proportional hazards regression identifying the most significant features
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
setwd('/data1/users/luym/project/ov_prediction/result/1.DNAseq')
write.table(list3,'coxDNAseq.txt',sep='\t')
setwd('/data1/users/luym/project/ov_prediction/result/1.DNAseq')
list3 <- read.table('coxDNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=TRUE)
}
```

Step2. Model Training and Predictions

```
# five-fold cross validation without PCA 
# LASSO and Cox * 100 times
library("glmnet")
library(Hmisc) 
library(randomForestSRC)
mat.100 <- matrix(nrow= 100, ncol=2)
for (i in 1:20) {
  nrFolds <- 5
  folds <- sample(rep_len(1:nrFolds, nrow(table3)))
  
  for (k in 1:nrFolds) {
    repeat{
      fold <- which(folds == k)
      data.train <- table3[-fold,]
      data.test <- table3[fold,]
      # setwd('/data1/users/luym/project/ov_prediction/result/1.DNAseq')
      setwd('/Users//zhangzhe/Project/ov_prognosis_drugprecision/result/repeat/1.DNAseq')
      
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
      if (length(listgene_lasso) >0) {
        break
      }
    }
    max.dev.index <- which.max(fit$dev.ratio)
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
    if (cindex >= 0.5) {
      mat.100[5*(i-1)+k,1] <- cindex
    } else {
      mat.100[5*(i-1)+k,1] <- 1-cindex
    }    
# random forest * 100 times    
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
      mat.100[5*(i-1)+k,2] <- 1-cindex
    }
  }
}
colnames(mat.100) <- c('lasso','rf.surv')
```

Step3. Save the results of calculated C-indexs

```
write.csv(mat.100, 'mat.100.5fold.csv')
```

## 2. Part2 tx.outcome.model

This part provided prognosis prediction for treatment outcome model. There are four folders (*1.DNAseq*, *2.RNAseq*, *3.miRNA*, *4.DNAmethy*), and in each of them, there are three files written by R and Python. The three scripts are used for creating files for figure3C. We will take *1.DNAseq* as an example for step-by-step instruction on running the code

The steps will be

```
Step1. Run R script *1.DNAseq.drug.r* for data collection and feature pre-selection. Create related files and save them
Step2. Run Python scripts *DNAseq.5fold.py* and *DNAseq.pca.5fold.py*. This step will create predicted results and save them
Step3. Back to run R script *1.DNAseq.drug.r* for C-index calculation and save the results.
```

## 3. Part3 pathway

This is pathway analyses. There is a single file written by R. Before this part, make sure you have the following files: *panther.slim.go.ids.txt*, *coxRNAseq.txt*, *RNASeq_pvalue.csv*, *coxDNAseq.txt*, *DNAseq.pvalue.csv*, *coxDNAmethy.txt*, and *DNAmethy_pvalue.csv*

## 4. Part4 plot

This part provided figure2-figure6 plotting code

### 4.1 

The R script *figure2.r* is for figure2. Before this, make sure you have the following files: *tcga.patient.txt*, *tcga.RNAseq.csv*, *icgc.RNAseq.csv*, *icgc.patient.txt*, *tgca.DNAmethylation.txt*, and *icgc.DNAmethylation.txt*

### 4.2 

The R script *figure3.r* is for figure3 which shows C-index of all the models. Before this step, you need first run all the scripts from part1 (*1.survival.model*) and part2 (*2.tx.outcome.model*) to create related files. Next, bind them together to create the following files: *survival.5fold.txt*, *survival.icgc.txt*, and *treatment.5fold.txt*. You can bind them manually by yourself. However, we have already uploaded this three files.
In detail,

4.2.1 The file *survival.5fold.txt* is used for figure 3A plotting. After running the scripts from part1 (*1.survival.model*), you will have four files named "mat.100.5fold.csv" and other four files named "mat.100.pca.5fold.csv" (from *1.DNAseq*, *2.RNAseq*, *3.miRNA*, *4.DNAmethy*), bind them together to generate the file *survival.5fold.txt*.

4.2.2 The file *survival.icgc.txt* is used for figure 3B plotting. After running the scripts from part1 (*1.survival.model*), you will have four files named "mat.100.icgc.csv" and other four files named "mat.100.pca.icgc.csv" (from *1.DNAseq*, *2.RNAseq*, *3.miRNA*, *4.DNAmethy*), bind them together to generate the file *survival.icgc.txt*.

4.2.3 The file *treatment.5fold.txt* is used for figure 3C plotting. After running the scripts from part2 (*2.tx.outcome.model*), you will have four files named "mat.100python.txoutcome.csv" (from *1.DNAseq*, *2.RNAseq*, *3.miRNA*, *4.DNAmethy*), bind them together to generate the file treatment.5fold.txt.

### 4.3 

The R script *figure4.r* is to plot figure4. Before this, make sure you have the following files: *tcga.patient.txt*, *tcga.RNAseq.csv*, *RNASeq_pvalue.csv*, *icgc.RNAseq.csv*, and *icgc.patient.txt*

### 4.4 

The R script *figure5.r* is for figure5. Just run it.

### 4.5 

The directory *figure6* including scripts and data for generating figure6 circos plot.
The circos plot have tutorials on http://circos.ca/documentation/tutorials/. Before the figure6 plot, we recommend you have read the tutorials and installed circos. For UNIX users, just run corresponding script named *circos.pathway.conf* (for Windows users, see in their tutorials)

```
../bin/circos -conf figure6/figure6E.dnamethy.sv.config/circos.pathway.conf 
```

The files you need to use are already provided in the corresponding folder and remember to change the file path on your local machine. The file *karyotype.dnamethy.txt* defines the axes, these are typically outer ring of BIOLOGICAL_PROCESS (bp), CELLULAR_COMPONENT (cc), MOLECULAR_FUNCTION (mf), and their size. The file *highlight.2level.txt* defines the position and color of parent terms of the significant GO terms (the inner ring), their labels are described by the file *labels.2level.txt*. The file *highlight.3level.txt* defines the position and color of the significant GO terms (the third ring), their labels are described by the file *labels.3level.txt*, and their average fold changes are described by the file *histogram.3level.txt*. The files *links.txt*, *links1.txt*, *links2.txt*, *links3.txt*, and *links4.txt* define the overlapping genes of significant pathways (the bezier curves). 

You can use the files above to generate figure6 directly, or use R scripts *circos.DNAseq.r*, *circos.RNAseq.r*, and *circos.DNAmethy.r* to recreate these files by yourself(the required files saved in folder named *data*).  You may notice that your created files are a little different to our provided files, that is because R package *org.Hs.eg.db* has been updated recently. 

