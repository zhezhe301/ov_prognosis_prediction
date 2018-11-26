# figure 4A
setwd('../')
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
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
write.table(list3,'coxRNAseq.txt',sep='\t')
list3 <-   read.table('coxRNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}
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
library("survival")
library("glmnet")
set.seed(1000)
nrFolds <- 5
folds <- sample(rep_len(1:nrFolds, nrow(trainset.aftercox)))
table.curve.all <- NA
for (k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- table3[-fold,]
  data.test <- table3[fold,]
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
  
  table.curve <- data.test[1:2]
  table.curve[,3] <- predict(coxph.model,newdata=data.test[,selectedXname],type="risk")
  table.curve.all <- rbind(table.curve,table.curve.all)
}
table.curve.all <- na.omit(table.curve.all)
for (i in 1:nrow(table.curve.all)) {
  if (table.curve.all[i,3]>=median(table.curve.all[,3])) {
    table.curve.all[i,4] <- 1
  } else {
    table.curve.all[i,4] <- 2
  }
}
setEPS()
postscript('figure4A.eps', width=9, height=6)
table.curve.all[,4] <- as.factor(table.curve.all[,4])
kmfit <- survfit(Surv(days,vital_status)~table.curve.all[,4], data=table.curve.all, conf.type = "log-log")
plot(kmfit, col=c('palevioletred','skyblue'))
survdiff(Surv(days,vital_status)~table.curve.all[,4], data=table.curve.all)
dev.off()

## figure 4B
set.seed(1000)
nrFolds <- 5
folds <- sample(rep_len(1:nrFolds, nrow(trainset.aftercox)))
table.curve.all <- NA
for (k in 1:nrFolds) {
  fold <- which(folds == k)
  data.train <- trainset.aftercox[-fold,]
  data.test <- trainset.aftercox[fold,]
  rf <- rfsrc( Surv(days, vital_status) ~ . , data = data.train, ntree = 200,nsplit = 2)
  rf.pred <- predict(rf,newdata = data.test)
  table.curve <- data.test[1:2]
  table.curve[,3] <- rf.pred$predicted
  table.curve.all <- rbind(table.curve,table.curve.all)
}
table.curve.all <- na.omit(table.curve.all)
for (i in 1:nrow(table.curve.all)) {
  if (table.curve.all[i,3]>=median(table.curve.all[,3])) {
    table.curve.all[i,4] <- 1
  } else {
    table.curve.all[i,4] <- 2
  }
}
setEPS()
postscript('figure4B.eps', width=9, height=6)
table.curve.all[,4] <- as.factor(table.curve.all[,4])
kmfit <- survfit(Surv(days,vital_status)~table.curve.all[,4], data=table.curve.all, conf.type = "log-log")
plot(kmfit, col=c('palevioletred','skyblue'),lwd=3)
legend(4000, .9, c("                    ","                    "), lty=c(1,1),lwd=2,
       col=c('palevioletred','skyblue'))
survdiff(Surv(days,vital_status)~table.curve.all[,4], data=table.curve.all)
dev.off()  

## figure 4C
library('plotly')
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- core[complete.cases(core[,21]),]
core <- core[complete.cases(core[,22]),]
core <- as.data.frame(core)
list1 <- read.csv('RNASeq_pvalue.csv',row.names = 1)
list1<- as.matrix(list1)
list2 <- sort(list1[,1])
list3 <- list2[1:1282]  ###<0.05
list3 <- list3[1:184]
RNAseqall <- read.csv('tcga.RNAseq.csv',row.names = 1)
RNAseq <- RNAseqall[which(rownames(RNAseqall)%in%rownames(core)), intersect(names(list3),colnames(RNAseqall))]
txoutcome <- core[rownames(RNAseq),23]
table <- cbind(txoutcome, RNAseq)
for (i in 1:length(colnames(table))) {
  table[which(is.na(table[,i])==TRUE),i] <- median(table[,i],na.rm=T)
}
table[,2:length(colnames(table))] <- log2(table[,2:length(colnames(table))] +1) 
table1 <- table[,-1]
table1 <- as.matrix(table1)
pca <- princomp(table1,cor = TRUE)
pca.table <- table1 %*% pca$loadings[,1:54]
pca.table <- cbind(txoutcome, pca.table)
pca.table <- as.data.frame(pca.table)
plot_ly(pca.table, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3,color = ~txoutcome, colors = c('#BF382A', '#0C4B8E'))%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

## figure 4D
total <- read.table('tcga.patient.txt',sep='\t',header=T,row.names=1) # 587 person
core <- total
core <- as.data.frame(core)
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
write.table(list3,'coxRNAseq.txt',sep='\t')
list3 <-   read.table('coxRNAseq.txt',sep = '\t')
table2 <- table[,rownames(list3)]
table2 <- cbind(table[,c(1,2)],table2)
table3 <- table2
for (i in 1:length(colnames(table3))) {
  table3[which(is.na(table3[,i])==TRUE),i] <- median(table3[,i],na.rm=T)
}
testset <- read.csv('icgc.RNAseq.csv')
rownames(testset) <- testset[,2]
testset <- testset[,-c(1,2)]
surv.info <- read.table('icgc.patient.txt',header=T,sep='\t',row.names = 1)
surv.info <- as.data.frame(surv.info)
vital_status <- surv.info[rownames(testset),1]
days <- surv.info[rownames(testset),6]
testcore <- cbind(days,vital_status,testset)
int.name <- intersect(colnames(testcore), rownames(list3))
table3 <- cbind(table3[,c(1,2)],table3[,int.name])
testcore <- cbind(testcore[,c(1,2)],testcore[,int.name])
table.train.pca <- table3[,-c(1,2)]
table.train.pca <- as.matrix(table.train.pca)
pca <- princomp(table.train.pca,cor = TRUE)
new.data.train <- table.train.pca %*% pca$loadings[,1:66]
new.data.train <- cbind(table3[,c(1,2)],new.data.train)
new.data.train <- new.data.train[-which(rownames(new.data.train)=='TCGA-25-2042'),]
new.data.train <- new.data.train[-which(rownames(new.data.train)=='TCGA-61-2102'),]
plot_ly(new.data.train, x = ~Comp.1, y = ~Comp.2, z = ~Comp.3,
        marker = list(color = ~days, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))%>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))