### All of the results were generated by the code from the previous steps in directory "1.survival.model", and "2.tx.outcome.model".

## figure3 A 
## four files named "mat.100.5fold.csv" and other four files named "mat.100.pca.5fold.csv" (from DNAseq, RNAseq, miRNAseq, and DNAmethy) were binded together to generate the file "survival.5fold.txt"
library(ggplot2)
setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction')
try <- read.table('survival.5fold.txt')
colnames(try) <- c('value','omic','method')
setwd('~/Project/ov_prognosis_drugprecision/plot/c-index/')
try$omic <- factor(try$omic,levels=c('DNA-seq','RNA-seq','miRNA-seq','DNA-methy'))
setEPS()
png('survival.5fold2.10.png',width=10,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
boxplot(value~as.factor(omic), data=try[1:800,], boxwex=0.15, at=1:4-0.3, subset= method=="Lasso",
        col='#fdae6185', notch=TRUE,outline=F, las=2, main="Omics", ylab="C-index", ylim= c(0.4,0.85),xlim=c(0.5, 4.5),xaxt="n",border= '#fdae61')
boxplot(value~omic, data=try[1:800,], boxwex=0.15, at=1:4+0.12, subset= method=="Random_Forest",
        col="#80cdc185", notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border="#80cdc1")

boxplot(value~omic, data=try[801:1600,], boxwex=0.15, at=1:4-0.12, subset= method=="Lasso+PCA",
        col="#d7191c85", notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border="#d7191c")
boxplot(value~omic, data=try[801:1600,], boxwex=0.15, at=1:4+0.3, subset= method=="Random_Forest+PCA",
        col='#01857185', notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border= '#018571')
abline(h=0.5, lty=5,col='red')
dev.off()

## figure3 B
## four files named "mat.100.icgc.csv" and other four files named "mat.100.pca.icgc.csv" (from DNAseq, RNAseq, miRNAseq, and DNAmethy) were binded together to generate the file "survival.icgc.txt"
setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction')
try <- read.table('survival.icgc.txt')
colnames(try) <- c('value','omic','method')
setwd('~/Project/ov_prognosis_drugprecision/plot/c-index/')
try$omic <- factor(try$omic,levels=c('DNA-seq','RNA-seq','miRNA-seq','DNA-methy'))
setEPS()
png('survival.icgc.10.png',width=10,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
boxplot(value~omic, data=try[1:800,], boxwex=0.15, at=1:4-0.3, subset= method=="Lasso",
        col='#fdae6185', notch=TRUE,outline=F, las=2, main="Omics", ylab="C-index", ylim= c(0.4,0.85),xlim=c(0.5, 4.5),border='#fdae61',xaxt="n")
boxplot(value~omic, data=try[1:800,], boxwex=0.15, at=1:4+0.12, subset= method=="Random_Forest",
        col="#80cdc185", notch=TRUE,outline=F, las=2, add=TRUE,border="#80cdc1",xaxt="n")

boxplot(value~omic, data=try[801:1600,], boxwex=0.15, at=1:4-0.12, subset= method=="Lasso+PCA",
        col="#d7191c85", notch=TRUE,outline=F, las=2, add=TRUE,border="#d7191c",xaxt="n")
boxplot(value~omic, data=try[801:1600,], boxwex=0.15, at=1:4+0.3, subset= method=="Random_Forest+PCA",
        col='#01857185', notch=TRUE,outline=F, las=2, add=TRUE,border='#018571',xaxt="n")

abline(h=0.5, lty=5,col='red')
dev.off()

## figure3 C
## four files named "mat.100python.txoutcome.csv" (from DNAseq, RNAseq, miRNAseq, and DNAmethy) were binded together to generate the file "treatment.5fold.txt"
setwd('~/Project/ov_prognosis_drugprecision/result/4.prediction')
try <- read.table('treatment.5fold.txt')
colnames(try) <- c('value','omic','method')
setwd('~/Project/ov_prognosis_drugprecision/plot/c-index/')
try$omic <- factor(try$omic,levels=c('DNA-seq','RNA-seq','miRNA-seq','DNA-methy'))
setEPS()
png('treatment.5fold2.10.png',width=15,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
boxplot(value~omic, data=try[1:1200,], boxwex=0.1, at=1:4-0.38, subset= method=="DT",
        col='#fdae6185', notch=TRUE,outline=F, las=2, main="Omics", ylab="C-index", ylim= c(0.4,0.9),xlim=c(0.5, 4.5),xaxt="n",border='#fdae61')
boxplot(value~omic, data=try[1:1200,], boxwex=0.1, at=1:4-0.08, subset= method=="Random_Forest",
        col="#80cdc185", notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border="#80cdc1" )
boxplot(value~omic, data=try[1:1200,], boxwex=0.1, at=1:4+0.23, subset= method=="DNN",
        col="#92c5de85", notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border="#92c5de" )
boxplot(value~omic, data=try[1201:2400,], boxwex=0.1, at=1:4-0.23, subset= method=="DT+PCA",
        col="#d7191c85", notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border="#d7191c")
boxplot(value~omic, data=try[1201:2400,], boxwex=0.1, at=1:4+0.08, subset= method=="Random_Forest+PCA",
        col='#01857185', notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border='#018571' )
boxplot(value~omic, data=try[1201:2400,], boxwex=0.1, at=1:4+0.38, subset= method=="DNN+PCA",
        col='#0571b085', notch=TRUE,outline=F, las=2, add=TRUE,xaxt="n",border='#0571b0')
abline(h=0.5, lty=5,col='red')
dev.off()
# legend( "topright",  inset=c(-0.2,0), title="Methods",
#         c("DT","DT+PCA","Random_Forest",'Random_Forest+PCA','DNN','DNN+PCA'),
#         fill=c('#fdae61',"#d7191c","#80cdc1",'#018571',"#92c5de",'#0571b0'))

