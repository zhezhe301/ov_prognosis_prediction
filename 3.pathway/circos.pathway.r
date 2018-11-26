#### Pathway analyses
library('org.Hs.eg.db')
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
eg2symbol <- as.list(x[mapped_genes])
xx <- as.list(org.Hs.egGO2EG)
xxx <- as.list(org.Hs.egGO2ALLEGS)
setwd('../')
table.go2gene <- read.table('panther.slim.go.ids.txt',stringsAsFactors = FALSE)
allgenome <- as.character(eg2symbol[as.data.frame(xx[table.go2gene[1,1]],stringsAsFactors= F)[,1]])
for (i in 1:length(table.go2gene[,1])) {
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    allgenome <- unique(union(allgenome, as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]) )) 
  }
}  
### RNA-seq
## survival model
total <- read.table('coxRNAseq.txt',header=T,row.names = 1)
p2.1 <- row.names(total)[which(total[,1]<0.05)]
mapgo2gene <- matrix(NA, ncol=6,nrow=length(table.go2gene[,1]))
rownames(mapgo2gene) <- table.go2gene[,1]
for (i in 1:length(table.go2gene[,1])) {
  b <- length(intersect(p2.1, allgenome)) 
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    mapgo2gene[i,1] <- length(intersect(p2.1,unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]))  ))
    mapgo2gene[i,2] <- b - length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) ))
    mapgo2gene[i,3] <- length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    mapgo2gene[i,4] <- length(allgenome) - length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    job <- matrix(c(mapgo2gene[i,1], mapgo2gene[i,2], mapgo2gene[i,3], mapgo2gene[i,4]),nrow=2)
    try(mapgo2gene[i,5] <- fisher.test(job)$p.value, silent=T) 
  }
}
mapgo2gene[,6] <- p.adjust(mapgo2gene[,5],method='fdr',n=length(mapgo2gene[,6]))
mapgo2gene2 <- mapgo2gene[order(mapgo2gene[,6]),]
colnames(mapgo2gene2) <- c('GeneinGenome','GenenotinGO','GennomeinGO','GenomenotinGO','Fisher','FDR')
## save the results (we have saved them in './4.plot/figure6/data')
write.csv(mapgo2gene2, '1RNAseq.survival.csv')

## treatment outcome model
total <- read.csv('RNASeq_pvalue.csv',row.names = 1)
p2.1 <- row.names(total)[which(total[,1]<0.05)]
mapgo2gene <- matrix(NA, ncol=6,nrow=length(table.go2gene[,1]))
# mapgo2gene <- as.data.frame(mapgo2gene)
rownames(mapgo2gene) <- table.go2gene[,1]
for (i in 1:length(table.go2gene[,1])) {
  b <- length(intersect(p2.1, allgenome)) 
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    mapgo2gene[i,1] <- length(intersect(p2.1,unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]))  ))
    mapgo2gene[i,2] <- b - length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) ))
    mapgo2gene[i,3] <- length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    mapgo2gene[i,4] <- length(allgenome) - length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    job <- matrix(c(mapgo2gene[i,1], mapgo2gene[i,2], mapgo2gene[i,3], mapgo2gene[i,4]),nrow=2)
    try(mapgo2gene[i,5] <- fisher.test(job)$p.value, silent=T) 
  }
}
mapgo2gene[,6] <- p.adjust(mapgo2gene[,5],method='fdr',n=length(mapgo2gene[,6]))
mapgo2gene2 <- mapgo2gene[order(mapgo2gene[,6]),]
colnames(mapgo2gene2) <- c('GeneinGenome','GenenotinGO','GennomeinGO','GenomenotinGO','Fisher','FDR')
## save the results (we have saved them in './4.plot/figure6/data')
write.csv(mapgo2gene2, '1RNAseq.tx.csv')

### DNAseq 
## survival model
total <- read.table('coxDNAseq.txt',header=T,row.names = 1)
p2.1 <- row.names(total)[which(total[,1]<0.05)]
mapgo2gene <- matrix(NA, ncol=6,nrow=length(table.go2gene[,1]))
rownames(mapgo2gene) <- table.go2gene[,1]
for (i in 1:length(table.go2gene[,1])) {
  b <- length(intersect(p2.1, allgenome)) 
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    mapgo2gene[i,1] <- length(intersect(p2.1,unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]))  ))
    mapgo2gene[i,2] <- b - length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) ))
    mapgo2gene[i,3] <- length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    mapgo2gene[i,4] <- length(allgenome) - length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    job <- matrix(c(mapgo2gene[i,1], mapgo2gene[i,2], mapgo2gene[i,3], mapgo2gene[i,4]),nrow=2)
    try(mapgo2gene[i,5] <- fisher.test(job)$p.value, silent=T) 
  }
}
mapgo2gene[,6] <- p.adjust(mapgo2gene[,5],method='fdr',n=length(mapgo2gene[,6]))
mapgo2gene2 <- mapgo2gene[order(mapgo2gene[,6]),]
colnames(mapgo2gene2) <- c('GeneinGenome','GenenotinGO','GennomeinGO','GenomenotinGO','Fisher','FDR')
## save the results (we have saved them in './4.plot/figure6/data')
write.csv(mapgo2gene2, '1DNAseq.survival.csv')

## treatment outcome model
total <- read.csv('DNAseq.pvalue.csv',row.names = 1)
p2.1 <- row.names(total)[which(total[,11]<0.05)]
mapgo2gene <- matrix(NA, ncol=6,nrow=length(table.go2gene[,1]))
rownames(mapgo2gene) <- table.go2gene[,1]
for (i in 1:length(table.go2gene[,1])) {
  b <- length(intersect(p2.1, allgenome)) 
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    mapgo2gene[i,1] <- length(intersect(p2.1,unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]))  ))
    mapgo2gene[i,2] <- b - length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) ))
    mapgo2gene[i,3] <- length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    mapgo2gene[i,4] <- length(allgenome) - length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    job <- matrix(c(mapgo2gene[i,1], mapgo2gene[i,2], mapgo2gene[i,3], mapgo2gene[i,4]),nrow=2)
    try(mapgo2gene[i,5] <- fisher.test(job)$p.value, silent=T) 
  }
}
mapgo2gene[,6] <- p.adjust(mapgo2gene[,5],method='fdr',n=length(mapgo2gene[,6]))
mapgo2gene2 <- mapgo2gene[order(mapgo2gene[,6]),]
colnames(mapgo2gene2) <- c('GeneinGenome','GenenotinGO','GennomeinGO','GenomenotinGO','Fisher','FDR')
## save the results (we have saved them in './4.plot/figure6/data')
write.csv(mapgo2gene2, '1DNAseq.tx.csv')


### DNAmethy 
## survival model
total <- read.table('coxDNAmethy.txt',header=T,row.names = 1)
p2.1 <- row.names(total)[which(total[,1]<0.05)]
mapgo2gene <- matrix(NA, ncol=6,nrow=length(table.go2gene[,1]))
rownames(mapgo2gene) <- table.go2gene[,1]
for (i in 1:length(table.go2gene[,1])) {
  b <- length(intersect(p2.1, allgenome)) 
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    mapgo2gene[i,1] <- length(intersect(p2.1,unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]))  ))
    mapgo2gene[i,2] <- b - length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) ))
    mapgo2gene[i,3] <- length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    mapgo2gene[i,4] <- length(allgenome) - length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    job <- matrix(c(mapgo2gene[i,1], mapgo2gene[i,2], mapgo2gene[i,3], mapgo2gene[i,4]),nrow=2)
    try(mapgo2gene[i,5] <- fisher.test(job)$p.value, silent=T) 
  }
}
mapgo2gene[,6] <- p.adjust(mapgo2gene[,5],method='fdr',n=length(mapgo2gene[,6]))
mapgo2gene2 <- mapgo2gene[order(mapgo2gene[,6]),]
colnames(mapgo2gene2) <- c('GeneinGenome','GenenotinGO','GennomeinGO','GenomenotinGO','Fisher','FDR')
## save the results (we have saved them in './4.plot/figure6/data')
write.csv(mapgo2gene2, '1DNAmethy.survival.csv')

## treatment outcome model
total <- read.csv('DNAmethy_pvalue.csv',row.names = 1)
p2.1 <- row.names(total)[which(total[,1]<0.05)]
mapgo2gene <- matrix(NA, ncol=6,nrow=length(table.go2gene[,1]))
rownames(mapgo2gene) <- table.go2gene[,1]
for (i in 1:length(table.go2gene[,1])) {
  b <- length(intersect(p2.1, allgenome)) 
  if (as.character(xx[table.go2gene[i,1]]) != "NULL") {
    mapgo2gene[i,1] <- length(intersect(p2.1,unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]]))  ))
    mapgo2gene[i,2] <- b - length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) ))
    mapgo2gene[i,3] <- length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    mapgo2gene[i,4] <- length(allgenome) - length(unique(as.character(eg2symbol[as.data.frame(xx[table.go2gene[i,1]],stringsAsFactors= F)[,1]])) )
    job <- matrix(c(mapgo2gene[i,1], mapgo2gene[i,2], mapgo2gene[i,3], mapgo2gene[i,4]),nrow=2)
    try(mapgo2gene[i,5] <- fisher.test(job)$p.value, silent=T) 
  }
}
mapgo2gene[,6] <- p.adjust(mapgo2gene[,5],method='fdr',n=length(mapgo2gene[,6]))
mapgo2gene2 <- mapgo2gene[order(mapgo2gene[,6]),]
colnames(mapgo2gene2) <- c('GeneinGenome','GenenotinGO','GennomeinGO','GenomenotinGO','Fisher','FDR')
## save the model (we have saved them in './4.plot/figure6/data')
write.csv(mapgo2gene2, '1DNAmethy.tx.csv')