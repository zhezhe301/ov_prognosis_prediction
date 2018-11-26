################ DNAseq ###############################################################
##### survival
library('org.Hs.eg.db')
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
eg2symbol <- as.list(x[mapped_genes])
xx <- as.list(org.Hs.egGO2EG)
#xxx <- as.list(org.Hs.egGO2ALLEGS)
setwd('../../../')
total <- read.table('coxDNAseq.txt',header=T,row.names = 1)
p2.1 <- row.names(total)[which(total[,1]<0.05)]
setwd('./4.plot/figure6/data')
table.go2gene <- read.csv('1DNAseq.survival.csv',stringsAsFactors= F)
table.gofamily <- read.csv('goid.new.csv',stringsAsFactors= F)
#table.allgo <- read.csv('goid2genes.csv')
table.go2gene <- table.go2gene[which(table.go2gene[,6]< 0.05),]
list.go <- table.go2gene[which(table.go2gene[,2] != 0),1]

#### BP ####
table <- matrix(NA, ncol = 11,nrow = length(list.go))
table[,1] <- list.go
for (i in 1:length(list.go)) {
  if (is.na(table.gofamily[which(table.gofamily[,2]==list.go[i]),3]) & length(rownames(as.data.frame(xx[table[i,1]],stringsAsFactors= F)))!=0) {
    table[i,2] <- length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table[i,1]],stringsAsFactors= F)[,1]])) ))
  }
}
table <- table[complete.cases(table[,2]),]     
for (i in 1:length(table[,1])){
  table[i,9] <- table.gofamily[which(table.gofamily[,2]==table[i,1]),7]
  table[i,5] <- table.gofamily[which(table.gofamily[,2]==table[i,1]),8]
  if (is.na(table[i,5])) {
    if (is.na(table[i,9])) {
      table[i,5] <- 'Other_bp'
    } else {
      table[i,5] <- 'GO:0008150'
    }
  }
}   
table <- table[order(table[,5]),]
table[1,3] <- 0
table[1,4] <- table[1,2]
for (i in 2:length(table[,1])) {
  table[i,3] <- as.numeric(table[i-1,4])
  table[i,4] <- as.numeric(table[i,2]) +  as.numeric(table[i,3])
} 
for (i in 1:length(table[,1])) {
  list.father <- which(table[,5]==table[i,5])
  table[i,6] <- sum(as.numeric(table[list.father,2]))
}
table[1,10] <- sum(as.numeric(table[,2]))
table[1,7] <- table[1,3]
table[1,8] <- table[1,6]
j <- 1
for (i in 1:length(table[,1])) {
  if (is.na(table[i,5])==F ) {
    if (!table[i,5] %in% table[(1:(i-1)),5]) {
      table[i,7] <- as.numeric(table[j,8])
      table[i,8] <- as.numeric(table[j,8]) + as.numeric(table[i,6])
      j <- i
    }
  }
}
table[,11] <- 'bp'
colnames(table) <- c('goid','count','start','end','father','f_count','f_start','f_end','grandpa','g_count','chromosome')
table.bp <- table
#setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
#write.csv(table, 'table.bp.csv')

table.bp.links <- matrix(NA, ncol=5,nrow = as.numeric(table[1,10]))
num <- 0
for (i in 1:length(table[,1])) {
  for (j in 1:as.numeric(table[i,2])) {
    table.bp.links[num+j,2] <-  intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table[i,1]],stringsAsFactors= F)[,1]])) )[j] 
    table.bp.links[num+j,4] <- num + j
    table.bp.links[num+j,3] <- num + j - 1
  }
  table.bp.links[((num+1):(num+j)),1] <- table[i,1]
  num <-  num + j
}
table.bp.links[,5] <- 'bp'
colnames(table.bp.links) <- c('goid','gene','start','end','chromosome')
#setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
# write.csv(table.bp.links, 'table.bp.links.csv')


### CC ###
table <- matrix(NA, ncol = 11,nrow = length(list.go))
table[,1] <- list.go
for (i in 1:length(list.go)) {
  if (is.na(table.gofamily[which(table.gofamily[,2]==list.go[i]),4])) {
    table[i,2] <- length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table[i,1]],stringsAsFactors= F)[,1]])) ))
  }
}
table <- as.data.frame(table,stringsAsFactors= F)
table <- table[complete.cases(table[,2]),]  
for (i in 1:length(table[,1])){
  table[i,9] <- table.gofamily[which(table.gofamily[,2]==table[i,1]),14]
  table[i,5] <- table.gofamily[which(table.gofamily[,2]==table[i,1]),15]
  if (is.na(table[i,5])) {
    if (is.na(table[i,9])) {
      table[i,5] <- 'Other_cc'
    } else {
      table[i,5] <- 'GO:0005575'
    }
  }
}  
table <- table[order(table[,5]),]
table[1,3] <- 0
table[1,4] <- table[1,2]
# for (i in 2:length(table[,1])) {
#   table[i,3] <- as.numeric(table[i-1,4])
#   table[i,4] <- as.numeric(table[i,2]) +  as.numeric(table[i,3])
# }    
for (i in 1:length(table[,1])) {
  list.father <- which(table[,5]==table[i,5])
  table[i,6] <- sum(as.numeric(table[list.father,2]))
}
table[1,10] <- sum(as.numeric(table[,2]),na.rm = T)
table[1,7] <- table[1,3]
table[1,8] <- table[1,6]
j <- 1
for (i in 1:length(table[,1])) {
  if (is.na(table[i,5])==F ) {
    if (!table[i,5] %in% table[(1:(i-1)),5]) {
      table[i,7] <- as.numeric(table[j,8])
      table[i,8] <- as.numeric(table[j,8]) + as.numeric(table[i,6])
      j <- i
    }
  }
}
table[,11] <- 'cc'
colnames(table) <- c('goid','count','start','end','father','f_count','f_start','f_end','grandpa','g_count','chromosome')
table.cc <- table
table.cc <- table.cc[complete.cases(table.cc[,1]),]
table.cc[1,9] <- table[1,5]
#setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
#write.csv(table, 'table.cc.csv')

table.cc.links <- matrix(NA, ncol=5,nrow = as.numeric(table[1,10]))
num <- 0
for (i in 1:length(table[,1])) {
  for (j in 1:as.numeric(table[i,2])) {
    table.cc.links[num+j,2] <-  intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table[i,1]],stringsAsFactors= F)[,1]])) )[j] 
    table.cc.links[num+j,4] <- num + j
    table.cc.links[num+j,3] <- num + j - 1
  }
  table.cc.links[((num+1):(num+j)),1] <- table[i,1]
  num <-  num + j
}
table.cc.links[,5] <- 'cc'
colnames(table.cc.links) <- c('goid','gene','start','end','chromosome')
#setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
# write.csv(table.cc.links, 'table.cc.links.csv')



### MF ###
table <- matrix(NA, ncol = 11,nrow = length(list.go))
table[,1] <- list.go
for (i in 1:length(list.go)) {
  if (is.na(table.gofamily[which(table.gofamily[,2]==list.go[i]),5])) {
    table[i,2] <- length(intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table[i,1]],stringsAsFactors= F)[,1]])) ))
  }
}
table <- table[complete.cases(table[,2]),]     
for (i in 1:length(table[,1])){
  table[i,9] <- table.gofamily[which(table.gofamily[,2]==table[i,1]),21]
  table[i,5] <- table.gofamily[which(table.gofamily[,2]==table[i,1]),22]
  if (is.na(table[i,5])) {
    if (is.na(table[i,9])) {
      table[i,5] <- 'Other_mf'
    } else {
      table[i,5] <- 'GO:0003674'
    }
  }
}   
table <- table[order(table[,5]),]
table[1,3] <- 0
table[1,4] <- table[1,2]
for (i in 2:length(table[,1])) {
  table[i,3] <- as.numeric(table[i-1,4])
  table[i,4] <- as.numeric(table[i,2]) +  as.numeric(table[i,3])
}    
for (i in 1:length(table[,1])) {
  list.father <- which(table[,5]==table[i,5])
  table[i,6] <- sum(as.numeric(table[list.father,2]))
}
table[1,10] <- sum(as.numeric(table[,2]))
table[1,7] <- table[1,3]
table[1,8] <- table[1,6]
j <- 1
for (i in 1:length(table[,1])) {
  if (is.na(table[i,5])==F ) {
    if (!table[i,5] %in% table[(1:(i-1)),5]) {
      table[i,7] <- as.numeric(table[j,8])
      table[i,8] <- as.numeric(table[j,8]) + as.numeric(table[i,6])
      j <- i
    }
  }
}
table[,11] <- 'mf'
colnames(table) <- c('goid','count','start','end','father','f_count','f_start','f_end','grandpa','g_count','chromosome')
table.mf <- table
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
# write.csv(table, 'table.mf.csv')

table.mf.links <- matrix(NA, ncol=5,nrow = as.numeric(table[1,10]))
num <- 0
for (i in 1:length(table[,1])) {
  for (j in 1:as.numeric(table[i,2])) {
    table.mf.links[num+j,2] <-  intersect(p2.1, unique(as.character(eg2symbol[as.data.frame(xx[table[i,1]],stringsAsFactors= F)[,1]])) )[j] 
    table.mf.links[num+j,4] <- num + j
    table.mf.links[num+j,3] <- num + j - 1
  }
  table.mf.links[((num+1):(num+j)),1] <- table[i,1]
  num <-  num + j
}
table.mf.links[,5] <- 'mf'
colnames(table.mf.links) <- c('goid','gene','start','end','chromosome')
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
# write.csv(table.mf.links, 'table.mf.links.csv')


### table #### 
table <- rbind(table.bp, table.cc, table.mf)
table<- as.matrix(table,stringsAsFactors= F)
karyotype <- matrix(NA, ncol =7 ,nrow=3)
karyotype[,1] <- 'chr'
karyotype[,2] <- '-'
karyotype[,3] <- c('bp','cc','mf')
karyotype[,4] <- c('BIOLOGICAL_PROCESS','CELLULAR_COMPONENT','MOLECULAR_FUNCTION')
karyotype[,5] <- 0
karyotype[1,6] <- table.bp[1,10] 
karyotype[2,6] <- table.cc[1,10]
karyotype[3,6] <- table.mf[1,10] 
karyotype[,7] <- c('try1','try2','try3')
setwd('../figure6A.dnaseq.sv.config')
write.table(karyotype,'karyotype.dnaseq.txt',row.names=F,col.names=F,quote= F,sep=' ')

highlight.2level <- matrix(NA, ncol =4 ,nrow =length(table[,1]))
for (i in 1:length(table[,1])) {
  if (is.na(table[i,7])==F) {
    highlight.2level[i,1] <- table[i,11]
    highlight.2level[i,2] <- table[i,7]
    highlight.2level[i,3] <- table[i,8]
    highlight.2level[i,4] <- paste('fill_color=chr',sample(1:20,1), sep='')
  }
}
highlight.2level <- na.omit(highlight.2level)
for (i in 1:length(highlight.2level[,1])) {
  if (i %% 2) {
    if (highlight.2level[i,1] == 'bp') {
      highlight.2level[i,4] <- 'fill_color=try1_a2'
    } else {
      if (highlight.2level[i,1] == 'cc') {
        highlight.2level[i,4] <- 'fill_color=try2_a2'
      } else {
        highlight.2level[i,4] <- 'fill_color=try3_a2'
      }
    }
  } else {
    if (highlight.2level[i,1] == 'bp') {
      highlight.2level[i,4] <- 'fill_color=try1_a3'
    } else {
      if (highlight.2level[i,1] == 'cc') {
        highlight.2level[i,4] <- 'fill_color=try2_a3'
      } else {
        highlight.2level[i,4] <- 'fill_color=try3_a3'
      }
    }
  }
}
write.table(highlight.2level,'highlight.2level.txt',row.names=F,col.names=F,quote= F,sep=' ')

highlight.3level <- matrix(NA, ncol =4 ,nrow =length(table[,1]))
for (i in 1:length(table[,1])) {
  highlight.3level[i,1] <- table[i,11]
  highlight.3level[i,2] <- table[i,3]
  highlight.3level[i,3] <- table[i,4]
  highlight.3level[i,4] <- paste('fill_color=chr',sample(1:20,1), sep='')
}
highlight.3level <- na.omit(highlight.3level)
for (i in 1:length(highlight.3level[,1])) {
  if (i %% 2) {
    if (highlight.3level[i,1] == 'bp') {
      highlight.3level[i,4] <- 'fill_color=try1_a2'
    } else {
      if (highlight.3level[i,1] == 'cc') {
        highlight.3level[i,4] <- 'fill_color=try2_a2'
      } else {
        highlight.3level[i,4] <- 'fill_color=try3_a2'
      }
    }
  } else {
    if (highlight.3level[i,1] == 'bp') {
      highlight.3level[i,4] <- 'fill_color=try1_a3'
    } else {
      if (highlight.3level[i,1] == 'cc') {
        highlight.3level[i,4] <- 'fill_color=try2_a3'
      } else {
        highlight.3level[i,4] <- 'fill_color=try3_a3'
      }
    }
  }
}
write.table(highlight.3level,'highlight.3level.txt',row.names=F,col.names=F,quote= F,sep=' ')



#### histogram ####
histogram <- matrix(NA, ncol=4, nrow =length(table[,1]))
for (i in 1:length(table[,1])) {
  histogram[i,1] <- table[i,11]
  histogram[i,2] <- table[i,3]
  histogram[i,3] <- table[i,4]
  histogram[i,4] <- -log(table.go2gene[which(table.go2gene[,1]== table[i,1]),6],10)
}
write.table(histogram,'histogram.3level.txt',row.names=F,col.names=F,quote= F,sep=' ')


### text labels ####
labels <- matrix(NA, ncol=4, nrow =length(table[,1]))
for (i in 1:length(table[,1])) {
  labels[i,1] <- table[i,11]
  labels[i,2] <- table[i,3]
  labels[i,3] <- table[i,4]
  labels[i,4] <- table[i,1]
}
write.table(labels,'labels.3level.txt',row.names=F,col.names=F,quote= F,sep=' ')

labels <- matrix(NA, ncol=4, nrow =length(table[,1]))
for (i in 1:length(table[,1])) {
  labels[i,1] <- table[i,11]
  labels[i,2] <- table[i,7]
  labels[i,3] <- table[i,8]
  labels[i,4] <- table.gofamily[which(table.gofamily[,2]==table[i,1] ),28] 
}
labels <- labels[complete.cases(labels[,2]),]
write.table(labels,'labels.2level.txt',row.names=F,col.names=F,quote= F,sep=' ')


### all links ###
table.alllinks <- rbind(table.bp.links, table.cc.links, table.mf.links)
k <- 1
for (i in 1:length(table.alllinks[,1])) {
  for (j in (i+1):length(table.alllinks[,1])) {
    if (j <= length(table.alllinks[,1])) {
      if (table.alllinks[i,2] == table.alllinks[j,2]) {
        k <- k +1
      }
    }
  }
}
table.links <- matrix(NA, nrow =k-1 , ncol = 6)
k <- 1
for (i in 1:length(table.alllinks[,1])) {
  for (j in (i+1):length(table.alllinks[,1])) {
    if (j <= length(table.alllinks[,1])) {
      if (table.alllinks[i,2] == table.alllinks[j,2]) {
        table.links[k,1] <- table.alllinks[i,5]
        table.links[k,2] <- table.alllinks[i,3]
        table.links[k,3] <- table.alllinks[i,4]
        table.links[k,4] <- table.alllinks[j,5]
        table.links[k,5] <- table.alllinks[j,3]
        table.links[k,6] <- table.alllinks[j,4]
        k <- k +1
      }
    }
  }
}
write.table(table.links,'links.txt',row.names=F,col.names=F,quote= F,sep=' ')

table.links1 <- table.links
for (i in 1:length(table.links1[,1])) {
  if (table.links1[i,1]!='bp') {
    table.links1[i,] <- NA
  } else {
    if (table.links1[i,4] != 'bp') {
      table.links1[i,] <- NA
    } 
  }
}
table.links1 <- na.omit(table.links1)
write.table(table.links1,'links1.txt',row.names=F,col.names=F,quote= F,sep=' ')

table.links2 <- table.links
for (i in 1:length(table.links2[,1])) {
  if (table.links2[i,1]!='cc') {
    table.links2[i,] <- NA
  } else {
    if (table.links2[i,4] != 'cc') {
      table.links2[i,] <- NA
    } 
  }
}
table.links2 <- na.omit(table.links2)
write.table(table.links2,'links2.txt',row.names=F,col.names=F,quote= F,sep=' ')

table.links3 <- table.links
for (i in 1:length(table.links3[,1])) {
  if (table.links3[i,1]!='mf') {
    table.links3[i,] <- NA
  } else {
    if (table.links3[i,4] != 'mf') {
      table.links3[i,] <- NA
    } 
  }
}
table.links3 <- na.omit(table.links3)
write.table(table.links3,'links3.txt',row.names=F,col.names=F,quote= F,sep=' ')

table.links4 <- table.links
for (i in 1:length(table.links4[,1])) {
  if (table.links4[i,1]==table.links4[i,4]) {
    table.links4[i,] <- NA
  }
}
table.links4 <- na.omit(table.links4)
write.table(table.links4,'links4.txt',row.names=F,col.names=F,quote= F,sep=' ')
