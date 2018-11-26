### figure 5A
# total <- read.table('coxDNAseq.txt',header=T,row.names = 1)
# p.sv1 <- row.names(total)
# total <- read.table('coxRNAseq.txt',header=T,row.names = 1)
# p.sv2 <- row.names(total)
# total <- read.table('coxDNAmethy.txt',header=T,row.names = 1)
# p.sv3 <- row.names(total)
# int <- intersect(intersect(p.sv1, p.sv2),p.sv3)
# ##### 
# #> intersect(p.sv1, p.sv2)
# #[1] "RIN2"   "FAM21A" "CXCL9" 
# #> intersect(p.sv1, p.sv3)
# #[1] "P4HB"    "ZNF702P"
# #> intersect(p.sv2, p.sv3)
# #[1] "EHD1"     "CACNA2D1"
# ################################################
# total <- read.csv('DNAseq.pvalue.csv',row.names = 1)
# p.tx1 <- row.names(total)[which(total[,11]<0.05)]
# total <- read.csv('RNASeq_pvalue.csv',row.names = 1)
# p.tx2 <- row.names(total)[which(total[,1]<0.05)]
# total <- read.csv('DNAmethy_pvalue.csv',row.names = 1)
# p.tx3 <- row.names(total)[which(total[,1]<0.05)]
# ####
# #> intersect(p.tx1, p.tx2)
# #[1] "ATP9B" "IRF6"  "PTPN5"
# #> intersect(p.tx1, p.tx3)
# #[1] "A1CF"      "ABHD1"     "ACP1"      "ADAMTS8"   "AKAP8"     "APOA1BP"  
# #[7] "ARL8B"     "ARPC1B"    "BRD9"      "C14orf166" "CABIN1"    "CCDC117"  
# #[13] "CERCAM"    "CHCHD1"    "CLTCL1"    "COL4A4"    "CRTAP"     "CXCL11"   
# #[19] "DNAJC21"   "DYRK4"     "ENTPD8"    "FAM162A"   "FBXO7"     "FRMPD1"   
# #[25] "GABARAP"   "GLRA2"     "GTPBP1"    "H3F3A"     "HNRNPA2B1" "JAG1"     
# #[31] "LAMP3"     "LEPROT"    "LRIG1"     "MAGOH"     "MCM5"      "MDC1"     
# #[37] "MOV10"     "MTHFS"     "MTO1"      "MYL6"      "NAA15"     "NDUFAF1"  
# #[43] "NIPBL"     "OMA1"      "PART1"     "PDE4A"     "PI4KA"     "PIGU"     
# #[49] "PLAA"      "POLE"      "PRKY"      "RAP2B"     "RREB1"     "SCGB3A2"  
# #[55] "SLC30A4"   "TAF4"      "TARS"      "TBCA"      "TNFRSF21"  "TNK2"     
# #[61] "TNPO3"     "TOP3B"     "TRIM27"    "TSPO"      "TSTD1"     "TWF2"     
# #[67] "UBE2E2"    "UBE2L3"    "USP8"      "ZFAND6"    "ZNF32"    
# #> intersect(p.tx2, p.tx3)
# #[1] "FANCA" "LEPR"  "NR4A2"

library('VennDiagram')
library('grid')
library('futile.logger')
n1 <- 892
n2 <- 1685
n3 <- 362
n12 <- 70
n13 <- 13
n23 <- 27
n123 <- 2
col.helas3 <- '#1f78b4'
col.huvec <- '#33a02c'
col.k562 <- '#fb9a99'
setEPS()
png('figure5A.png',width=6,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
draw.triple.venn(n1, n2, n3, n12, n23, n13, n123, category=c('DNA-seq','RNA-seq','DNA-methy'),col=c(col.helas3,col.huvec,col.k562),fill=c(col.helas3,col.huvec,col.k562),cat.dist=c(1,1,1)*0.26, lwd=5,cex=1.5,cat.cex=1.5,margin=0.1)
dev.off()


### figure 5B
n1 <- 43
n2 <- 1408
n3 <- 808
n12 <- 3
n13 <- 3
n23 <- 71
n123 <- 0
col.helas3 <- '#1f78b4'
col.huvec <- '#33a02c'
col.k562 <- '#fb9a99'
setEPS()
png('figure5B.png',width=6,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
draw.triple.venn(n1, n2, n3, n12,n23, n13,  n123, category=c('DNA-seq','RNA-seq','DNA-methy'),col=c(col.helas3,col.huvec,col.k562),fill=c(col.helas3,col.huvec,col.k562),cat.dist=c(1,1,1)*0.26, lwd=5,cex=1.5,cat.cex=1.5,margin=0.1)
dev.off()

### figure 5C
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.sv.config/')
# table <- read.table('labels.3level.txt',stringsAsFactors = F)
# lable1 <- table[,4]
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/rnaseq.sv.config/')
# table <- read.table('labels.3level.txt',stringsAsFactors = F)
# lable2 <- table[,4]
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnamethy.sv.config/')
# table <- read.table('labels.3level.txt',stringsAsFactors = F)
# lable3 <- table[,4]
# length(intersect(lable1,lable2))
library('VennDiagram')
library('grid')
library('futile.logger')
n1 <- 15
n2 <- 18
n3 <- 22
n12 <- 2
#  "GO:0005509"(calcium ion binding) "GO:0008009"(chemokine activity)
n13 <- 1
# "GO:0006351"(transcription, DNA-templated)
n23 <- 1
# "GO:0005743"(mitochondrial inner membrane)
n123 <- 0
col.helas3 <- '#1f78b4'
col.huvec <- '#33a02c'
col.k562 <- '#fb9a99'
setEPS()
png('figure5C.png',width=6,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
draw.triple.venn(n1, n2, n3, n12, n23, n13, n123, category=c('DNA-seq','RNA-seq','DNA-methy'),col=c(col.helas3,col.huvec,col.k562),fill=c(col.helas3,col.huvec,col.k562),cat.dist=c(1,1,1)*0.26, lwd=5,cex=1.5,cat.cex=1.5,margin=0.1)
dev.off()

### figure 5D
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnaseq.tx.config/')
# table <- read.table('labels.3level.txt',stringsAsFactors = F)
# lable1 <- table[,4]
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/rnaseq.tx.config/')
# table <- read.table('labels.3level.txt',stringsAsFactors = F)
# lable2 <- table[,4]
# setwd('~/Project/ov_prognosis_drugprecision/plot/circos/dnamethy.tx.config/')
# table <- read.table('labels.3level.txt',stringsAsFactors = F)
# lable3 <- table[,4]
# length(intersect(lable1,lable2))
library('VennDiagram')
library('grid')
library('futile.logger')
n1 <- 9
n2 <- 25
n3 <- 25
n12 <- 2
#  "GO:0006367"(transcription initiation from RNA polymerase II promoter) 
# "GO:0005737"(cytoplasm)
n13 <- 2
#  "GO:0006461"(protein complex assembly)
# "GO:0005737"(cytoplasm)
n23 <- 8
# [1] "GO:0007186" "GO:0016020" "GO:0005654" "GO:0005730" "GO:0005737" "GO:0004930"
# [7] "GO:0005515" "GO:0000166"
n123 <- 1
# "GO:0005737"(cytoplasm)
col.helas3 <- '#1f78b4'
col.huvec <- '#33a02c'
col.k562 <- '#fb9a99'
setEPS()
png('figure5D.png',width=6,height=6,fonts=c("serif", "Palatino"),units="in",res=300)
draw.triple.venn(n1, n2, n3, n12, n23, n13, n123, category=c('DNA-seq','RNA-seq','DNA-methy'),col=c(col.helas3,col.huvec,col.k562),fill=c(col.helas3,col.huvec,col.k562),cat.dist=c(1,1,1)*0.26, lwd=5,cex=1.5,cat.cex=1.5,margin=0.1)
dev.off()

