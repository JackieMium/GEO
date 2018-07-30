### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-15 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------

if(F){
  install.packages("devtools",
                   repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(devtools) 
  source("https://bioconductor.org/biocLite.R") 
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")  
  BiocInstaller::biocLite(c('airway','DESeq2','edgeR','limma'))
  
}

library(DESeq2)
library(edgeR)
library(limma)
library(airway)

if(F){
  library(airway)
  data(airway)
  exprSet=assay(airway)
  group_list=colData(airway)[,3]
  save(exprSet,group_list,file = 'airway_exprSet.Rdata')
}
rm(list = ls())
options(stringsAsFactors = F)
load(file = 'airway_exprSet.Rdata')
source('functions.R')
### ---------------
###
### Firstly run DEseq2 
###
### ---------------
suppressMessages(library(DESeq2)) 
group_list=c('untrt','trt' ,'untrt','trt' ,'untrt','trt','untrt','trt')
(colData <- data.frame(row.names=colnames(exprSet), 
                       group_list=group_list) )
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds)
res <- results(dds, 
               contrast=c("group_list","trt","untrt"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG = na.omit(DEG)
if(F){
  
  png("DESeq2_qc_dispersions.png", 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  rld <- rlogTransformation(dds)
  exprMatrix_rlog=assay(rld)
  
  x=apply(exprMatrix_rlog,1,mean)
  y=apply(exprMatrix_rlog,1,mad) 
  plot(x,y) 
  
  png("DESeq2_RAWvsNORM.png",height = 800,width = 800)
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(exprSet, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()
  
  
  
}
nrDEG=DEG
DEseq_DEG=nrDEG
nrDEG=DEseq_DEG[,c(2,6)]
colnames(nrDEG)=c('log2FoldChange','pvalue') 
draw_h_v(exprSet,nrDEG,'DEseq2')
source('functions.R')
### ---------------
###
### Then run edgeR 
###
### ---------------
library(edgeR)
d <- DGEList(counts=exprSet,group=factor(group_list))
keep <- rowSums(cpm(d)>1) >= 2
table(keep)
d <- d[keep, , keep.lib.sizes=FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples

design <- model.matrix(~0+factor(group_list))
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))
dge=d
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
# https://www.biostars.org/p/110861/
lrt <- glmLRT(fit,  contrast=c(-1,1)) 
nrDEG=topTags(lrt, n=nrow(dge))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_DEG =nrDEG 
nrDEG=edgeR_DEG[,c(1,5)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'edgeR')

source('functions.R')
### ---------------
###
### Then run edgeR 
###
### --------------- 
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design

dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('untrt-trt'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='untrt-trt', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
head(DEG_limma_voom)
nrDEG=DEG_limma_voom[,c(1,4)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'limma')


save(DEG_limma_voom,DEseq_DEG,edgeR_DEG,
     dds,exprSet,group_list,
     file = 'DEG_results.Rdata')


rm(list = ls())
load(file = 'DEG_results.Rdata')
source('functions.R')

nrDEG=DEG_limma_voom[,c(1,4)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'limma')

nrDEG=edgeR_DEG[,c(1,5)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'edgeR')


nrDEG=DEseq_DEG[,c(2,6)]
colnames(nrDEG)=c('log2FoldChange','pvalue') 
draw_h_v(exprSet,nrDEG,'DEseq2')




