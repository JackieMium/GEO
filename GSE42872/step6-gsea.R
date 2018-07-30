rm(list=ls())
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-09 20:11:07
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------


load(file='GSE42872_DEG.Rdata')
source('functions.R')
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(rownames(DEG), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
head(DEG)
DEG$SYMBOL = rownames(DEG)
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)
geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)

head(geneList)


## Molecular Signatures Database (MSigDb) 
d='~/biosoft/MSigDB/entrez/'
gmts=list.files(d,pattern = 'all')
gmts
if(T){
  msigdb <- lapply(gmts, function(gmtfile){
    geneset <- read.gmt(file.path(d,gmtfile)) 
    egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
    head(egmt)
    return(egmt)
  })
} 
save(msigdb,file='msigdb.Rdata')
 
load(file='msigdb.Rdata')
msigdb_df <- lapply(msigdb, function(x){
  cat(paste(dim(x@result)),'\n')
  x@result
})
gmts
df=do.call(rbind ,msigdb_df)
gseaplot(msigdb[[2]],'KEGG_CELL_CYCLE')  



