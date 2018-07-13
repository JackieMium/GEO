library('hugene10sttranscriptcluster.db')

if(FALSE){
  downGSE <- function(studyID = "GSE1009", destdir = ".") {
    
    library(GEOquery)
    eSet <- getGEO(studyID, destdir = destdir, getGPL = FALSE)
    
    exprSet = exprs(eSet[[1]])
    pdata = pData(eSet[[1]])
    
    write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
    write.csv(pdata, paste0(studyID, "_metadata.csv"))
    return(eSet)
    
  }
  downGSE('GSE42872')
}

load(file='output_data/GSE42872_raw_exprSet.Rdata')
exprSet=raw_exprSet

# explore annotation information
ids=toTable(hugene10sttranscriptclusterSYMBOL)
# unique probe-gene pairs
length(unique(ids$symbol))
# some multi-mapping examples
tail(sort(table(ids$symbol)))
# summary of mapping
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))

# how many probes are in the annotation file?
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
# select only probes that have an annotation entry
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)
# select probes that we have in expreSet in annotation file
ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
head(exprSet)

# annotate eSet and select max in case of multi-mapping
rm_Dup_anno <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}
new_exprSet <- rm_Dup_anno(exprSet,ids)

# not run since file is already there
save(new_exprSet,group_list, file='GSE42872_new_exprSet.Rdata')
load(file='GSE42872_new_exprSet.Rdata')

exprSet <- new_exprSet

# PCA and Clutering all-in-one function
# not run
if(TRUE){
  library(reshape2)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(exprSet))
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p)
  
  ## hclust
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10))
  png('hclust.png',res=120)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  
  ## PCA
  
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  png('pca.png',res=120)
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
  dev.off()
}
