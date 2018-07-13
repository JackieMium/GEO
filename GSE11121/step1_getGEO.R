if(FALSE){
  library(GEOquery)
  gset <- getGEO('GSE11121', destdir="raw_data/",
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  save(gset,file='output_data/GSE11121_eSet.Rdata')
}

load('output_data/GSE11121_eSet.Rdata')

b = gset[[1]]
raw_exprSet=exprs(b)
phe=pData(b)
# we will only use part of all cilnical information about these samples in this tutorial
colnames(phe)
phe=phe[,c(43:46,48)]
save(raw_exprSet,phe, file='output_data/GSE11121_raw_exprSet.Rdata')
