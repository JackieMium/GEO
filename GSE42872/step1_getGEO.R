if(FALSE){
  library(GEOquery)
  eSet <- getGEO('GSE42872', destdir="raw_data/",
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
  save(eSet,file='output_data/GSE42872_eSet.Rdata')
}
load('output_data/GSE42872_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
group_list=c(rep('control',3),rep('case',3))

# not run since files are already there
save(raw_exprSet,group_list,
     file='output_data/GSE42872_raw_exprSet.Rdata')
