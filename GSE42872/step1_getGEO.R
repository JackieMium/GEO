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

if(F){
  library(GEOquery)
  eSet <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE42872_eSet.Rdata')
}
load('GSE42872_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
group_list=c(rep('control',3),rep('case',3))
save(raw_exprSet,group_list,
     file='GSE42872_raw_exprSet.Rdata')
