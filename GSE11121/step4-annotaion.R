load(file='output_data/GSE11121_DEG.Rdata')
source('functions.R')
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(rownames(DEG), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df); head(DEG)
DEG$SYMBOL = rownames(DEG)
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)

gene_up= DEG[DEG$change == 'UP','ENTREZID'] 
gene_down=DEG[DEG$change == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )

data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = TRUE)

## KEGG pathway analysis
if(TRUE){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.05)
  head(kk.down)[,1:6]
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1

  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)

  ggsave(g_kegg,filename = 'output_plots/kegg_up_down.png')
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = 'output_plots/kegg_up_down_gsea.png')
}

### GO database analysis 

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)

if(TRUE){
  go_enrich_results <- lapply(g_list , function(gene) {
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Hs.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      
      print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = 'output_data/go_enrich_results.Rdata')
}


load(file = 'output_data/go_enrich_results.Rdata')

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    fpath <- ('output_plots/GO_dotplots/')
    png(paste0(fpath, fn),res=150,width = 1080)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}