# Best practice for mRNA microarray

 Note : Please **don't use it** if you are not the fan of our [biotrainee](http://www.bio-info-trainee.com/), Thanks.

### Install required packages  by the codes below:

```r
source("http://bioconductor.org/biocLite.R") 
install.packages('devtools')
BiocInstaller::biocLite("jmzeng1314/biotrainee")
library(biotrainee)
```

But if you are in China, you should use the codes below:

```r
install.packages("devtools", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
library(devtools) 
source("https://bioconductor.org/biocLite.R") 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")  
BiocInstaller::biocLite('org.Hs.eg.db')
install.packages("remotes",repos="https://mirror.lzu.edu.cn/CRAN/")
BiocInstaller::biocLite("jmzeng1314/biotrainee")
install.packages("pheatmap",repos="https://mirror.lzu.edu.cn/CRAN/")
```

It will install many other packages for you automately, such as : `ALL, CLL, pasilla, airway ,limma，DESeq2，clusterProfiler  ` , that's why it will take a long time to finish if all of these packages are not installed before in your computer. 

### Then run  step1 :

It always not very easy to download data if you are in China, so I also upload the   file `GSE42872_raw_exprSet.Rdata` , you can load it directly. 

```r
if(F){
  library(GEOquery)
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(gset,'GSE42872.gset.Rdata')
}
load('GSE42872_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
group_list=c(rep('control',3),rep('case',3))
save(raw_exprSet,group_list,
     file='GSE42872_raw_exprSet.Rdata')

```

### Then step2: 

Try to understand my codes, how did I filter the probes by the annotation of each microarry, and how I check the group information for the different samples in each experiment.

Including PCA and Cluster figures, as below:

![Cluster](http://www.bio-info-trainee.com/wp-content/uploads/2018/07/hclust.png)

![PCA](http://www.bio-info-trainee.com/wp-content/uploads/2018/07/pca.png)



Please ensure that you do run those codes by yourself !!!

### Then step3:

Normally we will do differential expression analysis for the microarray, and LIMMA is one of the best method, so I just use it. If the expression matrix(raw counts ) comes from mRNA-seq, you can also choose DESeq based on negative binomial (NB) distributions or baySeq and EBSeq.

Once DEG finished, we can choose top N genes for heatmap as below:

![heatmap](http://www.bio-info-trainee.com/wp-content/uploads/2018/07/DEG_top50_heatmap.png)

and volcano plot as below:

![](http://www.bio-info-trainee.com/wp-content/uploads/2018/07/volcano.png)

### Last step :

Annotation for the significantly changed genes, over-representation test or GSEA for GO/KEGG/biocarta/rectome/MsigDB and so on. 

![KEGG_GSEA](http://www.bio-info-trainee.com/wp-content/uploads/2018/07/kegg_up_down_gsea.png)

![KEGG-enrichment](http://www.bio-info-trainee.com/wp-content/uploads/2018/07/kegg_up_down.png)

### 最重要的是：

如果你觉得我的教程对你有帮助，请赞赏一杯咖啡哦！

如果你的赞赏超过了50元，请在扫描赞赏的同时留下你的邮箱地址，我会发送给你一个惊喜哦！

![](http://www.bio-info-trainee.com/wp-content/uploads/2016/09/jimmy-donate.jpg)

### 广告时间

关于我们

- 我的博客：生信菜鸟团 <http://www.bio-info-trainee.com/>
- 我们的论坛：生信技能树 <http://www.biotrainee.com/thread-1376-1-1.html>
- 我们的VIP社区：<https://vip.biotrainee.com/d/311->
- 我们的微信公众号: <https://mp.weixin.qq.com/s/egAnRfr3etccU_RsN-zIlg>
- 我们的知识星球: <https://t.zsxq.com/VjmQZNn>
- 我们的腾讯课堂： <https://biotree.ke.qq.com/>
- 请善用搜索功能：<http://weixin.sogou.com/>
- 发邮件向我反馈 jmzeng1314@163.com