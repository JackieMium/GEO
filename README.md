fork 自 Jimmy 的一个处理芯片数据的 Repo。
原 Repo 信息：

>
### ---------------
### Create: Jianming Zeng
### Date: 2018-07-09 20:11:07
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
### ---------------


**All credits goes to [jmzeng1314](https://github.com/jmzeng1314).**
————————————————————————————————————————————————————————————————————————————————

## Repo 内容概览：
- 流程涵盖芯片数据处理的常规流程，包括表达数据矩阵获取、数据读入和存储、芯片注释和 DEG 分析，DEG 结果的常见下游富集分析和简单可视化。
- 项目包含 GSE42872 和 GSE11122 两个不同芯片数据，放在单独文件夹下。
- 每个项目文件夹顶层为所有需要用到的 R 脚本，文件夹 `raw_data` 下存放原始数据，是 GEO 下载的表达矩阵数据；文件夹 `output_data` 存放中间数据或输出结果数据；`output_plots` 存放输出的图片。
- [GSE42872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) 来自 Affymetrix Human Gene 1.0 ST Array [transcript (gene) version] 芯片，6 个样本 3+3 分为两组。详细信息参考 GEO 页面介绍。
- [GSE11122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse11121) 来自 Affymetrix Human Genome U133A Array 芯片，200 个样本。详细信息参考 GEO 页面介绍。

## 使用指北
建议下载 Repo 后双击 `xxx.Rproj` 在 RStudio 中打开项目文件，然后依次打开打开 `stepNum.R` 这些文件就可以愉快地开始了。
注意：

- 每执行完一个脚本时建议重启 R，在 RStudio 中快捷键 `Ctrl + Shifr + F10` 即可；
- 由于原始数据、中间输出数据和图在 Repo 其实里已经有了，所以建议在执行到存储数据或存储图片时直接跳过，或者更好的做法是，更改文件名或者存储路径，这样你就可以把你输出的文件或者图片和 repo 里的做对比，看看你的输出（或者我的）有没有问题。



中国大陆网络环境建议切换镜像到国内并安装需要用到的包，举例:
```r
options(repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("pheatmap")

source("https://bioconductor.org/biocLite.R") 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
BiocInstaller::biocLite('org.Hs.eg.db')
```

## 部分结果图片

### 样本聚类和主成分分析：

![Cluster](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/hclust.png)

![PCA](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/pca.png)
### 差异基因火山图

![](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/volcano.png)


### 50 个差异基因的热图

![heatmap](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/DEG_top50_heatmap.png)

### KEGG 富集

![KEGG_GSEA](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/kegg_up_down_gsea.png)

![KEGG-enrichment](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/kegg_up_down.png)


### GO 富集

![GO-BP](https://github.com/JackieMium/GEO/blob/master/GSE42872/output_plots/GO_dotplots/dotplot_gene_diff_BP.png)