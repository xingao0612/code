# Perform quality control on single-cell data.
options(stringsAsFactors = F)
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(monocle)


rt=read.table(file = "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt",header=T,check.names=F,row.names = 1)

rt[1:5,1:5]

sce <- CreateSeuratObject(counts = rt,project = "seurat", min.cells = 3, min.features = 250,names.delim = "_",names.field = 1)

View(sce@meta.data)


###save(sce,file="sce.Rdata")

#load("sce.Rdata")

sce[["percent.mt"]]=PercentageFeatureSet(object = sce, pattern = "^MT-")

sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比

pdf("1.QC-VlnPlot.pdf",width = 8,height = 4.5)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), 
        ncol = 3)
dev.off()


plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("2.QC-FeatureScatter.pdf"),width = 18,height = 7)
CombinePlots(plots = list(plot1, plot2),legend = "right")
dev.off()
rm(plot1,plot2)


raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

pearplot_befor1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
pearplot_befor1
ggsave(filename = '3.QC_before1.pdf',plot = pearplot_befor1,he=7,wi=18)
#rm(sce)


pdf(file = "4.hist_of_mt.pdf",w=12,h=10)
hist(sce@meta.data$percent.mt,col = "red",breaks = 100)
dev.off()

pdf(file = "4.hist_of_nFeature_RNA.pdf",w=12,h=10)
hist(sce@meta.data$nFeature_RNA,col = "red",breaks = 100)
dev.off()

pdf(file = "4.hist_of_nFeature_RNA_截取.pdf",w=12,h=10)
hist(sce@meta.data$nFeature_RNA,col = "red",breaks = 100,xlim = c(0,2000))
dev.off()

pdf(file = "4.hist_of_nCount_RNA_RNA.pdf",w=12,h=10)
hist(sce@meta.data$nCount_RNA,col = "red",breaks = 100)
dev.off()


sce <-  subset(sce, 
               subset = 
                 nFeature_RNA > 500 & 
                 nCount_RNA > 1000 & 
                 #nCount_RNA < 20000 &
                 percent.mt < 15) 




clean_count <- table(sce@meta.data$orig.ident)

#比较一下过滤前后的细胞数差别
raw_count #过滤前
table(sce@meta.data$orig.ident)#过滤后

before <- raw_count
after <- table(sce@meta.data$orig.ident)
a <- rbind(before,after)

write.table(a,file = "5.过滤前后细胞数.txt",sep="\t",row.names=T,col.names=T)
rm(a)
rm(before)
rm(after)

summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$sample <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","sample")
counts[,1] <- as.numeric(counts[,1])
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = sample, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = "Set1") +
  theme(text=element_text(size=10),legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = '5.fit_cell_count.pdf',plot = fit_cell_count,width =13,height = 9)

pearplot_after1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                         pt.size = 0.1, 
                         ncol = 4)
pearplot_after1
ggsave(filename = '6.QC_after.pdf',plot = pearplot_after1,he=7,wi=18)

pearplot_befor1
pearplot_after1
qc_merge<- CombinePlots(plots = list(pearplot_befor1,pearplot_after1) , 
                        nrow=2, legend='none')
qc_merge
ggsave(filename = '6.qc_merge.pdf',plot = qc_merge,he=9,wi=18)


sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
#计算高变基因
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000)  #此处选择高变基因最大量#默认返回2000个基因
#
### 可视化前20个高变基因
top20 <- head(VariableFeatures(sce), 20)
plot1 <- VariableFeaturePlot(sce)  ##全部画出2000个基因且不显示特征基因
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)
ggsave(filename = '7.feat_plot2.pdf',plot = plot2,he=4,wi=6)
feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = '7.feat_20.pdf',plot = feat_20,he=7,wi=16)

scale.genes <-  rownames(sce)
sce <- ScaleData(sce, features = scale.genes)
#样本的分组
meta1<-data.frame(matrix(nrow=length(sce@meta.data$orig.ident), ncol=2))  
colnames(meta1)=c('Sample','Group1')
meta1$Sample=sce@meta.data$orig.ident
unique(meta1$Sample)

meta1[grep("KUL19-T",meta1$Sample),]$Group1="Core region"
meta1[grep("KUL21-T",meta1$Sample),]$Group1="Core region"
meta1[grep("KUL28-T",meta1$Sample),]$Group1="Core region"
meta1[grep("KUL30-T",meta1$Sample),]$Group1="Core region"
meta1[grep("KUL31-T",meta1$Sample),]$Group1="Core region"
meta1[grep("KUL01-T",meta1$Sample),]$Group1="Core region"
###########################################################

meta1[grep("KUL19-B",meta1$Sample),]$Group1="Border region"
meta1[grep("KUL21-B",meta1$Sample),]$Group1="Border region"
meta1[grep("KUL28-B",meta1$Sample),]$Group1="Border region"
meta1[grep("KUL30-B",meta1$Sample),]$Group1="Border region"
meta1[grep("KUL31-B",meta1$Sample),]$Group1="Border region"
meta1[grep("KUL01-B",meta1$Sample),]$Group1="Border region"
###########################################################
meta1[grep("KUL01-N",meta1$Sample),]$Group1="Normal mucosa"
meta1[grep("KUL19-N",meta1$Sample),]$Group1="Normal mucosa"
meta1[grep("KUL21-N",meta1$Sample),]$Group1="Normal mucosa"
meta1[grep("KUL28-N",meta1$Sample),]$Group1="Normal mucosa"
meta1[grep("KUL30-N",meta1$Sample),]$Group1="Normal mucosa"
meta1[grep("KUL31-N",meta1$Sample),]$Group1="Normal mucosa"



sce <- AddMetaData(sce, meta1$Sample,col.name = "Sample")
sce <- AddMetaData(sce, meta1$Group1,col.name = "Group1")
view(sce@meta.data)

meta2<-data.frame(matrix(nrow=length(sce@meta.data$orig.ident), ncol=2))  #先用细胞数建立一个空的数据框
colnames(meta2)=c('Sample','Group2')
meta2$Sample=sce@meta.data$orig.ident
unique(meta2$Sample)

meta2[grep("KUL19-T",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL21-T",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL28-T",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL30-T",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL31-T",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL01-T",meta2$Sample),]$Group2="Tumor"
###########################################################

meta2[grep("KUL19-B",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL21-B",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL28-B",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL30-B",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL31-B",meta2$Sample),]$Group2="Tumor"
meta2[grep("KUL01-B",meta2$Sample),]$Group2="Tumor"
###########################################################
meta2[grep("KUL01-N",meta2$Sample),]$Group2="Normal"
meta2[grep("KUL19-N",meta2$Sample),]$Group2="Normal"
meta2[grep("KUL21-N",meta2$Sample),]$Group2="Normal"
meta2[grep("KUL28-N",meta2$Sample),]$Group2="Normal"
meta2[grep("KUL30-N",meta2$Sample),]$Group2="Normal"
meta2[grep("KUL31-N",meta2$Sample),]$Group2="Normal"


sce <- AddMetaData(sce, meta2$Group2,col.name = "Group2")
view(sce@meta.data)


metadata <- sce@meta.data
metadata

write.csv(metadata,file = "metadata临时文件.csv")

metadata2 <- read.csv(file = "metadata临时文件.csv",row.names = 1)

metadata2 -> sce@meta.data
view(sce@meta.data)


save(sce,file = 'sce_qualified.RData')

###Quality control completed!

#############Dimensionality reduction analysis.###############################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
#加载sce
load('sce_qualified.RData')
#PCA降维，选择合适的拐点
sce <- RunPCA(sce, features = VariableFeatures(sce)) 
dimplot1 <- DimPlot(sce, reduction = "pca") 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca")  ##根据ElbowPlot确定PC数量
sc_pca <- dimplot1+elbowplot1
sc_pca
ggsave(filename = '1.sc_pca.pdf',plot = sc_pca,he=4.5,wi=9)


pbmc <- JackStraw(object = sce, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:20,reduction = "pca",ymax = 1)
ggsave(filename = "2.JackStrawPlot.pdf",w=10,h=9)



pdf("3.可视化前2个PC的top20个基因.pdf",h =15, w=15)
VizDimLoadings(sce, dims = 1:2, nfeatures = 20, reduction = "pca")
dev.off()

#前20个PC
pdf("3.前20个PC.pdf",h =15, w=15)
DimHeatmap(sce, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()


Dims <- 40
Resolution <- 0.3

write.table(Resolution,file = "4.Resolution.txt",col=F,row=F)

#########################################################
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#颜色
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#CCCCFF","#000000","#7B68EE","#9400D3","#A0522D","#800080","#D2B48C","#D2691E",
            "#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A",
            "#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513",
            "#DEB887")
length(table(sce@active.ident))
mycolor = allcolour[1:length(table(sce@active.ident))]


cluster.frequency.table <- sce@meta.data %>%
  dplyr::count(seurat_clusters) %>%
  dplyr::mutate(freq = n / sum(n)*100) %>%
  ungroup()%>%as.data.frame()

cluster.frequency.table

pie(cluster.frequency.table$n, labels=round(cluster.frequency.table$freq,2),radius=1.0, main = "Percentage of Cluster", col=mycolor)   
legend("right",legend=unique(cluster.frequency.table$seurat_clusters),bty="n",fill=mycolor)


cluster.frequency.sample=data.frame()
for (i in as.character(unique(sce@meta.data$Group1))){
  data1<-sce@meta.data[which(sce@meta.data$Group1==i),]
  dat1 <- data1 %>%
    dplyr::group_by(Group1) %>%
    dplyr::count(seurat_clusters) %>%
    dplyr::mutate(freq = n / sum(n)*100) %>%
    ungroup()%>%as.data.frame()
  cluster.frequency.sample=rbind(cluster.frequency.sample,dat1)
}
head(cluster.frequency.sample)
cluster.freq.sample<-tidyr::spread(data=cluster.frequency.sample[,c("Group1","seurat_clusters","freq")],
                                   key=Group1, value=freq)
cluster.freq.sample[is.na(cluster.freq.sample)]<-0
head(cluster.freq.sample)


cluster.freq<-ggplot(data=cluster.frequency.sample, mapping=aes(x=Group1,y=freq,fill=seurat_clusters))+
  geom_bar(stat='identity',width=0.9)+coord_polar(theta="y",start = 0)+
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())+
  scale_fill_manual(values=mycolor)
cluster.freq
pdf('5.cluster_freq.pdf',he=7,wi=9)
cluster.freq
dev.off()
write.csv(cluster.frequency.sample,file ="5.cluster.frequency.csv")

### UMAP
sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")

#可视化
sc_umap = DimPlot(sce,cols=mycolor,
                  reduction="umap",
                  #reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_umap
ggsave('6.sc_umap_cluster.pdf',sc_umap,he=7,wi=7)


sc_umap_group1 = DimPlot(sce,cols=mycolor,group.by='Sample',
                         reduction="umap",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 


sc_umap_group1
ggsave('6.sc_umap_sample.pdf',sc_umap_group1,he=7,wi=7)

sc_umap_group2 = DimPlot(sce,cols=mycolor,group.by='Group1',
                         reduction="umap",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 


sc_umap_group2
ggsave('6.sc_umap_group.pdf',sc_umap_group2,he=7,wi=7)


pdf("6.CellCluster-umapPlot_SamGroupPC.pdf",width = 30,height = 6)
DimPlot(object = sce, cols=mycolor,
        split.by ="Sample", 
        pt.size=0.5,reduction = "umap")
dev.off()


pdf("6.CellCluster-UmapPlot_splitbyGroup.pdf",width = 13,height = 6)
DimPlot(object = sce, 
        split.by ="Group1", cols=mycolor,
        pt.size=0.5,reduction = "umap")
dev.off()


##########################################################################

###tsne 降维
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)#1000次迭代


#可视化

sc_tsne = DimPlot(sce,cols=mycolor,
                  #reduction="umap",
                  reduction="tsne",
                  label = "F", 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 
sc_tsne
ggsave('7.sc_tsne_cluster.pdf',sc_tsne,he=3,wi=4)




sc_tsne_group1 = DimPlot(sce,cols=mycolor,group.by='Sample',
                         reduction="tsne",
                         label = "T", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 


sc_tsne_group1
ggsave('7.sc_tsne_sample.pdf',sc_tsne_group1,he=3,wi=4)

sc_tsne_group2 = DimPlot(sce,cols=mycolor,group.by='Group1',
                         reduction="tsne",
                         label = "F", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  ) 


sc_tsne_group2
ggsave('7.sc_tsne_group.pdf',sc_tsne_group2,he=3,wi=4.5)


pdf("7.CellCluster-TSNEPlot_SamGroupPC.pdf",width = 30,height = 6)
DimPlot(object = sce, cols=mycolor,
        split.by ="Sample", 
        pt.size=0.5,reduction = "tsne")
dev.off()


pdf("7.CellCluster-TSNEPlot_splitbyGroup.pdf",width = 10,height = 3.7)
DimPlot(object = sce, cols=mycolor,
        split.by ="Group1", 
        pt.size=0.5,reduction = "tsne")
dev.off()



table(sce@meta.data$Sample)

b <- table(sce@meta.data$Sample,
           sce@meta.data$seurat_clusters)
b


write.csv(b, file = "7.cell.num.cluster.csv")


meta_data <- sce@meta.data
plot_data <- data.frame(table(meta_data$orig.ident,meta_data$seurat_clusters))
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)

pdf("8.PC_Cell_Chart.pdf",width = 5,height = 4)

ggplot(plot_data,aes(x = Var1,y = Percentage,fill = Var2)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + labs(fill = "Cluster")

dev.off()
rm(meta_data)


#marker基因的筛选################


#寻找差异基因时的差异倍数
Logfc = 0.3
#差异基因时最小的表达比例，默认是wilcox秩和检验
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, 
                              min.pct = Minpct,only.pos = T, test.use = "wilcox")
#保存全部Marker的结果便于手工EXCEL筛选
write.table(sce.markers,
            file="8.total_marker_genes_tsnePC.txt",
            sep="\t",quote = F,row.names = F)
length(unique(sce.markers$gene))#查看一下筛选前的基因数量


sce.markers["8.pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene)) #查看一下筛选后的基因数
head(sce.markers)

write.table(sce.markers,'8.scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')


marker.sig <- sce.markers %>% 
  mutate(Ratio = round(pct.1/pct.2,3)) %>%
  filter(p_val_adj <= 0.05) 
  


for(cluster_id in unique(marker.sig$cluster)){
  # cluster.markers <- FindMarkers(experiment.aggregate, ident.1 = cluster, min.pct = 0.3)
  # cluster.markers <- as.data.frame(cluster.markers) %>% 
  #   mutate(Gene = rownames(cluster.markers))
  cl4.genes <- marker.sig %>% 
    filter(cluster == cluster_id) %>%
    arrange(desc(avg_log2FC))
  cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),4),"gene"]
  
  #VlnPlot
  pvn <- VlnPlot(sce, features = cl4.genes,ncol = 2)
  pdf(paste0("9.MarkerGene-VlnPlot_cluster",cluster_id,"_tsne_","PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
  
  #Feather plot 
  pvn <- FeaturePlot(sce,features=cl4.genes,ncol = 2)
  pdf(paste0("9.MarkerGene-FeaturePlot_cluster",cluster_id,"_tsne_","PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
  
  #RidgePlot
  pvn<-RidgePlot(sce, features = cl4.genes, ncol = 2)
  pdf(paste0("9.MarkerGene-RidgePlot_cluster",cluster_id,"_tsne_","PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
}
rm(cl4.genes,cluster_id,pvn)


### 选择前5个marker基因
Top5 <- sce.markers %>% 
  group_by(cluster) %>% 
  slice_max(n =5, order_by = avg_log2FC)  
Top5 <- unique(Top5$gene)

sc_marker_dotplot <- DotPlot(object = sce, 
                             features = Top5,
                             cols=c("blue", "red"),
                             scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) 


#top-marker基因dotplot
sc_marker_dotplot
ggsave(filename = '10.sc_marker_dotplot.pdf',
       plot = sc_marker_dotplot,
       height = 9,width = 25)

#dotplot换个颜色
pdf("10.MarkerGene-DotPlot_all_cluster_tsne_PC.pdf",width = 20,height = 6)
DotPlot(sce, features = Top5)+
  RotatedAxis()
dev.off()



#热图展示
library(viridisLite)
sc_marker_heatmap<- DoHeatmap(object = sce,
                              features = Top5,
                              group.colors = mycolor,
                              label = F) + 
  ggtitle("Top 5 Marker Genes") + 
  theme(plot.title = element_text(hjust = 0.5)) 
sc_marker_heatmap
ggsave(filename = '10.sc_marker_heatmap.pdf',
       plot = sc_marker_heatmap,
       width = 12,height = 12)


##保存数据后续分析使用
save(sce,file = 'sce2.RData')

################Oxidative gene set scoring analysis.#################################################################################################

# BiocManager::install("GSVA")
rm(list = ls())
library(Seurat)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
#library(IOBR)
##创建gmt文件转list函数
# gmt2list <- function(gmtfile){
#   sets <- as.list(read_lines(gmtfile))
#   for(i in 1:length(sets)){
#     tmp = str_split(sets[[i]], '\t')
#     n = length(tmp[[1]])
#     names(sets)[i] = tmp[[1]][1]
#     sets[[i]] = tmp[[1]][3:n]
#     rm(tmp, n)
#   }
#   return(sets)
# }
#读取基因集数据库
# s.sets = gmt2list("h.all.v7.4.symbols.gmt")
load("sce3_已手工注释.RData")

gmtfile <- "h.all.v7.5.1.symbols.gmt"

hallmark <- read.gmt(gmtfile)

#hallmark$term <- gsub('HALLMARK_','',hallmark$term)

hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

DimPlot(sce)
#提取单细胞表达矩阵
expr <- as.matrix(sce@assays$RNA@counts)
meta <- sce@meta.data

##ssGSEA分析
es.matrix = gsva(expr, 
                 hallmark.list, 
                 kcdf="Poisson",
                 method="ssgsea", 
                 abs.ranking=T ,parallel.sz=3)


##结果可视化
library(pheatmap)
library(patchwork)

#绘制热图

pheatmap(es.matrix, show_rownames=1, show_colnames=0, annotation_col=meta,
         fontsize_row=5, width=15, height=12)


#挑选感兴趣的基因集（通路）绘制featureplot
##这里也可以将TCGA数据库得到的signature评分映射上
es <- data.frame(t(es.matrix),stringsAsFactors=F)
sce <- AddMetaData(sce, es)
head(sce@meta.data)

save(sce,file = "sce.GSEA_score_Hallmark.Rdata")
#load("sce_GSEA_score.Rdata")
write.csv(sce@meta.data,file = "meta.data.GSEA.results.csv")

pathway="HALLMARK_PI3K_AKT_MTOR_SIGNALING"

FeaturePlot(sce,features = pathway,label = T,reduction = "tsne")+
  scale_color_gradient2(high = "red",low = "blue")

ggsave(paste0(pathway,"1.pdf"),w=6,h=5)


FeaturePlot(sce,features = pathway,label = T,reduction = "tsne")+
  scale_color_gradient2(high = "red",low = "blue")

ggsave(paste0(pathway,"2.pdf"),w=6,h=5)

##分组作图
FeaturePlot(sce,features = pathway,label = F,reduction = "tsne", split.by = "Group1",)
ggsave(paste0(pathway,"3.pdf"),w=12,h=5)


VlnPlot(sce,features = pathway,group.by = "Group1")
ggsave(paste0(pathway,"4.pdf"),w=6,h=5)

VlnPlot(sce,features = pathway,group.by = "Group2")
ggsave(paste0(pathway,"5.pdf"),w=6,h=5)

RidgePlot(sce,features = pathway)
ggsave(paste0(pathway,"6.pdf"),w=9,h=5)
##结束

########cellchat analysis########################################################################################################

rm(list=ls())
#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) ###set the compute core
options(future.globals.maxSize = 60000 * 1024^2)
getwd()
library(stringr)

#memory.limit()
#BiocManager::install("igraph",force = T)

#setwd("~/Documents_PC/scRNA-seq/Data/")
####prepare the data for cellchat analysis####
load("./sce3_cell.anno.RData")
#View(sce@meta.data)

table(Idents(sce))

sce<-sce[,Idents(sce)!="Doublets"]

sce@meta.data[1:5,]

table(sce$Group1)
levels(Idents(sce))

normaldata<-sce[,sce$Group1=="Border region"]
tumordata<-sce[,sce$Group1=="Core region"]

normaldata<-normaldata[,1:8368]  
table(Idents(normaldata))

tumordata<-tumordata[,1:7340]
table(Idents(tumordata))
###perform the cellchat analysis##

#devtools::install_github("sqjin/CellChat")
library(CellChat)

####normal组构建cellchart对象
normal.input <- GetAssayData(normaldata, assay = "RNA", slot = "data") # normalized data matrix

labels <- factor(normaldata$CellType_anno,levels=levels(Idents(sce)))
labels

meta <- data.frame(group = labels, row.names = rownames(normaldata@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_")
meta$group <- factor(meta$group)

cellchat_normal <- createCellChat(object = normal.input, meta = meta, group.by = "group")
saveRDS(cellchat_normal,file="normalcellchat_1.rds")

###tumor组代码同上
tumor.input <- GetAssayData(tumordata, assay = "RNA", slot = "data") # normalized data matrix
labels <- factor(tumordata$CellType_anno,levels=levels(Idents(sce)))
meta <- data.frame(group = labels, row.names = rownames(tumordata@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_")
meta$group <- factor(meta$group)
cellchat_tumor <- createCellChat(object = tumor.input, meta = meta, group.by = "group")
saveRDS(cellchat_tumor,file="tumorcellchat_1.rds")



####prepare the deconvolution analysis gene signature####


rm(list=ls()) 

#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(stringr)
#devtools::install_github("sqjin/CellChat")
library(CellChat)
#setwd("~/Documents_PC/scRNA-seq/Data")
normal_cellchat<-readRDS(file="normalcellchat_1.rds") ##as tumor border
tumor_cellchat<-readRDS(file="normalcellchat_1.rds")  ##as tumor core

normal_cellchat@DB <- CellChatDB.human  
gc()
normal_cellchat <- subsetData(normal_cellchat)  
# subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 2) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore

##atumor border
normal_cellchat <- identifyOverExpressedGenes(normal_cellchat) #找细胞分组差异高表达基因
normal_cellchat <- identifyOverExpressedInteractions(normal_cellchat)    
normal_cellchat <- projectData(normal_cellchat, PPI.human)  ###需要网络
normal_cellchat <- computeCommunProb(normal_cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
normal_cellchat <- filterCommunication(normal_cellchat, min.cells = 1)
normal_cellchat <- computeCommunProbPathway(normal_cellchat)
normal_cellchat <- netAnalysis_computeCentrality(normal_cellchat, slot.name = "netP")


saveRDS(normal_cellchat,file="normalcellchat_1.rds")

#normal_cellchat <- readRDS("normalcellchat_1.rds")

normal_cellchat <- aggregateNet(normal_cellchat)

##接下来可视化看结果
groupSize <- as.numeric(table(normal_cellchat@idents))
groupSize
table(normal_cellchat@idents)


par(mfrow = c(1,2), xpd=TRUE)  ##
netVisual_circle(normal_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(normal_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off() 

####提取配受体数量
mat <- normal_cellchat@net$count

mat <- normal_cellchat@net$weight

####下面针对其中一个细胞亚型分析其对其他细胞的interaction的Ligand-receptor数量
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[1, ] <- mat[1, ]  ##先查看mat2每一行是什么细胞决定修改参数
netVisual_circle(mat2, vertex.weight = groupSize,arrow.size = 0.2, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[1])


###查看具体的信号通路
df.net <- subsetCommunication(normal_cellchat)   ##分泌型。表面型等等
#df.net_tumor <- subsetCommunication(tumor_cellchat)


levels(df.net$source)
unique(df.net$pathway_name)##便于查看有哪些通路可以用
pathways.show <- c("PDGF")   ##选择自己感兴趣的通路

#table(df.net$pathway_name)
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#dev.new()
seq(1,4)

vertex.receiver = seq(1,4) # a numeric vector. 选择1到4作为recevor
netVisual_aggregate(normal_cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")

par(mfrow = c(1,1), xpd=TRUE)  
netVisual_aggregate(normal_cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "chord")

#?netVisual_aggregate
#par(mfrow=c(1,1))

#netVisual_aggregate(normal_cellchat, signaling = pathways.show,vertex.receiver = vertex.receiver, layout = "circle")
netVisual_aggregate(normal_cellchat, signaling = pathways.show, layout = "circle")




netAnalysis_contribution(normal_cellchat, signaling = pathways.show)


####get specific  cell 
levels(normal_cellchat@idents)

netVisual_bubble(normal_cellchat, sources.use = c(1), targets.use = c(2:10), remove.isolate = FALSE)###epithelial
netVisual_bubble(normal_cellchat, sources.use = c(2), targets.use = c(8,9), remove.isolate = FALSE)##TNK
netVisual_bubble(normal_cellchat, sources.use = c(5), targets.use = c(1:4,6:10), remove.isolate = FALSE)##myeloid cell
#saveRDS(normal_cellchat,file="./Cellchat/normalcellchat_analysis.rds")



######tumor core#####
##肿瘤组的一样操作

###We have run these in the workstation###

tumor_cellchat<-readRDS(file="tumorcellchat_1.rds")
tumor_cellchat@DB <- CellChatDB.human

#gc()
tumor_cellchat <- subsetData(tumor_cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel

tumor_cellchat <- identifyOverExpressedGenes(tumor_cellchat)
tumor_cellchat <- identifyOverExpressedInteractions(tumor_cellchat)
tumor_cellchat <- projectData(tumor_cellchat, PPI.human)
tumor_cellchat <- computeCommunProb(tumor_cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
tumor_cellchat <- filterCommunication(tumor_cellchat, min.cells = 1)
tumor_cellchat <- computeCommunProbPathway(tumor_cellchat)
tumor_cellchat <- netAnalysis_computeCentrality(tumor_cellchat, slot.name = "netP") 


saveRDS(tumor_cellchat,file="tumorcellchat_1.rds")
#tumor_cellchat<-readRDS("tumorcellchat_1.rds")



tumor_cellchat <- aggregateNet(tumor_cellchat)

#par(mfrow = c(2,2), xpd=TRUE)


groupSize <- as.numeric(table(tumor_cellchat@idents))
groupSize
table(tumor_cellchat@idents)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(tumor_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(tumor_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(tumor_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(normal_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(tumor_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(normal_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


##########Fisher's test (two-sided) to calculate the P value for the difference in cell number in tumour core and border#############################################################################################
library(ggplot2)

load('sce3_cell_anno.RData')

cell_num<-as.data.frame(table(sce@active.ident,sce$Group1))
colnames(cell_num)=c('clusters','type','cell_num') 
cell_num=tidyr::spread(cell_num,type,cell_num)
cell_num

write.csv(cell_num,file = "3.cell_num.csv")

cell_num1=data.frame()
for (i in 1:nrow(cell_num)){
  print(i)
  pval=fisher.test(matrix(c(cell_num[i,3],cell_num[i,2],
                            sum(cell_num[,3])-cell_num[i,3],
                            sum(cell_num[,2])-cell_num[i,2]),
                          nrow = 2,ncol = 2)
                   ,alternative = "two.sided")$p.value
 
  cell_num1[i,'C_celltype']=cell_num[i,3]
  cell_num1[i,'C_no_celltype']=sum(cell_num[,3])-cell_num[i,3]
  cell_num1[i,'B_celltype']=cell_num[i,2]
  cell_num1[i,'B_no_celltype']=sum(cell_num[,2])-cell_num[i,2]
  cell_num1[i,'p.val']=pval
  cell_num1[i,'fc']=(cell_num[i,3]/(sum(cell_num[,3])-cell_num[i,3]))/(cell_num[i,2]/(sum(cell_num[,2])-cell_num[i,2]))
}

cell_num1$cell_name=cell_num$clusters
 # 使用FDR校正矫正P值
cell_num1$p.adj=p.adjust(cell_num1$p.val, method = "FDR")


cell_num1

a=subset(cell_num1, p.adj<0.05)
b=subset(cell_num1, p.val<0.05)


