
#Microbial α diversities###################################
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)

df <- read.csv("TCGA.CRC.micro.count.csv", row.names = 1,check.names = F)

df=as.data.frame(t(df))

Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)#spe.rich =sobs

index <- as.data.frame(cbind(Shannon, Simpson, Richness))
tdf <- t(df)
tdf<-ceiling(as.data.frame(t(df)))
#计算obs，chao，ace指数
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]#统一行名
#将obs，chao，ace指数与前面指数计算结果进行合并
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
#计算Pielou及覆盖度
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)

write.table(cbind(sample=c(rownames(index)),index),'diversity.index.txt', row.names = F, sep = '\t', quote = F)

##差异性计算及绘图
#读读入数据及分组文件
index <- read.delim('diversity.index.txt', header = T, row.names = 1)
##figure:take shannon for example
index$samples <- rownames(index)
#读入分组文件
groups <- read.delim('group.txt',header = T, stringsAsFactors = F)

head(groups)
#colnames(groups)[1:3] <- c('samples','group')#改列名

df2 <- merge(index,groups,by.x ='samples', by.y = "ID")

#绘图
#--Shannon

#s设置比较组
group=levels(factor(df2$cluster))
df2$cluster=factor(df2$cluster, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#Shannon
my_colors <- c("#4995C6", "#EE4431")

ggboxplot(df2, x="cluster", y="Shannon", fill = "cluster", 
         xlab="OSRGcluster", ylab="Shannon Index",palette = my_colors,
         #add="jitter",
         add.params = list(size=0.5),
         title = "")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",label="p.signif") + theme(legend.position = "none") +labs(title = NULL)

ggsave(filename = "Shannon Index.pdf", w=3.5,h=3)


#--Simpson
#Simpson

ggboxplot(df2, x="cluster", y="Simpson", fill = "cluster", 
          xlab="OSRGcluster", ylab="Simpson Index",palette = my_colors,
          #add="jitter",
          add.params = list(size=0.5),
          title = "")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",label="p.signif") + theme(legend.position = "none") +labs(title = NULL)

ggsave(filename = "Simpson Index.pdf", w=3.5,h=3)



#Richness
#Richness

ggboxplot(df2, x="cluster", y="Richness", fill = "cluster", 
          xlab="OSRGcluster", ylab="Richness Index",palette = my_colors,
          #add="jitter",
          add.params = list(size=0.5),
          title = "")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",label="p.signif") + theme(legend.position = "none") +labs(title = NULL)

ggsave(filename = "Richness Index.pdf", w=3.5,h=3)



#--Pielou
#Pielou

ggboxplot(df2, x="cluster", y="Pielou", fill = "cluster", 
          xlab="OSRGcluster", ylab="Pielou Index",palette = my_colors,
          #add="jitter",
          add.params = list(size=0.5),
          title = "")+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",label="p.signif") + theme(legend.position = "none") +labs(title = NULL)

ggsave(filename = "Pielou Index.pdf", w=3.5,h=3)

##########.Bray-Curtis################################################################################################################

otu <- read.table("TCGA.CRC.micro.count.txt",row.names = 1,header = T,sep = "\t")
#otu <- t(otu)

#计算Bray-Crutis距离。
library(vegan)
bray <- vegdist(otu,method = "bray", binary=F)

bray <- as.matrix(bray)
write.table(bray,"bray-crutis.txt",sep = "\t")
library(pheatmap)
library(RColorBrewer)
pheatmap(bray,color = colorRampPalette(brewer.pal(7,"RdYlBu"))(100))

dune.env=read.table("group.txt", header = T, row.names = 1)
dune_pcoa <- cmdscale(bray, k=3, eig=T)

dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)

colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)

samesample=intersect(rownames(dune_pcoa_points), rownames(dune.env))
dune_pcoa_points=dune_pcoa_points[samesample,]

dune_pcoa_result <- cbind(dune_pcoa_points, dune.env)

head(dune_pcoa_result)
#画图
library(ggplot2)

ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=cluster)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=1
  ) + stat_ellipse(level=0.6) +
  theme_classic()

ggsave(filename = "PCoA.pdf",w=5,h=3)


library(ggalt)
ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=cluster, group = cluster)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=5) + 
  geom_encircle(aes(fill=cluster), alpha = 0.1, show.legend = F) +
  theme_classic() + coord_fixed(1)

# 基于bray-curtis距离进行计算
set.seed(123)

otu2=otu[rownames(dune.env),]

dune.div <- adonis2(otu2 ~ dune.env$cluster, data = dune.env, permutations = 999, method="bray")

dune.div
#把统计检验结果加到PcOA的图上。
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,3), "; P-value: ", dune.div$`Pr(>F)`)

# install.packages("ggalt")
library(ggalt)
ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=cluster, group = cluster)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title=dune_adonis) +
  geom_point(size=1) + 
  geom_encircle(aes(fill=cluster), alpha = 0.1, show.legend = F) +
  theme_classic() + coord_fixed(1)

my_colors <- c("#4995C6", "#EE4431")
ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=cluster)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title=dune_adonis) +
  geom_point(size=1
  ) + stat_ellipse(level=0.6) +
  theme_classic()+
  scale_color_manual(values = my_colors)  

ggsave(filename = "PCoA_with_pvalue.pdf",w=4.5,h=3)


############Differential analysis of microbial abundance between different OSRGclusters####################################################################################################

rt=read.csv("TCGA.CRC.micro.count.csv", row.names = 1)
rt=as.data.frame(t(rt))
head(rt)

pattern <- "g__\\w+\\b" 
matching_names <- grep(pattern, rownames(rt), value = TRUE, ignore.case = TRUE)

rt2=rt[matching_names,]

rt2$OTU=rownames(rt2)

library(stringr)
spli<- str_split_fixed(rt2$OTU, "g__", 2)
spli <- spli[,-1] #删除冗余列
rt2$OTU=spli

which(duplicated(rt2$OTU))
rownames(rt2)=rt2$OTU
rt2=rt2[,-790]

write.csv(rt2, "TCGA.CRC.micro.single.names.csv")

library(DESeq2)
library(ggplot2)
library(ggrepel)

group <- read.table(file="group.txt",sep="\t",
                    header=T,check.names=FALSE,row.names=1 )

group$cluster<-factor(group$cluster,levels = c('B','A')) #B比A

rt3=rt2[,rownames(group)]

rt3=rt3+1
#构建DESeqDataSet对象
df_dds <- DESeqDataSetFromMatrix(countData = rt3, colData = group, design = ~cluster)
#差异分析
df_dds_res <- DESeq(df_dds)
suppressMessages(df_dds_res)


df_res <- results(df_dds_res)
# 根据p-value进行重新排序
df_res = df_res[order(df_res$pvalue),]
df_res #查看结果
summary(df_res)
#合并数据
df <-  merge(as.data.frame(df_res),
             as.data.frame(counts(df_dds_res,normalize=TRUE)),
             by="row.names",sort=FALSE)

df1<-df[!is.na(df$padj),]
#数据分类——根据其中log2FoldChange、padj指标对OTU进行分类
df1$SG<-as.factor(ifelse(df1$padj<0.05&abs(df1$log2FoldChange)>=1,"Y","N"))

df1$label<-ifelse(df1$padj<0.05&abs(df1$log2FoldChange)>=2,"Y","N")
df1$label<-ifelse(df1$label == 'Y', as.character(df1$Row.names), '')

df1$SG <- factor(df1$SG, levels = c('Y', 'N'), labels = c('differences','No difference'))

p <- ggplot(df1, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = SG),alpha=0.6, size=2)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks=seq(-10,10, 2))+
  geom_vline(xintercept = c(-2, 2), lty=3,color = 'black', lwd=0.5) + 
  geom_hline(yintercept = -log10(0.05), lty=3,color = 'black', lwd=0.5) +
  scale_color_manual(values = c( 'red','grey'))+
  labs(title="Cluster B vs. Cluster A",
       x = 'log2 fold change',
       y = '-log10 pvalue')+
  geom_text_repel(aes(x = log2FoldChange,
                      y = -log10(padj),          
                      label=label),                       
                  max.overlaps = 1000,
                  size=3,
                  box.padding=unit(0.8,'lines'),
                  point.padding=unit(0.8, 'lines'),
                  segment.color='black',
                  show.legend=FALSE)
p
ggsave(filename = "volcano.pdf", w=7.5,h=5)

write.csv(df1, "all.results.csv", row.names = F)
















