
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
#setwd("D:\\biowolf\\m6aLnc\\26.model")      
rt=read.table("uniSigExp_TCGA.txt", header=T, sep="\t", check.names=F, row.names=1)     #��ȡ�����ļ�
rt$futime[rt$futime<=0]=0.003  

##ȫ�����������������
inTrain<-createDataPartition(y=rt[,3], p=0.8, list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train), train)
  testOut=cbind(id=row.names(test), test)

#���������
 write.table(trainOut,file=paste0("./",sam.name,"/���13.train.data.txt"),sep="\t",quote=F,row.names=F)
 write.table(testOut,file=paste0("./",sam.name,"/���13.test.data.txt"),sep="\t",quote=F,row.names=F)

#��ȡ�ļ������������ļ���������
rt=read.table(paste0("./",sam.name,"/���13.train.data.txt"), header=T, sep="\t", check.names=F, row.names=1)
test=read.table(paste0("./",sam.name,"/���13.test.data.txt"), header=T, sep="\t", check.names=F, row.names=1)
sameGene=intersect(colnames(rt)[3:ncol(rt)], colnames(test)[3:ncol(test)])
rt=rt[,c("futime","fustat",sameGene)]
#rt$futime=rt$futime/365
#��ȡGEO���ݿ��ļ������������ļ���������
GEO1=read.table("GSE17536.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
sameGene=intersect(colnames(rt)[3:ncol(rt)], colnames(GEO1)[3:ncol(GEO1)])
GEO1=GEO1[,c("futime","fustat",sameGene)]
GEO1$futime=GEO1$futime/365 
#GEO[,3:ncol(GEO)]=log2(GEO[,3:ncol(GEO)]+1) 

#��ȡGEO���ݿ��ļ������������ļ���������
GEO2=read.table("GSE39582.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
sameGene2=intersect(colnames(rt)[3:ncol(rt)], colnames(GEO2)[3:ncol(GEO2)])
GEO2=GEO2[,c("futime","fustat",sameGene2)]
GEO2$futime=GEO2$futime/365 
#GEO[,3:ncol(GEO)]=log2(GEO[,3:ncol(GEO)]+1) 

#��ȡȫ�������ϲ��������ļ������������ļ���������
allsample=read.table("All.sample.diff&time.txt", header=T, sep="\t", check.names=F, row.names=1)
sameGene3=intersect(colnames(rt)[3:ncol(rt)], colnames(allsample)[3:ncol(allsample)])
allsample=allsample[,c("futime","fustat",sameGene3)]
allsample$futime=allsample$futime/365  
#GEO[,3:ncol(GEO)]=log2(GEO[,3:ncol(GEO)]+1) 


#��TCGA��ѵ���������ݹ���ģ��
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime, rt$fustat))
fit=glmnet(x, y, family="cox", maxit=1000)
#lasso�ع�ͼ��
pdf(paste0("./",sam.name,"/���13.lasso.lambda.pdf"))
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
#������֤ͼ��
cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
pdf(paste0("./",sam.name,"/���13.lasso.cvfit.pdf"))
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

#�����ػ���ϵ��
coef=coef(fit, s=cvfit$lambda.1se)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
n <- length(geneCoef[,1])

write.table(geneCoef, file=paste0("./",sam.name,"/���13.lasso.geneCoef_",n,".txt"), sep="\t", quote=F, row.names=F)


#���train�����ֵ
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab1=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab1),outTab1),file=paste0("./",sam.name,"/���13.risk_TCGA.train.txt"),sep="\t",quote=F,row.names=F)
trainRiskOut <- cbind(id=rownames(outTab1),outTab1)

#���test�����ֵ
#rt=read.table(icgcFile, header=T, sep="\t", check.names=F, row.names=1)
#rt$futime=rt$futime/365
#test[,3:ncol(test)]=log2(test[,3:ncol(test)]+1)  #ȡlog
testFinalGeneExp=test[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab2=cbind(test[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab2),outTab2),file=paste0("./",sam.name,"/���13.risk_TCGA.Test.txt"),sep="\t",quote=F,row.names=F)
testRiskOut <- cbind(id=rownames(outTab2),outTab2)

#���ȫ�����ߵķ�������
allRiskOut=rbind(trainRiskOut, testRiskOut)
write.table(allRiskOut,file=paste0("./",sam.name,"/���13.allRisk_TCGA.txt"),sep="\t",quote=F,row.names=F)

#���GEO�����ֵ
#rt=read.table(icgcFile, header=T, sep="\t", check.names=F, row.names=1)
#rt$futime=rt$futime/365
#rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
GEOFinalGeneExp=GEO1[,lassoGene]
GEOScore=apply(GEOFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(GEOScore>median(trainScore),"high","low"))
outTab=cbind(GEO1[,outCol],riskScore=as.vector(GEOScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file=paste0("./",sam.name,"/���13.risk.GSE17536.txt"),sep="\t",quote=F,row.names=F)

##################################################
#rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
GEOFinalGeneExp=GEO2[,lassoGene]
GEOScore=apply(GEOFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(GEOScore>median(trainScore),"high","low"))
outTab=cbind(GEO2[,outCol],riskScore=as.vector(GEOScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file=paste0("./",sam.name,"/���13.risk.GSE39582.txt"),sep="\t",quote=F,row.names=F)

##ȫ����������
#rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
GEOFinalGeneExp=allsample[,lassoGene]
GEOScore=apply(GEOFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(GEOScore>median(trainScore),"high","low"))
outTab=cbind(allsample[,outCol],riskScore=as.vector(GEOScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file=paste0("./",sam.name,"/���13.risk.ALL.txt"),sep="\t",quote=F,row.names=F)


