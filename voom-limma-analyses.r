###################################################
# Analysis differential gene/exon expression from Chickpea RNA-seq data from Gary R et al. 
# (Sci Rep. 2016 Jan 13;6:19228. doi: 10.1038/srep19228) under salinity stress and control 
# with voom-limma, alignment with Subread and count summarise with FeatureCount

# Zhi-Ping Feng
# 2016
###################################################
library(limma)
library(edgeR)
library(reshape)
library(gridBase)
library(grid)

count.raw<-read.table("IndianSamplee-counts-ICC4598-v2.txt",header=TRUE,sep='\t')
exon.counts<-count.raw[,7:22]
exon.counts.exp<-exon.counts[rowSums(exon.counts)!=0,]
dim(exon.counts.exp)
#[1] 121851    16

gene.ID<-unlist(lapply(strsplit(rownames(exon.counts.exp),".exon"),"[",1))
tmp<-data.frame(gene.ID,exon.counts.exp)

## merge counts from multiple exons
gene.counts<-cast(melt(tmp,id="gene.ID"),gene.ID ~ variable,sum)
dim(gene.counts)
gene.idx<-rowSums(cpm(gene.counts)>1) >= 4
table(gene.idx) 
#FALSE  TRUE 
#6378 17662 
g.count.f<-gene.counts[gene.idx,]
dim(g.count.f)

save(g.count.f,file="Indina-JM-geneCounts.f.RData")
save(exon.counts.exp,file="Indina-JM-exonCounts.f.RData")

# Differential expression analysis for g.count.f
# similar to exon.counts.exp as well.

group0<-c(rep("C",2),rep("S",2))
design0<-matrix(0, 4, 2)
rownames(design0)<-group0
colnames(design0)<-c("interscept","S2C")
design0[,1]<-c(rep(1,4))
design0[,2]<-c(rep(0,2),rep(1,2))

for (i in c(2,6,10,14)) {

data<-g.count.f[,c(i:(i+3))]
d<-DGEList(counts=data,group=group0)
nf<-calcNormFactors(d$counts,method="TMM")  
### change “none” to “TMM for TMM normalization
y1<-voom(d$counts,plot=F,design0,lib.size=colSums(d$counts)*nf)

## with sample quality weight
#y2<-voomWithQualityWeights(d$counts,plot=F,replace.weights=T,design0,lib.size=colSums(d$counts)*nf) 
#aw <- arrayWeights(y2$E, design = design0)
#voom(d$counts,plot=T,design0,weight=aw,lib.size=colSums(d$counts)*nf)
#plotMDS(y2,labels=colnames(g.count.f)[(i:(i+3))],col=cols[(i-1):(i+2)],main="weighted @TMM")

#names(aw)<-colnames(g.count.f)[(i:(i+3))]
#barplot(aw,col=cols[(i-1):(i+2)],las=2,main="weight")
#abline(h = 1, col = 2, lty = 2)

fit <- lmFit(y1,design0)
# with wight
# fit <- lmFit(y2,design0)

fit2 <- eBayes(fit,trend=T,robust=T)
options(digits=3)

top<-topTable(fit2,coef=2,n=nrow(fit2),genelist =gene.ID)
print (summary(decideTests(fit2,adjust.method="BH",p.value=0.05,lfc=1.5)))

cat(paste("FDR<0.05 ",length(rownames(top[top$adj.P.Val<0.05 & abs(top$logFC)>1.5,]))),"\n")
cat(paste("P.value<0.05 & |logFC|>1.5 ",length(rownames(top[top$P.Value<0.05 & abs(top$logFC)>1.5,]))),"\n")
n<-length(rownames(top[top$adj.P.Val<0.05 & abs(top$logFC)>1.5,]))
write.csv(head(top,n=n),file=filenames.exon[j])
write.csv(ann.new[ann.new$Gene.ID%in%head(top,n=n)$ID,],file=paste("annotation",filenames.exon[j],sep=""))
#head(top,n=n)$ID
j=j+1
}
