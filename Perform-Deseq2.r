source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library('DESeq2')

c<-cbind(c(1:6),c(7:12))

par(mfrow=c(1,2),mar=c(7,5,5,3))
for (i in 1:2){
  
countData<-as.matrix(d1.add[,c[,i]])

ddsHTSeq<-DESeqDataSetFromMatrix(countData,DataFrame(condition),~condition)

dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)
results(dds)
summary(res)
sum(res$padj < 0.01, na.rm=TRUE)
deseq<-res[complete.cases(res),]
de<-deseq[deseq$padj< 0.01,]
if(i==1){
  write.csv(res,file="deseq283.csv")
  smoothScatter(log2(deseq$baseMean+1),deseq$log2FoldChange,main="(a)",cex=0.4,col="blue",pch=16) 
  points(log2(de$baseMean+1),de$log2FoldChange,cex=0.4,col="magenta",pch=16) 
  abline(h=0,lwd=2,col="red")
  res1<-res
  }

else if(i==2) {
  write.csv(res,file="deseq8261.csv")
 # plot(log2(res$baseMean),res$log2FoldChange,ylim=c(-2,2),main="(b)",cex=0.4,col="blue",pch=16 )
  smoothScatter(log2(deseq$baseMean+1),deseq$log2FoldChange,main="(b)",cex=0.4,col="blue",pch=16) 
  points(log2(de$baseMean+1),de$log2FoldChange,cex=0.4,col="magenta",pch=16) 
  abline(h=0,lwd=2,col="red")
  res2<-res
}
}


plot(res1$baseMean+1, -log10(res1$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3),main="(a) ICC283")
plot(res2$baseMean+1, -log10(res2$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3),main="(b)  ICC8261")
