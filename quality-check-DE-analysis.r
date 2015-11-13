##############################################################################
#  Pipeline for data quality control and differential expression analyses
#  based on RNA-seq 
#  Supplementary codes to the manuscript "Avoiding pitfalls in differential 
#  gene expression analyses based on Plant RNA-seq "
#
#  Dr Zhi-Ping Feng  Nov 2015
############################################################################## 

## install limma and edgeR if you run this first time
source("http://bioconductor.org/biocLite.R")
biocLite()

library(limma)
library(edgeR)

## set a working directory
pwd<-getwd()
setwd(pwd)

## load example data
lname<-load("Chickpea-roots-30Days.RData",sep='')
lname
## check the data size
dim(Data)
# [1] 52755    18

## check the sizes of libraries
lib.size<-colSums(Data)
range(lib.size)
# [1]  1040689 17113372


## rename the data
d1<-Data

## define color for plots
###  for 18 read libraries
cols<-c(rep("darkgrey",3),rep("red",6),rep("darkgreen",5),rep("blue",4))

### for 12 merged libraries
col3<-cols[c(1:3,4,6,8,10,12,14,15,17:18)]

### merge the technical replicates
d1.add<-cbind(d1[,1:3],d1[,4]+d1[,5],d1[,6]+d1[,7],d1[,8]+d1[,9],d1[,10]+d1[,11],d1[,12]+
              d1[,13],d1[,14],d1[,15]+d1[,16],d1[,17],d1[,18])

colnames(d1.add)<-c(colnames(d1)[1:3],paste(colnames(d1)[c(4,6,8,10,12)],"m",sep=""),
                    colnames(d1)[14],paste(colnames(d1)[15],"m",sep=""),
                    colnames(d1)[17:18])

## Figure S1 barplots show the library size
par(mfrow=c(1,2),mar=c(9,5,5,1))
barplot(colSums(d1),main="(a)",col=cols,las=2)
barplot(colSums(d1.add),main="(b)",col=col3,las=2)


## correlations for each pair of the technical replicates

m<-cbind(c(4,5),c(6,7),c(8,9),c(10,11),c(12,13),c(15,16))
par(mfrow=c(2,3),mar=c(6,5,5,3))

for (i in 1:6){
  cc<-cor(log2(pmax(as.matrix(d1[,m[1,i]]),1)),log2(pmax(as.matrix(d1[,m[2,i]]),1)),method = "spearman")
  
  plot(log2(pmax(as.matrix(d1[,m[1,i]]),1)),log2(pmax(as.matrix(d1[,m[2,i]]),1)),
       xlab=paste(colnames(d1)[m[1,i]]),ylab=paste(colnames(d1)[m[2,i]]),
       
       pch=16,col="blue",cex=0.6,main=paste("Spearman cc=",format(cc,digits=2)))
}

## MA-plot for each pair of the technical replicates
for (i in 1:6) {
  x<-(log2(pmax(as.matrix(d1[,m[1,i]]),1))+log2(pmax(as.matrix(d1[,m[2,i]]),1)))/2
  y<-(log2(pmax(as.matrix(d1[,m[1,i]]),1))-log2(pmax(as.matrix(d1[,m[2,i]]),1)))
  
  plot(x,y,
       main=paste(colnames(d1) [m[1,i]],"vs",colnames(d1)[m[2,i]]),
                  xlab="A",ylab="M" ,pch=16,col="blue",cex=0.6)
       
       points(x[abs(y)>=5],y[abs(y)>=5],col="red")
       abline(h=0,col="red")
}

## correlations for biological replicates
m<-cbind(c(1,2),c(1,3),c(4,5),c(4,6),c(7,8),c(7,9),c(10,11),c(10,12))
par(mfrow=c(2,4),mar=c(6,5,5,1))
for (i in c(1:8)){
  cc<-cor(log2(pmax(as.matrix(d1.add[,m[1,i]]),1)),log2(pmax(as.matrix(d1.add[,m[2,i]]),1)),method = "spearman")
  plot(log2(pmax(as.matrix(d1.add[,m[1,i]]),1)),log2(pmax(as.matrix(d1.add[,m[2,i]]),1)),
       xlab=paste("",colnames(d1.add)[m[1,i]]),ylab=paste("",colnames(d1.add)[m[2,i]]),
       pch=16,col="darkgreen",cex=0.6,main=paste("Spearman cc:",format(cc,digits=2)))
}


## MA-plots for biological replicates
par(mfrow=c(2,4),mar=c(6,5,5,1))
for (i in c(1:8)){
  x<-(log2(pmax(as.matrix(d1.add[m[1,i]]),1))+log2(pmax(as.matrix(d1.add[,m[2,i]]),1)))/2
  y<-(log2(pmax(as.matrix(d1.add[,m[1,i]]),1))-log2(pmax(as.matrix(d1.add[,m[2,i]]),1)))
  
  plot(x,y,
       main=paste(colnames(d1.add)[m[1,i]],"vs",colnames(d1.add)[m[2,i]]),
       xlab="A",ylab="M" ,
       pch=16,col="darkgreen",cex=0.6)
  points(x[abs(y)>=5],y[abs(y)>=5],col="red")
  abline(h=0,col="red")
}


## remove the genes weakly expressed

d1.idx<-rowSums(cpm(d1)>1) >= 3
    table(d1.idx) 
                
    table(d1.idx)/dim(d1)[1] 
    d1.add.idx<-rowSums(cpm(d1.add)>1) >= 3 
    table(d1.add.idx) 
    
    d1.f<-d1[d1.idx,]
    d1.add.f<-d1.add[d1.add.idx,]
    
## check the dimension of the data    
    dim(d1.f)
#   [1] 17291    18
    dim(d1.add.f)
#   [1] 16683    12    
    
## check the expression after filtering the lowly expressed genes    
### boxplots 
    par(mfrow=c(1,2),mar=c(9,5,5,1))
    I<-log2(pmax(as.matrix(d1.f),1))
    boxplot(I,main="(a)",las=2,col=cols)
    I<-log2(pmax(as.matrix(d1.add.f),1))
    boxplot(I,main="(b)",las=2,col=col3)
### RLE plots   
    par(mfrow=c(1,2),mar=c(7,5,5,3))
    I<-log2(pmax(as.matrix(d1.f),1))
    makeRLEplot(I,main="(a)",col=cols)
    I<-log2(pmax(as.matrix(d1.add.f),1))
    makeRLEplot(I,main="(b)",col=col3)
 
### prepare data for differential expression analysis       
    lab1<-colnames(d1)  
    lab2<-as.factor(paste(unlist(lapply(strsplit(colnames(d1.add.f),"_"),"[",1)),"_",unlist(lapply(strsplit(colnames(d1.add.f),"_"),"[",2)),sep=""))
    
##   sample types
    lab2
# [1] 30C_283  30C_283  30C_283  30S_283  30S_283  30S_283 
# [7] 30C_8261 30C_8261 30C_8261 30S_8261 30S_8261 30S_8261


### define a design matrix (factor of interest), indicates which RNA samples 
### want to been applied. Each row of the design matrix corresponds to a 
### library of the experiment and each column corresponds to a coefficient 
### that is used to describe the RNA sources in the experiment
    
    group0<-as.factor(lab2)
    design0<-matrix(0, 12, 5)
    rownames(design0)<-lab2
    design0[,1]<-c(rep(1,12))
    
    for (i in 1:4){
      design0[rownames(design0)%in%levels(group0)[i],i+1]<-1
    }
    colnames(design0)<-c("interscept",levels(group0))
 
 ## check desing0
    design0
#         interscept 30C_283 30C_8261 30S_283 30S_8261
#30C_283           1       1        0       0        0
#30C_283           1       1        0       0        0
#30C_283           1       1        0       0        0
#30S_283           1       0        0       1        0
#30S_283           1       0        0       1        0
#30S_283           1       0        0       1        0
#30C_8261          1       0        1       0        0
#30C_8261          1       0        1       0        0
#30C_8261          1       0        1       0        0
#30S_8261          1       0        0       0        1
#30S_8261          1       0        0       0        1
#30S_8261          1       0        0       0        1
    
## normalization and estimate mean variance trend
### no weight
    d<-DGEList(counts=d1.add.f,group=group0)
    nf<-calcNormFactors(d$counts,method="none")  
    ### change “none” to “TMM for TMM normalization
    y1<-voom(d$counts,plot=T,design0[,-1],lib.size=colSums(d$counts)*nf)
    
## Figure 1 boxplot and  RLE plot after cpm normalization
    
    par(mfrow=c(1,2),mar=c(7,5,5,3))
    ### boxplots
    boxplot(y1$E,main="(a)  Boxplot",col=col3,las=2)
    ### RLE plots
    makeRLEplot(y1$E,main="(b)  RLE plot",col=col3)
 
  ## Figure S8 boxplot and  RLE plot without outliers after cpm normalization
    
    par(mfrow=c(1,2),mar=c(7,5,5,3))
    ### boxplots
    boxplot(y1$E[,-c(6,8)],main="(a)  Boxplot",col=col3[-c(6,8)],las=2)
    ### RLE plots
    makeRLEplot(y1$E[,-c(6,8)],main="(b)  RLE plot",col=col3[-c(6,8)])
 
    d0<-DGEList(counts=d1.add.f[,-c(6,8)],group=group0[-c(6,8)])
	nf<-calcNormFactors(d0$counts,method="none")  
	y0<-voom(d0$counts,plot=T,design0[-c(6,8),-1],lib.size=colSums(d0$counts)*nf)
	
	#### define a function for PCA plot
    PCAplot<-function(I,main=main,col=col,lab=lab){
      pca<-prcomp(t(I),cor=T,scale. = T)
      ev<-pca$sdev^2
      ev.pc<-as.numeric(format(100*ev/sum(ev),digits=2))
      xlim<-range(pca$x[,1])*2
      ylim<-range(pca$x[,2])*2
      plot(pca$x[,1:2],main=main,type='n',xlim=xlim,ylim=ylim,
           xlab=paste("PC1: ",ev.pc[1],"%",sep=''),ylab=paste("PC2: ",ev.pc[2],"%",sep=''),cex.lab=1.5,col=col)
      text(pca$x[,1:2],labels=lab,col=col,cex=0.8,xlim=xlim,ylim=ylim)
      return(ev.pc[1:10])
    }
    
## draw a PCA plot
	PCAplot.all(y0$E,main="",col=col3[-c(6,8)],lab=colnames(d1.add.f)[-c(6,8)])

    
### Undereight the observations from variable samples
    d<-DGEList(counts=d1.add.f,group=group0)
    y2<-voomWithQualityWeights(d$counts,plot=F,replace.weights=F,design0[,-2],lib.size=colSums(d$counts)*nf) 
    
## Figure 2. 
### PCA plots before and after weight the observations from variable samples
## plot PCA    
    PCAplot.all(y1$E,main="(c)  PCA plot",col=c(rep("black",3),col3[4:12]),lab=colnames(d1.add.f))
    PCAplot.all(y2,main="(d)  PCA plot",col=c(rep("black",3),col3[4:12]),lab=colnames(d1.add.f))
    
    
## Figure S10. scatter plots of the mean-variance trend and weight of variable samples
    d<-DGEList(counts=d1.add.f,group=group0)
    nf<-calcNormFactors(d$counts,method="none")  
###### perform voom-cpm normalization; change “none” to “TMM” for TMM normalization
    par(mfrow=c(1,2),mar=c(7,5,5,3))
     y1<-voom(d$counts,plot=T,design0[,-1],lib.size=colSums(d$counts)*nf)
     y2<-voomWithQualityWeights(d$counts,plot=F,replace.weights=T,design0[,-1],lib.size=colSums(d$counts)*nf) 
     voom(d$counts,plot=T,design0[,-1],weight=aw,lib.size=colSums(d$counts)*nf)

## Figure S11     
     aw <- arrayWeights(y2$E, design = design0[,-1])
     names(aw)<-colnames(d1.add.f)
     barplot(aw,col=col3,las=2)
     abline(h = 1, col = 2, lty = 2)

#####################################################################################################     
## Differential gene expression analysis  by voom-limma without underweighting variable libraries
#####################################################################################################
 
  d<-DGEList(counts=d1.add.f,group=group0)
  nf<-calcNormFactors(d$counts,method="TMM")  
  y1<-voom(d$counts,plot=F,design0[,-1],lib.size=colSums(d$counts)*nf)
  for (i in c(2,3)){
    fit <- lmFit(y1,design0[,-i])
    fit2 <- eBayes(fit)
    options(digits=3)
    j=i+1
    top<-topTable(fit2,coef=j,n=nrow(fit2))
  # cat(paste("FDR<0.1 & |logFC|>3 & up=",length(rownames(top[top$adj.P.Val<0.1 & top$logFC>3,]))),"\n")
  # cat(paste("FDR<0.1 & |logFC|>3 & down=",length(rownames(top[top$adj.P.Val<0.1 & top$logFC<(-3),]))),"\n")
  	cat(paste("FDR<0.1 & |logFC|>3 & up=",length(rownames(top[top$adj.P.Val<0.1,]))),"\n")
  	cat(paste("FDR<0.1 & |logFC|>3 & down=",length(rownames(top[top$adj.P.Val<0.1,]))),"\n")
  }

#####################################################################################################
#   Differential gene expression analysis with weight observations from variable samples 
#####################################################################################################

  par(mfrow=c(2,2),mar=c(5,4,3,1))
  d<-DGEList(counts=d1.add.f,group=group0)
  nf<-calcNormFactors(d$counts,method="TMM")  ## or change "TMM" to "none", the results are very similar
  y2<-voomWithQualityWeights(d$counts,plot=F,replace.weights=T,design0[,-1],lib.size=colSums(d$counts)*nf) 
  par(mfrow=c(1,2),mar=c(7,4,3,1))
  for (i in c(2,3)){
    fit <- lmFit(y2,design0[,-i])
    fit2 <- eBayes(fit)
   	options(digits=3)
   	j=i+1
    top<-topTable(fit2,coef=j,n=nrow(fit2))
#  cat(paste("FDR<0.1 & |logFC|>3 & up=",length(rownames(top[top$adj.P.Val<0.1 & top$logFC>3,]))),"\n")
#  cat(paste("FDR<0.1 & |logFC|>3 & down=",length(rownames(top[top$adj.P.Val<0.1 & top$logFC<(-3),]))),"\n")
    cat(paste("FDR<0.1 & |logFC|>3 & up=",length(rownames(top[top$adj.P.Val<0.1,]))),"\n")
    cat(paste("FDR<0.1 & |logFC|>3 & down=",length(rownames(top[top$adj.P.Val<0.1,]))),"\n")

  	write.csv(top,file=paste("DERseults-",colnames(design0)[i+2],"-vs-",colnames(design0)[i],sep=""))
  	
   if(i==2){
    	j=i+1
    	smoothScatter(fit2$Amean,fit2$coefficients[,j],cex=0.4,ylab="logFC",xlab="Mean expression",
       	main="(a)  30S_283 vs 30C_283")
    	abline(h=0,col='red',lty=1,lwd=2)
    	de1.stat<-rownames(top[top$adj.P.Val<0.1 & abs(top$logFC)>3,])
  	#   de1.stat<-rownames(top[top$adj.P.Val<0.1 ,])
    	points(fit2$Amean[names(fit2$Amean)%in%de1.stat],fit2$coefficients[names(fit2$Amean)%in%de1.stat,3],
           col="red",cex=0.4,pch=16)      
 # 		points(fit2$Amean[names(fit2$Amean)%in%DE.283$X],fit2$coefficients[names(fit2$Amean)%in%DE.283$X,3],
 #      col="magenta",cex=0.4,pch=16)
  
  }
  if (i==3){  
      	j=i+1
      	smoothScatter(fit2$Amean,fit2$coefficients[,j],cex=0.4,ylab="logFC",xlab="Mean expression",
           main="(b)  30S_8261 vs 30C_8261")
      	abline(h=0,col='red',lty=1,lwd=2)
     	de1.stat<-rownames(top[top$adj.P.Val<0.1 & abs(top$logFC)>3,])
    #  	de1.stat<-rownames(top[top$adj.P.Val<0.1 ,])
     	points(fit2$Amean[names(fit2$Amean)%in%de1.stat],fit2$coefficients[names(fit2$Amean)%in%de1.stat,j],
             col="red",cex=0.4,pch=16)      
  #   	points(fit2$Amean[names(fit2$Amean)%in%DE.8261$Gene],fit2$coefficients[names(fit2$Amean)%in%DE.8261$Gene,j],
  #         col="magenta",cex=0.4,pch=16)
      
  }}



