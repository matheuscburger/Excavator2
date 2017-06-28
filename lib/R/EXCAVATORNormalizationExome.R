vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))

options("scipen" = 20)

### Setting input paths for RC, mappability, GC-content and target ###
ProgramFolder <- split.vars[1]
PathInVec <- split.vars[2]
ExpName <- split.vars[3]
TargetName <- split.vars[4]
Assembly <- split.vars[5]

TargetFolder<-file.path(ProgramFolder,"data","targets",Assembly,TargetName)

### Setting paths for RC, mappability, GC-content e target ###
PathGC <- file.path(TargetFolder,"GCC")
PathMAP <- file.path(TargetFolder,"MAP")
GCFiles <- list.files(PathGC)
MAPFiles <- list.files(PathMAP)

### Loading target chromosomes ###
TargetChrom <- file.path(TargetFolder,paste(TargetName,"_chromosome.txt",sep=""))
CHR<-readLines(con = TargetChrom, n = -1L, ok = TRUE, warn = TRUE,encoding = "unknown")
unique.chrom<-strsplit(CHR," ")[[1]]

Path2ExomeRC <- file.path(ProgramFolder,"lib","R","LibraryExomeRC.R")
source(Path2ExomeRC)



#### Loading RC Data ####
PathRC <- file.path(PathInVec,"RC")
RCFiles <- list.files(PathRC)
RCTL <- loadRC(PathRC,RCFiles,unique.chrom)

#### Setting output file and folders ##
library(Hmisc)
PathNorm <- file.path(PathInVec,"RCNorm")
PathImages <- file.path(PathInVec,"Images")
FileOut <- file.path(PathNorm,paste(ExpName,".NRC.RData",sep=""))



### Loading target, mappability and GC content Files ###
TargetOut <- loadTarget(TargetFolder,unique.chrom,TargetName)
MyTarget <- TargetOut[[1]]
GCContentT <- TargetOut[[2]]*100
MAPT <- TargetOut[[3]]*100


chrom <- as.character(MyTarget[,1])
start <- as.integer(MyTarget[,2])
end <- as.integer(MyTarget[,3])
GeneName <- as.character(MyTarget[,4])
Class <- as.character(MyTarget[,5])
L <- end-start
Position <- (start+end)/2


RCTMatrixL <- t(t(RCTL)/L)
GCContentL <- GCContentT

indInTarget<-which(Class == "IN")
indOutTarget<-which(Class == "OUT")

GCContentLIn <-  GCContentL[indInTarget]
GCContentLOut <-  GCContentL[indOutTarget]
MAPTIn<-MAPT[indInTarget]
MAPTOut<-MAPT[indOutTarget]

### Normalization for GC Content, mappability and exon size if InTarget###

RCTL <- RCTMatrixL



### Exon size normalization only if InTarget###
step <- 5
RCLNormListIn <- CorrectSize(RCTL[indInTarget],L[indInTarget],step)
RCLNormIn <- RCLNormListIn$RCNorm


### Mappability normalization ###
step <- 5
RCMAPNormListIn <- CorrectMAP(RCLNormIn,MAPTIn,step)
RCMAPNormIn <- RCMAPNormListIn$RCNorm

RCMAPNormListOut <- CorrectMAP(RCTL[indOutTarget],MAPTOut,step)
RCMAPNormOut <- RCMAPNormListOut$RCNorm


### GC-content Normalization ###
step <- 5
RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContentLIn,step)
RCGCNormIn <- RCGCNormListIn$RCNorm

RCGCNormListOut <- CorrectGC(RCMAPNormOut,GCContentLOut,step)
RCGCNormOut <- RCGCNormListOut$RCNorm


### Images generation ###
### Quantiles for exon size only for InTarget ###


step <- 5
QuantileLList <- QuantileLength(RCTL[indInTarget],L[indInTarget],step)
MedianVecL <- QuantileLList$Median
QuantileVecL <- QuantileLList$Quantile

SizeBias <- file.path(PathImages,"InSizeBias.pdf")
pdf(SizeBias,height=10,width=15)
par(mfrow=c(2,1))
yup <- median(RCTL[indInTarget],na.rm=T)*3
stepseq <- seq(0,max(L[indInTarget]),by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Exon Length (bp)"
labely <- "WMRC"
errbar(xstep,MedianVecL,QuantileVecL[,2],QuantileVecL[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,xlim=c(0,160),ylim=c(0,yup),lty=3,cex.axis=1.5,cex.lab=1.5)
#title("Exon Length Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),cex=2)
axis(2,cex=2)


### Quantile of normalized exon length only for InTarget###
step <- 5
QuantileLListNIn <- QuantileLength(RCLNormIn,L[indInTarget],step)
MedianVecLNIn <- QuantileLListNIn$Median
QuantileVecLNIn <- QuantileLListNIn$Quantile

#SizeBiasNormIn <- file.path(PathImages,"SizeBiasNorm.pdf")
#pdf(SizeBiasNormIn,height=10,width=15)
yup <- median(RCLNormIn,na.rm=T)*2
stepseq <- seq(0,max(L[indInTarget]),by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Exon Length (bp)"
labely <- "WMRC"
errbar(xstep,MedianVecLNIn,QuantileVecLNIn[,2],QuantileVecLNIn[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,xlim=c(0,160),ylim=c(0,yup),lty=3,cex.axis=1.5,cex.lab=1.5)
title("Normalized Read Count Exon Length Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
#ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),cex=2)
axis(2,cex=2)
dev.off()



### Quantiles for mappability ###
step <- 5
## IN
QuantileMapListIn <- QuantileMAP(RCTL[indInTarget],MAPTIn,step)
MedianVecMAPIn <- QuantileMapListIn$Median
QuantileVecMAPIn <- QuantileMapListIn$Quantile

MAPBiasIn <- file.path(PathImages,"InMAPBias.pdf")
pdf(MAPBiasIn,height=10,width=15)
par(mfrow=c(2,1))
yup <- max(QuantileVecMAPIn[,1]) + max(QuantileVecMAPIn[,2])
#yup <- median(RCTL[indInTarget],na.rm=T)*4
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Mappability"
labely <- "WMRC"
errbar(xstep,MedianVecMAPIn,QuantileVecMAPIn[,2],QuantileVecMAPIn[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.7,cex.lab=1.7)
#title("Mappability Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),font.axis=2,cex=2)
axis(2,cex=2)

### Quantile of normalized mappability ###
step <- 5
#IN
QuantileMapListNIn <- QuantileMAP(RCMAPNormIn,MAPTIn,step)
MedianVecMAPNIn <- QuantileMapListNIn$Median
QuantileVecMAPNIn <- QuantileMapListNIn$Quantile
#MAPBiasNormIn <- file.path(PathImages,"InMAPBiasNorm.pdf")
#pdf(MAPBiasNormIn,height=10,width=15)
yup <- median(RCMAPNormIn,na.rm=T)*3
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Mappability"
labely <- "Normalized WMRC"
errbar(xstep,MedianVecMAPNIn,QuantileVecMAPNIn[,2],QuantileVecMAPNIn[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.7,cex.lab=1.7)
#title("Normalized Read Count Mappability Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),cex=2)
axis(2,cex=2)
dev.off()

### Quantiles for mappability ###
## OUT
QuantileMapListOut <- QuantileMAP(RCTL[indOutTarget],MAPTOut,step)
MedianVecMAPOut <- QuantileMapListOut$Median
QuantileVecMAPOut <- QuantileMapListOut$Quantile

MAPBiasOut <- file.path(PathImages,"OutMAPBias.pdf")
pdf(MAPBiasOut,height=10,width=15)
par(mfrow=c(2,1))
par(mar=c(5.1,5.1,4.1,2.1))
yup <- max(QuantileVecMAPOut[,1]) + max(QuantileVecMAPOut[,2])
#yup <- median(RCTL[indOutTarget],na.rm=T)*6
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Mappability"
labely <- "WMRC"
errbar(xstep,MedianVecMAPOut,QuantileVecMAPOut[,2],QuantileVecMAPOut[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.7,cex.lab=1.7)
#title("Mappability Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),cex=2)
axis(2,cex=2)
#dev.off()

### Quantile of normalized mapping ###
step <- 5
#OUT
QuantileMapListNOut <- QuantileMAP(RCMAPNormOut,MAPTOut,step)
MedianVecMAPNOut <- QuantileMapListNOut$Median
QuantileVecMAPNOut <- QuantileMapListNOut$Quantile
#MAPBiasNormOut <- file.path(PathImages,"OutMAPBiasNorm.pdf")
#pdf(MAPBiasNormOut,height=10,width=15)
yup <- max(MedianVecMAPNOut) + max(QuantileVecMAPNOut[,2])
#yup <- median(RCMAPNormOut,na.rm=T)*3
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "Mappability"
labely <- "Normalized WMRC"
errbar(xstep,MedianVecMAPNOut,QuantileVecMAPNOut[,2],QuantileVecMAPNOut[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.7,cex.lab=1.7)
#title("Normalized Read Count Mappability Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
dev.off()


#### Quantiles for GC-content ###
step <- 5
#IN
QuantileGCListIn <- QuantileGC(RCTL[indInTarget],GCContentLIn,step)
MedianVecGCIn <- QuantileGCListIn$Median
QuantileVecGCIn <- QuantileGCListIn$Quantile

GCBiasIn <- file.path(PathImages,"InGCBias.pdf")
pdf(GCBiasIn,height=10,width=15)
par(mfrow=c(2,1))
yup <- max(MedianVecGCIn) + max(QuantileVecGCIn[,2])
#yup <- median(RCTL[indInTarget],na.rm=T)*4
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "GC-Content (%)"
labely <- "WMRC"
errbar(xstep,MedianVecGCIn,QuantileVecGCIn[,2],QuantileVecGCIn[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.7,cex.lab=1.7)
#title("GC-Content Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),cex=2)
axis(2,cex=2)
#dev.off()


#### Quantile of normalized GC-content ###
#IN
step <- 5
QuantileGCListNIn <- QuantileGC(RCGCNormIn,GCContentLIn,step)
MedianVecGCNIn <- QuantileGCListNIn$Median
QuantileVecGCNIn <- QuantileGCListNIn$Quantile
#GCBiasNormIn <- file.path(PathImages,"InGCBiasNorm.pdf")
#pdf(GCBiasNormIn,height=10,width=15)
yup <- median(RCGCNormIn,na.rm=T)*3
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "GC-Content (%)"
labely <- "Normalized WMRC"
errbar(xstep,MedianVecGCNIn,QuantileVecGCNIn[,2],QuantileVecGCNIn[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.7,cex.lab=1.7)
#title("Normalized Read Count GC-Content Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"),cex=2)
axis(2,cex=2)
dev.off()

#### Quantiles for GC-content ###
#OUT
QuantileGCListOut <- QuantileGC(RCTL[indOutTarget],GCContentLOut,step)
MedianVecGCOut <- QuantileGCListOut$Median
QuantileVecGCOut <- QuantileGCListOut$Quantile
GCBiasOut <- file.path(PathImages,"OutGCBias.pdf")
pdf(GCBiasOut,height=10,width=15)
par(mfrow=c(2,1))
yup <- max(MedianVecGCOut) + max(QuantileVecGCOut[,2])
#yup <- median(RCTL[indOutTarget],na.rm=T)*3
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "GC-Content (%)"
labely <- "WMRC"
errbar(xstep,MedianVecGCOut,QuantileVecGCOut[,2],QuantileVecGCOut[,1],axes=F,frame.plot=TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.2,cex.lab=1.3)
#title("GC-Content Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
#dev.off()


#OUT
step <- 5
QuantileGCListNOut <- QuantileGC(RCGCNormOut,GCContentLOut,step)
MedianVecGCNOut <- QuantileGCListNOut$Median
QuantileVecGCNOut <- QuantileGCListNOut$Quantile
#GCBiasNormOut <- file.path(PathImages,"OutGCBiasNorm.pdf")
#pdf(GCBiasNormOut,height=10,width=15)
yup <- max(MedianVecGCNOut) + max(QuantileVecGCNOut[,2])
#yup <- median(RCGCNormOut,na.rm=T)*3
stepseq <- seq(0,100,by=step)
xstep <- c(1:(length(stepseq)-1))
labelx <- "GC-Content (%)"
labely <- "Normalized WMRC"
errbar(xstep,MedianVecGCNOut,QuantileVecGCNOut[,2],QuantileVecGCNOut[,1],axes=F,frame.plot = TRUE,xlab=labelx,ylab=labely,add=F,ylim=c(0,yup),lty=3,cex.axis=1.2,cex.lab=1.3)
#title("Normalized Read Count GC-Content Bias",cex.main=2)
#text((length(stepseq)-3),2000,"a",cex=3)
ax.x <- stepseq[2:length(stepseq)]
axis(1,xstep,labels=formatC(ax.x, format="fg"))
axis(2)
dev.off()

##### Saving normalized RC file###

RCGCNormIn[which(RCGCNormIn==0)] <- min(RCGCNormIn[which(RCGCNormIn!=0)])
RCGCNormOut[which(RCGCNormOut==0)] <- min(RCGCNormOut[which(RCGCNormOut!=0)])

### Combining read count for In and Out target


RCNorm <- rep(NA,length(RCTL))
RCNorm[indInTarget]<- RCGCNormIn
RCNorm[indOutTarget]<- RCGCNormOut
MatrixNorm <- cbind(chrom,Position,start,end,GeneName,RCNorm,Class)
FileOut <- file.path(PathNorm,paste(ExpName,".NRC.RData",sep=""))
save(MatrixNorm,file=FileOut)













