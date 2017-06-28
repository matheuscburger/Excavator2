################ Correzione dei RC dal GC content secondo Yoon ####################
CorrectGC<-function(RC,GCContent,step)
{


stepseq<-seq(0,100,by=step)

MasterMedian<-median(RC,na.rm=T)
MedianGC<-rep(0,length(stepseq)-1)
RCNormMedian<-RC
for (i in 1:(length(stepseq)-1))
{
if (i==1)
{
ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1])
}
if (i!=1)
{
ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])
}
if (length(ind)>0)
{
m<-median(RC[ind],na.rm=T)
if (m>0)
{
MedianGC[i]<-m
RCNormMedian[ind]<-RC[ind]*MasterMedian/m
}
}
}
RCNormList<-list()
RCNormList$Median<-MedianGC
RCNormList$StepGC<-stepseq[1:(length(stepseq)-1)]
RCNormList$RCNorm<-RCNormMedian
RCNormList
}

################ Correzione dei RC Per il Bin Size ####################
CorrectSize<-function(RC,L,step)
{


stepseq<-seq(0,max(L),by=step)

MasterMedian<-median(RC,na.rm=T)
MedianL<-rep(0,length(stepseq)-1)
RCNormMedian<-RC
for (i in 1:(length(stepseq)-1))
{
if (i==1)
{
ind<-which(L>=stepseq[i] & L<=stepseq[i+1])
}
if (i!=1)
{
ind<-which(L>stepseq[i] & L<=stepseq[i+1])
}
if (length(ind)>0)
{
m<-median(RC[ind],na.rm=T)
if (m>0)
{
MedianL[i]<-m
RCNormMedian[ind]<-RC[ind]*MasterMedian/m
}
}
}
RCNormList<-list()
RCNormList$Median<-MedianL
RCNormList$RCNorm<-RCNormMedian
RCNormList
}

################ Correzione dei RC dalla Mappability ####################
CorrectMAP<-function(RC,MAPContent,step)
{


stepseq<-seq(0,100,by=step)

MasterMedian<-median(RC,na.rm=T)
MedianMAP<-rep(0,length(stepseq)-1)
RCNormMedian<-RC
for (i in 1:(length(stepseq)-1))
{
if (i==1)
{
ind<-which(MAPContent>=stepseq[i] & MAPContent<=stepseq[i+1])
}
if (i!=1)
{
ind<-which(MAPContent>stepseq[i] & MAPContent<=stepseq[i+1])
}

if (length(ind)>0)
{
m<-median(RC[ind],na.rm=T)
if (m>0)
{
MedianMAP[i]<-m
RCNormMedian[ind]<-RC[ind]*MasterMedian/m
}
}
}
RCNormList<-list()
RCNormList$Median<-MedianMAP
RCNormList$StepMAP<-stepseq[1:(length(stepseq)-1)]
RCNormList$RCNorm<-RCNormMedian
RCNormList
}

################ Calcolo dei percentili dei RC rispetto alla mappability ####################
QuantileMAP<-function(RCTL,MAP,step)
{
stepseq<-seq(0,100,by=step)
QuantileVecMAP<-c()
MedianVecMAP<-c()

for (i in 1:(length(stepseq)-1))
{
if (i==1)
{
ind<-which(MAP>=stepseq[i] & MAP<=stepseq[i+1])
}
if (i!=1)
{
ind<-which(MAP>stepseq[i] & MAP<=stepseq[i+1])
}
if (length(ind)>0)
{
QuantileVecMAP<-rbind(QuantileVecMAP,quantile(RCTL[ind],c(0.1,0.9)))
MedianVecMAP<-c(MedianVecMAP,median(RCTL[ind],na.rm=T))
}
if (length(ind)==0)
{
MedianVecMAP<-c(MedianVecMAP,0)
QuantileVecMAP<-rbind(QuantileVecMAP,c(0,0))
}
}

Outlist<-list()
Outlist$Median<-MedianVecMAP
Outlist$Quantile<-QuantileVecMAP
Outlist
}

################ Calcolo dei percentili dei RC rispetto al GC-content ####################
QuantileGC<-function(RCTL,GCContent,step)
{
stepseq<-seq(0,100,by=step)


QuantileVec<-c()

MedianVec<-c()

for (i in 1:(length(stepseq)-1))
{
if (i==1)
{
ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1])
}
if (i!=1)
{
ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])
}

if (length(ind)>0)
{
QuantileVec<-rbind(QuantileVec,quantile(RCTL[ind],c(0.1,0.9),na.rm=T))
MedianVec<-c(MedianVec,median(RCTL[ind],na.rm=T))

}
if (length(ind)==0)
{
MedianVec<-c(MedianVec,0)
QuantileVec<-rbind(QuantileVec,c(0,0))
}
}
Outlist<-list()
Outlist$Median<-MedianVec
Outlist$Quantile<-QuantileVec
Outlist
}

################ Calcolo dei percentili dei RC rispetto alla Length ####################
QuantileLength<-function(RCTL,L,step)
{
stepseq<-seq(0,max(L),by=step)
QuantileVecL<-c()
MedianVecL<-c()

for (i in 1:(length(stepseq)-1))
{
if (i==1)
{
ind<-which(L>=stepseq[i] & L<=stepseq[i+1])
}
if (i!=1)
{
ind<-which(L>stepseq[i] & L<=stepseq[i+1])
}
if (length(ind)>0)
{
QuantileVecL<-rbind(QuantileVecL,quantile(RCTL[ind],c(0.1,0.9)))
MedianVecL<-c(MedianVecL,median(RCTL[ind],na.rm=T))
}
if (length(ind)==0)
{
MedianVecL<-c(MedianVecL,0)
QuantileVecL<-rbind(QuantileVecL,c(0,0))
}
}

Outlist<-list()
Outlist$Median<-MedianVecL
Outlist$Quantile<-QuantileVecL
Outlist
}

#### Funzione per caricare i file di RC in un unico vettore #####
loadRC<-function(PathRC,RCFiles,unique.chrom)
{


RCT<-c()
stringSupp<-paste(".RC.",unique.chrom,".RData",sep="")

for (i in 1:length(unique.chrom))
{
  indRC<-grep(stringSupp[i],RCFiles,fixed=T)
  FileRC<-RCFiles[indRC]
  load(file.path(PathRC,FileRC))
  RCT<-c(RCT,RC)
}
RCT
}



###### Funzione per caricare Mappabilità e GC-Content del Target ####
loadTarget<-function(PathTarget,unique.chrom,TargetName)
{


fileTarget<-file.path(PathTarget,paste(TargetName,".RData",sep=""))

load(fileTarget)
PathGC<-file.path(PathTarget,"GCC")
PathMAP<-file.path(PathTarget,"MAP")

GCFiles<-list.files(PathGC)
MAPFiles<-list.files(PathMAP)
#GCFilesSupp<-paste(".",GCFiles,sep="")
#MAPFilesSupp<-paste(".",GCFiles,sep="")

unique.chromGC<-paste("GCC.",unique.chrom,".RData",sep="")
unique.chromMAP<-paste("Map.",unique.chrom,".RData",sep="")

GCContentT<-c()
MAPT<-c()

for (i in 1:length(unique.chrom))
{
  indGC<-which(GCFiles==unique.chromGC[i])
  indMAP<-which(MAPFiles==unique.chromMAP[i])
  FileGC<-GCFiles[indMAP]
  FileMAP<-MAPFiles[indGC]
  load(file.path(PathGC,FileGC))
  GCContentT<-c(GCContentT,GCContent)
  load(file.path(PathMAP,FileMAP))
  MAPT<-c(MAPT,MapMed)

}
Result<-list()
Result[[1]]<-MyTarget
Result[[2]]<-GCContentT
Result[[3]]<-MAPT
Result
}




