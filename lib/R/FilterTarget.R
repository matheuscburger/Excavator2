vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))
BedIn <- split.vars[1]
ProgramFolder <- split.vars[2]
target.name <- split.vars[3]
assembly <- split.vars[4]
step <- as.numeric(split.vars[5])
flank<-200

options("scipen"=20)

if (assembly=="hg19"){
  CoordIn <- paste(ProgramFolder,"data/support",assembly,"ChromosomeCoordinate_HG19.txt",sep="/")
  FileGap <- paste(ProgramFolder,"data/support",assembly,"GapHg19.UCSC.txt",sep="/") 
}

if (assembly=="hg38"){
  CoordIn <- paste(ProgramFolder,"data/support",assembly,"ChromosomeCoordinate_HG38.txt",sep="/")
  FileGap <- paste(ProgramFolder,"data/support",assembly,"GapHg38.UCSC.txt",sep="/") 
}




CoordTable<-read.table(CoordIn,sep="\t",quote="\"",fill=T,header=F)
ChrCoord<-as.character(CoordTable[,1])
StartCoord<-as.numeric(CoordTable[,2])
EndCoord<-as.numeric(CoordTable[,3])



BedTable<-read.table(BedIn,sep="\t",quote="\"",fill=T,header=F)

Chr<-as.character(BedTable[,1])
Start<-as.numeric(BedTable[,2])
End<-as.numeric(BedTable[,3])

if (nchar(Chr[1])<3)
{
  ChrVec<-c(1:22,"X")
  ChrCoord<-c(1:22,"X")
}
if (nchar(Chr[1])>3)
{
  ChrVec<-paste("chr",c(1:22,"X"),sep="")
  ChrCoord<-paste("chr",c(1:22,"X"),sep="")
}


emptyS<-c()
emptyE<-c()
emptyChr<-c()

stepP<-step+2*flank
Totalempty<-c()
TotalTarget<-c()
for (i in 1:length(ChrVec))
{
  indCoordC<-which(ChrCoord==ChrVec[i])
  StartCoordC<-StartCoord[indCoordC]
  EndCoordC<-EndCoord[indCoordC]
  indC<-which(Chr==ChrVec[i]) 
  StartC<-Start[indC]
  EndC<-End[indC]
  emptyS<-c(StartCoordC,EndC[1:(length(EndC))])
  emptyE<-c(StartC[1:length(StartC)],EndCoordC)
  
  emptySize<-emptyE-emptyS
  indF<-which(emptySize>=stepP)
  StartCF<-emptyS[indF]
  EndCF<-emptyE[indF]
  emptySizeF<-emptySize[indF]
  
  TotalTarget<-rbind(TotalTarget,cbind(Chr[indC],StartC,EndC,rep("IN",length(indC))))
  
  TotalemptyC<-c()
  for (k in 1:length(indF))
  {
    numstep<-floor(emptySizeF[k]/step)
    seqstep<-seq(0,numstep*step,by=step)
    emptystep<-StartCF[k]+seqstep+flank
    emptystepS<-emptystep[1:(length(emptystep)-1)]+1
    emptystepE<-emptystep[2:length(emptystep)]
    TotalemptyC<-rbind(TotalemptyC,cbind(rep(ChrVec[i],length(emptystep)-1),emptystepS,emptystepE,rep("OUT",length(emptystep)-1)))
  }
  
  
  
  Totalempty<-rbind(Totalempty,TotalemptyC)
  
}

FullTarget<-rbind(TotalTarget,Totalempty)
ChrFull<-as.character(FullTarget[,1])
StartFull<-as.numeric(FullTarget[,2])
EndFull<-as.numeric(FullTarget[,3])
StatusFull<-as.character(FullTarget[,4])

FinalTarget<-c()
for (i in 1:length(ChrVec))
{
  indC<-which(ChrFull==ChrVec[i]) 
  StartFullC<-StartFull[indC]
  EndFullC<-EndFull[indC]
  StatusFullC<-StatusFull[indC]
  indS<-sort(StartFullC,index.return=T)$ix
  FinalTarget<-rbind(FinalTarget,cbind(ChrFull[indC],StartFullC[indS],EndFullC[indS],StatusFullC[indS]))
  
}

### Filtro Gapmeri e Gap##
GapTable<-read.table(FileGap,sep="\t",header=F,skip=1,fill=TRUE,quote = "")



  if (nchar(Chr[1])>3)
{
  GapChrIn<-as.character(GapTable[,2])
}else{
  GapChrIn<-sub("chr","",as.character(GapTable[,2]))
} 

indAlt <- grep("_",GapChrIn) 
GapChr <- GapChrIn[-indAlt]
GapStart<-as.numeric(GapTable[-indAlt,3])
GapEnd<-as.numeric(GapTable[-indAlt,4])

FinalTargetF<-c()
for (c in 1:length(ChrVec))
{
  indGap<-which(GapChr==ChrVec[c])
  GapPos<-cbind(GapStart[indGap],GapEnd[indGap])
  indC<-which(FinalTarget[,1]==ChrVec[c])
  FinalTargetC<-cbind(FinalTarget[indC,],rep(0,length(indC)))
  for (i in 1:nrow(FinalTargetC)) 
  {
    if(FinalTargetC[i,5]!=1)
    { 
      for (j in 1:nrow(GapPos))
      {
        if(length(intersect(findInterval(as.vector(FinalTargetC[i,c(2,3)]),as.vector(GapPos[j,c(1,2)])),1))!=0)
        {FinalTargetC[i,5] <- 1}
      }
    }
  }
  indF <- which(FinalTargetC[,5]==0)
  FinalTargetF<-rbind(FinalTargetF,FinalTargetC[indF,c(1,2,3,4)])
}

#new format
a<-paste("a",seq(1,nrow(FinalTargetF)),sep="")
FinalTargetF<-cbind(FinalTargetF[,c(1:3)],a,FinalTargetF[,4])








MyTarget <- FinalTargetF
MyChr <- unique(FinalTargetF[,1])


dir.create(file.path(ProgramFolder,"data","targets",assembly))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name,"GCC"))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name,"MAP"))
dir.create(file.path(ProgramFolder,"data","targets",assembly,target.name,"FRB"))
write.table(data.frame(rbind(MyChr)),file.path(ProgramFolder,"data","targets",assembly,target.name,paste(target.name,"_chromosome.txt",sep="")),col.names=F,row.names=F,quote=F)
save(MyTarget,file=file.path(ProgramFolder,"data","targets",assembly,target.name,paste(target.name,".RData",sep="")))
write.table(MyTarget,file.path(ProgramFolder,"data","targets",assembly,target.name,paste("Filtered.txt",sep="")),col.names=F,row.names=F,sep="\t",quote=F)

