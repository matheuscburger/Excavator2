vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


###  Setting input paths for normalized read count and experimental design ###
ProgramFolder <- split.vars[1]
ExperimentalFile <- split.vars[2]
OutputFolder <- split.vars[3]
TargetName <- split.vars[4]
Assembly <- split.vars[5]
Processors <- as.numeric(split.vars[6])
mode <- split.vars[7]


PathOut<-file.path(OutputFolder,".tmp")
ExperimentalFileSplit<-unlist(strsplit(ExperimentalFile,"/"))
ExperimentalFileLabel<-ExperimentalFileSplit[length(ExperimentalFileSplit)]
SettingLabel<-paste(" --assembly ",Assembly," --output ",OutputFolder," --target ",TargetName," --mode ",mode,sep="")
### Load and set experimental design ###
ExperimentalTable <- read.table(ExperimentalFile,sep=" ",quote="",header=F)
LabelName <- as.character(ExperimentalTable[,1])
PathInVec <- as.character(ExperimentalTable[,2])
ExpName <- as.character(ExperimentalTable[,3])


if (mode=="pooling")
{
  indT <- grep("T",LabelName)
  ExperimentalTableFinal<-ExperimentalTable[indT,]
}


if (mode=="paired")
{
  indT <- grep("T",LabelName)
  indC <- grep("C",LabelName)
  ExperimentalTableT<-ExperimentalTable[indT,]
  ExperimentalTableC<-ExperimentalTable[indC,]
  LabelNameT<-LabelName[indT]
  LabelNameC<-LabelName[indC]
  
  
  NumT<-unlist(strsplit(LabelNameT,"T"))[seq(2,length(indT)*2,by=2)]
  NumC<-unlist(strsplit(LabelNameC,"C"))[seq(2,length(indC)*2,by=2)]
  indTS<-sort(NumT,index.return=T,decreasing = FALSE)$ix
  indCS<-sort(NumC,index.return=T,decreasing = FALSE)$ix
  ExperimentalTableTS<-rbind(ExperimentalTableT[indTS,])
  ExperimentalTableCS<-rbind(ExperimentalTableC[indCS,])
  indseqT<-seq(1,length(indTS)*2,by=2)
  indseqC<-seq(2,length(indCS)*2,by=2)
  ExperimentalTableFinal<-ExperimentalTable
  
  ExperimentalTableFinal[indseqT,]<-ExperimentalTableTS
  ExperimentalTableFinal[indseqC,]<-ExperimentalTableCS

}


if (mode=="paired")
{
  
NExp<-nrow(ExperimentalTableFinal)/2
if (Processors>NExp)
{
  Processors <- NExp
}
}

if (mode!="paired")
{
  
  NExp<-nrow(ExperimentalTableFinal)
  if (Processors>NExp)
  {
    Processors <- NExp
  }
}




if (mode!="paired")
{
  Q <- NExp%/%Processors
  R <- NExp%%Processors
  ExpPart<-c(rep(Q+1,times=R),rep(Q,times=Processors-R))
  
ShVector<-c()
StartInd<-1
for (i in 1:length(ExpPart))
{
  EndInd<-StartInd+ExpPart[i]-1
  ExperimentalTableSplit<-ExperimentalTableFinal[c(StartInd:EndInd),]
  FileOut<-file.path(PathOut,paste(ExperimentalFileLabel,"_",i,".txt",sep=""))
  write.table(ExperimentalTableSplit,FileOut,col.names=F,row.names=F,sep=" ",quote=F)
  ShVector<-c(ShVector,paste("perl ",file.path(ProgramFolder,"lib","perl","AnalyzePerla.pl "),FileOut,SettingLabel,sep=""))
  StartInd<-EndInd+1
  }

}

if (mode=="paired")
{
  Q <- NExp%/%Processors
  R <- NExp%%Processors
  ExpPart<-c(rep(Q+1,times=R),rep(Q,times=Processors-R))*2
  ShVector<-c()
  StartInd<-1
  for (i in 1:length(ExpPart))
  {
    EndInd<-StartInd+ExpPart[i]-1
    ExperimentalTableSplit<-ExperimentalTableFinal[c(StartInd:EndInd),]
    FileOut<-file.path(PathOut,paste(ExperimentalFileLabel,"_",i,".txt",sep=""))
    write.table(ExperimentalTableSplit,FileOut,col.names=F,row.names=F,sep=" ",quote=F)
    ShVector<-c(ShVector,paste("perl ",file.path(ProgramFolder,"lib","perl","AnalyzePerla.pl "),FileOut,SettingLabel,sep=""))
    StartInd<-EndInd+1
  }
  
}

FileShOut<-file.path(PathOut,"ParallelAnalyzePerla.sh")
write.table(ShVector,FileShOut,col.names=F,row.names=F,sep="\t",quote=F)

