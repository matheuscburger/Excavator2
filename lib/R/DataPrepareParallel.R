vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


###  Setting input paths for normalized read count and experimental design ###
ProgramFolder <- split.vars[1]
ExperimentalFile <- split.vars[2]
labeltmp <- split.vars[3]
TargetName <- split.vars[4]
Assembly <- split.vars[5]
Processors <- as.numeric(split.vars[6])



PathOut<-file.path(ProgramFolder,labeltmp)
ExperimentalFileSplit<-unlist(strsplit(ExperimentalFile,"/"))
ExperimentalFileLabel<-ExperimentalFileSplit[length(ExperimentalFileSplit)]
SettingLabel<-paste(" --assembly ",Assembly," --target ",TargetName,sep="")
### Load and set experimental design ###
ExperimentalTable <- read.table(ExperimentalFile,sep=" ",quote="",header=F)
LabelName <- as.character(ExperimentalTable[,1])
PathInVec <- as.character(ExperimentalTable[,2])
ExpName <- as.character(ExperimentalTable[,3])

NExp<-nrow(ExperimentalTable)
if (Processors>NExp)
{
  Processors <- NExp
}


Q <- NExp%/%Processors
R <- NExp%%Processors
ExpPart<-c(rep(Q+1,times=R),rep(Q,times=Processors-R))




ShVector<-c()
StartInd<-1
for (i in 1:length(ExpPart))
{
  EndInd<-StartInd+ExpPart[i]-1
  ExperimentalTableSplit<-ExperimentalTable[c(StartInd:EndInd),]
  FileOut<-file.path(PathOut,paste(ExperimentalFileLabel,"_",i,".txt",sep=""))
  write.table(ExperimentalTableSplit,FileOut,col.names=F,row.names=F,sep=" ",quote=F)
  ShVector<-c(ShVector,paste("perl ",file.path(ProgramFolder,"lib","perl","ReadPerla.pl "),FileOut,SettingLabel,sep=""))
  StartInd<-EndInd+1
}


FileShOut<-file.path(PathOut,"ParallelReadPerla.sh")
write.table(ShVector,FileShOut,col.names=F,row.names=F,sep="\t",quote=F)

