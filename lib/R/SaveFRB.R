vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))
chrsel <- split.vars[1]
FRBPath <- split.vars[2]


File_In<-file.path(FRBPath,"FRB.txt")
FRBMat<-read.table(File_In,sep="\n",quote="\"",fill=T,header=F)
FRBMat<-as.character(FRBMat[,1])
indCoord<-seq(1,length(FRBMat),by=2)
indSeq<-seq(2,length(FRBMat),by=2)

FRBData<-cbind(unlist(strsplit(unlist(strsplit(FRBMat[indCoord],":"))[seq(2,length(indCoord)*2,by=2)],"-"))[seq(1,length(indCoord)*2,by=2)],FRBMat[indSeq])

save(FRBData,file=file.path(FRBPath,paste("FRB.",chrsel,".RData",sep="")))
