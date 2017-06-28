FileIn<-"/home/ngs1/AnalisiCNV/TargetInput/ProvaGC1.txt"
vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))
chrsel <- split.vars[1]
GCCPath <- split.vars[2]


FileIn<-file.path(GCCPath,"GCC.txt")


GCContent <- scan(FileIn)



save(GCContent,file=file.path(GCCPath,paste("GCC.",chrsel,".RData",sep="")))
