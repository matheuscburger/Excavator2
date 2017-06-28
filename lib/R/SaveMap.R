vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))
chrsel <- split.vars[1]
MAPPath <- split.vars[2]


FileIn<-file.path(MAPPath,"Mapout.txt")


MappingTable <- read.table(FileIn,sep="\t",quote="",header=F)

MapMed<-MappingTable[,6]


save(MapMed,file=file.path(MAPPath,paste("Map.",chrsel,".RData",sep="")))
