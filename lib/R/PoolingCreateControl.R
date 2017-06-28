vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


##  Setting input paths for normalized read count and experimental design ###
ProgramFolder <- split.vars[1]
DataFolder <- split.vars[2]
TargetFolder <- split.vars[3]
ExperimentalFile <- split.vars[4]
ExperimentalDesign <- split.vars[5]
TargetName <- split.vars[6]



### Load and set experimental design ###
ExperimentalTable <- read.table(ExperimentalFile,sep=" ",quote="",header=F)
LabelName <- as.character(ExperimentalTable[,1])
PathInVec <- as.character(ExperimentalTable[,2])
ExpName <- as.character(ExperimentalTable[,3])

### Loading target chromosomes ###
TargetChrom <- file.path(TargetFolder,paste(TargetName,"_chromosome.txt",sep=""))
CHR<-readLines(con = TargetChrom, n = -1L, ok = TRUE, warn = TRUE,encoding = "unknown")
unique.chrom<-strsplit(CHR," ")[[1]]

Path2ExomeRC <- file.path(ProgramFolder,"lib","R","LibraryExomeRC.R")
source(Path2ExomeRC)


### Create the vector for the experimental design ###
if (ExperimentalDesign=="pooling")
{
  indC <- grep("C",LabelName)
  PathInVecC<-PathInVec[indC]
  ExpNameC<-ExpName[indC]
  ### Create the RC matrix with all the Experiments ### 
  
  PathRC <- file.path(PathInVecC[1],"RCNorm")
  FileIn<-file.path(PathRC,paste(ExpNameC[1],".NRC.RData",sep=""))
  load(FileIn)
  
  MetaData<-MatrixNorm[,1:5]
  Class<-MatrixNorm[,7]
  RCNorm<-as.numeric(MatrixNorm[,6])
  if (length(PathInVecC)>1)
  {
  for (i in 2:length(PathInVecC))
  {
    PathRC <- file.path(PathInVecC[i],"RCNorm")
    FileIn<-file.path(PathRC,paste(ExpNameC[i],".NRC.RData",sep=""))
    load(FileIn)
    RCNorm<-RCNorm+as.numeric(MatrixNorm[,6])
    }
  RCNorm<-RCNorm/length(PathInVecC)
  }
  MatrixNorm<-cbind(MetaData,RCNorm,Class)
}

dir.create(file.path(DataFolder,"Control"))
dir.create(file.path(DataFolder,"Control","RCNorm"))

FileOut<-file.path(DataFolder,"Control","RCNorm","Control.NRC.RData")
save(MatrixNorm,file=FileOut)




