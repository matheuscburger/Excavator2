vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


RC_Folder <- split.vars[1]
GC_Folder <- split.vars[2]
R_Target_Path <- split.vars[3]
R_Norm_Path <- split.vars[4]
Sample_Name <- split.vars[5]

RCFiles <- list.files(RC_Folder)
GCFiles <- list.files(GC_Folder)

unique.chrom <- paste("chr",c(seq(1,22,by=1),"X","Y"),".",sep="")

GCContentT <- c()
RCT <- c()

for (i in 1:length(unique.chrom))
{
	FileRC <- RCFiles[grep(unique.chrom[i],RCFiles,fixed=T)]
	FileGC <- GCFiles[grep(unique.chrom[i],GCFiles,fixed=T)]

	load(file.path(RC_Folder,FileRC))
	RCT <- c(RCT,RC)

	load(file.path(GC_Folder,FileGC))
	GCContentT <- c(GCContentT,GCContent)
}

load(R_Target_Path)
chrom <- as.character(MyTarget[,1])
start <- as.integer(MyTarget[,2])
end <- as.integer(MyTarget[,3])
L <- end-start+2
Pos<-(start+end)/2

RCTL <- RCT/L
GCContentL <- (GCContentT/L)*100

source(R_Norm_Path)

step <- 5
RCLNormGCList <- CorrectGC(RCTL,GCContentL,step)
RCLNormGC <- RCLNormGCList$DocNorm
MatrixFinal <- cbind(chrom,Pos,L,RCLNormGC)

save(MatrixFinal,file=file.path(RC_Folder,paste(Sample_Name,".RC.GC.RData",sep="")))
