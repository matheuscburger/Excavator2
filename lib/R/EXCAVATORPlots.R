vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


###  Setting input paths for normalized read count and experimental design ###
DataFolder <- split.vars[1]
ExperimentalFile <- split.vars[2]


### Load and set experimental design ###
ExperimentalTable <- read.table(ExperimentalFile,sep=" ",quote="",header=F)
LabelName <- as.character(ExperimentalTable[,1])
ExpName <- as.character(ExperimentalTable[,3])
  

### Create the vector for the experimental design ###
indT <- grep("T",LabelName)
ExpTest <- c(ExpName[indT])



for (zz in 1:length(ExpTest))
{
FileSeg <- file.path(DataFolder,"Results",ExpTest[zz],paste("HSLMResults_",ExpTest[zz],".txt",sep=""))
FileCall <- file.path(DataFolder,"Results",ExpTest[zz],paste("FastCallResults_",ExpTest[zz],".txt",sep=""))
PathOut <- file.path(DataFolder,"Plots",ExpTest[zz])

############   ##################

y <- read.table(FileSeg,sep="\t",quote="\"",fill=T,header=T)

log2RSeq <- as.numeric(as.character(y[,3]))
chrSeq <- as.character(y[,1])
SegSeq <- as.numeric(as.character(y[,4]))
PositionSeq <- as.numeric(as.character(y[,2]))

UniqueChr <- unique(chrSeq)

z <- read.table(FileCall,sep="\t",quote="\"",fill=T,header=T)

chrCall <- as.character(z[,1])
StartCall <- as.numeric(as.character(z[,2]))
EndCall <- as.numeric(as.character(z[,3]))
Call <- as.numeric(as.character(z[,5]))

for (i in 1:length(UniqueChr))
{
	indSeq <- which(chrSeq==UniqueChr[i])
	PositionSeqC <- PositionSeq[indSeq]
	log2RSeqC <- log2RSeq[indSeq]
	SegSeqC <- SegSeq[indSeq]
	
	FileName <- paste("PlotResults_", UniqueChr[i],".pdf",sep="")
	
	FilePlot <- file.path(PathOut, FileName)
	pdf(FilePlot,height=10,width=15)
	par(mfrow=c(2,1))
	plot(PositionSeqC,log2RSeqC,ylim=c(-3,3),main=UniqueChr[i],pch=16,cex=0.5,xlab="Position",ylab="log2ratio")
	lines(PositionSeqC, SegSeqC,lwd=2,col="red")
	abline(h=0,lty=2,lwd=1,col="green")
	
	plot(PositionSeqC,log2RSeqC,ylim=c(-3,3),type="l",main=UniqueChr[i],lwd=0.5,col="grey",xlab="Position",ylab="log2ratio")
	
	indCall <- which(chrCall==UniqueChr[i])
	if (length(indCall)!=0)
	{
		StartCallC <- StartCall[indCall]
		EndCallC <- EndCall[indCall]
		CallC <- Call[indCall]
		for (j in 1:(length(indCall)))
		{
			if (CallC[j]==1)
			{
				rect(xleft=StartCallC[j], ybottom=0, xright=EndCallC[j], ytop=3, density = NA, angle = 45, col = "red",border="red")
			}
			if (CallC[j]==2)
			{
				rect(xleft=StartCallC[j], ybottom=0, xright=EndCallC[j], ytop=3, density = NA, angle = 45, col = "darkred",border="darkred")
			}
			if (CallC[j]==-1)
			{
				rect(xleft=StartCallC[j], ybottom=-3, xright=EndCallC[j], ytop=0, density = NA, angle = 45, col = "green",border="green")
			}
			if (CallC[j]==-2)
			{
				rect(xleft=StartCallC[j], ybottom=-3, xright=EndCallC[j], ytop=0, density = NA, angle = 45, col = "darkgreen",border="darkgreen")
			}
		}
	}
	
	lines(PositionSeqC,log2RSeqC,ylim=c(-3,3),main=UniqueChr[i],lwd=0.5,col="grey")
	abline(h=0,lwd=1,col="black")
	dev.off()
}
}
