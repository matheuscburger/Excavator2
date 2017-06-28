######## Starting Condition Fucntion ########

StartCond<-function(mdata,thru,thrd)

{
muvec<- c(-3,-1,0,0.58,1)
sdvec<-c(0.01,0.01,0.01,0.01,0.01)
lvec<-c(-50,-1.5,-thrd,thru,0.9)
uvec<-c(-1.5,-thrd,thru,0.9,50)

for (i in 1:5)
{
u<-uvec[i]
l<-lvec[i]
ind<-which(mdata<=u & mdata>=l)
if (length(ind)==0)
{
muvec[i]<-muvec[i]
sdvec[i]<-sdvec[i]
}
if (length(ind)==1)
{
#muvec[i]<-mdata[ind]
muvec[i]<-muvec[i]
sdvec[i]<-sdvec[i]
}
if (length(ind)>1)
{
#muvec[i]<-mean(mdata[ind])
muvec[i]<-muvec[i]
sdvec[i]<-sd(mdata[ind])
}
}
sdvec[which(sdvec<0.001)]<-0.001
Result<-list()
Result$muvec<-muvec
Result$sdvec<-sdvec
Result
}

################### EM Algorithm ##############

EMFastCall<-function(mdata,thru,thrd)
  
{
  
  StartPar<-StartCond(mdata,thru,thrd)
  muvec<-StartPar$muvec
  sdvec<-StartPar$sdvec
  prior<-c(0.05,0.1,0.7,0.1,0.05)
  lvec<-c(-20,-1.3,-thrd,thru,0.9)
  uvec<-c(-1.3,-thrd,thru,0.9,20)
  LikeliNew<-sum(PosteriorP(mdata,muvec,sdvec,prior)*prior)
  
  threshold<-1e-05
  for (i in 1:1000)
  {
    muvecold<-muvec
    LikeliOld<-LikeliNew
    taux<-EStep(mdata,muvec,sdvec,prior,lvec,uvec)
    MResult<-MStep(mdata,taux,muvec,sdvec)
    muvec<-MResult$mu
    sdvec<-MResult$sdev
    prior<-MResult$prior
    LikeliNew<-sum(PosteriorP(mdata,muvec,sdvec,prior)*prior)
    if (abs(LikeliNew-LikeliOld)<threshold)
    {	
      break
    }
    
  }
  Results<-list()
  Results$muvec<-muvec
  Results$sdvec<-sdvec
  Results$prior<-prior
  Results$iter<-i
  Results$bound<-cbind(lvec,uvec)
  Results
}




######## GFCT Function ##########
gfct<-function(x,moy,sdev,l,u)
{
	(dnorm(x,mean=moy,sd=sdev)*(x<=u)*(x>=l))/(pnorm(u,mean=moy,sd=sdev)-pnorm(l,mean=moy,sd=sdev))
}


######## Posterior Probability Function ##########
PosteriorP<-function(mdata,muvec,sdvec,prior)
{
		Ncount<-length(muvec)
		Posterior <- matrix(0,nrow=length(mdata),ncol=length(muvec))
		deno <- 0
		normaldataMat<-c()

		for (j in 1:Ncount)
		{
			moy<-muvec[j]
			sdev<-sdvec[j]
			normaldataMat<-cbind(normaldataMat,dnorm(mdata,mean=moy,sd=sdev))
		}
		tauMat<-c()
		for (j in 1:Ncount)
		{
			tauMat<-cbind(tauMat,c(prior[j]*(t(normaldataMat[,j]))))		}
			
		deno <- rowSums(tauMat)
		Posterior<- tauMat/deno
		Posterior
}



####### Expectation Step #############
EStep<-function(mdata,muvec,sdvec,prior,lvec,uvec)
{
	Ncount<-length(muvec)
	taux <- matrix(0,nrow=length(mdata),ncol=Ncount)
	deno <- 0
	normaldataMat<-c()
	for (g in 1:Ncount)
	{
		mu<-muvec[g]
		sdev<-sdvec[g]
		l<-lvec[g]
		u<-uvec[g]
		if (sum((mdata<=u)*(mdata>=l))!=0)
		{
		normaldataVec<-c(gfct(t(mdata),mu,sdev,l,u))
		normaldataVec[which(normaldataVec==Inf)]<-100
		normaldataMat<-cbind(normaldataMat,normaldataVec)
		} else
		{
		normaldataMat<-cbind(normaldataMat,rep(0,length(mdata)))
		}
	}
	tauMat<-c()
	for (j in 1:Ncount)
	{
		tauMat<-cbind(tauMat,c(prior[j]*(t(normaldataMat[,j]))))
	}
	deno <- rowSums(tauMat)
  ind0<-which(deno==0)
  if (length(ind0)!=0)
  {
    for (k in 1:length(ind0))
    {
      indmin<-which.min(abs(muvec-mdata[ind0[k]]))
      tauMat[ind0[k],indmin]<-1
    }
    deno <- rowSums(tauMat)
  }
	taux<- tauMat/deno
taux
}

########### Maximization Step ###################
MStep<-function(mdata,taux,muvec,sdvec)
{
	N <- ncol(taux)
	moy<-c()
	sdev<-c()
	ptmp <- colSums(taux)
	pnew<-c()
	for (j in 1:N)
		{
		if (ptmp[j]!=0)
		{
			moy[j] <- muvec[j]
			#moy[j] <- (taux[,j]%*%mdata)/ptmp[j]
			if (sqrt((taux[,j]%*%((mdata-moy[j])^2))/ptmp[j])< 1e-100)
			{ 
				sdev[j]<-sdvec[j]
			} else
			{
			sdev[j]<- sqrt((taux[,j]%*%((mdata-moy[j])^2))/ptmp[j])
			}
			pnew[j]<-ptmp[j]/length(mdata)
		} else
		{
		moy[j] <- muvec[j]
		sdev[j]<- sdvec[j]
		pnew[j]<-1e-06
		}
		}
      pnew<-pnew/sum(pnew)
	ParamResult<-list()
	ParamResult$mu<-moy
	ParamResult$sdev<-sdev
	ParamResult$prior<-pnew
ParamResult
}




####### Format Data for FastCall ########
MakeData<-function(TotalTable, infoPos.StartEnd)
{
	if (infoPos.StartEnd) {
	MetaTable<-TotalTable[,1:4]
	} else
	{
	MetaTable<-TotalTable[,1:3]
	}
	if (infoPos.StartEnd) {
	NumericTable<-TotalTable[,5:ncol(TotalTable)]
	} else
	{
	NumericTable<-TotalTable[,4:ncol(TotalTable)]
	}
	NExp<-ncol(NumericTable)/2
	TableSeg<-as.matrix(NumericTable[,(NExp+1):(NExp*2)])
  
	SummaryData<-c()
	for (i in 1:NExp)
	{
	  segdata<-as.numeric(TableSeg[,i])
	  NData<-NumericTable[,i]
	  startseg<-c(1,(1+which(diff(segdata)!=0)))
	  endseg<-c(which(diff(segdata)!=0),length(segdata))
	  sdvec<-c()
	  for (j in 1:length(startseg))
	  {
	    sdvec<-c(sdvec,sd(NData[startseg[j]:endseg[j]]))
	  }
	  SummaryData<-rbind(SummaryData,cbind(rep.int(i,length(startseg)),startseg,endseg,segdata[startseg],sdvec))
	}

	DataList<-list()
	DataList$SummaryData<-SummaryData
	DataList$MetaTable<-MetaTable
	DataList
}


############# Make output data #####
DataRecomp<-function(SummaryData,MetaData,out)
{
	ExpClassTotal<-c()
	ExpProbTotal<-c()
	NumExp<-length(unique(SummaryData[,1]))
	for (i in 1:NumExp)
	{
		indExp<-which(SummaryData[,1]==i)
		StartEndMat<-rbind(SummaryData[indExp,2:3])
		OutClass<-out[indExp,1]
		OutProb<-out[indExp,2]
		ExpProb<-rep(0,StartEndMat[length(indExp),2])
		ExpClass<-rep(0,StartEndMat[length(indExp),2])
		for (j in 1:length(indExp))
		{
			ExpProb[StartEndMat[j,1]:StartEndMat[j,2]]<-OutProb[j]
			ExpClass[StartEndMat[j,1]:StartEndMat[j,2]]<-OutClass[j]
		}
	ExpClassTotal<-cbind(ExpClassTotal,ExpClass)
	ExpProbTotal<-cbind(ExpProbTotal,ExpProb)
	}
TotalTableFastCall<-cbind(MetaData,ExpClassTotal,ExpProbTotal)
TotalTableFastCall
}



### Draw legends (by CGHweb team) ###
vcbar <- function(zlim, colV) 
 {
  n <- length(colV)
  plot(NA, xlim = c(0,1), ylim = c(1,n), axes = FALSE, xlab = "", ylab = "",
       main = "", xaxs = "i", yaxs = "i")
  abline(h = 1:n, col = colV, lwd = 1)
  axis(4, at = round(seq(1,n,length.out=5)), las = 2, cex.axis = 0.6,
       label = formatC(seq(zlim[1], zlim[2], length.out = 5),
         format = "f", digits = 2))
  box()
  mtext("Log-Ratio", side = 3, line = 0.5, cex = 0.7, font = 2,
        at = 0.9)
}

gllegend <-
function(Class.names)  {
if (length(Class.names) > 2)
{
  plot(NA, xlab = "", ylab = "", main = "", axes = FALSE,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  legend(x = "center", col = c("red","orange", "green","green4"), pch = 15, cex = 0.9,
         pt.cex = 1.2, bg = "white", legend = Class.names, bty = "n")
} else
{
 plot(NA, xlab = "", ylab = "", main = "", axes = FALSE,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  legend(x = "center", col = c("orange", "green"), pch = 15, cex = 0.8,
         pt.cex = 1.2, bg = "white", legend = Class.names, bty = "n")

}
}

seglegend <-
function(Method)  {
  plot(NA, xlab = "", ylab = "", main = "", axes = FALSE,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  legend(x = "center", col = c("grey","black"), pch=c(19,1), lty=1,lwd=c(0,1.5), cex = 0.9,
         pt.lwd = c(0.9,0), bg = "white", legend = c("data",Method), bty = "n")
}


############# Assegnazione delle Label Statistiche ##########
LabelAss<-function(P0,mdata)
{
	CallResults<-rep(0,length(mdata))
	ProbResults<-c()
	
	indcall<-max.col(P0)
	CallResults[indcall==1]<- -2
	CallResults[indcall==2]<- -1
	CallResults[indcall==3]<- 0
	CallResults[indcall==4]<- 1
	CallResults[indcall==5]<- 2

	for (i in 1:length(mdata))
	{
	ProbResults<-c(ProbResults,P0[i,indcall[i]])
	}

Results<-c()
Results<-cbind(CallResults,ProbResults)
Results
}

###### Funzione per Salvare i risultati di FastCall in formato VCF 4.0 ####
VCFWindowCreate<-function(Assembly,DataFolder,ExpLabelOut,TargetFolder,SummaryData,MetaData,out,DataFilt)
{
  StringFormat<-'##fileformat=VCFv4.0'
  StringDate<-paste("##fileDate=",format(Sys.time(), "%Y%d%m"),sep="")
  StringSource<-'##source=EXCAVATOR2v1.0'
  StringReference<-paste("##reference=",Assembly,sep="")
  StringAssembly<-'##assembly=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta'
  StringINFO1<-'##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
  StringINFO2<-'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  StringINFO3<-'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
  StringINFO4<-'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
  StringFORMAT1<-'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  StringFORMAT2<-'##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">'
  StringFORMAT3<-'##FORMAT=<ID=CNF,Number=1,Type=Float,Description="Copy number genotype fraction for imprecise events">'
  StringFORMAT4<-'##FORMAT=<ID=FCL,Number=1,Type=Float,Description="FastCall inference Label">'
  StringFORMAT5<-'##FORMAT=<ID=FCP,Number=1,Type=Float,Description="FastCall Posterior Probability">'
  StringFORMAT6<-'##FORMAT=<ID=L2R,Number=1,Type=Float,Description="Normalized log2Ratio value">'
  StringALT1<-'##ALT=<ID=CNV,Description="Copy number variable region">'
  HeaderMat<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
  Header<-c(StringFormat,StringDate,StringSource,StringReference,StringAssembly,StringINFO1,StringINFO2,StringINFO3,StringINFO4,StringALT1,StringFORMAT1,StringFORMAT2,StringFORMAT3,StringFORMAT4,StringFORMAT5,StringFORMAT6)
  FormatField<-"GT:CN:CNF:FCL:FCP:L2R"
  HeadMatSample<-c(HeaderMat,ExpLabelOut)
  
  
  indSig<-which(out[,1]!=0)
  

  
  
  if (length(indSig)!=0)
  {
    
    
    outSig<-rbind(out[indSig,])
    P0Sig<-rbind(P0[indSig,])
    SummaryDataSig<-rbind(SummaryData[indSig,])
    CNFSig<-round(2*(2^SummaryDataSig[,4]),digit=2)
    CNSig<-round(CNFSig)
    
    Ntot<-sum((SummaryDataSig[,3]+1)-SummaryDataSig[,2])
    MatSig<-matrix(0,nrow=Ntot,ncol=8)

    indS<-1
    for (j in 1:length(indSig))
    {
      
    indP<-c(SummaryDataSig[j,2]:SummaryDataSig[j,3])
    indE<-indS+length(indP)-1
    MatSig[c(indS:indE),]<-cbind(MetaData[indP,1],MetaData[indP,3],MetaData[indP,4],DataFilt[indP],rep(CNFSig[j],length(indP)),rep(CNSig[j],length(indP)),rep(outSig[j,1],length(indP)),rep(round(outSig[j,2],digit=2),length(indP)))
    indS<-indE+1
    }
    
    
    ChromoSig<-MatSig[,1]
    StartSig<-MatSig[,2]
    EndSig<-MatSig[,3]
    l2rSig<-MatSig[,4]
    CNFSig<-MatSig[,5]
    CNSig<-MatSig[,6]
    CallSig<-MatSig[,7]
    ProbSig<-MatSig[,8]
    
    
    
    indEnde<-grep("e",EndSig)
    EndSige<-EndSig[indEnde]
    indStarte<-grep("e",StartSig)
    StartSige<-StartSig[indStarte]
    
    options(scipen = 10)
    
    StartSig[indStarte]<-as.character(as.numeric(StartSige))
    EndSig[indEnde]<-as.character(as.numeric(EndSige))
    L<-as.character(as.numeric(EndSig)-as.numeric(StartSig)+1)
    
    MatOut<-matrix(NA,ncol=10,nrow=length(ChromoSig))
    ChromoSigU<-unique(ChromoSig)
    
    indStart<-1
    for (i in 1:length(ChromoSigU))
    {
      
      FRBIn<-file.path(TargetFolder,"FRB",paste("FRB.",ChromoSigU[i],".RData",sep=""))
      load(FRBIn)
      FRBPos<-as.character(FRBData[,1])
      FRBRef<-as.character(FRBData[,2])
      
      indC<-which(ChromoSig==ChromoSigU[i])
      ChromoSigC<-ChromoSig[indC]
      EndSigC<-EndSig[indC]
      StartSigC<-StartSig[indC]
      CNSigC<-CNSig[indC]
      CNFSigC<-CNFSig[indC]
      CallSigC<-CallSig[indC]
      ProbSigC<-ProbSig[indC]
      l2rSigC<-l2rSig[indC]
      LC<-L[indC]
      RefC<-FRBRef[which(!is.na(match(FRBPos,StartSigC)))]
      AltC<-rep("<CNV>",length(RefC))
      indEnd<-indStart+length(indC)-1
      InfoField<-paste("IMPRECISE;SVTYPE=CNV;END=",EndSigC,";SVLEN=",LC,";",sep="")
      GenoField<-paste("1/1:",CNSigC,":",CNFSigC,":",CallSigC,":",ProbSigC,":",l2rSigC,sep="")
      MatOut[c(indStart:indEnd),]<-cbind(ChromoSigC,StartSigC,".",RefC,AltC,".","PASS",InfoField,FormatField,GenoField)
      indStart<-indEnd+1
    }
  }
  if (length(indSig)==0)
  {
    MatOut<-c()
  }
  MatOut<-rbind(HeadMatSample,MatOut)
  
  FileOut<-file.path(DataFolder,"Results",ExpLabelOut,paste("EXCAVATORWindowCall_",ExpLabelOut,".vcf",sep=""))
  
  zz <- file(FileOut, "w")
  cat(Header,file=zz,sep = "\n")
  write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
  close(zz)
  
}



###### Funzione per Salvare i risultati di FastCall per ogni Regione del target in formato VCF 4.0 ####
VCFRegionCreate<-function(Assembly,DataFolder,ExpLabelOut,TargetFolder,SummaryData,MetaData,out)
{
  StringFormat<-'##fileformat=VCFv4.0'
  StringDate<-paste("##fileDate=",format(Sys.time(), "%Y%d%m"),sep="")
  StringSource<-'##source=EXCAVATOR2v1.0'
  StringReference<-paste("##reference=",Assembly,sep="")
  StringAssembly<-'##assembly=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta'
  StringINFO1<-'##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
  StringINFO2<-'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  StringINFO3<-'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
  StringINFO4<-'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
  StringFORMAT1<-'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  StringFORMAT2<-'##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">'
  StringFORMAT3<-'##FORMAT=<ID=CNF,Number=1,Type=Float,Description="Copy number genotype fraction for imprecise events">'
  StringFORMAT4<-'##FORMAT=<ID=FCL,Number=1,Type=Float,Description="Label Inferred by FastCall algorithm">'
  StringFORMAT5<-'##FORMAT=<ID=FCP,Number=1,Type=Float,Description="Posterior Probability inferred by FastCall algorithm">'
  StringALT1<-'##ALT=<ID=CNV,Description="Copy number variable region">'
  HeaderMat<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
  Header<-c(StringFormat,StringDate,StringSource,StringReference,StringAssembly,StringINFO1,StringINFO2,StringINFO3,StringINFO4,StringALT1,StringFORMAT1,StringFORMAT2,StringFORMAT3,StringFORMAT4,StringFORMAT5)
  FormatField<-"GT:CN:CNF:FCL:FCP"
  HeadMatSample<-c(HeaderMat,ExpLabelOut)
  
  
  indSig<-which(out[,1]!=0)
  
  
  
  
  if (length(indSig)!=0)
  {
    
    
    outSig<-rbind(out[indSig,])
    P0Sig<-rbind(P0[indSig,])
    SummaryDataSig<-rbind(SummaryData[indSig,])
    CNFSig<-round(2*(2^SummaryDataSig[,4]),digit=2)
    CNSig<-round(CNFSig)
    
    
    MatSig<-matrix(0,nrow=length(indSig),ncol=7)
    
    
    for (j in 1:length(indSig))
    {
      
      indS<-SummaryDataSig[j,2]
      indE<-SummaryDataSig[j,3]
      MatSig[j,]<-cbind(MetaData[indS,1],MetaData[indS,3],MetaData[indE,4],CNFSig[j],CNSig[j],outSig[j,1],round(outSig[j,2],digit=2))
    }
    
    
    ChromoSig<-MatSig[,1]
    StartSig<-MatSig[,2]
    EndSig<-MatSig[,3]
    CNFSig<-MatSig[,4]
    CNSig<-MatSig[,5]
    CallSig<-MatSig[,6]
    ProbSig<-MatSig[,7]
    
    
    
    indEnde<-grep("e",EndSig)
    EndSige<-EndSig[indEnde]
    indStarte<-grep("e",StartSig)
    StartSige<-StartSig[indStarte]
    
    options(scipen = 10)
    
    StartSig[indStarte]<-as.character(as.numeric(StartSige))
    EndSig[indEnde]<-as.character(as.numeric(EndSige))
    L<-as.character(as.numeric(EndSig)-as.numeric(StartSig)+1)
    
    MatOut<-matrix(NA,ncol=10,nrow=length(ChromoSig))
    ChromoSigU<-unique(ChromoSig)
    
    indStart<-1
    for (i in 1:length(ChromoSigU))
    {
      
      FRBIn<-file.path(TargetFolder,"FRB",paste("FRB.",ChromoSigU[i],".RData",sep=""))
      load(FRBIn)
      FRBPos<-as.character(FRBData[,1])
      FRBRef<-as.character(FRBData[,2])
      
      indC<-which(ChromoSig==ChromoSigU[i])
      ChromoSigC<-ChromoSig[indC]
      EndSigC<-EndSig[indC]
      StartSigC<-StartSig[indC]
      CNSigC<-CNSig[indC]
      CNFSigC<-CNFSig[indC]
      CallSigC<-CallSig[indC]
      ProbSigC<-ProbSig[indC]
      LC<-L[indC]
      RefC<-FRBRef[which(!is.na(match(FRBPos,StartSigC)))]
      AltC<-rep("<CNV>",length(RefC))
      indEnd<-indStart+length(indC)-1
      InfoField<-paste("IMPRECISE;SVTYPE=CNV;END=",EndSigC,";SVLEN=",LC,";",sep="")
      GenoField<-paste("1/1:",CNSigC,":",CNFSigC,":",CallSigC,":",ProbSigC,sep="")
      MatOut[c(indStart:indEnd),]<-cbind(ChromoSigC,StartSigC,".",RefC,AltC,".","PASS",InfoField,FormatField,GenoField)
      indStart<-indEnd+1
    }
  }
  if (length(indSig)==0)
  {
    MatOut<-c()
  }
  MatOut<-rbind(HeadMatSample,MatOut)
  
  FileOut<-file.path(DataFolder,"Results",ExpLabelOut,paste("EXCAVATORRegionCall_",ExpLabelOut,".vcf",sep=""))
  
  zz <- file(FileOut, "w")
  cat(Header,file=zz,sep = "\n")
  write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
  close(zz)
  
}








