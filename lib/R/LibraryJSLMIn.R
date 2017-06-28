###################### Starting Parameter ####################
ParamEstSeq <- function(DataMatrix,omega)
{
	T=ncol(DataMatrix)
	NExp<-nrow(DataMatrix)
	sigmax<-c()
	mi<-c()
	
	for (i in 1:NExp)
	{
		thru<-as.numeric(quantile(DataMatrix[i,],0.99))
		thrd<-as.numeric(quantile(DataMatrix[i,],0.01))
		mi[i]<-0
		sigmax[i]<-var(DataMatrix[i,which(DataMatrix[i,]<=thru & DataMatrix[i,]>=thrd)])
	}
	
	smu<-sqrt(omega*sigmax)
	sepsilon<-sqrt((1-omega)*sigmax)
	Results<-list()
	Results$mi<-mi
	Results$smu<-smu
	Results$sepsilon<-sepsilon
	
	Results
}


###################### MUK Estimation ####################
MukEst <- function(DataMatrix,mw)
{
	NExp<-dim(DataMatrix)[1]
	if (NExp==1)
	{
		muk<-rbind(seq(-1,1,by=0.1))
	}
	if (NExp>1)
	{
		DataMatrix[which(DataMatrix>1)]<-1
		DataMatrix[which(DataMatrix< -1)]<- -1
		
		DataMatrixA<-c()
		for (i in 1:NExp)
		{
			DataMatrixA<-rbind(DataMatrixA,SMA(DataMatrix[i,], n=mw))
		}
		
		DataMatrixA<-DataMatrixA[,2:length(DataMatrixA[1,])]
		
		
		binsize=0.2
		binVec<-c(seq(-1,-0.2,by=binsize),0,seq(0.2,1,by=binsize))
		binmean<-c(seq(-0.9,-0.3,by=binsize),0,0,seq(0.3,0.9,by=binsize))
		
		
		DataQuant<-DataMatrixA
		
		for (i in 1:(length(binVec)-1))
		{
			DataQuant[which(DataMatrixA>binVec[i] & DataMatrixA<=binVec[i+1])]<-binmean[i]
		}
		
		muk<-unique(DataQuant,  MARGIN = 2)
		muk<-muk[,-1]
	}
	muk
}


################# SLM Inomogeneo #########
JointSegIn <- function(DataMatrix,muk,mi,smu,sepsilon,Pos,omega,eta,stepeta)
{
####  Calcolo il vettore delle covariate dipendenti dalla distanza ###
#stepeta<-200000
	CovPos<-diff(Pos)
	CovPosNorm<-CovPos/stepeta
	etavec<-eta+((1-eta)*exp(log(eta)/CovPosNorm))
	
	
### Calcolo i parametri del modello SLM
	NCov<-length(etavec)
	K0<-ncol(muk)
	etav<-log(rep(1,K0)*(1/K0))
	T=ncol(DataMatrix)
	NExp<-nrow(DataMatrix)
	
	
	
	
####  Calcolo le Matrici di Emissione e Transizione ####
	P<-matrix(data=0,nrow=K0,ncol=(K0*NCov))
	G<-matrix(data=0,nrow=K0,ncol=K0)
	emission<-matrix(data=0,nrow=K0,ncol=T)
	out<-.Fortran("transemisi",as.vector(muk),as.vector(mi),as.double(etavec),as.integer(NCov),as.matrix(DataMatrix),as.integer(K0),as.integer(NExp),as.vector(smu),as.vector(sepsilon),as.integer(T),as.matrix(G),as.matrix(P),as.matrix(emission))
	
	P<-out[[12]]
	emission<-out[[13]]
	
	
	
	
####### Algoritmo di Viterbi #########
	psi<-matrix(data=0,nrow=K0,ncol=T)
	path<-c(as.integer(rep(0,T)))
	out2<-.Fortran("bioviterbii",as.vector(etav),as.matrix(P),as.matrix(emission),as.integer(T),as.integer(K0),as.vector(path),as.matrix(psi))
	s<-out2[[6]]
	
	
	sortResult <- SortState(s)
	TotalPredBreak<-sortResult[[3]]
	TotalPredBreak
}

################# Funzione per riordinare gli stati della segmentazione #########
SortState <- function(s)
{
	l<-1
	seg<-c()
	brek<-c()
	t<-1
	for (k in 1:(length(s)-1))
	{
		if (s[k]!=s[k+1])
		{
			brek[t]<-k
			t<-t+1
			if (length(which(seg==s[k]))==0)
			{
				seg[l]<-s[k]
				l<-l+1
			}
		}
	}
	brek<-c(0,brek,length(s))
	if (length(which(seg==s[length(s)]))==0)
	{
		seg<-c(seg,s[length(s)])
	}
	
	s0<-c()
	for (k in 1:length(seg))
	{
		s0[which(s==seg[k])]<-k
	}
	
	SortResult<-list()
	SortResult[[1]]<-s0
	SortResult[[2]]<-seg
	SortResult[[3]]<-brek
	SortResult
	
}

################# Funzione calcolare i valori medi dei segmenti e le regioni alterate dai dati di breakpoint #########
SegResults <- function(DataSeq,TotalPredBreak)
{
	
	regioncall<-c()
	TotalPred<-c()
	NExp<-nrow(DataSeq)
	for (j in 1:NExp)
	{
		s<-rep(0,ncol(DataSeq))
		for (i in 1:(length(TotalPredBreak)-1))
		{
			s[(TotalPredBreak[i]+1):TotalPredBreak[i+1]]<-median(DataSeq[j,(TotalPredBreak[i]+1):TotalPredBreak[i+1]])
		}
		TotalPred<-rbind(TotalPred,s)
	}
	
	Result<-TotalPred
	Result
}


###### Funzione per filtrare i piccoli shift #####
FilterSeg <- function(TotalPredBreak,FW)
{
  controllength<-diff(TotalPredBreak)
	
	indF<-which(controllength<=FW)
	if (length(indF)!=0)
	{
    if (indF[1]==1)
    {
      indF[1]<-2
      indF<-unique(indF)
      TotalPredBreak1 <- TotalPredBreak[-(indF)]
    }
    if (indF[1]!=1)
    {
      TotalPredBreak1 <- TotalPredBreak[-(indF)]
    }
	}
	if (length(indF)==0)
	{
		TotalPredBreak1<-TotalPredBreak
	}
	TotalPredBreak1
}



###### Funzioni per lo smooth e la rimozione degli outlier ###
trimmed.variance <- function(genomdat, trim=0.025)
{
    n <- length(genomdat)
    n.keep <- round((1-2*trim)*(n-1))
    inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))
}


inflfact <- function(trim)
{
    a <- qnorm(1-trim)
    x <- seq(-a,a,length=10001)
    x1 <- (x[-10001] + x[-1])/2
    1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
}

filterOut<- function(genomdat, smooth.region,outlier.SD.scale,smooth.SD.scale,trim)
{
	
	genomdatOut<-genomdat
	trimmed.SD <- sqrt(trimmed.variance(genomdat, trim))
	oSD <- outlier.SD.scale*trimmed.SD
	sSD <- smooth.SD.scale*trimmed.SD
	k <- smooth.region
	nbhd<-c(-k:-1,1:k)
	n<-length(genomdat)
	for (i in 1:n)
	{
		xi <- genomdat[i]
		nbhd <- 1+nbhd
		xnbhd <- genomdat[nbhd[nbhd>0 & nbhd <=n]]
		if (xi > max(xnbhd) + oSD)
		{
			xi <- median(c(xi,xnbhd)) + sSD
		}
		
		if (xi < min(xnbhd) - oSD)
		{
			xi <- median(c(xi,xnbhd)) - sSD
		}
		genomdatOut[i]<-xi
	}
	genomdatOut
}




