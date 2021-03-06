       SUBROUTINE BIOVITERBII(ETAV,P,EMISSION,T,KTILDE,PATH,PSI)
       
	   IMPLICIT NONE
       INTEGER T,KTILDE,IND,PATH(T),PSI(KTILDE,T),I,J,K,COUNTP
       DOUBLE PRECISION NORM,NORM1,NUMMAX
       DOUBLE PRECISION ETAV(KTILDE),EMISSION(KTILDE,T)
       DOUBLE PRECISION P(KTILDE,KTILDE),DELTA(KTILDE,T)
       DOUBLE PRECISION NDELTA(KTILDE,T),PDELTA(KTILDE,T)
       
       
       DO 202 I=1,KTILDE
          DELTA(I,1)=ETAV(I)+EMISSION(I,1)
          PSI(I,1)=0
 202   CONTINUE

       COUNTP=0
       DO 203 I=2,T
          DO 213 J=1,KTILDE
             NUMMAX=0.0
             NUMMAX=DELTA(1,I-1)+P(1,(J+COUNTP))
             IND=1
             DO 223 K=2,KTILDE
                IF((DELTA(K,I-1)+P(K,(J+COUNTP))).GT.NUMMAX) THEN
                   NUMMAX=(DELTA(K,I-1)+P(K,(J+COUNTP)))
                   IND=K
                ENDIF
 223         CONTINUE
             PSI(J,I)=IND
             DELTA(J,I)=NUMMAX+EMISSION(J,I)
 213      CONTINUE
          COUNTP=COUNTP+KTILDE
 203   CONTINUE
       
       NUMMAX=0.0
       NUMMAX=DELTA(1,T)
       IND=1
       DO 253 K=2,KTILDE
          IF(DELTA(K,T).GT.NUMMAX) THEN
             NUMMAX=DELTA(K,T)
             IND=K
           ENDIF
 253   CONTINUE
       
       PATH(T)=IND
       DO 263 K=T-1,1,-1
          PATH(K)=PSI(PATH(K+1),K+1)
 263   CONTINUE

       RETURN
       END 
	



       SUBROUTINE TRANSEMISI(MUK,MI,ETA,NCOV,TOTALSEQ,KTILDE,
     c NUMSEQ,SMU,SEPSILON,T,G,P,EMISSION)

	   IMPLICIT NONE

       INTEGER T,KTILDE,NUMSEQ,I,J,K,COUNTP,NCOV
       DOUBLE PRECISION NORM,GSUP,ETA(NCOV),PI,ELNSUM
       DOUBLE PRECISION MI(NUMSEQ),P(KTILDE,KTILDE*NCOV)
       DOUBLE PRECISION GVECT(KTILDE)
       DOUBLE PRECISION G(KTILDE,KTILDE)
       DOUBLE PRECISION EMISSION(KTILDE,T)
       DOUBLE PRECISION MUK(NUMSEQ,KTILDE)
       DOUBLE PRECISION TOTALSEQ(NUMSEQ,T)
       DOUBLE PRECISION SMU(NUMSEQ),SEPSILON(NUMSEQ)
       PARAMETER(PI=3.14159265358979)
       
      
       DO 627 I=1,KTILDE
             GSUP=0
             DO 617 J=1,NUMSEQ
                GSUP=GSUP+(-(((MUK(J,I)-MI(J))**2)/(2*SMU(J)**2)))
 617         CONTINUE
             GVECT(I)=GSUP
 627   CONTINUE
       
       NORM=GVECT(1)
       IF (KTILDE.GT.1) THEN
       DO 607 J=2,KTILDE
          NORM=ELNSUM(NORM,GVECT(J))
 607   CONTINUE
       ENDIF

       DO 700 I=1,KTILDE
          DO 702 J=1,KTILDE
             G(I,J)=GVECT(J)-NORM
 702      CONTINUE
 700   CONTINUE
              
       COUNTP=0
       DO 710 I=1,NCOV
           DO 720 J=1,KTILDE
               DO 750 K=1,KTILDE
                   IF (J.EQ.K) THEN
                      P(J,(K+COUNTP))=ELNSUM(LOG(1-ETA(I)),(LOG(ETA(I))+
     c                G(J,K)))
                      ELSE
                      P(J,(K+COUNTP))=LOG(ETA(I))+G(J,K)
                   ENDIF
 750           CONTINUE
 720       CONTINUE
           COUNTP=COUNTP+KTILDE
 710   CONTINUE


       DO 730 J=1,KTILDE
          DO 740 K=1,T
             DO 741 I=1,NUMSEQ
              EMISSION(J,K)=EMISSION(J,K)+LOG(1/
     c        (SQRT(2*PI)*SEPSILON(I)))+(-0.5*((TOTALSEQ(I,K)-
     c        MUK(I,J))/(SEPSILON(I)))**2)
 741         CONTINUE
 740      CONTINUE
 730   CONTINUE
       RETURN
       END 
	






      DOUBLE PRECISION FUNCTION ELNSUM(X,Y)
	  
	  IMPLICIT NONE
      DOUBLE PRECISION X,Y
      IF (X.GT.Y) THEN
         ELNSUM=X+LOG(1+EXP(Y-X))
      ELSE
         ELNSUM=Y+LOG(1+EXP(X-Y))
      ENDIF

      RETURN
      END

