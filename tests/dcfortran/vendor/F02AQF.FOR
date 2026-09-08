      SUBROUTINEF02AQF(N,LOW,UPP,MACHEP,H,IH,VECS,IVECS,WR,WI,CNT,IFAIL)
C    MARK 2 RELEASE. NAG COPYRIGHT 1972   
C     MARK 3 REVISED. 
C     MARK 4 REVISED. ER AM3-29 
C 
C        HQR2         
C 
C        1ST DECEMBER 1971      
C 
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(IH,N),VECS(IVECS,N),WR(N),WI(N),CNT(N)
      INTEGER CNT,UPP,EN,EN2,UPP1  
      REAL*8 NORM,MACHEP
      LOGICAL NOTLAS
      ISAVE= IFAIL    
C ISOLATED ROOTS       
      IF (LOW.LE.1) GOTO 1      
      J=LOW-1         
      DO 2 I=1,J      
      WR(I)=H(I,I)    
      WI(I)=0.0       
      CNT(I)=0        
    2 CONTINUE        
    1 IF (UPP.GE.N) GOTO 3      
      J=UPP+1         
      DO 4 I=J,N      
      WR(I)=H(I,I)    
      WI(I)=0.0       
      CNT(I)=0        
    4 CONTINUE        
    3 EN=UPP
      T=0.0 
    5 IF (EN.LT.LOW) GOTO 100   
      ITS=0 
      NA=EN-1         
C   LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT         
    6 IF (LOW+1.GT.EN) GOTO 8   
      LOW1=LOW+1      
      DO 7 LL=LOW1,EN 
      L=EN+LOW1-LL    
      IF(DABS(H(L,L-1)).LE.(MACHEP *(DABS(H(L-1,L-1))+DABS(H(L,L)))))
     2 GOTO 9
    7 CONTINUE        
    8 L=LOW 
    9 X=H(EN,EN)      
      IF (L.EQ.EN) GOTO 30      
      Y=H(NA,NA)      
      W=H(EN,NA)*H(NA,EN)       
      IF (L.EQ.NA) GOTO 31      
      IF (ITS.EQ.30) GOTO 101   
C   FORM SHIFT 
      IF ((ITS.NE.10).AND.(ITS.NE.20)) GOTO 10      
      T=T+X 
      IF (LOW.GT.EN) GOTO 80    
      DO 11 I=LOW,EN  
      H(I,I)=H(I,I)-X 
   11 CONTINUE        
C***********DABS
   80 S=DABS(H(EN,NA))+DABS(H(NA,EN-2))     
      X=0.75*S        
      Y=X   
      W=-0.4375*S**2  
   10 ITS=ITS+1       
C   LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS         
      IF (L.GT.EN-2) GOTO 13    
      EN2=EN-2        
      DO 12 MM=L,EN2  
      M=L+EN2-MM      
      Z=H(M,M)        
      R=X-Z 
      S=Y-Z 
      P=(R*S-W)/H(M+1,M)+H(M,M+1)         
      Q=H(M+1,M+1)-Z-R-S        
      R=H(M+2,M+1)    
C*******DABS
      S=DABS(P)+DABS(Q)+DABS(R)    
      P=P/S 
      Q=Q/S 
      R=R/S 
      IF (M.EQ.L) GOTO13        
C*******DABS- ALTERED CORRECTLY?
      IF((DABS(H(M,M-1))*(DABS(Q)+DABS(R))).LE.(MACHEP *DABS(P)
     2 *(DABS(H(M-1,M-1))+DABS(Z)+DABS(H(M+1,M+1))))) GOTO 13       
   12 CONTINUE        
   13 M2=M+2
      IF (M2.GT.EN) GOTO 81     
      DO 14 I=M2,EN   
      H(I,I-2)=0.0    
   14 CONTINUE        
   81 M3=M+3
      IF (M3.GT.EN) GOTO 82     
      DO 15 I=M3,EN   
      H(I,I-3)=0.0    
   15 CONTINUE        
   82 IF (M.GT.NA) GOTO 83      
      DO 16 K=M,NA    
      NOTLAS =K.NE.NA 
      IF (K.EQ.M) GOTO 17       
      P=H(K,K-1)      
      Q=H(K+1,K-1)    
      R=0   
      IF (NOTLAS ) R=H(K+2,K-1) 
C*******DABS
      X=DABS(P)+DABS(Q)+DABS(R)    
      IF (X.EQ.0.0) GOTO 16     
      P=P/X 
      Q=Q/X 
      R=R/X 
C*******DSQRT
   17 S=DSQRT(P**2+Q**2+R**2)    
      IF (P.LT.0.0) S=-S        
      IF (K.NE.M) GOTO 19       
      IF (L.NE.M) H(K,K-1)=-H(K,K-1)
      GOTO 20         
   19 H(K,K-1)=-S*X   
   20 P=P+S 
      X=P/S 
      Y=Q/S 
      Z=R/S 
      Q=Q/P 
      R=R/P 
C   ROW MODIFICATION     
      DO 21 J=K,N     
      P=H(K,J)+Q*H(K+1,J)       
      IF (.NOT.NOTLAS ) GOTO 22 
      P=P+R*H(K+2,J)  
      H(K+2,J)=H(K+2,J)-P*Z     
   22 H(K+1,J)=H(K+1,J)-P*Y     
      H(K,J)=H(K,J)-P*X         
   21 CONTINUE        
      J=EN  
      IF (K+3.LT.EN) J=K+3      
C   COLUMN MODIFICATION  
      DO 23 I=1,J     
      P=X*H(I,K)+Y*H(I,K+1)     
      IF (.NOT.NOTLAS ) GOTO 24 
      P=P+Z*H(I,K+2)  
      H(I,K+2)=H(I,K+2)-P*R     
   24 H(I,K+1)=H(I,K+1)-P*Q     
      H(I,K)=H(I,K)-P 
   23 CONTINUE        
C   ACCUMULATE TRANSFORMATIONS     
      IF (LOW.GT.UPP) GOTO 16   
      DO 25 I=LOW,UPP 
      P=X*VECS(I,K)+Y*VECS(I,K+1)         
      IF (.NOT.NOTLAS ) GOTO 27 
      P=P+Z*VECS(I,K+2)         
      VECS(I,K+2)=VECS(I,K+2)-P*R         
   27 VECS(I,K+1)=VECS(I,K+1)-P*Q         
      VECS(I,K)=VECS(I,K)-P     
   25 CONTINUE        
   16 CONTINUE        
   83 GOTO 6
C   ONE ROOT FOUND       
   30 WR(EN)=X+T      
      H(EN,EN)=WR(EN) 
      WI(EN)=0.0      
      CNT(EN)=ITS     
      EN =NA
      GOTO 5
C   TWO ROOTS FOUND      
   31 P=(Y-X)/2.0     
      Q=P**2+W        
C******DABS,DSQRT
      Z=DSQRT(DABS(Q))  
      X=X+T 
      H(EN,EN)=X      
      H(NA,NA)=Y+T    
      CNT(EN)=-ITS    
      CNT(NA)=ITS     
      IF (Q.LT.0.0) GOTO 28     
C   REAL PAIR  
      IF (P.LT.0.0) Z=P-Z       
      IF (P.GT.0.0) Z=P+Z       
      WR(NA)=X+Z      
      WR(EN)=WR(NA)   
      IF(Z.NE.0.0) WR(EN)=X-W/Z 
      WI(NA)=0.0      
      WI(EN)=0.0      
      X=H(EN,NA)      
C*****DSQRT
      R=DSQRT(X**2+Z**2)         
      P=X/R 
      Q=Z/R 
C   ROW MODIFICATION     
      DO 29 J=NA,N    
      Z=H(NA,J)       
      H(NA,J)=Q*Z+P*H(EN,J)     
      H(EN,J)=Q*H(EN,J)-P*Z     
   29 CONTINUE        
C   COLUMN MODIFICATION  
      DO 32 I=1,EN    
      Z=H(I,NA)       
      H(I,NA)=Q*Z+P*H(I,EN)     
      H(I,EN)=Q*H(I,EN)-P*Z     
   32 CONTINUE        
C   ACCUMULATE TRANSFORMATIONS     
      DO 33 I=LOW,UPP 
      Z=VECS(I,NA)    
      VECS(I,NA)=Q*Z+P*VECS(I,EN)         
      VECS(I,EN)=Q*VECS(I,EN)-P*Z         
   33 CONTINUE        
      GOTO 34         
C   COMPLEX PAIR         
   28 WR(NA)=X+P      
      WR(EN)=X+P      
      WI(NA)=Z        
      WI(EN)=-Z       
   34 EN=EN-2         
      GOTO 5
C   ALL ROOTS FOUND NOW BACKSUBSTITUTE       
  100 NORM=0.0        
      K=1   
      DO 35 I=1,N     
      DO 36 J=K,N     
C*****DABS
   36 NORM=NORM+DABS(H(I,J))     
      K=I   
   35 CONTINUE        
      NORM=NORM*MACHEP
C   BACKSUBSTITUTION     
      DO 37 KK=1,N    
      EN=N+1-KK
      P=WR(EN)        
      Q=WI(EN)        
      NA=EN-1         
      IF (Q.NE.0.0) GOTO 38     
C   REAL VECTOR
      M=EN  
      H(EN,EN)=1.0    
      IF (NA.LT.1) GOTO 37      
      DO 39 II=1,NA   
      I=NA+1-II       
      W=H(I,I)-P      
      R=H(I,EN)*H(EN,EN)        
      IF (M.GT.NA) GOTO 40      
      DO 41 J=M,NA    
   41 R=R+H(I,J)*H(J,EN)        
   40 IF (WI(I).GE.0.0) GOTO 42 
      Z=W   
      S=R   
      GOTO 45         
   42 M=I   
      IF (WI(I).GT.0.0) GOTO 44 
      IF(W.NE.0.0) GOTO 70      
C****DABS
      IF (DABS(R).LT.10.0*NORM) GOTO 71
      H(I,EN)=-R      
      I1=I+1
      DO 72 J=I1,EN   
      H(J,EN)=H(J,EN)*NORM
   72 CONTINUE        
      GOTO 45 
   70 H(I,EN)=-R/W    
      GOTO 45 
   71 H(I,EN)=-R/NORM 
      GOTO 45         
C   SOLVE REAL EQUATIONS 
   44 X=H(I,I+1)      
      Y=H(I+1,I)      
      Q=(WR(I)-P)**2+WI(I)**2   
      T=(X*S-Z*R)/Q   
      H(I,EN)=T       
      H(I+1,EN)=(-S-Y*T)/Z      
      IF (DABS(X).GT.DABS(Z)) H(I+1,EN)=(-R-W*T)/X    
C*******DABS
   45 CONTINUE        
   39 CONTINUE        
C  END REAL VECTOR      
      GOTO 37         
   38 IF (Q.GT.0.0) GOTO 37     
C  COMPLEX VECTOR ASSOCIATED WITH LAMBDA=P-I*Q        
      M=NA  
C*******DABS
      IF (DABS(H(EN,NA)).LE.DABS(H(NA,EN))) GOTO 48   
      H(NA,NA)=-(H(EN,EN)-P)/H(EN,NA)     
      H(NA,EN)=-Q/H(EN,NA)      
      GOTO 49         
   48 ZZ=-H(NA,EN)    
      CALL A02ACF(ZZ,0.0,H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))      
   49 H(EN,NA)=1.0    
      H(EN,EN)=0.0    
      IF (NA.LT.2) GOTO 37      
      NA1=NA-1        
      DO 99 II=1,NA1  
      I=1+NA1-II      
      W=H(I,I)-P      
      RA=H(I,EN)      
      SA=0.0
      IF (M.GT.NA) GOTO 52      
      DO 53 J=M,NA    
      RA=RA+H(I,J)*H(J,NA)      
      SA=SA+H(I,J)*H(J,EN)      
   53 CONTINUE        
   52 IF (WI(I).GE.0.0) GOTO 54 
      Z=W   
      R=RA  
      S=SA  
      GOTO 99         
   54 M=I   
      IF (WI(I).EQ.0.0) GOTO 56 
C   SOLVE COMPLEX EQUATIONS        
      X=H(I,I+1)      
      Y=H(I+1,I)      
      VR=(WR(I)-P)**2+WI(I)**2-Q**2       
      VI=(WR(I)-P)*2.0*Q        
C******DABS
      IF ((VR.EQ.0.0).AND.(VI.EQ.0.0)) VR=MACHEP *NORM*(DABS(W)+
     2DABS(Q)+DABS(X)+DABS(Y)+DABS(Z))       
      CALL A02ACF(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA,VR,VI,H(I,NA),H(I,EN))    
C*******DABS
      IF (DABS(X).LE.DABS(Z)+DABS(Q)) GOTO 57
      H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X         
      H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X         
      GOTO 99         
   57 CALL A02ACF(-R-Y*H(I,NA),-S-Y*H(I,EN),Z,Q,H(I+1,NA),H(I+1,EN))    
      GOTO 99         
   56 CALL A02ACF(-RA,-SA,W,Q,H(I,NA),H(I,EN))      
   99 CONTINUE        
C  END COMPLEX VECTOR   
   37 CONTINUE        
C   END BACKSUBSTITUTION 
C   VECTORS OF ISOLATED ROOTS      
   51 IF (LOW.LE.1) GOTO 59     
      LOW1=LOW-1      
      DO 58 I=1,LOW1  
      I1=I+1
      DO 58 J=I1,N    
   58 VECS(I,J)=H(I,J)
   59 IF (UPP.GE.N) GOTO 60     
      UPP1=UPP+1      
      DO 61 I=UPP1,N  
      I1=I+1
      DO 61 J=I1,N    
   61 VECS(I,J)=H(I,J)
C   MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C   VECTORS OF ORIGINAL FULL MATRIX
   60 DO 62 JJ=LOW,N  
      J=LOW+N-JJ      
      M=UPP 
      IF (J.LE.UPP) M=J         
      L=J-1 
      IF (WI(J).GE.0.0) GOTO 63 
      DO 64 I=LOW,UPP 
      Y=0.0 
      Z=0.00
      IF (LOW.GT.M) GOTO 66     
      DO 65 K=LOW,M   
      Y=Y+VECS(I,K)*H(K,L)      
   65 Z=Z+VECS(I,K)*H(K,J)      
   66 VECS(I,L)=Y     
      VECS(I,J)=Z     
   64 CONTINUE        
      GOTO 62         
   63 IF (WI(J).NE.0.0) GOTO 62 
      DO 68 I=LOW,UPP 
      Z=0.0 
      IF (LOW.GT.M) GOTO 68     
      DO 69 K=LOW,M   
   69 Z=Z+VECS(I,K)*H(K,J)      
   68 VECS(I,J)=Z     
   62 CONTINUE        
      IFAIL=0         
      RETURN
101   CONTINUE		!IFAIL NOT ZERO-FAULT SOMEWHERE
      RETURN
      END
