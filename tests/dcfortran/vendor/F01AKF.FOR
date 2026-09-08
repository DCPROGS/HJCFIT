      SUBROUTINE F01AKF(N,K,L,A,IA,INT)   
C    MARK 2 RELEASE. NAG COPYRIGHT 1972   
C     MARK 4 REVISED. 
C 
C        DIRHES       
C 
C        1ST AUGUST 1971        
C 
C      DOUBLE PRECISION YD       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C ABS->DABS AND DBLE,SNGL CALLS REMOVED
      DIMENSION A(IA,N),INT(N)  
      K1=K+1
      IF (K1.GT.L) RETURN       
      DO 1 J=K1,L     
      M=J   
      X=0.0 
      DO 2 I=J,L      
      IF(DABS(A(I,J-1)).LE.DABS(X)) GO TO 2 
      X=A(I,J-1)      
      M=I   
    2 CONTINUE        
      INT(J)=M        
      IF(M.EQ.J) GO TO 3        
      DO 4 I=K,N      
      Y=A(M,I)        
      A(M,I)=A(J,I)   
    4 A(J,I)=Y        
      DO 5 I=1,L      
      Y=A(I,M)        
      A(I,M)=A(I,J)   
    5 A(I,J)=Y        
    3 IF(X.EQ.0.0) GO TO 6      
      IF(J.EQ.L) GO TO 6        
      J1=J+1
      DO 7 I=J1,L     
    7 A(I,J-1)=A(I,J-1)/X       
    6 DO 8 I=1,L      
      Y=A(I,J)        
      YD=Y		!INSTEAD OF NEXT LINE IN DP VERSION
C     YD=DBLE(Y)      
      IF(X.EQ.0.0) GO TO 9      
      IF(J.EQ.L) GO TO 9        
      J1=J+1
      DO 10 M=J1,L    
C  10 YD=YD+DBLE(A(I,M))*A(M,J-1)
10    YD=YD+A(I,M)*A(M,J-1)	!REPLACES PREV LINE IN DP VERSION
    9 IJ=J  
      IF(I.LE.J) IJ=I-1         
      YD=-YD
      IF((K+1).GT.IJ) GO TO 12  
      DO 11 M=K1,IJ   
C  11 YD=YD+DBLE(A(I,M-1))*A(M,J)         
11    YD=YD+A(I,M-1)*A(M,J)	!DITTO
C  12 Y=SNGL(YD)
12    Y=YD			!DITTO
    8 A(I,J)=-Y       
    1 CONTINUE        
      RETURN
      END   
