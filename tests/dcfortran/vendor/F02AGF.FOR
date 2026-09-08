      SUBROUTINE F02AGF(A,IA,N,RR,RI,VR,IVR,VI,IVI,INT,IFAIL) 
C    MARK 2 RELEASE. NAG COPYRIGHT 1972   
C     MARK 3 REVISED. 
C     MARK 4 REVISED. ER *M1A-29
C 
C     EIGENVALUES AND EIGENVECTORS OF REAL UNSYMMETRIC MATRIX 
C 
C        1ST AUGUST 1971        
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 MACHEP,MAX 
C 4 ABS->DABS. SQRT->DSQRT
      DIMENSION A(IA,N),VR(IVR,N),VI(IVI,N),RR(N),RI(N), INT(N)
      ISAVE = IFAIL   
      IFAIL = 1       
      MACHEP=2.0**(-47)
      CALL F01AKF(N,1,N,A,IA,INT)         
      CALL F01APF(N,1,N,INT,A,IA,VR,IVR)  
      CALL F02AQF(N,1,N,MACHEP,A,IA,VR,IVR,RR,RI,INT,IFAIL)   
      IF (IFAIL.EQ.0) GO TO 100
      RETURN
  100 DO 1 I = 1,N    
      IF(RI(I).EQ.0.0) GOTO 2   
      IF (RI(I).GT.0.0) GOTO 3  
      DO 4 J=1,N      
      VR(J,I)=VR(J,I-1)         
      VI(J,I)=-VI(J,I-1)        
    4 CONTINUE        
      GOTO 1
    2 DO 5 J=1,N      
      VI(J,I)=0.0     
    5 CONTINUE        
      GOTO 1
    3 DO 6 J=1,N      
      VI(J,I)=VR(J,I+1)         
    6 CONTINUE        
    1 CONTINUE        
      DO 7 I=1,N      
      SUM=0.0         
      MAX=0.0         
      DO 12 J=1,N     
      IF (DABS(VR(J,I)).LE.MAX) GOTO 13    
      MAX=DABS(VR(J,I))
   13 IF (DABS(VI(J,I)).LE.MAX) GOTO 12    
      MAX=DABS(VI(J,I))
   12 CONTINUE        
      DO 11 J=1,N     
      VR(J,I)=VR(J,I)/MAX       
      VI(J,I)=VI(J,I)/MAX       
   11 CONTINUE        
      MAX=0.0         
      DO 8 J=1,N      
      TERM=VR(J,I)**2+VI(J,I)**2
      SUM=SUM+TERM    
      IF (TERM.LE.MAX) GOTO 9   
      MAX=TERM        
      C=VR(J,I)       
      D=-VI(J,I)      
    9 CONTINUE        
    8 CONTINUE        
      SUM=SUM*(C**2+D**2)       
      SUM=DSQRT(SUM)   
      DO 10 J=1,N     
      TERM=VR(J,I)    
      VR(J,I)=(VR(J,I)*C-VI(J,I)*D)/SUM   
      VI(J,I)=(D*TERM+C*VI(J,I))/SUM      
   10 CONTINUE        
    7 CONTINUE        
      RETURN
      END   
