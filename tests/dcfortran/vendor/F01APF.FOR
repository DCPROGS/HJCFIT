      SUBROUTINE F01APF (N, LOW, IUPP, INT, H, IH, V, IV)     
C    MARK 2 RELEASE. NAG COPYRIGHT 1972   
C     MARK 4 REVISED. 
C 
C        DIRTRANS     
C 
C        1ST AUGUST 1971        
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(IV, N), H(IH, N), INT(N)
      DO 1 I = 1, N   
      DO 2 J = 1,N    
         V(I, J) = 0.0
    2 CONTINUE        
         V(I, I) = 1.0
    1 CONTINUE        
         LOW1 = LOW + 1         
      IF (LOW1.GT.IUPP) RETURN  
      DO 3 II = LOW1, IUPP      
         I = LOW1 + IUPP - II   
         I1 = I - 1   
      IF (LOW1.GT.I1) GOTO 6    
      DO 4 J = LOW1, I1         
         V(I, J) = H(I, J-1)    
    4 CONTINUE        
    6 M=INT(I)        
      IF( M.EQ.I )  GOTO 3      
      DO 5 J = LOW1, IUPP       
         X = V(M, J)  
         V(M, J) = V(I, J)      
         V(I, J) = X  
    5 CONTINUE        
    3 CONTINUE        
      RETURN
      END   
