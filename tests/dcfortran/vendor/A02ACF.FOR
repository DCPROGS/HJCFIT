      SUBROUTINE A02ACF(XXR, XXI, YYR, YYI, ZR, ZI) 
C    MARK 2 RELEASE. NAG COPYRIGHT 1972   
C 
C        CDIV         
C 
C        1ST AUGUST 1971        
C 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C TWO ABS CALLS CHANGED TO DABS
      XR=XXR
      XI=XXI
      YR=YYR
      YI=YYI
      IF(DABS(YR).LE.DABS(YI)) GO TO 1      
      H=YI/YR         
      YR=H*YI+YR      
      ZR=(XR+H*XI)/YR 
      ZI=(XI-H*XR)/YR 
      GO TO 2         
    1 H=YR/YI         
      YI=H*YR+YI      
      ZR=(H*XR+XI)/YI 
      ZI=(H*XI-XR)/YI 
    2 RETURN
      END
