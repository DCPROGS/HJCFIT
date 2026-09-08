	subroutine MATSCL2(phi,QAF,col,kA,kF,scal,KP,KQ,KC)
C Variant of MATSCL for non-square matrices
C CONVERTS MATRIX Q (kA X kF) TO A SCALAR BY PREMULTIPLYING
C BY VECTOR PHI ( 1 X kA) AND POSTMULTIPLYING BY VECTOR col (kF X 1).
c RESULT IN SCAL.
C	KP= DECLARED DIM OF PHI(1,KP)
C	KQ= DECLARED DIM OF Q(KQ,KQ)
C	KE= DECLARED DIM OF col(KE,1)
C
	REAL*8 phi(1,KP),QAF(KQ,KQ),col(KC,1),scal,zero
	REAL*8 T
C
	zero=0.0d0
	scal=zero
	do 1 i=1,kA
	T=zero
	do 2 j=1,kF
	T=T+QAF(i,j)*col(j,1)
2	continue
	scal=scal+phi(1,i)*T
1	continue
	RETURN
	END


