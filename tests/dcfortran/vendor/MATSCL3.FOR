	subroutine MATSCL3(phi,Amat,m,end,k,SCAL,kp,kq,ke)
C
c MATSCL3 is version of MATSCL in which matrix is 3D, Amat(i,j,m)
c and m is specified as a parameter. Top left ka*ka elements used
c
C CONVERTS MATRIX TO A SCALAR BY PREMULTIPLYING
C BY VECTOR PHI ( 1 X k) AND POSTMULTIPLYING BY VECTOR
C END (k X 1). RESULT IN SCAL.
C	KP= DECLARED DIM OF PHI(1,KP)
C	KQ= DECLARED DIM OF Q(KQ,KQ)
C	KE= DECLARED DIM OF END(KE,1)
C
	REAL*8 PHI(1,KP),Amat(kq,kq,kq),END(KE,1),SCAL,ZERO
	REAL*8 T
C
	ZERO=0.0D0
	SCAL=ZERO
	DO 1 I=1,K
	T=ZERO
	DO 2 J=1,K
	T=T+amat(i,j,m)*end(j,1)	!ACCUM ITH ELEMENT OF amat*END
c	print 3,I,J,amat(I,J,m),END(J,1)
c3	FORMAT(2I5,2G13.6)
2	CONTINUE
	SCAL=SCAL+PHI(1,I)*T
c	print 4,I,PHI(1,I),SCAL
c4	FORMAT(' I,PHI(I),S= ',I8,2G13.6)
1	CONTINUE
	RETURN
	END


