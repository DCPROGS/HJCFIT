	subroutine RAC3d(ROW,COL,amat,I1,I2,J1,J2,MAX,W,W1,
     &	kdr,kdc,kda1,kda2,kda3,kdw)
C RAC3d is version of RAC3 that allows declared dimensions of AMAT
c to be unequal, and specifies dimension of w and w1 too.
C RAC3 is version of RAC1 that uses AMAT rather than EM.EN
c
C SUBROUTINE TO CALC COEFFICIENTS W1(M) (REAL*8) AS
C ROW*A(M)*COL.
C RAC1 IS MORE ELABORATE VERSION OF RAC() THAT IS DESIGNED TO USE
C A SUBSECTION ONLY OF A(M), VIZ:
C   ROWS I1 TO I2 OF A(M) =NR ROWS, NR=I2-I1+1
C   COLS J1 TO J2 OF A(M) =NC COLUMNS, NC=J2-J1+1
C   MAX=SIZE OF A(M)=NUMBER OF W1() VALUES PRODUCED
C ASSUMES THAT ROW() IS 1XNR (I.E. USES 1ST NR ELEMENTS OF ROW())
C ASSUMES THAT COL() IS 1XNC (I.E. USES 1ST NC ELEMENTS OF COL())
C 	THE A(I,J,M) ARE SUPPLIED AS EM,EN.
C ALSO OUTPUTS REAL*4 VALUES IN W(M)
C   DECLARED DIMENSIONS IN CALLING PROG:
C	KDR FOR ROW(1,NR)
C	KDC FOR COL(NC,1)
C	KDA FOR EM(NA,NA),EN(NA,NA)
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION ROW(1,KDR),COL(KDC,1),AMAT(kda1,kda2,kda3),W1(kdw)
	REAL*4 W(kdw)
C
	nR=i2-i1+1
	nC=j2-j1+1
C
	do m=1,max
	   w1(m)=0.0d0
	   do L=1,nR
		do n=1,nC
		   L1=L+i1-1	!=i1 to i2
		   n1=n+j1-1	!=j1 to j2
		   W1(M)=W1(M)+row(1,L)*amat(L1,n1,m)*col(n,1)
		enddo
	   enddo
	   w(m)=SNGL(w1(m))
	enddo

c	do 62 m=1,max
c	W1(M)=0.0D0
c	do 63 L=1,nR
c	do 63 n=1,nC
c	L1=L+i1-1	!=i1 to i2
c	n1=n+j1-1	!=j1 to j2
c63	W1(M)=W1(M)+row(1,L)*amat(L1,n1,m)*col(n,1)
c	W(m)=SNGL(W1(m))
c62	CONTINUE
c
	RETURN
	END
