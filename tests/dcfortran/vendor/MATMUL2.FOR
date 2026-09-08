	subroutine MATMUL2(A,B,R,na,n,nb,scalar,
     & ka1,ka2,kb1,kb2,kr1,kr2,ndscale)
c
c Modif of MATMUL, 05/07/04 09:54am, designed to be more resistant to
c floating point overflows
c Maximum for real*8 is 10^308 = exp(709)
c NB makes no sense to scale each element of result differently! Would need
c to return whole array of scale factors, ndscale(i,j).

c Returns with ndscale=0 if OK
c If ndscale>0 then whole if result in R is scaled down by a factor
c of 10^(-20*ndscale) (both A1 and B1 scaled down ndscale times by factors of
c of 1.d-10
c
c Old notes
C***MODIF 7-OCT-81 SO THAT (1)OUTPUT ARRAY CAN BE SAME AS INPUT
C AND (2) RESULT IS MULTIPLIED BY SCALAR
C DOUBLE PRECISION SUBROUTINE FOR MATRIX MULTIPLICATION WITH
C MATRICES IN NORMAL 2-DIMENSIONAL FORM. CAN EQUALLY USE FOR
C VECTOR*MATRIX OR VECTOR*VECTOR IF VECTORS DECLARED AS B(1,J) ETC
C	CALCS R=A*B WHERE A IS (NA X N),B IS (N X NB) SO R IS (NA X NB).
C KA1,KB1,KR1 ETC=DIMENSIONS OF A,B,R IN CALLING PROG
c Examples  of calls to MATMUL etc
c	call MATMUL(GDB,GBD,HDD,kD,kB,kD,one,
c     & km,km,km,km,km,km)
c	call MATMUL(phi,GED,row2,1,kE,kD,one,1,km,km,km,1,km)
c	call MATMUL(QBD,end1,col1,kB,kD,1,one,km,km,km,1,km,1)
c	call VECMUL(row2,unit,den,kD,km,km)	!postmult by uD
C
	real*8 A(ka1,ka2),B(kb1,kb2),R(kr1,kr2)
	allocatable::temp,A1,B1
	real*8 temp(:,:),A1(:,:),B1(:,:)
	real*8 ZERO,SCALAR
	real*8 ai,bj,x,d	!for overflow check
C
c corrected 05/11/05 03:35pm 
c	ALLOCATE(temp(na,nb),A1(na,nb),B1(na,nb))
	zero=0.0d0
	nd=max(na,nb,n)
	ALLOCATE(temp(nd,nd),A1(nd,nd),B1(nd,nd))
	do i=1,nd
	   do j=1,nd
		temp(i,j)=zero
		A1(i,j)=zero
		B1(i,j)=zero
	   enddo
	enddo
c
	d=1.0d0	!scale fac
	ndscale=0
1	continue
c Define A1, B1 as internal (possibly svcaled down) versions of A,B so
c the input arrays are not changed
	do i=1,nA
	   do j=1,nB
		A1(i,j)=d*A(i,j)
		B1(i,j)=d*B(i,j)
	   enddo
	enddo
c
	do i=1,nA
	   do j=1,nB
		temp(i,j)=zero
		do m=1,n
c overflow check (NB could get bj=0. but in this case overflow not possible)
		   ai=dabs(A1(i,m))
		   bj=dabs(B1(m,j))
		   if(ai.gt.bj) then	!make bj the bigger
			x=ai
			ai=bj
			bj=x
		   endif
		   if(bj.gt.1.d-50) then
			if(ai.gt.1.d308/bj) then	!mult below would overflow
			   ndscale=ndscale+1
			   d=d*1.d-10
			   goto 1	!until OK
			endif
		   endif
c now the multiplication
		   temp(i,j)=temp(i,j) + A1(i,m)*B1(m,j)
		enddo
	   enddo
	enddo
c
	do i=1,nA
	   do j=1,nB
		R(i,j)=TEMP(i,j)*scalar
	   enddo
	enddo
C
	DEALLOCATE(temp,A1,B1)
	RETURN
	END


