	subroutine MATINV2(array,k,kd1,recip,kd2,calcdet,det,ndscale)
c Alternative matrix inversion routine based on LU decomposition
c (Numerical Recipes 2nd ed, p 40).  Uses same LUDCMP routine as used
c by DETERM2.  Call same as for original matinv routine, with addition
c of calcdet=true if determinant to be returned too,
c  determinant=det*10(**ndscale)
c ndscale=number of factors of 10^10 by which determinant has to scaled
c to prevent overflow/underflow
c Needs real*8 version of LUDCMP, as at present in determ2.for file -now separate
c
c Maximum for real*8 is 10^308 = exp(709)
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 ARRAY(kd1,kd1),RECIP(kd2,kd2),det
	logical calcdet
	allocatable::a1,index
	real*8 a1(:,:)
	integer*4 index(:)
c=	COMMON/determ/det
c
	if(k.eq.1) then
	   recip(1,1)=1.d0/array(1,1)
	   det=array(1,1)
	   ndscale=0
	   RETURN
	endif
c
	ALLOCATE(a1(k,k),index(k))
c
c Copy input to a1 so not destroyed, and make unit matrix in recip()
	do i=1,k
	   do j=1,k
		a1(i,j)=array(i,j)
		recip(i,j)=0.0d0
	   enddo
	   recip(i,i)=1.d0
	enddo
c
	call LUDCMPd(a1,k,k,index,d,nerr)	!returns d=-1 or +1
	if(nerr.eq.1) then
	   det=0.0d0
	   RETURN
	endif
c
c Calc determinant if requested
	ndscale=0
	if(calcdet) then
c=       do i=1,n
       do i=1,k
c overflow check
1	   di=dlog(dabs(D))
	   ai=dlog(dabs(A1(i,i)))
c	   if(di+ai.gt.308.d0) then	!mult below would overflow
	   if(di+ai.gt.700.d0) then	!mult below would overflow
		ndscale=ndscale+1
		d=d*1.d-10
		goto 1	!until OK
	   endif
c Underflow check
2	   di=dlog(dabs(D))
	   ai=dlog(dabs(A1(i,i)))
c	   if(di+ai.lt.-308.d0) then	!mult below would overflow
	   if(di+ai.lt.-700.d0) then	!mult below would overflow
		ndscale=ndscale-1
		d=d*1.d10
		goto 2	!until OK
	   endif
c The actual loop
	   D=D*A1(i,i)
	 enddo
	 DET=d
	endif
c
c  Calculate inverse
	do j=1,k
	   call LUBKSBd(a1,k,k,index,recip(1,j))	!r(1,j) is address of jth col
	enddo
c calc determimant too?
c
	deallocate(A1,index)
c
	RETURN
	end


