	subroutine EQOC_HJC(Q,peq,k,km1,km)
c Version of eqocc2 for HJCFIT
c To solve peq*Q=0 by Hawkes method
c
c Modif 02/05/02 12:26pm
c NB in HJCASYMP this routine is use to solve rW(s)=0 , ru=1 (HJC 1992, p394), and
c in thus case I think  negative entries may be allowable, so report error
c only if sum is near zero, not if it is negative (must not be negative
c for occupancies of course, so separate routine would be an idea)
c
c See TEQ.FOR for tests.
	real*8 Q(km1,km1),peq(km)
	allocatable::X,XT,Z,uT
	real*8 X(:,:),XT(:,:),Z(:,:),uT(:,:)
	real*8 det
	COMMON/determ/det			!for matinv
	common/perr/nerr2		!for errors in eqoc_hjc
c
	real*8 one,s
c
	k1=k+1
	one=1.0d0
	nerr2=0	!error flag
	ALLOCATE(X(k1,k1),XT(k1,k1),Z(k1,k1),uT(1,k1))
c
	do j=1,k
	   uT(1,j)=one		!define
	enddo
c
	do i=1,k
	   do j=1,k
		X(i,j)=Q(i,j)
	   enddo
	enddo
c
	do i=1,k
        X(i,k1)=one 	!add column of ones
	enddo
c
	do i=1,k
	   do j=1,k1
		XT(j,i)=X(i,j)
	   enddo
	enddo
c
      call MATMUL(X,XT,Z,k,k1,k,one,
     &	k1,k1,k1,k1,k1,k1)
c
c      call MATINV(Z,k,k1,Z,k1)		!inv(Z) in Z
      call MATINV2(Z,k,k1,Z,k1,.false.,det,ndscale)		!inv(Z) in Z
c
      call MATMUL(uT,Z,peq,1,k,k,one,
     &	1,k1,k1,k1,1,km)
c Check sum to 1
	s=0.0d0
	do i=1,k
	   s=s+peq(i)
	enddo
c NB in HJCASYMP this routine is use to solve rW(s)=0 , ru=1 (HJC 1992, p394), and
c in thus case I think  negative entries may be allowable, so report error
c only if sum is near zero, not if it is negative (must not be negative
c for occupancies of course, so separate routine would be an idea)
c	if(s.lt.1.d-100) then
	if(dabs(s).lt.1.d-100) then
	   nerr2=1
	else
	   do i=1,k
		peq(i)=peq(i)/s
	   enddo
	endif
c
	DEALLOCATE(X,XT,Z,uT)
	RETURN
	end

