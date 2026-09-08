	subroutine EQOC_red(Q,peq,k,km1,km)
c Version of eqocc2 for HJCFIT
c To solve peq*Q=0 by reduced matrix method, as alternative to 'adding
c column of ones.  Latter method, as usually used, can be very ill-conditioned
c when SS' has some entries = 1 but other entries ~10^10 Hawkes method.
c
c May be useful in HJCASYMP where eqocc routine is use to solve
c rW(s)=0 , ru=1 (HJC 1992, p394), and
c In thus case I think  negative entries are allowable, so report error
c only if sum is near zero, not if it is negative (must not be negative
c for occupancies of course, so separate routine would be an idea)
c
c See TEQ.FOR for tests.
	real*8 Q(km1,km1),peq(km)
	allocatable::R,invR,row
	real*8 R(:,:),invR(:,:),row(:,:)
	allocatable::unit
	real*8 unit(:,:)
	real*8 det
	COMMON/determ/det			!for matinv
	common/perr/nerr2		!for errors in eqoc_hjc
c
	real*8 one,s
c
c reduced matrix, R, is (k-1) x (k-1)
	k1=k-1
	one=1.0d0
	nerr2=0	!error flag
	ALLOCATE(R(k1,k1),invR(k1,k1),row(1,k))
	ALLOCATE(unit(k1,k1))
c
	do j=1,k-1
	   row(1,j)=Q(k,j)		!define last row
	enddo
c
	do i=1,k-1
	   do j=1,k-1
		R(i,j)=Q(i,j)-Q(k,j)	!subract last row, omit last col
	   enddo
	enddo
c
      call MATINV2(R,k1,k1,invR,k1,.false.,det,ndscale)	!inv(R) in invR
c
c	Check inverse
      call MATMUL(R,invR,unit,k1,k1,k1,one,
     &	k1,k1,k1,k1,k1,k1)
	do i=1,k1
	   if(dabs(unit(i,i)-one).gt.1.d-6) nerr2=2	!diagonal bad
	   if(dabs(unit(i,i)-one).gt.1.d-2) nerr2=4	!diagonal very bad
	   do j=1,k1
		if(i.ne.j) then
		   if(dabs(unit(i,j)).gt.1.d-6) nerr2=3	!off-diagonal bad
		   if(dabs(unit(i,j)).gt.1.d-2) nerr2=4	!off-diagonal very bad
		endif
	   enddo
	enddo
c end of check
c An even better check.
c p'=-row*inv(R)
      call MATMUL(row,invR,peq,1,k1,k1,-one,
     &	1,k,k1,k1,1,km)
c Define peq(k) so sum to 1
	s=0.0d0
	do i=1,k-1
	   s=s+peq(i)
	enddo
	peq(k)=one - s
c NB in HJCASYMP this routine is use to solve rW(s)=0 , ru=1 (HJC 1992, p394), and
c in thus case negative entries are allowable, so report error
c only if sum is near zero, not if it is negative (must not be negative
c for occupancies of course)
	if(dabs(s).lt.1.d-100) then
	   nerr2=1
	endif
c An even better check, peq*Q=0? (now done in checkrW.for)
c	call MATMUL(peq,Q,row,1,k,k,one,
c     &	1,km,km1,km1,1,k)
cc result should be all zeros in row(1,j)
c	do j=1,k
c	   if(dabs(row(1,j)).gt.1.d-6) nerr2=5
c	enddo
c	print 407, (row(1,j),j=1,k)
c407	format(' rW = ',7(g13.6,2x))
c
	DEALLOCATE(R,row,invR)
	DEALLOCATE(unit)
	RETURN
	end

