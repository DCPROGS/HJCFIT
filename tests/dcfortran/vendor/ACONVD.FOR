	subroutine ACONVd(amp,eigen,ncomp,tres,discprt,km)
c To renormalise areas that are defined from t=tres to infinity, to range
c t=0 to infinity, and scale to sum to 1.
c Expects EIGEN to contain negative values
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 eigen(km),amp(km)
	allocatable::area
	real*8 area(:)
	real*8 dexp1
	logical discprt
c
	allocate(area(ncomp))
	atot=0.0d0
	do m=1,ncomp
	   rate=-eigen(m)		!positive rate constant
	   if(dabs(rate).ge.1.d-30) then
		area(m)=amp(m)/rate
		area(m)=area(m)*dexp1(rate*tres)	!scale area
		atot=atot+area(m)
	   else
c		area(m)=-1.d6		!error signal
		area(m)=0.d0
	   endif
	enddo
c
      print 22
      if(discprt) write(8,22)
22	format(
     &' Areas for asymptotic pdf renormalised for t=0 to infinity',
     &' (and sum=1)',/,
     &' (so areas can be compared with ideal pdf)')
	do 23 m=1,ncomp
	  area(m)=area(m)/atot
        print 24,m,area(m)
        if(discprt) write(8,24) m,area(m)
24	  format(i3,4x,g13.6)
23	continue
      print 2
      if(discprt) write(8,2)
2	format(/)
	deallocate(area)
	RETURN
	end


