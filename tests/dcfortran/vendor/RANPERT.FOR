	subroutine RANPERT(theta,thmin,kfit,perfac)
c For HJCFIT (called in HJCLIK). To put into theta() a version of the 'previous best'
c parameters into which a random perturbation has been added, to prevent looping
c The perturbation is in the range +/- perfac, e.g. perfac=0.1 for
c a perturbation that gives a value in the range 90% to 110% of thmin(i)
c
	real*8 theta(200),thmin(200),perfac
	real*8 u2,drandom
c
	u2=drandom()			!uniform(0,1)
c
	do i=1,kfit
	   u2=drandom()			!uniform(0,1)
	   theta(i)=thmin(i)*(u2*2.d0*perfac + (1.d0-perfac))
	enddo
	RETURN
	end

