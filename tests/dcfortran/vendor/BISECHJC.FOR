	subroutine BISEChjc(FUNC,X1,X2,Y,XOUT,YOUT,epsx,epsy,nsmax,
     & NERR,NDISP,OMIT)
c Double precision version of BISEC0
c
c Special version for HJCFIT 02/13/06 09:35am.
c  There can be a bad numerical problem in finding asymptotic roots.
c For this want to find root, s, such that det[W(s)]=0.  The problem
c arises when det[W(s)] is so sensitive to s that even when x1 and x2
c are identical in all sig figs, the determinant can be far from zero.
c When epsx is used to control convergence, the routine may not get
c as close as it can, because, although s will be very precise, det[W(s)]
c may not be as close to zero as is achievable.
c Two things can be done to help.
c (1) even when using epsx for convergence (so have fixed number of steps)
c   check that yout is no longer changing at end, and if it is, then increase
c   the number of steps (up to maximum of 100)
c (2) Keep record of yerr=yout-y.  If an earlier iteration happens to give
c	a smaller yerr than the final iteration, then use the best one
c	for output.
c (3) nsmax=max number of steps added to arguments
c
c Basic bisection subroutine in which all values that are needed
c by FUNC are in common, and FUNC is a function subroutine
c Modif 02/07/03 10:23am so x1, x2 not altered by call
c Modif. so  nerr=0 OK
c		 nerr=1 if f(x1) and f(x2) both less than Y
c	       nerr=2 if f(x1) and f(x2) both greater than Y
c and nerr positive if f(x) increases with x
c     nerr negative if f(x) decreases with x
c	ABORT with F2: if aborted returns with nerr=3
c  NDISP =>1 prints every NDISPth iteration
c  NDISP=0 no print 'except F2 to abort'
c  NDISP<0 no print at all
c
C*** VERSION OF BISECT WHICH,IF OMIT SET TO TRUE, OMITS THE FIRST
C TWO FUNCTION EVALUATIONS AT X1,X2 SO MUST ENSURE BEFORE CALL
C THAT FUNC(X1) < Y AND THAT X1,X2 REALLY INCLUDE ROOT
c EPSx=accuracy for X (ignored if neg- i.e. Nstep not calc so loops
c	until EPSy is satisfied; max Nstep set to 1000)
c EPSY=accuracy for func(X)-Y (ignored if neg)
c
C SOLVES FUNC(X)=Y BY BISECTION GIVEN Y AND GUESSES X1,AND
C X2 THAT INCLUDE SOLUTION.
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 FUNC
	logical OMIT
	logical use_expx
	integer*2 ktype
	logical KBHIT		!Utility lib
	character ch,getch
	external FUNC
c
	NDSAV=0
	x1sav=x1
	x2sav=x2
	use_expx=epsy.lt.0.0d0
	if(ndisp.ge.0) then
	   print 121
121	   format(' Bisection may be aborted with F2')
	   print 12
12	   FORMAT( ' Step #     X(out)      Y(out)',/)
	endif
c
	nerr=0
	IF(OMIT) GOTO 1
	fx1=func(x1)
	fx2=func(x2)
	if(ndisp.gt.0) then
	   print 25,x1,fx1
	   print 25,x2,fx2
25	   format(' f(',g13.6,')= ',g13.6)
	endif
c
c check for 'same sign' error and set nerr
c Modif. so nerr=1 if f(x1) and f(x2) both less than Y
c	    nerr=2 if f(x1) and f(x2) both greater than Y
c and nerr positive if f(x) increases with x
c     nerr negative if f(x) decreases with x
	if((fx1-y)*(fx2-y).gt.0.0d0) then
98	   if(fx1.lt.y.and.fx2.lt.y) nerr=1
  	   if(fx1.gt.y.and.fx2.gt.y) nerr=2
	   if(x2-x1.gt.1.d-30) then
		if((fx2-fx1)/(x2-x1).lt.0.0d0) nerr=-nerr
	   endif
	   if(ndisp.gt.0) then
		print 24,x1,fx1,x2,fx2
24		format(
     &	  ' BISECD: Guesses don''t include root:',
     &	  2(' f(',g13.6,') = ',g13.6))
	   endif
	   goto 9
	endif
c
c swap x1, x2 if necessary
	if(fx1.gt.y) then
	   temp=x1
	   x1=x2
	   x2=temp
	endif
1	continue
c
	if(epsx.gt.0.) then
	   nstep=int4(dlog(dabs(X1-X2)/epsx)/dlog(2.0d0)+0.5d0)
	   if(nstep.gt.nsmax) nstep=nsmax
	endif
C
c
	ns=0
cc	DO 2 NS=1,NSTEP
22	continue		!return here for next iteration until ns=nstep (if epsx used)
	if(KBHIT()) then
	   ch=GETCH(ktype)
	   if(ktype.eq.0.and.ichar(ch).eq.60) then
		nerr=3	!abort
		goto 9
	   endif
	endif
	ns=ns+1
	xout=0.5d0*(x1+x2)
	Yout=func(Xout)
	if(ndisp.gt.0.and.((ns-ndsav).ge.ndisp)) then
	   print 4,NS,XOUT,YOUT
4	   FORMAT(I5,3X,G13.6,2X,G13.6)
	   NDSAV=NS
	endif
	yerr=yout-y
	if(yout.lt.y) then
	   x1=xout
	else if(yout.gt.y) then
	   x2=xout
	endif
	if(ns.eq.1) then
	   yerrmin=dabs(yerr)	!define initial value
	else if(dabs(yerr).lt.yerrmin) then
	   yerrmin=dabs(yerr)
	   ybest=yout
	   xbest=xout
	endif
c
	if(yout.eq.y) goto 9		!done
c
	if(use_expx) then 	!use epsy to control convergence
c but add bit to increase number of steps if y-error still falling
c when last step reached

	   if(dabs(yerr-yerrlast).gt.1.d-10.and.ns.eq.nstep) then
		nstep=nsmax
	   endif
	   yerrlast=yerr
	   if(ns.lt.nstep) goto 22		! otherwise use nstep
	else
	   if(dabs(yerr).gt.epsy) goto 22
c	   goto 9	!done
	endif
c Normal end.
c
c check that an earlier iteration not better
c
9	continue
	yout=ybest
	xout=xbest
	x1=x1sav	!restore
	x2=x2sav
	RETURN
c
	END



