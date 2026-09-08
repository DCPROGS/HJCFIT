	subroutine BISECd(FUNC,X1,X2,Y,XOUT,YOUT,EPSx,epsy,
     & NERR,NDISP,OMIT)
c Double precision version of BISEC0
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
c	real*8 a1,a2
c	real*4  a3
	logical OMIT
	integer*2 ktype
	logical KBHIT		!Utility lib
	character ch,getch
	external FUNC
c
	NDSAV=0
	x1sav=x1
	x2sav=x2
	if(ndisp.ge.0) print 121
121	format(' Bisection may be aborted with F2')
	if(ndisp.gt.0) print 12
12	FORMAT( ' Step #     X(out)      Y(out)',/)
	NERR=0
	IF(OMIT) GOTO 1
	fx1=func(x1)
	if(ndisp.gt.0) print 25,x1,fx1
	fx2=func(x2)
	if(ndisp.gt.0) print 25,x2,fx2
25	format(' f(',g13.6,')= ',g13.6)
	IF((FX1-Y)*(FX2-Y).gt.0.0d0) GOTO 98	!SAME SIGN=ERROR
	IF(FX1.LT.Y) GOTO 1
	P=X1
	X1=X2
	X2=P
1	CONTINUE
c
	Nstep=1000	!max if epsy used to control convergence
	if(epsx.gt.0.) then
	   nstep=int4(dlog(dabs(X1-X2)/epsx)/dlog(2.0d0)+0.5d0)
	endif
C
	ns=0
c	DO 2 NS=1,NSTEP
22	continue		!return here for next iteration
	if(KBHIT()) then
		ch=GETCH(ktype)
		if(ktype.eq.0.and.ichar(ch).eq.60) goto 999	!abort
	endif
	ns=ns+1
	XOUT=0.5d0*(X1+X2)
	Yout=func(Xout)
	if(ndisp.lt.1.or.((ns-ndsav).lt.ndisp)) goto 10
	print 4,NS,XOUT,YOUT
4	FORMAT(I5,3X,G13.6,2X,G13.6)
	NDSAV=NS
10	continue
	erry=yout-y
	IF(YOUT-Y) 5,6,7
5	X1=XOUT
	GOTO 2
6	GOTO 99
7	X2=XOUT
c2	CONTINUE
2	if(epsy.gt.0.d0) goto 21 	!use epsy to control convergence
	if(ns.lt.nstep) goto 22		! otherwise use nstep
	goto 99
21	if(dabs(erry).gt.epsy) goto 22
	goto 99
C
c Set NERR
c Modif. so nerr=1 if f(x1) and f(x2) both less than Y
c	    nerr=2 if f(x1) and f(x2) both greater than Y
c and nerr positive if f(x) increases with x
c     nerr negative if f(x) decreases with x
98	if(fx1.lt.y.and.fx2.lt.y) nerr=1
  	if(fx1.gt.y.and.fx2.gt.y) nerr=2
	if(x2-x1.gt.1.d-30) then
	   if((fx2-fx1)/(x2-x1).lt.0.0d0) nerr=-nerr
	endif
	if(ndisp.gt.0) then
	   print 24,x1,fx1,x2,fx2
24	   format(
     &	' BISECD: Guesses don''t include root:',
     &	2(' f(',g13.6,') = ',g13.6))
	endif
c
99	x1=x1sav	!restore
	x2=x2sav
	RETURN
c after aborting
999	nerr=3
	x1=x1sav	!restore
	x2=x2sav
	RETURN
	END



