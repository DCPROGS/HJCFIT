	real*8 function DETWA(s)
c To calculate det(W(S)) for HJC distributions.  Also returns, via COMMON
c the matrix WA(s) via common
c DET(W(s)) for OPEN times
c For open times W(s) is kA*kA matrix
c For shut times W(s) is kF*kF matrix
c
c Modified 07/10/03 02:21pm to use determ2 and to treat overflow and underflow
c better -this is increasing problem in mechs with many states. Determ2 returns
c ndscale = number of factors of 1.d-10 that have been applied to keep det
c from overflowing (negative ndscale = factors of 1.d10 to prevent underflow)
c common/det4/ndscale,idsign
c added so ndscale, idsign available elsewhere
c idsign=1 if det>0, idsign=-1 if det<0, idsign=0 if det=0
c
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c Keep commons /detw1 and /detw3 so make QAF etc fixed size (use PARAMETER so
c easy to change max kA)
	PARAMETER(kAx=60) 	!max number of open states fir following
	real*8 QAF(kAx,100),QFA(100,kAx),QAA(kAx,kAx),QFF(100,100)
	real*8 EXPQA(kAx,kAx),EXPQF(100,100)
	real*8 WA(kAx,kAx),WF(100,100)
c=	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA   !for DETWA,DETWF,dARSdS -now param
c=	common/detw3/WA,WF			  !for DETWA,DETWF -now params
c=	real*8 QAF(10,10),QFA(10,10),QAA(10,10),QFF(10,10)
c=	real*8 EXPQF(10,10),EXPQA(10,10),WA(10,10),WF(10,10)
c Local
	real*8 W(100,100),Q1(100,100),det
c	integer lwork(100),mwork(100)	!not needed for determ2
	real*8 dexp1
	common/KBLK/kA,kB,kC,kD
	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF
	common/detw2/tres,km,nerr				   !for DETWA,DETWF
	common/detw3/WA,WF					  !for DETWA,DETWF
	common/det4/ndscale,idsign
c
	kF=kB+kC+kD
	if(km.ne.100) then
	   print *,' km is not equal to 100 in DETWA'
	   STOP
	endif
c exp[-(sI-QFF)tres]= exp(-s*tres)*exp(QFF*tres)=exp(-s*tres)*EXPQF
c Calc I-exp[-(sI-QFF)tres] in this loop, in W. Also get (sI-QFF) in Q1 in
c same loop
	one=1.0d0
	s1=-s*tres
	nerr=0
	if(dabs(s1).gt.709.d0) then
	   nerr=4
	   RETURN
	endif
	est=dexp1(s1)
	do i=1,kF
	   do j=1,kF
		if(i.ne.j) then
	 	   W(i,j)=-est*EXPQF(i,j)
	 	   Q1(i,j)=-QFF(i,j)
		else
		   W(i,j)=one-est*EXPQF(i,j)
		   Q1(i,j)=s-QFF(i,j)
		endif
	   enddo
	enddo
c
c Calc inv(sI-QFF) in Q1
      call MATINV2(Q1,kF,km,Q1,km,.false.,det,ndscale)
c Postmult by I-exp[-(sI-QFF)tres]
      call MATMUL(Q1,W,W,kF,kF,kF,one,
     & km,km,km,km,km,km)
c Postmult by QFA
      call MATMUL(W,QFA,W,kF,kF,kA,one,
     & km,km,km,kAx,km,km)
c Premult by QAF
	call MATMUL(QAF,W,W,kA,kF,kA,one,
     & kAx,km,km,km,km,km)
c and add QAA to complete H(s), then calc W(s)=sI-H(s)
	do i=1,kA
	   do j=1,kA
		if(i.ne.j) then
		   W(i,j)=-(QAA(i,j)+W(i,j))
		else
		   W(i,j)=s-(QAA(i,j)+W(i,j))
		endif
	   enddo
	enddo
c Copy to WA for return via COMMON
	do i=1,kA
	   do j=1,kA
		WA(i,j)=W(i,j)
	   enddo
	enddo
c Now calc the determinant of W(s)
c DETV altered to prevent overflow: value of DET
c which is output must be mult by 10**ndexp to get correct determinant
c	call ARRAYD(2,kA,kA,km,km,W,W)	!convert to 1-D storage
cc	call DETV(W,kA,det,lwork,mwork,ndexp)	!DETV is within EQOCCD
c	call DETERM1(W,kA,det,lwork,mwork,ndexp,nerr)
c	DETWA=det*(10.d0**ndexp)
c Try determinant routine from numerical recipes
c	subroutine DETERM2(A,n,ndim,det,ndscale)
c ndscale=number of factors of 10^-10 that det was multiplied by prevent
c overflow/underflow
c For most calcs need only the sign of det, or sign of slope
c dsign(a,b) transfers sign of b to sign of a
c dsign(1.d0,det) gives +1 or -1
c	subroutine DETERM2(A,n,ndim,det,ndscale)
	ndim=100	!declared size of W
	call DETERM2(W,kA,ndim,det,ndscale)
	detsign=dsign(1.d0,det)		!+1 or -1 (real*8) det=0.d0 gives 1 not 0
	idsign=int4(detsign)
	if(det.eq.0.d0) idsign=0
c Correct det for scaling (and leave with ndscale=0) if possible
c Otherwise return with scaled det and ndscale not zero
	if(ndscale.ne.0) then
c	   dn=dlog(10.d0**(10*ndscale))	!arg can overflow!
	   dn=dfloat(10*ndscale)*dlog(10.d0)	!this form OK
	   dd=dlog(dabs(det))
	   if(dn+dd.lt.308.d0) then	!det can be corrected before leaving
		det=det*10.d0**(10*ndscale)
		ndscale=0	!DET has correct value now
	   endif
	endif
c
	DETWA=det
c
	RETURN
	end




	
