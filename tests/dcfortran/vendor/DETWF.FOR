	real*8 function DETWF(s)
c To calculate det(W(S)) for HJC distributions.  Also returns, via COMMON
c the matrix WF(s) via common
c DET(W(s)) for SHUT times
c For open times W(s) is kA*kA matrix
c For shut times W(s) is kF*kF matrix
c
c NB it should never matter if W(s) is scaled, because it find zeros
c we need only the sign of det[W(s)] for bisection, and to know
c the point at which the determinant is zero -scaling of either W(s) or
c its determinant by any constant factor should not affect either
c of these.
c
c Modified 05/07/04 03:12pm to use MATMUL2 which scales to prevent
c overflows during multiplication -returns a product that is scaled
c by a factor of 10^(-20*ndscale1), but if ndscale1 is not zero, need
c to ensure that final value of W(s) is scaled accordingly.
c
c W(s)=sI -H(s)
c where, for open times,
c H(s)=QAA + QAF*(inv(sI-QFF)*(I - exp(-(sI-QFF)*tres))*QFA
c and for shut times
c H(s)=QFF + QFA*(inv(sI-QAA)*(I - exp(-(sI-QAA)*tres))*QAF
c term  1     2       3              4                   5
c This has form 1 + 2*3*4*5
c so accumulate total scal factor for 2*3*4*5 and scale 1=QFF by same
c factor before adding it. That gives H(s)' = H(s)*10^(-nds) where nds=total
c scale factor, and then get W(s)'=10^(-nds)sI - H(s)'
c
c Modified 07/30/03 09:23am to look at properties of H(s) with view to
c using Ball's method to find initial guesses for bisection of  roots (see
c also HJCroots.mcd) -now removed again
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
	PARAMETER(kAx=60) 	!max number of open states for following
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
c for Ball method
c	real*8 H(100,100),HFmat(100,100,100),eigHFF(100)
c	real*4 s4a,s4b
c	character*1 ans
c	logical debug,caplock
c	logical discprt,first
c	common/dp/discprt
c
	common/KBLK/kA,kB,kC,kD
	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF
	common/detw2/tres,km,nerr			!for HJCLIK,DETWA,DETWF
	common/detw3/WA,WF					   !for DETWA,DETWF
	common/det4/ndscale,idsign
	real*8 dexp1
c
	kF=kB+kC+kD
	if(km.ne.100) then
	   print *,' km is not equal to 100 in DETWF'
	   STOP
	endif
c
c exp[-(sI-QAA)tres]= exp(-s*tres)*exp(QAA*tres)=exp(-s*tres)*EXPQA
c Calc I-exp[-(sI-QAA)tres] in this loop, in W. Also get (sI-QAA) in Q1 in
c same loop
	one=1.0d0
	nerr=0
	nds=0		!total scale factor for W(s) = 10^(-nds)
	s1=-s*tres
	if(dabs(s1).gt.709.d0) then
	   nerr=4
	   RETURN
	endif
	est=dexp1(s1)
	do i=1,kA
	   do j=1,kA
		if(i.ne.j) then
		   W(i,j)=-est*EXPQA(i,j)
		   Q1(i,j)=-QAA(i,j)
		else
		W(i,j)=one-est*EXPQA(i,j)
		Q1(i,j)=s-QAA(i,j)
		endif
	   enddo
	enddo
c
c Calc inv(sI-QAA) in Q1
      call MATINV2(Q1,kA,km,Q1,km,.false.,det,ndscale)
c scale factor for inverse is 10^(-10*ndscale)
	nds=nds+10*ndscale
c
c Try to avoid floating point overflows by using MATMUL2 here
c Should return ndscale1=0 -the product is scaled down by factor
c of 10^(-20*ndscale1) but input arrays not changed (so altered
c to make output array always different name from input)
c
c Postmult by I-exp[-(sI-QAA)tres]
c      call MATMUL2(Q1,W,W,kA,kA,kA,one,
      call MATMUL2(Q1,W,WF,kA,kA,kA,one,
     & km,km,km,km,km,km,ndscale1)
	nds=nds + 20*ndscale1	!scale fac for output in WF
c Postmult by QAF
c      call MATMUL2(W,QAF,W,kA,kA,kF,one,
      call MATMUL2(WF,QAF,W,kA,kA,kF,one,
     & km,km,kAx,km,km,km,ndscale2)
	nds=nds + 20*ndscale2	!scale fac for output in W
c Premult by QFA
	call MATMUL2(QFA,W,WF,kF,kA,kF,one,
     & km,kAx,km,km,km,km,ndscale3)
	nds=nds + 20*ndscale3	!scale fac for output in WF
c Copy scaled WF to W
	do i=1,kF
	   do j=1,kF
		W(i,j)=WF(i,j)
	   enddo
	enddo

c	if(ndscale1.ne.0) then
c	   print 7
c7	   format(
c     &' ERROR: numerical overflow in matrix multiplication in DETWF')
c	endif
c and add QFF to complete H(s), then calc W(s)=sI-H(s)
c Must scale QFF by same factor as W before adding
	if(nds.eq.0) then
	   sfac=one
	else if(nds.le.308) then
	   sfac=10.d0**(-dfloat(nds))
	else
	   nerr=5
	endif
c
	do i=1,kF
	   do j=1,kF
		if(i.ne.j) then
		   W(i,j)=-sfac*QFF(i,j) - W(i,j)
		else
		   W(i,j)=sfac*(s-QFF(i,j)) - W(i,j)
		endif
	   enddo
	enddo
c Calculate H(s) explicitly for Ball method
c	do i=1,kF
c	   do j=1,kF
c	      H(i,j)=QFF(i,j)+W(i,j)
c	   enddo
c	enddo
cc  then calc W(s)=sI-H(s)
c	do i=1,kF
c	   do j=1,kF
c		if(i.ne.j) then
c		   W(i,j)=-H(i,j)
c		else
c		   W(i,j)=s-H(i,j)
c		endif
c	   enddo
c	enddo
c Copy to WF for return via COMMON
	do i=1,kF
	   do j=1,kF
		WF(i,j)=W(i,j)
	   enddo
	enddo
c Now calc the determinant of W(s)
c
c DETV altered to prevent overflow: value of DET
c which is output must be mult by 10**ndexp to get correct determinant
c
c 01/21/02 05:11pm Still get occasional overflow in DETV
c e.g. with 5x5 matrix which is near diagonal and elements all
c around 1.d174. Can do rough check for this by looking at
c product of diagonals first
c 07/10/03 12:55pm can remove this bit because DETERM2 should never
c overflow!
c	sumdiag=0.0d0
c	do i=1,kF
c	   sumdiag=sumdiag+dlog10(dabs(W(i,i)))
c	enddo
c	if(sumdiag.gt.234.d0) then
c	   nerr=10		!overflow of det() likely
c	   RETURN
c	endif
c
c	call ARRAYD(2,kF,kF,km,km,W,W)	!convert to 1-D storage
cc	call DETV(W,kF,det,lwork,mwork,ndexp)	!DETV is in EQOCCD
c	call DETERM1(W,kF,det,lwork,mwork,ndexp,nerr)
c	print 4,det,ndexp
c4	format(' det, ndexp = ',g13.6,i10)
c	DETWF=det*(10.d0**ndexp)
c Try determinant routine from numerical recipes
c	subroutine DETERM2(A,n,ndim,det,ndscale)
c ndscale=number of factors of 10^-10 that det was multiplied by prevent
c overflow/underflow
c For most calcs need only the sign of det, or sign of slope
c dsign(a,b) transfers sign of b to sign of a
c dsign(1.d0,det) gives +1 or -1
c	subroutine DETERM2(A,n,ndim,det,ndscale)
	ndim=100	!declared size of W
	call DETERM2(W,kF,ndim,det,ndscale)
	detsign=dsign(1.d0,det)		!+1 or -1 (real*8)
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
	DETWF=det
c
	RETURN
	end
