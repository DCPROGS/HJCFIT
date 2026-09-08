	subroutine GfuncF(s,ng)
c Subroutine to return Frank Ball's g(s)=ng=number of eigenvalues of H(s) that
c  are equal to or less than s, for SHUT times (GFUNCA is for open times).
c  H(s) is calculated as in DETWF.
c Input = s
c Output = ng
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c Keep commons /detw1 and /detw3 so make QAF etc fixed size (use PARAMETER so
c easy to change max kA)
	PARAMETER(kAx=60) 	!max number of open states for following
	real*8 QAF(kAx,100),QFA(100,kAx),QAA(kAx,kAx),QFF(100,100)
	real*8 EXPQA(kAx,kAx),EXPQF(100,100)
	common/KBLK/kA,kB,kC,kD
	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF
	common/detw2/tres,km,nerr			!for HJCLIK,DETWA,DETWF
	real*8 dexp1
c
	real*8 Q1(100,100),det
	real*8 H(100,100),HFmat(100,100,100),eigHFF(100)
	logical debug,caplock,discprt
	common/dp/discprt
c
	debug=caplock()
c
	kF=kB+kC+kD
	if(km.ne.100) then
	   print *,' km is not equal to 100 in GfuncF'
	   STOP
	endif
c
c exp[-(sI-QAA)tres]= exp(-s*tres)*exp(QAA*tres)=exp(-s*tres)*EXPQA
c Calc I-exp[-(sI-QAA)tres] in this loop, in H. Also get (sI-QAA) in Q1 in
c same loop
	one=1.0d0
	nerr=0
	s1=-s*tres
	if(dabs(s1).gt.709.d0) then
	   nerr=4
	   RETURN
	endif
	est=dexp1(s1)
	do i=1,kA
	   do j=1,kA
		if(i.ne.j) then
c		   W(i,j)=-est*EXPQA(i,j)
		   H(i,j)=-est*EXPQA(i,j)
		   Q1(i,j)=-QAA(i,j)
		else
c		   W(i,j)=one-est*EXPQA(i,j)
		   H(i,j)=one-est*EXPQA(i,j)
		   Q1(i,j)=s-QAA(i,j)
		endif
	   enddo
	enddo
c
c Calc inv(sI-QAA) in Q1
      call MATINV2(Q1,kA,km,Q1,km,.false.,det,ndscale)
c Postmult by I-exp[-(sI-QAA)tres]
c      call MATMUL(Q1,W,W,kA,kA,kA,one,
      call MATMUL(Q1,H,H,kA,kA,kA,one,
     & km,km,km,km,km,km)
c Postmult by QAF
c      call MATMUL(W,QAF,W,kA,kA,kF,one,
      call MATMUL(H,QAF,H,kA,kA,kF,one,
     & km,km,kAx,km,km,km)
c Premult by QFA
c	call MATMUL(QFA,W,W,kF,kA,kF,one,
	call MATMUL(QFA,H,H,kF,kA,kF,one,
     & km,kAx,km,km,km,km)
c Calculate H(s) explicitly for Ball method
	do i=1,kF
	   do j=1,kF
	      H(i,j)=QFF(i,j)+H(i,j)
	   enddo
	enddo
c Calculate eigenvalues of H(s)
	call QMAT5(H,HFmat,kF,eigHFF,ibad,km,km,km)
c   Find g(s) = number of eigenvalues of H(s) that are equal to or less than s
	ng=0
	do i=1,kF
	   if(eigHFF(i).le.s) ng=ng+1
	enddo
	if(debug) then
	   print 3,(eigHFF(i),i=1,kF),s,ng
	   if(discprt) write(8,3) (eigHFF(i),i=1,kF),s,ng
3	   format(' s = ',g13.6,' eig(H) = ',/,5(4(g13.6,3x)),/,
     &	' s = ',g13.6,' number of eigenvalues of H(s).le.s = ',i3)
	endif
	RETURN
	end

