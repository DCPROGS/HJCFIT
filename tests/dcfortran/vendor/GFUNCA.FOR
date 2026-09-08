	subroutine GFUNCA(s,ng)
c Subroutine to return Frank Ball's g(s)=number of eigenvalues of H(s) that
c  are equal to or less than s, for open times (GFUNCF is for shut times).
c  H(s) is calculated as in DETWA.
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c Keep commons /detw1 and /detw3 so make QAF etc fixed size (use PARAMETER so
c easy to change max kA)
	PARAMETER(kAx=60) 	!max number of open states fir following
	real*8 QAF(kAx,100),QFA(100,kAx),QAA(kAx,kAx),QFF(100,100)
	real*8 EXPQA(kAx,kAx),EXPQF(100,100)
	real*8 dexp1
	common/KBLK/kA,kB,kC,kD
	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF
	common/detw2/tres,km,nerr				   !for DETWA,DETWF
c
	real*8 Q1(100,100),det
	real*8 H(100,100),HAmat(100,100,100),eigHAA(100)
	logical debug,caplock,discprt
	common/dp/discprt
c
	debug=caplock()
c
	kF=kB+kC+kD
	if(km.ne.100) then
	   print *,' km is not equal to 100 in GfuncA'
	   STOP
	endif
c
c exp[-(sI-QFF)tres]= exp(-s*tres)*exp(QFF*tres)=exp(-s*tres)*EXPQF
c Calc I-exp[-(sI-QFF)tres] in this loop, in H. Also get (sI-QFF) in Q1 in
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
	 	   H(i,j)=-est*EXPQF(i,j)
	 	   Q1(i,j)=-QFF(i,j)
		else
		   H(i,j)=one-est*EXPQF(i,j)
		   Q1(i,j)=s-QFF(i,j)
		endif
	   enddo
	enddo
c
c Calc inv(sI-QFF) in Q1
      call MATINV2(Q1,kF,km,Q1,km,.false.,det,ndscale)
c Postmult by I-exp[-(sI-QFF)tres]
      call MATMUL(Q1,H,H,kF,kF,kF,one,
     & km,km,km,km,km,km)
c Postmult by QFA
      call MATMUL(H,QFA,H,kF,kF,kA,one,
     & km,km,km,kAx,km,km)
c Premult by QAF
	call MATMUL(QAF,H,H,kA,kF,kA,one,
     & kAx,km,km,km,km,km)
c Add QAA to complete H(s)
c Calculate H(s) explicitly for Ball method
	do i=1,kA
	   do j=1,kA
	      H(i,j)=QAA(i,j)+H(i,j)
	   enddo
	enddo
c Calculate eigenvalues of H(s)
	call QMAT5(H,HAmat,kA,eigHAA,ibad,km,km,km)
c   Find g(s) = number of eigenvalues of H(s) that are equal to or less than s
	ng=0
	do i=1,kA
	   if(eigHAA(i).le.s) ng=ng+1
	enddo
	if(debug) then
	   print 3,(eigHAA(i),i=1,kA),s,ng
	   if(discprt) write(8,3) (eigHAA(i),i=1,kA),s,ng
3	   format(' s = ',g13.6,' eig(HAA) = ',/,5(4(g13.6,3x)),/,
     &	' s = ',g13.6,' number of eigenvalues of H(s).le.s = ',i3)
	endif
c
	RETURN
	end

