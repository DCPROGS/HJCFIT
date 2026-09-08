	subroutine CeGAF(Q1,GAF,GFA,kA,kF,EXPQF,
     & kq1,kq2,kAF1,kAF2,kFA1,kFA2,ke1,ke2)
	real*8 EXPQF(ke1,ke2),one,det
	real*8 GAF(kAF1,kAF2),GFA(kFA1,kFA2)
	real*8 Q1(kq1,kq2)
c
c To calculate eGAF in HJCDIST (or eGFA when A,F interchanged in call)
c for calc of initial vector. Result returned in Q1; other args are inputs.
c Modif 12/01/03 09:49am by using matinv2 for inversion
c Modified 10/09/01 09:05am so ALL declared dimensions included in call
	one=1.0d0
c Get I-EXPQF in Q1
	do i=1,kF
	   do j=1,kF
		if(i.eq.j) then
		Q1(i,j)=one-EXPQF(i,j)
		else
		   Q1(i,j)=-EXPQF(i,j)
		endif
	   enddo
	enddo
c premult by GAF
      call MATMUL(GAF,Q1,Q1,kA,kF,kF,one,
     & kAF1,kAF2,kq1,kq1,kq1,kq1)
c postmult by GFA
      call MATMUL(Q1,GFA,Q1,kA,kF,kA,one,
     & kq1,kq1,kFA1,kFA2,kq1,kq1)
c Subtract this from unit matrix
	do i=1,kA
	   do j=1,kA
		if(i.eq.j) then
		   Q1(i,j)=one-Q1(i,j)
		else
		   Q1(i,j)=-Q1(i,j)
		endif
	   enddo
	enddo
c Invert it
      call MATINV2(Q1,kA,kq1,Q1,kq1,.false.,det,ndscale)
c Postmult by GAF
      call MATMUL(Q1,GAF,Q1,kA,kA,kF,one,
     & kq1,kq2,kAF1,kAF2,kq1,kq1)
c Postmult by exp(QFF*tres)
      call MATMUL(Q1,EXPQF,Q1,kA,kF,kF,one,
     & kq1,kq2,ke1,ke2,kq1,kq1)
c
	RETURN
	end

