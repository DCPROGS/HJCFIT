	subroutine PDFopen(QM,phio,area,tau,ncomp,km)
c PDF for all openings
c Altered so outputs tau (ms), area (single prec) and ncomp (as req in PDFIT)
c Outputs eigen,w1 as req for printing of PDFOUT1 etc
c Part of family of subroutines to do single channel calcs via
c set of discrete subroutines that can be called as req.
c   There is a problem in designating submatrices in that notation
c is somewat different in burst and cluster programs, eg F or T for all
c shut states. Here just work directly from kA,kB,kC,kD
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*4 tau(km),area(km)
	real*8 QM(km,km),phio(1,km)
	allocatable::W1,UCOL,col1,eigen,Q1,amat
	real*8  Q1(:,:),amat(:,:,:)
	real*8  W1(:),UCOL(:,:),col1(:,:),eigen(:)
	integer AA
c	integer AA,AB,AC,AD,BA,BB,BC,CB,CA,EE,AF,FA,FF,FD
c	integer BD,CD,DF,DA,TT,AT,TA,AH,BH,HA,HB,GG
c	common/subblk/AA,AB,AC,AD,BA,BB,BC,CB,CA,EE,AF,FA,FF,FD,
c     & BD,CD,DF,DA,TT,AT,TA,AH,BH,HA,HB,GG
c	common/calblk/kE,kG,kF,kT,kH,k,zero,one
	common/KBLK/KA,KB,KC,KD
c
c	kE=kA+kB
c	kG=kA+kB+kC
c	kF=kB+kC
c	kH=kB+kC
c	kT=kB+kC+kD
c	k=kA+kB+kC+kD
	AA=11
	one=1.0d0
	ALLOCATE(Q1(kA,kA),amat(kA,kA,kA))
	ALLOCATE(W1(kA),UCOL(kA,1),col1(kA,1),eigen(kA))
	do i=1,kA
	   ucol(i,1)=one
	enddo
c
	call SUBMAT(QM,AA,Q1,km,km,kA,kA)	!QAA in Q1
	call MATMUL(Q1,ucol,col1,kA,kA,1,-one,
     &	kA,kA,kA,1,kA,1)  !-QAA*uA in col1
	call QMAT5(Q1,Amat,kA,eigen,IBAD,kA,kA,kA)
	do m=1,kA		!calc w1(m)
	   call MATSCL3(phio,Amat,m,col1,kA,w1(m),km,kA,kA)
	enddo
c
	ncomp=kA
	mbad=0
	do m=1,ncomp
	  x=dabs(eigen(m))
	  if(x.lt.1.0d-35) then
		tau(m)=1.e35
		area(m)=1.0
		mbad=m
	  else
		tau(m)=1000./sngl(x)
		area(m)=sngl(w1(m)/x)
	  endif
	enddo
c area set to 1.0 for bad component, so set others to zero
	if(mbad.ne.0) then
	   do m=1,ncomp
		if(m.ne.mbad) area(m)=0.0		!area(mbad)=1 (above)
	   enddo
	endif
c
	DEALLOCATE(Q1,amat,W1,ucol,col1,eigen)
	RETURN
	end


