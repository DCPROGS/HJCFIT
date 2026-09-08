	subroutine PDFshut(QM,phio,area,tau,ncomp,km)
c Altered so outputs tau (ms), area (single prec) and ncomp (as req in PDFIT)
c PDF for all shuttings=phio*GAT*exp(QTT*t)*(-QTT)*uT
c Outputs eigen,w1 as req for printing of PDFOUT1 etc
c
c Part of family of subroutines to do single channel calcs via
c set of discrete subroutines that can be called as req.
c   There is a problem in designating submatrices in that notation
c is somewat different in burst and cluster programs, eg F or T for all
c shut states. Here just work directly from kA,kB,kC,kD
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*4 tau(km),area(km)
	real*8 QM(km,km),phio(1,km)
	allocatable::W1,UCOL,col1,row1,eigen,Q1,amat
	real*8  Q1(:,:),amat(:,:,:)
	real*8  W1(:),UCOL(:,:),col1(:,:),row1(:,:),eigen(:)
c	real*4 tau(10),area(10)
c	dimension W1(10),ucol(10,1),col1(10,1),eigen(10),phio(1,10)
c	dimension QM(10,10),Q1(10,10),amat(10,10,10),row1(1,10)
	integer AT,TT
	common/KBLK/KA,KB,KC,KD
c
	if(kD.ne.0) then
	   print 1
1	   format(' WARNING: kD NOT ZERO IN PDFSHUT')
	endif
c
	kT=kB+kC+kD
	kF=kT
	AT=17
	TT=77
	one=1.0d0
	ALLOCATE(Q1(kT,kT),amat(kT,kT,kT))
	ALLOCATE(W1(kT),UCOL(kT,1),col1(kT,1),row1(1,kT),eigen(kT))
	do i=1,kT
	   ucol(i,1)=one
	enddo
c Initial vector=phio*GAT
	call GMAT1(AT,QM,Q1,km,kT)	!GAT in Q1
      call MATMUL(phio,Q1,row1,1,kA,kT,one,
     &  1,kT,kT,kT,1,kT)  !phio*GAT=row1
c Final vector
	call SUBMAT(QM,TT,Q1,km,km,kT,kT)	!QTT in Q1
	call MATMUL(Q1,ucol,col1,kT,kT,1,-one,
     &	kT,kT,kT,1,kT,1)  !-QTT*uT in col1
c Calc pdf
	call QMAT5(Q1,Amat,kT,eigen,IBAD,kT,kT,kT)
	do m=1,kT		!calc w1(m)
	   call MATSCL3(row1,Amat,m,col1,kT,w1(m),kT,kT,kT)
	enddo
c
	ncomp=kT
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
	DEALLOCATE(Q1,amat,W1,ucol,col1,row1,eigen)
	RETURN
	end





