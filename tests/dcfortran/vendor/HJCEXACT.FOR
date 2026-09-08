	subroutine HJCEXACT(QD,k,kA,kF,EXPQF,QEXPQF,EXPQA,QEXPQA,QAF,QFA,
     & Z00A,Z10A,Z11A,Z00F,Z10F,Z11F,eigen,kAm,kFm,kAx,ibad,eig1)
c
c To calc matrices for exact HJC distributions, for HJCLIK/HJCFIT (based on
c calcs in \fortran\calc\HJCDIST)
c (NB km is in common)
c Have already got (calc in HJCASYMP) QAF*exp(QFF*tres) in QEXPQF (kA*kF)
c Have already got (calc in HJCASYMP) QFA*exp(QAA*tres) in QEXPQA (kF*kA)
c where EXPQF=exp(QFF*tres), and EXPQA=exp(QAA*tres)
c To calculate likelihood for times in 1st deadtime need matrix, Z00(m) say,
c	Z00(m)=Ci00(m)*QAF*EXPQF
c for t in 2nd deadtime need
c	Z10(m)=Ci10(m)*QAF*EXPQF
c	Z11(m)=Ci11(m)*QAF*EXPQF
c and similarly for shut
c Array sizes
c Z00A(kA,kF,k),Z00F(kF,kA,k),Z10A(kA,kF,k),Z10A(kA,kF,k)
c Z11A(kA,kF,k),Z11F(kF,kA,k)
c  -now make the local ones allocatable
c D(kA,kA,k),D(kF,kF,k)	!later,
c Ci10(kA,kA,k)
c Ci10(kF,kF,k)	!later
c Ci11(kA,kA,k)
c Ci11(kF,kF,k)	!later
c Ci00(kF,kF,k)	!not used
	real*8 tres,one
	real*8 QD(100,100)
	real*8 eigen(km),eig1
c	real*8 QAF(10,10),QFA(10,10),EXPQA(10,10),EXPQF(10,10)
	real*8 QAF(kAx,100),QFA(100,kAx)
	real*8 EXPQA(kAx,kAx),EXPQF(100,100)
c	real*8 Z00A(10,10,10),Z10A(10,10,10),Z11A(10,10,10)
c	real*8 Z00F(10,10,10),Z10F(10,10,10),Z11F(10,10,10)
c	real*8 QEXPQA(10,10),QEXPQF(10,10)
	real*8 Z00A(kAm,kFm,km),Z10A(kAm,kFm,km),Z11A(kAm,kFm,km)
	real*8 Z00F(kFm,kAm,km),Z10F(kFm,kAm,km),Z11F(kFm,kAm,km)
	real*8 QEXPQA(kFm,kAM),QEXPQF(kAm,kFm)
c Local arrays
c	real*8 Q1(10,10),Q2(10,10),Q3(10,10)
c	real*8 QMAT(10,10,10),D(10,10,10),eigen(10)
c	real*8 Ci10(10,10,10),Ci11(10,10,10)
	allocatable:: Q1,Q2,Q3
	real*8 Q1(:,:),Q2(:,:),Q3(:,:)
	allocatable:: QMAT,D
	real*8 QMAT(:,:,:),D(:,:,:)
	allocatable:: Ci10,ci11
	real*8 Ci10(:,:,:),Ci11(:,:,:)
	logical caplock,debug
	logical discprt
	common/detw2/tres,km,nerr			!for HJCLIK,DETWA,DETWF
	common/dp/discprt
	common/deb/idebug,idebug1
c
	debug()=caplock()
c
c NB here and elsewhere, these could be made k=number of states, rather than
c km=100
	ALLOCATE(Q1(k,k),Q2(k,k),Q3(k,k))
	ALLOCATE(QMAT(k,k,k),D(k,k,k))
	ALLOCATE(Ci10(k,k,k),Ci11(k,k,k))
	one=1.d0
c
c EXACT DISTRIBUTIONS
c NB The index i is used to denote components (1,..,k) in the names Ci10
c etc, because this notation is used in the papers. But in the program, the
c index m (and m1) are used for this purpose, so that i,j can be used for
c the elements of these arrays as usual.  Note that m (and m1) go from 1 to k.
c Expand whole of Q
c Get eigenvalues of -Q [as in HJC papers; not eigenvalues of Q which leads
c to sign problems in calculation of Ci10(m) and g10(m), at least for
c shut states -odd number of states? -took days to discover this!]
	do i=1,k
	   do j=1,k
	    	Q1(i,j)=-QD(i,j)
	   enddo
	enddo
	call QMAT5(Q1,Qmat,k,eigen,ibad,k,k,km)
	if(dabs(eigen(1)).gt.1.d-15) then
	   eig1=eigen(1)	!for printout
	   eigen(1)=0.0d0
	endif
c  Test pairs of eigenvalues (as used below) to check whether any are equal
c (this will give divide by zero below if true)
	do m=1,k
	   do m1=1,k
		if(m.ne.m1) then
		 if(dabs(eigen(m1)-eigen(m)).lt.1.d-15) then
		   print 1,m,eigen(m),m1,eigen(m1)
	         if(discprt) write(8,1) m,eigen(m),m1,eigen(m1)
1		   format(' Two ''equal'' eigenvalues of Q in HJCEXACT',/,
     &		' eig(',i3,') = ',g13.6,' eig(',i3,') = ',g13.6)
		 endif
		endif
	   enddo
	enddo
c
	if(debug()) then
	   call ATYPd(QD,'Q Matrix',k,k,km,km)
	   call ATYPd(EXPQF,' EXPQF  ',kF,kF,km,km)
	   call ATYPd(QEXPQF,'QEXPQF  ',kA,kF,kAm,kFm)
         print 110
         if(discprt) write(8,110)
110	   format(' A matrices for whole -Q matrix')
	   do m=1,k
		call ATYPD3(qmat,m,'A matrix',k,k,k,k,k)
	   enddo
	endif
c  (1) OPEN TIME
c get Z00A(m)=C00(m)*QAF*EXPQ=[A(m)]AA*QAF*EXPQ (kA*kA x kA*kF = kA*kF)
	do m=1,k
	   do i=1,kA
		do j=1,kF
		   Z00A(i,j,m)=0.0d0
		   do n=1,kA
			Z00A(i,j,m)=Z00A(i,j,m)+QMAT(i,n,m)*QEXPQF(n,j)
		   enddo
		enddo
	   enddo
	enddo
	if(debug()) then
	   do m=1,k
		call ATYPD3(Z00A,m,' Z00A   ',kA,kF,kAm,kFm,km)
	   enddo
	endif
c
c NOW 2nd deadtime
c For f1(t) first calculate the k matrices D(m)=[A(m)]AF*expQF*QFA (kA*kA),
c and hence the k values of g10A(m) and g11A(m).
	call MATMUL(expQF,QFA,Q1,kF,kF,kA,one,   !expQF*QFA in Q1
     & km,km,km,kAx,k,k)
	do m=1,k
c       need version of SUBMAT3 that allows all array dimensions to be specified
	   call SUBMAT3d(QMAT,m,15,Q2,k,k,k,k,k)   !15=AF section of QMAT(m) in Q2
	   call MATMUL(Q2,Q1,Q3,kA,kF,kA,one,   !D(m) in Q3
     &    k,k,k,k,k,k)
	   do i=1,kA
		do j=1,kA
		   D(i,j,m)=Q3(i,j)			!copy to D
		enddo
	   enddo
	   if(debug()) call ATYPD3(D,m,'D matrix',kA,kA,k,k,k)
	enddo
c
c D(i,j,m) now contains the D matrix for component m (kA*kA)
c Now do summation to get Ci10() array from the D() and the Ci00 which are
c the AA section of QMAT().
	do m=1,k
	   do i=1,kA
		do j=1,kA
		   Ci10(i,j,m)=0.0d0       !zero Ci10
		enddo
	   enddo
	enddo
c
	do m=1,k
	   do m1=1,k
c need [D(m)*Ci00(m1) + D(m1)*Ci00(m)]/(eigen(m1)-eigen(m))
c where D(m)=D(i,j,m) and Ci00(m)=[Qmat(i,j,m)]AA
		if(m.ne.m1) then
c===matmul3 assumes all arrays same dimensions -D, qmat and Q2 are all dim=k here
		   call MATMUL3(D,QMAT,Q2,kA,kA,kA,m,m1,1,kA,1,kA,
     &	     1,kA,1,kA,k)
		   call MATMUL3(D,QMAT,Q3,kA,kA,kA,m1,m,1,kA,1,kA,
     &	     1,kA,1,kA,k)
		   do i=1,kA
			do j=1,kA
			   Ci10(i,j,m)=Ci10(i,j,m)+(Q2(i,j)+Q3(i,j))/
     &			(eigen(m1)-eigen(m))
			enddo
		   enddo
		endif
	   enddo	!end of m1 loop
	   if(debug()) call ATYPD3(Ci10,m,' Ci10   ',kA,kA,k,k,k)
	enddo	!end of m loop
c
c and get Ci11(m)=D(m)*Ci00(m)
	do m=1,k
	   call MATMUL3(D,QMAT,Q2,kA,kA,kA,m,m,1,kA,1,kA,
     &  1,kA,1,kA,k)
c     &  1,kA,1,kA,km)
	   do i=1,kA
		do j=1,kA
		   Ci11(i,j,m)=Q2(i,j)
		enddo
	   enddo
	   if(debug()) call ATYPD3(Ci11,m,' Ci11   ',kA,kA,k,k,k)
c Multiply Ci10 and Ci11 by QEXPQF to get Z10A and Z11A
	   do i=1,kA
		do j=1,kF
		   Z10A(i,j,m)=0.0d0
		   Z11A(i,j,m)=0.0d0
		   do n=1,kA
			Z10A(i,j,m)=Z10A(i,j,m)+Ci10(i,n,m)*QEXPQF(n,j)
			Z11A(i,j,m)=Z11A(i,j,m)+Ci11(i,n,m)*QEXPQF(n,j)
		   enddo
		enddo
	   enddo
	enddo    !end of m loop
c
c End of exact open time calcs
c EXACT SHUT TIME CALCS
c Get Z00F(m)=C00(m)*QFA*EXPQA=[A(m)]FF*QFA*EXPQA (kF*kF x kF*kA = kF*kA)
	do m=1,k
	   do i=1,kF
		i1=i+kA			!=kA+1 to k, ie F states in QMAT
		do j=1,kA
		   Z00F(i,j,m)=0.0d0
		   do n=1,kF
			n1=n+kA		!=kA+1 to k, ie F states in QMAT
			Z00F(i,j,m)=Z00F(i,j,m)+QMAT(i1,n1,m)*QEXPQA(n,j)
		   enddo
		enddo
	   enddo
	enddo
c
c For f1(t) first calculate the k matrices D(m)=[A(m)]FA*expQA*QAF (kF*kF),
c and hence the k values of g10F(m) and g11F(m).
	call MATMUL(expQA,QAF,Q1,kA,kA,kF,one,   !expQA*QAF in Q1 (kA*kF)
     & kAx,kAx, kAx,km, k,k)
	do m=1,k
c       Use version of SUBMAT3 that allows all array dimensions to be specified
	   call SUBMAT3d(QMAT,m,51,Q2,k,k,k,k,k)   !51=FA section of QMAT(m) in Q2
	   call MATMUL(Q2,Q1,Q3,kF,kA,kF,one,
     &     k,k, k,k, k,k)   !D(m) in Q3
	   do i=1,kF
		do j=1,kF
		   D(i,j,m)=Q3(i,j)			!copy to D
		enddo
	   enddo
	   if(debug()) call ATYPD3(D,m,'D matrix',kF,kF,k,k,k)
	enddo
c All D() now calc
c D(i,j,m) now contains the D matrix for component m (kF*kF)
c
c Now do summation to get Ci10() array from the D() and the Ci00 which are
c the FF section of QMAT().
	do m=1,k
	   do i=1,kF
		do j=1,kF
		   Ci10(i,j,m)=0.0d0       !zero Ci10
		enddo
	   enddo
	enddo
c
	do m=1,k
	   do m1=1,k
c need [D(m)*Ci00(m1) + D(m1)*Ci00(m)]/(eigen(m1)-eigen(m))
c where D(m)=D(i,j,m) and Ci00(m)=[Qmat(i,j,m)]FF
c NB NEED TO USE ONLY THE FF SECTION OF QMAT(k*k) with D (kF*kF)
c===matmul3 assumes all arrays same dimensions -D, qmat and Q2 are all dim=k here
		if(m.ne.m1) then
		   call MATMUL3(D,QMAT,Q2,kF,kF,kF,m,m1,1,kF,1,kF,
     &	     kA+1,k,kA+1,k,k)
		   call MATMUL3(D,QMAT,Q3,kF,kF,kF,m1,m,1,kF,1,kF,
     &	      kA+1,k,kA+1,k,k)
		   do i=1,kF
			do j=1,kF
	 		   Ci10(i,j,m)=Ci10(i,j,m)+(Q2(i,j)+Q3(i,j))/
     &			(eigen(m1)-eigen(m))
			enddo
		   enddo
		endif
	    enddo	!end of m1 loop
	   if(debug()) call ATYPD3(Ci10,m,' Ci10   ',kF,kF,k,k,k)
	enddo		!end of m loop
c
c and get Ci11(m)=D(m)*Ci00(m)
	do m=1,k
	   call MATMUL3(D,QMAT,Q2,kF,kF,kF,m,m,1,kF,1,kF,
     &      kA+1,k,kA+1,k,k)
	   do i=1,kF
		do j=1,kF
149		   Ci11(i,j,m)=Q2(i,j)
		enddo
	   enddo
	   if(debug()) call ATYPD3(Ci11,m,' Ci11   ',kF,kF,k,k,k)
c Multiply Ci10 and Ci11 by QEXPQA to get Z10F and Z11F
	   do i=1,kF
		do j=1,kA
		   Z10F(i,j,m)=0.0d0
		   Z11F(i,j,m)=0.0d0
		   do n=1,kF
			Z10F(i,j,m)=Z10F(i,j,m)+Ci10(i,n,m)*QEXPQA(n,j)
			Z11F(i,j,m)=Z11F(i,j,m)+Ci11(i,n,m)*QEXPQA(n,j)
		   enddo
		enddo
	   enddo
	enddo          !end of m loop
C End of exact shut time calcs
	DEALLOCATE(Q1,Q2,Q3,QMAT,D,Ci10,Ci11)
	RETURN
	end
