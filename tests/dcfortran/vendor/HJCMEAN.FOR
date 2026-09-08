	subroutine HJCMEAN(QD,tres,hmopen,hmshut,
     &   k,kA,kF,kAm,kFm,km)
c ===all these args are defined in hjcdisp for call
c	subroutine HJCMEAN(g00,g10,g11,eigen,tres,k,km),===
c
c To calculate exact mean open or shut time from HJC pdf
c Called in HJCDISP

cc g00 (etc) in call is g00A for open times, g00F for shut times, for exact pdf
cc amp() tau() are for asymptotic pdf (calc in HJCDISP) -as in hjcmean.mcd
c NO -much better, and exact method, is to use dARSdS for open times (as
c in MOPADJC for conditional mean open times, but simpler!).
c And write equivalant dFRSdS for shut times.
c
c Modif 08/10/04 11:03am, so that QEXPQA, QEXPQF, phiA, phiF are all
c recalculated here, locally, from current QD (if conc changed, QD may not
c correspond to global values of these matrices).
c They are removed from arguments.
c Also make special versions of dARSdS and dARFdS that are local to
c this subroutine and have QAA, EXPQA atc as arguments, to avoid
c corruption of the global values in common/detw1/
c
c Output = overall mean open and shut times in real*8 hmopen, hmshut (msec)
c
c For mean open
c  E(topen) = phiA*{[-d((e)GAF(s)/ds]s=0}*uF
c  		= tres + phiA*{[-d(A)R(s)/ds]s=0}*QAF*exp(QFF*tres)*uF
c and analagously for for E(tshut)
c
c EXPQF=exp(QFF*tres), and EXPQA=exp(QAA*tres) already calc in HJCEXACT
c (and values are in common/detw1/ from last call to HJCLIK)
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 QD(km,km)
	real*8 ucol(100,1),col1(100,1)
	allocatable:: QAF,QFA,QAA,QFF
	real*8 QAF(:,:),QFA(:,:),QAA(:,:),QFF(:,:)
	allocatable::GAF,GFA
	real*8 GAF(:,:),GFA(:,:)
	allocatable::invQAA,invQFF
	real*8 invQAA(:,:),invQFF(:,:)

C=	common/phiblk/phiF(1,100),phiA(1,100)	!now local
	allocatable::DARS,DFRS,Q1,Q2
c
	allocatable::QEXPQA,QEXPQF,phiF,phiA		!local
	real*8 DARS(:,:),DFRS(:,:),Q1(:,:),Q2(:,:),phiF(:,:),phiA(:,:)
	real*8 QEXPQA(:,:),QEXPQF(:,:)
	allocatable EXPQA,EXPQF
	real*8 EXPQA(:,:),EXPQF(:,:)
c==	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF,dARSdS
	allocatable::FMAT,EIGFF,AMAT,EIGAA
	real*8 FMAT(:,:,:), EIGFF(:),AMAT(:,:,:), EIGAA(:)
	INTEGER AF,FA,AA,FF
	logical deb
c
	deb=.false.
	AA=11
	AF=15		!5=CODE FOR F(=B+C) SECTION
	FA=51
	FF=55
	one=1.d0
	do i=1,km
	   ucol(i,1)=one
	enddo
	ALLOCATE(DARS(kAm,kAm),DFRS(kFm,kFm))
	ALLOCATE(QEXPQA(kFm,kAM),QEXPQF(kAm,kFm),
     &	phiF(1,100),phiA(1,100))
	ALLOCATE(EXPQA(kAm,kAm),EXPQF(kFm,kFm))	!differnt from dARSdS
	kXm=kAm+kFm
	ALLOCATE(Q1(kXm,kXm),Q2(kXm,kXm))
	ALLOCATE(FMAT(kFm,kFm,kFm),EIGFF(kFm),
     &	AMAT(kAm,kAm,kAm),EIGAA(kAm))
	ALLOCATE(GAF(kAm,kFm),GFA(kFm,kAm))
	ALLOCATE(invQAA(kAm,kAm),invQFF(kFm,kFm))
	ALLOCATE(QAF(kAm,kFm),QFA(kFm,kAm),QAA(kAm,kAm),QFF(kFm,kFm))
c

c First recalculate phiA, phiF from current QD, as in hjcasymp
c Calculate GAF and GFA -do explicitly (GMATn gets complicated when arrays
c all have different declared sizes)
	call SUBMAT(QD,AA,QAA,km,km,kAm,kAm)
	call SUBMAT(QD,AF,QAF,km,km,kAm,kFm)
	call MATINV2(QAA,kA,kAm,invQAA,kAm,.false.,det,ndscale)
	call MATMUL(invQAA,QAF,GAF,kA,kA,kF,-one,
     &	kAm,kAm,kAm,kFm,kAm,kFm)

	call SUBMAT(QD,FF,QFF,km,km,kFm,kFm)
	call SUBMAT(QD,FA,QFA,km,km,kFm,kAm)
	call MATINV2(QFF,kF,kFm,invQFF,kFm,.false.,det,ndscale)
	call MATMUL(invQFF,QFA,GFA,kF,kF,kA,-one,
     &	kFm,kFm,kFm,kAm,kFm,kAm)
c Expand QFF
	call QMAT5(QFF,Fmat,kF,eigFF,ibad,kFm,kFm,kFm)
c Calc exp(QFF*tres) in EXPQF
	do i=1,kF
	   do j=1,kF
		EXPQF(i,j)=0.d0
		do m=1,kF
		   EXPQF(i,j)=EXPQF(i,j) + FMAT(i,j,m)*dexp1(eigFF(m)*tres)
		enddo
	   enddo
	enddo
c Expand QAA
	call QMAT5(QAA,Amat,kA,eigAA,ibad,kAm,kAm,kAm)
c Calc exp(QAA*tres) in EXPQA
	do i=1,kA
	   do j=1,kA
		EXPQA(i,j)=0.d0
		do m=1,kA
		   EXPQA(i,j)=EXPQA(i,j) + AMAT(i,j,m)*dexp1(eigAA(m)*tres)
		enddo
	   enddo
	enddo
c Calc INITIAL VECTOR: first calc eGAF*eGFA
c Use subroutine to get eGAF in Q1
	call CeGAF(Q1,GAF,GFA,kA,kF,EXPQF,
     & kXm,kXm,kAM,kFm,kFm,kAm,kFm,kFm)
	if(deb) call ATYPD(Q1,' eGAF   ' ,kA,kF,kXm,kXm)
c Now get eGFA by repeating with A and F interchanged
c Use subroutine to get eGFA in Q2
	call CeGAF(Q2,GFA,GAF,kF,kA,EXPQA,
     & kXm,kXm,kFm,kAm,kAm,kFm,kAm,kAm)
	if(deb) call ATYPD(Q2,' eGFA   ' ,kF,kA,kXm,kXm)
	call HJCphi(Q1,Q2,kA,kF,phiA,kXm,kXm,km)
c
c Repeat for phiF
	call HJCphi(Q2,Q1,kF,kA,phiF,kXm,kXm,km)
c
c Second recalculate QEXPQA,QEXPQF from current QD
c Get QAF*EXPQF in QEXPQF: this is kA*kF
      call MATMUL(QAF,EXPQF,QEXPQF,kA,kF,kF,one,
     &   kAm,kFm,kFm,kFm,kAm,kFm)
c Get QFA*EXPQA in QEXPQA -kF*kA
      call MATMUL(QFA,EXPQA,QEXPQA,kF,kA,kA,one,
     &   kFm,kAm,kAm,kAm,kFm,kAm)
c
c Now ready to calculate means
c Mean open
c local version of darsds, see below
	call dARSdS1(QD,tres,DARS,QAA,QAF,QFF,QFA,EXPQF,EXPQA,
     &	kA,kF,kAm,kFm,km)
c Already have QAF*exp(QFF*tres) in QEXPQF (kA*kF)
      call MATMUL(DARS,QEXPQF,Q1,kA,kA,kF,one,kAm,kAm,kAm,kFm,kXm,kXm)
c Need Q1*uF -keep in col1 (kA*1)
      call MATMUL(Q1,ucol,col1,kA,kF,1,one,kXm,kXm,km,1,km,1)
c premultiply my phiA
	call VECMUL(phiA,col1,hmopen,kA,km,km)
c Add tres
	hmopen=tres + hmopen
c
c Now same for mean shut
	call dFRSdS1(QD,tres,DFRS,QAA,QAF,QFF,QFA,EXPQF,EXPQA,
     & kA,kF,kAm,kFm,km)
c Already have QFA*exp(QAA*tres) in QEXPQA (kF*kA)
      call MATMUL(DFRS,QEXPQA,Q1,kF,kF,kA,one,kFm,kFm,kFm,kAm,kXm,kXm)
c Need Q1*uA -keep in col1 (kF*1)
      call MATMUL(Q1,ucol,col1,kF,kA,1,one,kXm,kXm,km,1,km,1)
c premultiply my phiF
	call VECMUL(phiF,col1,hmshut,kF,km,km)
c Add tres
	hmshut=tres + hmshut
c
c Convert to milliseconds before return
	hmopen=1.d3*hmopen
	hmshut=1.d3*hmshut
c
	DEALLOCATE(DARS,DFRS,Q1,Q2,QEXPQA,QEXPQF,phiF,phiA)
	DEALLOCATE(EXPQA,EXPQF)
	DEALLOCATE(FMAT,EIGFF,AMAT,EIGAA)
	DEALLOCATE(GAF,GFA,invQAA,invQFF)
	DEALLOCATE(QAF,QFA,QAA,QFF)
	return
	end

	subroutine dARSdS1(QD,tres,DARS,QAA,QAF,QFF,QFA,EXPQF,EXPQA,
     & kA,kF,kAm,kFm,km)

c To evaluate [-dAR(s)/ds] at s=0 for means etc of HJC distributions.
c [GAF,SFF,GFA do not need to be args]
c SFF=I-exp(QFF*tres).
c Result is output as the kAxkA matrix, DARS
c
c DARSDS1 is local version within HJCMEAN that has common/detw1/ removed
c and values passed as arguments (dimensions not same as in dARSdS)
c
c First evaluate [dVA(s)/dS]s=0
c NOT: Q1=-inv(QAA)*GAF*SFF*GFA - GAF*SFF*inv(QFF)*GFA + tres*GAF*exp(QFF*tres)*GFA
c CORRECTION FOR SIGN OF LAST TERM (AGH 26-May-92) -should be
c Q1=-inv(QAA)*GAF*SFF*GFA - GAF*SFF*inv(QFF)*GFA - tres*GAF*exp(QFF*tres)*GFA
c Then DARS=inv(VA)*(QAA**-2) - inv(VA)*Q1*inv(VA)*inv(QAA)
c = inv(VA)[inv(QAA) - Q1*inv(VA)]*inv(QAA)
c where VA=I - GAF*SFF*GFA
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	PARAMETER(kAx=60)		!arrays changed so not used now
	real*8 QD(km,km),one,det
c=	real*8 QAA(10,10),QFF(10,10),QAF(10,10),QFA(10,10),SFF(10,10)
c=	real*8 EXPQF(10,10),EXPQA(10,10)
c=	real*8 QAF(kAx,100),QFA(100,kAx),QAA(kAx,kAx),QFF(100,100)
c=	real*8 EXPQA(kAx,kAx),EXPQF(100,100)
	real*8 QAF(kAm,kFm),QFA(kFm,kAm),QAA(kAm,kAm),QFF(kFm,kFm)
	real*8 EXPQA(kAm,kAm),EXPQF(kFm,kFm)
	real*8 DARS(kAm,kAm)
c===	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF,dARSdS
c All the rest are local apart from QD
c Should they be allocated? See TMEM.FOr and TMEMA -uses much less memory
c if allocated

	ALLOCATABLE::Q1,Q2,Q3,GAF,GFA,INVQAA,INVQFF,SFF
	real*8 Q1(:,:),Q2(:,:),Q3(:,:)
	real*8 GAF(:,:),GFA(:,:),INVQAA(:,:),INVQFF(:,:),SFF(:,:)
	integer AA,AF,FA,FF

c
c EXPQF=exp(QFF*tres), and EXPQA=exp(QAA*tres) already calc in HJCEXACT
c (and values are in common/detw1/ from last call to HJCLIK)
	AA=11
	AF=15	!5=CODE FOR F(=B+C) SECTION
	FA=51
	FF=55
	one=1.0d0
	ALLOCATE(Q1(km,km),Q2(km,km),Q3(km,km),
     & GAF(kAm,kFm),GFA(kFm,kAm),INVQAA(kAm,kAm),INVQFF(kFm,kFm),
     & SFF(kFm,kFm))
c
c Calculate GAF and GFA -do explicitly (GMATn gets complicated when arrays
c all have different declared sizes)
c all GMAT2 replaced by GMAT4
c	subroutine GMAT4(AB,QM,GAB,invQAA,km1,km2,kg1,kg2,ki1,ki2)
c	call GMAT4(FA,QD,GFA,INVQFF,km,km,kFm,kAm,kFm,kFm)
	call SUBMAT(QD,AA,QAA,km,km,kAm,kAm)
	call MATINV2(QAA,kA,kAm,invQAA,kAm,.false.,det,ndscale)
	call SUBMAT(QD,AF,Q1,km,km,km,km)
	call MATMUL(invQAA,Q1,GAF,kA,kA,kF,-one,
     &	kAm,kAm,km,km,kAm,kFm)
c	call ATYPD(GAF,'   GAF= ',kA,kF,kAm,kFm)
	call SUBMAT(QD,FF,QFF,km,km,kFm,kFm)
	call MATINV2(QFF,kF,kFm,invQFF,kFm,.false.,det,ndscale)
	call SUBMAT(QD,FA,Q1,km,km,km,km)
	call MATMUL(invQFF,Q1,GFA,kF,kF,kA,-one,
     &	kFm,kFm,km,km,kFm,kAm)
c	call ATYPD(GFA,'   GFA= ',kF,kA,kFm,kAm)
c check
	call MATMUL(QAA,invQAA,Q1,kA,kA,kA,one,
     & kAm,kAm,kAm,kAm,kAm,kAm)
c	call ATYPD(Q1,'  unit= ',kA,kA,km,km)
	call MATMUL(QFF,invQFF,Q1,kF,kF,kF,one,
     & kFm,kFm,kFm,kFm,km,km)
c	call ATYPD(Q1,'  unit= ',kF,kF,km,km)
c
c	call ATYPD(GFA,'   GFA= ',kF,kA,kFm,kAm)
c calc SFF=I-EXPQF
	do i=1,kF
	   do j=1,kF
		SFF(i,j)=-EXPQF(i,j)
		if(i.eq.j) SFF(i,j)=SFF(i,j)+one
	   enddo
	enddo
c	call ATYPD(SFF,'   SFF= ',kF,kF,kFm,kFm)
c Last term of Q1 above = tres*GAF*exp(QFF*tres)*GFA  in Q1
	call MATMUL(expQF,GFA,Q1,kF,kF,kA,tres,	  !tres*expQF*GFA in Q1
     & kFm,kFm,kFm,kAm,km,km)
c	call ATYPD(Q1,'  Q1-1= ',kF,kA,km,km)
c Correct sign of 3rd term (see above)
	call MATMUL(GAF,Q1,Q1,kA,kF,kA,-one,  !-tres*GAF*expQF*GFA in Q1
     & kAm,kFm,km,km,km,km)
c	call ATYPD(Q1,'  Q1-2= ',kA,kA,km,km)
c 2nd term of Q1 above = GAF*SFF*inv(QFF)*GFA in Q2
	call MATMUL(GAF,SFF,Q2,kA,kF,kF,one,	  !GAF*SFF in Q2
     & kAm,kFm,kFm,kFm,km,km)
c	call ATYPD(Q2,'  Q2-1= ',kA,kF,km,km)
c	call ATYPD(invQFF,' invQFF ',kF,kF,kFm,kFm)

	call MATMUL(Q2,invQFF,Q2,kA,kF,kF,one,	  !GAF*SFF*inv(QFF) in Q2
     & km,km,kFm,kFm,km,km)
c	call ATYPD(Q2,'  Q2-2= ',kA,kF,km,km)
	call MATMUL(Q2,GFA,Q2,kA,kF,kA,-one,	  !-GAF*SFF*inv(QFF)*GFA in Q2
     & km,km,kFm,kAm,km,km)
c	call ATYPD(Q2,'   Q2-3 ',kA,kA,km,km)
c first term = -inv(QAA)*GAF*SFF*GFA in Q3
	call MATMUL(invQAA,GAF,Q3,kA,kA,kF,one,	  !inv(QAA)*GAF in Q3
     & kAm,kAm,kAm,kFm,km,km)
c	call ATYPD(Q3,'   Q3-1 ',kA,kF,km,km)
	call MATMUL(Q3,SFF,Q3,kA,kF,kF,one,		   !inv(QAA)*GAF*SFF in Q3
     & km,km,kFm,kFm,km,km)
c	call ATYPD(Q3,'   Q3-2 ',kA,kF,km,km)
	call MATMUL(Q3,GFA,Q3,kA,kF,kA,-one,	   !-inv(QAA)*GAF*SFF*GFA in Q3
     & km,km,kFm,kAm,km,km)
c	call ATYPD(Q3,'   Q3-3 ',kA,kA,km,km)
c Add Q1+Q2+Q3; result in Q1
      call MATADD(Q1,Q2,Q1,kA,kA,
     & km,km,km,km,km,km)
      call MATADD(Q1,Q3,Q1,kA,kA,
     & km,km,km,km,km,km)
c	call ATYPD(Q1,'   Q1-3 ',kA,kA,km,km)
c Calc VA=I - GAF*SFF*GFA in Q2
	call MATMUL(GAF,SFF,Q2,kA,kF,kF,one,
     & kAm,kFm,kFm,kFm,km,km)
c	call ATYPD(Q2,'   Q2-4 ',kA,kF,km,km)
	call MATMUL(Q2,GFA,Q2,kA,kF,kA,one,
     & km,km,kFm,kAm,km,km)
c	call ATYPD(Q2,'   Q2-5 ',kA,kA,km,km)
      call MATIMIN(Q2,Q2,kA,km)	!(I-Q2) in Q2
c	call ATYPD(Q2,'   Q2-6 ',kA,kA,km,km)
      call MATINV2(Q2,kA,km,Q2,km,.false.,det,ndscale)	!inv(VA) in Q2
c	call ATYPD(Q2,'   Q2-7 ',kA,kA,km,km)
c Calc DARS=inv(VA)*(QAA**-2) - inv(VA)*Q1*inv(VA)*inv(QAA)
c = inv(VA)[inv(QAA) - Q1*inv(VA)]*inv(QAA)
	call MATMUL(Q1,Q2,Q3,kA,kA,kA,-one,
     & km,km,km,km,km,km)                       !-Q1*inv(VAA) in Q3
c	call ATYPD(Q3,'   Q3-4 ',kA,kA,km,km)
      call MATADD(invQAA,Q3,Q3,kA,kA,
     & kAm,kAm,km,km,km,km)                     !inv(QAA)-Q1*inv(VAA) in Q3
c	call ATYPD(Q3,'   Q3-5 ',kA,kA,km,km)
	call MATMUL(Q2,Q3,Q3,kA,kA,kA,one,
     & km,km,km,km,km,km)                       !inv(VA)*Q3 in Q3
c	call ATYPD(Q3,'   Q3-6 ',kA,kA,km,km)
	call MATMUL(Q3,invQAA,DARS,kA,kA,kA,one,
     & km,km,kAm,kAm,kAm,kAm)                   !inv(VA)*Q3*invQAA in DARS
c	call ATYPD(DARS,'  DARS= ',kA,kA,kAm,kAm)
c
	DEALLOCATE(Q1,Q2,Q3,GAF,GFA,INVQAA,INVQFF)
	RETURN
	end


	subroutine dFRSdS1(QD,tres,DFRS,QAA,QAF,QFF,QFA,EXPQF,EXPQA,
     & kA,kF,kAm,kFm,km)
c To evaluate [-dFR(s)/ds] at s=0 for means etc of HJC shut time distributions.
c (Modified from dARSdS)
c [GAF,SFF,GFA do not need to be args]
c
c DFRSDS1 is local verion within HJCMEAN that has common/detw1/ removed
c and values passed as arguments (dimensions not same as in dARSdS)
c
c SAA=I-exp(QAA*tres).
c Result is output as the kFxkF matrix, DFRS
c
c Rest of notes are from dARSdS, so interchange A, F in them.
c
c First evaluate [dVA(s)/dS]s=0
c NOT: Q1=-inv(QAA)*GAF*SFF*GFA - GAF*SFF*inv(QFF)*GFA + tres*GAF*exp(QFF*tres)*GFA
c CORRECTION FOR SIGN OF LAST TERM (AGH 26-May-92) -should be
c Q1=-inv(QAA)*GAF*SFF*GFA - GAF*SFF*inv(QFF)*GFA - tres*GAF*exp(QFF*tres)*GFA
c Then DARS=inv(VA)*(QAA**-2) - inv(VA)*Q1*inv(VA)*inv(QAA)
c = inv(VA)[inv(QAA) - Q1*inv(VA)]*inv(QAA)
c where VA=I - GAF*SFF*GFA
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	PARAMETER(kAx=60)		!arrays changed so not used now
	real*8 QD(km,km),one,det
c=	real*8 QAA(10,10),QFF(10,10),QAF(10,10),QFA(10,10),SFF(10,10)
c=	real*8 EXPQF(10,10),EXPQA(10,10)
c=	real*8 QAF(kAx,100),QFA(100,kAx),QAA(kAx,kAx),QFF(100,100)
c=	real*8 EXPQA(kAx,kAx),EXPQF(100,100)
	real*8 QAF(kAm,kFm),QFA(kFm,kAm),QAA(kAm,kAm),QFF(kFm,kFm)
	real*8 EXPQA(kAm,kAm),EXPQF(kFm,kFm)
c==	real*8 DARS(kAm,kAm)
	real*8 DFRS(kFm,kFm)
c===	common/detw1/QAA,QAF,QFF,QFA,EXPQF,EXPQA		  !for DETWA,DETWF,dARSdS
c All the rest are local apart from QD
c Should they be allocated? See TMEM.FOr and TMEMA -uses much less memory
c if allocated

	ALLOCATABLE::Q1,Q2,Q3,GAF,GFA,INVQAA,INVQFF,SAA
	real*8 Q1(:,:),Q2(:,:),Q3(:,:)
	real*8 GAF(:,:),GFA(:,:),INVQAA(:,:),INVQFF(:,:),SAA(:,:)
	integer AA,AF,FA,FF

c
c EXPQF=exp(QFF*tres), and EXPQA=exp(QAA*tres) already calc in HJCEXACT
c (and values are in common/detw1/ from last call to HJCLIK)
	AA=11
	AF=15	!5=CODE FOR F(=B+C) SECTION
	FA=51
	FF=55
	one=1.0d0
	ALLOCATE(Q1(km,km),Q2(km,km),Q3(km,km),
     & GAF(kAm,kFm),GFA(kFm,kAm),INVQAA(kAm,kAm),INVQFF(kFm,kFm),
     & SAA(kAm,kAm))
c
c Calculate GAF and GFA -do explicitly (GMATn gets complicated when arrays
c all have different declared sizes)
c all GMAT2 replaced by GMAT4
c	subroutine GMAT4(AB,QM,GAB,invQAA,km1,km2,kg1,kg2,ki1,ki2)
c	call GMAT4(FA,QD,GFA,INVQFF,km,km,kFm,kAm,kFm,kFm)
c Calc GFA
	call SUBMAT(QD,FF,QFF,km,km,kFm,kFm)
	call MATINV2(QFF,kF,kFm,invQFF,kFm,.false.,det,ndscale)
	call SUBMAT(QD,FA,Q1,km,km,km,km)
	call MATMUL(invQFF,Q1,GFA,kF,kF,kA,-one,
     &	kFm,kFm,km,km,kFm,kAm)
c	call ATYPD(GFA,'   GFA= ',kF,kA,kFm,kAm)
c calc GAF
	call SUBMAT(QD,AA,QAA,km,km,kAm,kAm)
	call MATINV2(QAA,kA,kAm,invQAA,kAm,.false.,det,ndscale)
	call SUBMAT(QD,AF,Q1,km,km,km,km)
	call MATMUL(invQAA,Q1,GAF,kA,kA,kF,-one,
     &	kAm,kAm,km,km,kAm,kFm)
c	call ATYPD(GAF,'   GAF= ',kA,kF,kAm,kFm)
c
c	call MATMUL(QAA,invQAA,Q1,kA,kA,kA,one,
c     & kAx,kAx,kAm,kAm,km,km)
c	call ATYPD(Q1,'  unit= ',kA,kA,km,km)
c	call MATMUL(QFF,invQFF,Q1,kF,kF,kF,one,
c     & km,km,kFm,kFm,km,km)
c	call ATYPD(Q1,'  unit= ',kF,kF,km,km)
c
c	call ATYPD(GFA,'   GFA= ',kF,kA,kFm,kAm)
c calc SAA=I-EXPQA
	do i=1,kA
	   do j=1,kA
		SAA(i,j)=-EXPQA(i,j)
		if(i.eq.j) SAA(i,j)=SAA(i,j)+one
	   enddo
	enddo
c	call ATYPD(SAA,'   SAA= ',kA,kA,kAm,kAm)
	call MATMUL(expQA,GAF,Q1,kA,kA,kF,tres,	  !tres*expQA*GAF in Q1
     & kAm,kAm,kAm,kFm,km,km)
c	call ATYPD(Q1,'  Q1-1= ',kA,kF,km,km)
c Correct sign of 3rd term (see above)
	call MATMUL(GFA,Q1,Q1,kF,kA,kF,-one,  !-tres*GFA*expQA*GAF in Q1
     & kFm,kAm,km,km,km,km)
c	call ATYPD(Q1,'  Q1-2= ',kF,kF,km,km)
c 2nd term of Q1 above = GFA*SAA*inv(QAA)*GAF in Q2
	call MATMUL(GFA,SAA,Q2,kF,kA,kA,one,	  !GFA*SAA in Q2
     & kFm,kAm,kAm,kAm,km,km)
c	call ATYPD(Q2,'  Q2-1= ',kF,kA,km,km)
c	call ATYPD(invQAA,' invQAA ',kA,kA,kAm,kAm)

	call MATMUL(Q2,invQAA,Q2,kF,kA,kA,one,	  !GFA*SAA*inv(QAA) in Q2
     & km,km,kAm,kAm,km,km)
c	call ATYPD(Q2,'  Q2-2= ',kF,kA,km,km)
	call MATMUL(Q2,GAF,Q2,kF,kA,kF,-one,	  !-GFA*SAA*inv(QAA)*GAF in Q2
     & km,km,kAm,kFm,km,km)
c	call ATYPD(Q2,'   Q2-3 ',kF,kF,km,km)
c first term = -inv(QFF)*GFA*SAA*GAF in Q3
	call MATMUL(invQFF,GFA,Q3,kF,kF,kA,one,	  !inv(QFF)*GFA in Q3
     & kFm,kFm,kFm,kAm,km,km)
c	call ATYPD(Q3,'   Q3-1 ',kF,kA,km,km)
	call MATMUL(Q3,SAA,Q3,kF,kA,kA,one,		   !inv(QFF)*GFA*SAA in Q3
     & km,km,kAm,kAm,km,km)
c	call ATYPD(Q3,'   Q3-2 ',kF,kA,km,km)
	call MATMUL(Q3,GAF,Q3,kF,kA,kF,-one,	   !-inv(QFF)*GFA*SAA*GAF in Q3
     & km,km,kAm,kFm,km,km)
c	call ATYPD(Q3,'   Q3-3 ',kF,kF,km,km)
c Add Q1+Q2+Q3; result in Q1
      call MATADD(Q1,Q2,Q1,kF,kF,
     & km,km,km,km,km,km)
      call MATADD(Q1,Q3,Q1,kF,kF,
     & km,km,km,km,km,km)
c	call ATYPD(Q1,'   Q1-3 ',kF,kF,km,km)
c Calc VF=I - GFA*SAA*GAF in Q2
	call MATMUL(GFA,SAA,Q2,kF,kA,kA,one,
     & kFm,kAm,kAm,kAm,km,km)
c	call ATYPD(Q2,'   Q2-4 ',kF,kA,km,km)
	call MATMUL(Q2,GAF,Q2,kF,kA,kF,one,
     & km,km,kAm,kFm,km,km)
c	call ATYPD(Q2,'   Q2-5 ',kF,kF,km,km)
      call MATIMIN(Q2,Q2,kF,km)	!(I-Q2) in Q2
c	call ATYPD(Q2,'   Q2-6 ',kF,kF,km,km)
      call MATINV2(Q2,kF,km,Q2,km,.false.,det,ndscale)	!inv(VA) in Q2
c	call ATYPD(Q2,'   Q2-7 ',kF,kF,km,km)
c Calc DFRS=inv(VF)*(QFF**-2) - inv(VF)*Q1*inv(VF)*inv(QFF)
c = inv(VF)[inv(QFF) - Q1*inv(VF)]*inv(QFF)
	call MATMUL(Q1,Q2,Q3,kF,kF,kF,-one,
     & km,km,km,km,km,km)                       !-Q1*inv(VFF) in Q3
c	call ATYPD(Q3,'   Q3-4 ',kF,kF,km,km)
      call MATADD(invQFF,Q3,Q3,kF,kF,
     & kFm,kFm,km,km,km,km)                     !inv(QFF)-Q1*inv(VFF) in Q3
c	call ATYPD(Q3,'   Q3-5 ',kF,kF,km,km)
	call MATMUL(Q2,Q3,Q3,kF,kF,kF,one,
     & km,km,km,km,km,km)                       !inv(VF)*Q3 in Q3
c	call ATYPD(Q3,'   Q3-6 ',kF,kF,km,km)
	call MATMUL(Q3,invQFF,DFRS,kF,kF,kF,one,
     & km,km,kFm,kFm,kFm,kFm)                   !inv(VF)*Q3*invQFF in DFRS
c	call ATYPD(DARS,'  DFRS= ',kF,kF,kFm,kFm)
c
	DEALLOCATE(Q1,Q2,Q3,GAF,GFA,INVQAA,INVQFF)
	RETURN
	end



