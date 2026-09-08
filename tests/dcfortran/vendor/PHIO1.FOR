	subroutine PHIo1(QM,Peq,phi,km)
c Initial vector for openings
c Part of family of subroutines to do single channel calcs via
c set of discrete subroutines that can be called as req.
c   There is a problem in designating submatrices in that notation
c is somewat different in burst and cluster programs, eg F or T for all
c shut states. Here just work directly from kA,kB,kC,kD
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	dimension QM(100,100),Q1(100,100),Peq(100),phi(1,100)
	integer TA
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
c
	kT=kB+kC+kD
c	k=kA+kB+kC+kD
	TA=71
	zero=0.0d0
	one=1.0d0
	call SUBMAT(QM,TA,Q1,km,km,km,km)		!QTA in Q1
	SUM=ZERO
	DO 22 J=1,KA		!PREMULT BY pT(inf) TO GET PHI(o)
	phi(1,j)=ZERO
	DO 23 L=1,kT
23	phi(1,j)=phi(1,j)+Peq(L+kA)*Q1(L,j)
22	SUM=SUM+phi(1,j)
c
c During fitting can get 0/0 here so if this occurs arbitrarily set phi(1)=1
	if(sum.gt.1.0d-35) then
	   do 9 j=1,kA
	   phi(1,j)=phi(1,j)/sum	!NORMALISE PHI
9	   continue
	else
	   phi(1,1)=one
	   if(kA.gt.1) then
		do 10 j=2,kA
10		phi(1,j)=zero
	   endif
	endif
c
	RETURN
	end

