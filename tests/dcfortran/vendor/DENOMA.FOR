	subroutine DENOMA(m,s,QAF,QFF,QFA,EXPQF,rowvecA,colvecA,
     & tres,kA,kF,denA,kAx,kAm,kFm,km)
c To calculate normalising denominator for areas in HJC OPEN TIME distributions.
c denA=r(m)*W'(s(m))*c(m)
c 09/15/01 06:16pm Extra array dimensions added to call for 100 state version
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	real*8 QAF(10,10),QFA(10,10),QFF(10,10)
c	real*8 EXPQF(10,10)
c 	real*8 rowvecA(10,10),colvecA(10,10)
	real*8 QAF(kAx,100),QFA(100,kAx),QFF(100,100)
	real*8 EXPQF(100,100)
	real*8 rowvecA(kAm,kAm),colvecA(kAm,kAm)
	allocatable::SFF,Q1,Q2,Q3
	real*8 SFF(:,:),Q1(:,:),Q2(:,:),Q3(:,:)   !local
	real*8 row(1,100),col(100,1)				!local
	real*8 dexp1,det
	logical debug,caplock,discprt
	common/dp/discprt
c
	debug()=caplock()
c
c exp[-(sI-QFF)tres]= exp(-s*tres)*exp(QFF*tres)=exp(-s*tres)*EXPQF
c Calc SFF(s)= I-exp[-(sI-QFF)tres] in this loop, in SFF.
c Also get (sI-QFF) in Q1 in same loop
c Modif 05/02/03 03:31pm for case where kA > kF -in this case
c some of following arrays may not be big enough if all
c declared as (kFm,kFm)
	kmm=kFm
	if(kAm.gt.kFm) kmm=kAm		!max(kAm,kFm)
	one=1.0d0
	ALLOCATE(SFF(kFm,kFm),Q1(kmm,kmm),Q2(kmm,kmm),Q3(kmm,kmm))
c
	est=dexp1(-s*tres)
	do 1 i=1,kF
	do 1 j=1,kF
	 if(i.ne.j) then
	   SFF(i,j)=-est*EXPQF(i,j)
	   Q1(i,j)=-QFF(i,j)
	 else
	   SFF(i,j)=one-est*EXPQF(i,j)
	   Q1(i,j)=s-QFF(i,j)
	 endif
1	continue
c
c Calc inv(sI-QFF) in Q1
      call MATINV2(Q1,kF,kmm,Q1,kmm,.false.,det,ndscale)
c
      call MATMUL(SFF,Q1,Q2,kF,kF,kF,one,
     & kFm,kFm,kmm,kmm,kmm,kmm)		!SFF(s)*inv(sI-QFF)
c
	do 2 i=1,kF
	do 2 j=1,kF
	 if(i.ne.j) then
	   Q3(i,j)=-tres*SFF(i,j)
	 else
	   Q3(i,j)=tres*(one-SFF(i,j))	!Q3=tres*(I-SFF(s))
	 endif
2	continue
c subtract
	do 3 i=1,kF
	do 3 j=1,kF
	Q2(i,j)=Q2(i,j)-Q3(i,j)			!Q2=Q2-Q3
3	continue
c
      call MATMUL(QAF,Q2,Q2,kA,kF,kF,one,
     & kAx,km,kmm,kmm,kmm,kmm)			!QAF*ditto in Q2
      call MATMUL(Q2,Q1,Q2,kA,kF,kF,one,
     & kmm,kmm,kmm,kmm,kmm,kmm)                 !QAF*Q2*inv(sI-QFF)
      call MATMUL(Q2,QFA,Q2,kA,kF,kA,one,
     & kmm,kmm,km,kAx,kmm,kmm)			!QAF*Q2*inv(sI-QFF)*QFA in Q2
	do 4 i=1,kA
	do 4 j=1,kA
	 if(i.ne.j) then
	   Q2(i,j)=Q2(i,j)			!I+ditto=W'(s) in Q2
	 else
	   Q2(i,j)=one+Q2(i,j)
	 endif
4	continue
c Q2 now contains W'[s(m)]
	do 5 j=1,kA
5	 row(1,j)=rowvecA(m,j)
	do 6 i=1,kA
6	 col(i,1)=colvecA(i,m)
	if(debug()) then
	   call ATYPD(Q2,'dW(s)/ds' ,kA,kA,kFm,kFm)
         call VECMUL(row,col,scal,kA,km,km)
         print 7,m,s,scal
         if(discprt) write(8,7) m,s,scal
7	   format(' m, root(m), row*col = ',i4,2g13.6)
	endif
      call MATMUL(Q2,col,col,kA,kA,1,one,kmm,kmm,km,1,km,1)
      call VECMUL(row,col,denA,kA,km,km)	!denA=denominator
	DEALLOCATE(SFF,Q1,Q2,Q3)
	RETURN
	end

