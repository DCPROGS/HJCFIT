	subroutine DENOMF(m,s,QFA,QAA,QAF,EXPQA,rowvecF,colvecF,
     & tres,kA,kF,denF,kAx,kAm,kFm,km)
c To calculate normalising denominator for areas in HJC SHUT TIME distributions.
c denA=r(m)*W'(s(m))*c(m)
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	real*8 QFA(10,10),QAF(10,10),QAA(10,10)
c	real*8 EXPQA(10,10)
c	real*8 rowvecF(10,10),colvecF(10,10)
	real*8 QAF(kAx,100),QFA(100,kAx),QAA(kAx,kAx)
	real*8 EXPQA(kAx,kAx)
	real*8 rowvecF(kFm,kFm),colvecF(kFm,kFm)
	allocatable::SAA,Q1,Q2,Q3
	real*8 SAA(:,:),Q1(:,:),Q2(:,:),Q3(:,:)   !local
	real*8 row(1,100),col(100,1)				!local
	real*8 dexp1,det
	logical debug,caplock,discprt
	common/dp/discprt
c
	debug()=caplock()
c
c exp[-(sI-QAA)tres]= exp(-s*tres)*exp(QAA*tres)=exp(-s*tres)*EXPQA
c Calc SAA(s)= I-exp[-(sI-QAA)tres] in this loop, in SAA.
c Also get (sI-QAA) in Q1 in same loop
c Modif 05/02/03 03:31pm for case where kA > kF -in this case
c some of following arrays may not be big enough if all
c declared as (kFm,kFm)
	kmm=kFm
	if(kAm.gt.kFm) kmm=kAm		!max(kAm,kFm)
	one=1.0d0
c	ALLOCATE(SAA(kAm,kAm),Q1(kFm,kFm),Q2(kFm,kFm),Q3(kFm,kFm))
	ALLOCATE(SAA(kAm,kAm),Q1(kmm,kmm),Q2(kmm,kmm),Q3(kmm,kmm))
c
	est=dexp1(-s*tres)
	do 1 i=1,kA
	do 1 j=1,kA
	 if(i.ne.j) then
	   SAA(i,j)=-est*EXPQA(i,j)
	   Q1(i,j)=-QAA(i,j)
	 else
	   SAA(i,j)=one-est*EXPQA(i,j)
	   Q1(i,j)=s-QAA(i,j)
	 endif
1	continue
c
c Calc inv(sI-QAA) in Q1
      call MATINV2(Q1,kA,kmm,Q1,kmm,.false.,det,ndscale)
c
      call MATMUL(SAA,Q1,Q2,kA,kA,kA,one,
     & kAm,kAm,kmm,kmm,kmm,kmm)		!SAA(s)*inv(sI-QAA)
c
	do 2 i=1,kA
	do 2 j=1,kA
	 if(i.ne.j) then
	   Q3(i,j)=-tres*SAA(i,j)
	 else
	   Q3(i,j)=tres*(one-SAA(i,j))	!Q3=tres*(I-SAA(s))
	 endif
2	continue
c subtract
	do 3 i=1,kA
	do 3 j=1,kA
	Q2(i,j)=Q2(i,j)-Q3(i,j)			!Q2=Q2-Q3
3	continue
c
      call MATMUL(QFA,Q2,Q2,kF,kA,kA,one,
     & km,kAx,kmm,kmm,kmm,kmm)			!QFA*ditto in Q2
      call MATMUL(Q2,Q1,Q2,kF,kA,kA,one,
     & kmm,kmm,kmm,kmm,kmm,kmm)                 !QFA*Q2*inv(sI-QFF)
      call MATMUL(Q2,QAF,Q2,kF,kA,kF,one,
     & kmm,kmm,kAx,km,kmm,kmm)			!QFA*Q2*inv(sI-QAA)*QAF in Q2
	do 4 i=1,kF
	do 4 j=1,kF
	 if(i.ne.j) then
	   Q2(i,j)=Q2(i,j)			!I+ditto=W'(s) in Q2
	 else
	   Q2(i,j)=one+Q2(i,j)
	 endif
4	continue
c Q2 now contains W'[s(m)]
	do 5 j=1,kF
5	 row(1,j)=rowvecF(m,j)
	do 6 i=1,kF
6	 col(i,1)=colvecF(i,m)
	if(debug()) then
	   call ATYPD(Q2,'dW(s)/ds',kF,kF,kFm,kFm)
         call VECMUL(row,col,scal,kF,km,km)
         print 7,m,s,scal
         if(discprt) write(8,7) m,s,scal
7	   format(' m, root(m), row*col = ',i4,2g13.6)
	endif
      call MATMUL(Q2,col,col,kF,kF,1,one,kmm,kmm,km,1,km,1)
      call VECMUL(row,col,denF,kF,km,km)	!denF=denominator
	DEALLOCATE(SAA,Q1,Q2,Q3)
	RETURN
	end

