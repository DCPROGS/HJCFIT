	subroutine HJCphi(eGAF,eGFA,kA,kF,phiA,kAm,kFm,km)
c
c To solve phi(I-eGAF*eGFA)=0
c
c=	real*8 eGAF(10,10),eGFA(10,10),phiA(1,10),pdum(10)
	real*8 eGAF(kAM,kFm),eGFA(kFm,kAm),phiA(1,100),pdum(100)
	real*8 Q3(100,100)
c
	if(kA.eq.1) then
	   phiA(1,1)=1.d0
	   RETURN
	endif
c
c eGAF*eGFA in Q3
      call MATMUL(eGAF,eGFA,Q3,kA,kF,kA,1.0d0,
     & kAm,kFm,kFm,kAm,km,km)
c Subtract from unit
	do i=1,kA
	   do j=1,kA
		if(i.eq.j) then
		   Q3(i,j)=1.0d0-Q3(i,j)
		else
		   Q3(i,j)=-Q3(i,j)
		endif
	   enddo
	enddo
c Solve for phi (Q3 is kAxkA)
c	call EQOCC1(Q3,pdum,kA,km)
	call EQOC_HJC(Q3,pdum,kA,km,km)
	do j=1,kA
	   phiA(1,j)=pdum(j)
	enddo
c
	RETURN
	end
