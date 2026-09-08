	subroutine MATMUL3(A,B,R,NA,N,NB,m1,m2,ia1,ia2,ja1,ja2,
     & ib1,ib2,jb1,jb2,km)
c MATMUL3 is modif of matmul, to multiply two 3-D arrays, ie
c to evaluate in the 2-D result R(i,j) the product A(m1)*B(m2) where
c the arrays are A(i,j,m) and B(i,j,m) (01/20/92 04:04pm)
c Fixed (01/23/92 10:37am) to use rows ia1 to ia2, and cols ja1 to ja2 of A,
c and similarly for B.
c NB Output CANNOT be same as input (dimensions different!)
C
	real*8 A(km,km,km),B(km,km,km),R(km,km)
	real*8 ZERO
c
	if((nA.ne.ia2-ia1+1).or.
     & (n.ne.ja2-ja1+1).or.
     & (n.ne.ib2-ib1+1).or.
     & (nB.ne.jb2-jb1+1)) then
	   call BELL(2)
	   print 20
20	   format(' ERROR IN ARGUMENTS OF MATMUL3')
	   STOP
	endif
c
	zero=0.0d0
	do 1 i=1,nA
	   i1=i+ia1-1	!i1=ia1 to ia2
	   do 2 j=1,nB
		j2=j+jb1-1	!j1=jb1 to jb2
		R(i,j)=zero
		do 3 k=1,n
		   k1=k+ja1-1	!k1=ja1 to ja2
		   k2=k+ib1-1	!k2=ib1 to ib2
		   R(i,j)=R(i,j)+A(i1,k1,m1)*B(k2,j2,m2)
3		continue
2	   continue
1	continue
c
	RETURN
	END


