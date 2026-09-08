	subroutine MATADD(A,B,R,ki,kj,
     & KA1,KA2,KB1,KB2,KR1,KR2)
c Matrix addition R = A + B (double prec,dimension)
C R can be same as A or B
	real*8 A(KA1,KA2),B(KB1,KB2),R(KR1,KR2)
c
	do 1 i=1,ki
	do 2 j=1,kj
	R(i,j)=A(i,j) + B(i,j)
2	continue
1	continue
	RETURN
	END

