	subroutine MATIMIN(Q,Q1,kA,kdim)
c Calculates Q1=(I-Q) for the top left kA x kA part
c Q and Q1 can be the same
	real*8 Q(kdim,kdim),Q1(kdim,kdim)
c
	do 1 i=1,kA
	do 1 j=1,kA
	Q1(i,j)=-Q(i,j)
1	if(i.eq.j) Q1(i,j)=1.0d0 + Q1(i,j)
c
	RETURN
	end



