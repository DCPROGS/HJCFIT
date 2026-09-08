	subroutine MATTRAN1(Q,QT,k1,k2,km1,km2)
c Get transpose of Q (top left k1*k2) in QT (which is thus k2*k1)
	allocatable::Q0
	real*8 Q0(:,:)
	real*8 Q(km1,km1),QT(km2,km2)
c
c So it will work if Q and QT are same in the call, get result in Q0 at first
	km=max(km1,km2)
	ALLOCATE(Q0(km,km))
	do 1 i=1,k1
	do 1 j=1,k2
1	Q0(i,j)=Q(j,i)
c
	do 2 i=1,k1
	do 2 j=1,k2
2	QT(i,j)=Q0(i,j)

	DEALLOCATE(Q0)
	RETURN
	end
