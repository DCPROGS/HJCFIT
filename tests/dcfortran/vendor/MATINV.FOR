	subroutine MATINV(ARRAY,K,KD1,RECIP,KD2)
c Modified 10/04/93 08:26pm to making lwork,kwork allocatable
c so no restrictions on array size
c 06/17/90 08:26pm Modified by inclusion of COMMON/determ/det
c so determinant can be accessed by progs that declare this common
C DOUBLE PRECISION SUBROUTINE TO INVERT MATRIX WITHOUT DESTROYING
C IT (THOUGH ARRAY AND RECIP CAN BE SAME IF DESIRED).
C***MODIF (1-DEC-81) SO MINV NOT CALLED IF K=1.
C I/O IN USUAL 2-DIMENSIONAL FORM.
C	ARRAY=INPUT MATRIX
C	K=DIMENSION - INVERTS TOP LEFT K X K ELEMENTS
C	KD1=DECLARED DIMENSION OF ARRAY IN CALLING PROG
C	KD2=DECLARED DIMENSION OF RECIP IN CALLING PROG
C	RECIP=OUTPUT-INVERSE OF ARRAY
C
	real*8 ARRAY(KD1,KD1),RECIP(KD2,KD2),DET
c	integer LWORK(20),MWORK(20)
	integer LWORK(:),MWORK(:)
	allocatable:: Lwork,mwork
	COMMON/determ/det
C
	ALLOCATE(lwork(k),mwork(k))
	if(k.eq.1) then
	   recip(1,1)=1.0d0/array(1,1)
	   det=array(1,1)
	   RETURN
	endif
C COPY INPUT
	DO 1 I=1,K
	DO 1 J=1,K
1	RECIP(I,J)=ARRAY(I,J)
C CONVERT TO 1-D FORM
	call ARRAYD(2,K,K,KD2,KD2,RECIP,RECIP)
C INVERT ARRAY GIVING RECIP
	call MINVD(RECIP,K,DET,LWORK,MWORK,ndexp)
c Check for overflow (but can take log(det) only if non-zero)
	if(dabs(det).gt.1.d-300) then
	 if(dfloat(ndexp).le.(308.d0-dlog10(dabs(det)))) then
	   det=det*(10.d0**ndexp)		!check for overflow before mult?
	 else
	   print 3
3	   format(' *** Determinant > 10^308 in MATINV',/)
	  det=1.d308
	 endif
	endif
	if(dabs(det).lt.1.d-10) then
	   print 2,det
2	   FORMAT( ' ***Determinant = ',g13.6,' in MATINV',/)
	endif
C CONVERT BACK TO TWO DIMENSIONAL FORM
	call ARRAYD(1,K,K,KD2,KD2,RECIP,RECIP)
c
	DEALLOCATE(lwork,mwork)
	RETURN
	END


