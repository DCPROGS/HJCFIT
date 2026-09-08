	subroutine GMAT1(AB,QM,GAB,km,kr)
C GMAT1 is like GMAT but superfluous array dimensions removed from
c call as in GMAT2 (from which it differs only in that inv(QAA) not returned)
c
c Modif 07/21/99 07:04am to allow any array size
c Modif 12/19/91 06:59pm so get print out of determinant AND 'AB' if determ
c is small in MATINV
c
C SUBROUTINE TO CALCULATE PROB TRANSITION MATRICES ,GAB ETC, FROM Q
C	GAB=-INV(QAA)*QAB (UP TO 8 X 8 IN SIZE)
C	AB=INTEGER CODE FOR REQ MATRIX
C	QM=DOUBLE PREC Q
C	km=DIM OF QM IN MAIN PROG (SQUARE)
C	KR=DIM OF GAB  IN MAIN PROG (SQUARE)
C
	allocatable::temp
	real*8 temp(:,:)
	real*8 QM(km,km),GAB(kr,kr)
	real*8 one
	real*8 det
	character ch1*1,ch2*2
	logical discprt
	INTEGER AB
	COMMON/KBLK/KA,KB,KC,KD
	COMMON/determ/det			!from MATINV
	common/dp/discprt
C
	kdim=km
	allocate(temp(kdim,kdim))
	one=1.0D0
C FIRST SEPARATE THE DIGITS OF 'AB' AS IN SUBMAT
	I=AB/10	!LEFT HAND DIGIT
	J=MOD(AB,10)	!RIGHT HAND DIGIT
C RESULT WILL BE OF SIZE (K1 X K2) WHERE
	IF(I.EQ.1) K1=KA
	IF(I.EQ.2) K1=KB
	IF(I.EQ.3) K1=KC
	IF(I.EQ.4) K1=KD
	IF(I.EQ.5) K1=KB+KC	!=F
	IF(I.EQ.6) K1=KA+KB	!=E
	IF(I.EQ.7) K1=KB+KC+KD	!=T
	if(i.eq.8) k1=kC+kD
C
	IF(J.EQ.1) K2=KA
	IF(J.EQ.2) K2=KB
	IF(J.EQ.3) K2=KC
	IF(J.EQ.4) K2=KD
	IF(J.EQ.5) K2=KB+KC
	IF(J.EQ.6) K2=KA+KB
	IF(J.EQ.7) K2=KB+KC+KD
	if(j.eq.8) k2=kC+kD
C
C GET SUBSCRIPT FOR QAA ETC
	II=I*10+I
C GET QII IN GAB
	call SUBMAT(QM,II,GAB,km,km,kr,kr)
C INVERT IT GIVING INV(QII) IN GAB.  Q11 IS SIZE (K1 X K1)
	call MATINV(GAB,K1,kr,GAB,kr)
	if(dabs(DET).lt.1.d-10) then
	   ch1=char(i+64)       !1-4 =A-D
	   ch2=char(j+64)
	   if(i.eq.5) ch1='F'
	   if(j.eq.5) ch2='F'
	   if(i.eq.6) ch1='E'
	   if(j.eq.6) ch2='E'
	   if(i.eq.7) ch1='T'
	   if(j.eq.7) ch2='T'
	   print 2,ch1//ch2,det
	   if(discprt) write(8,2) ch1//ch2,det
2	   format(' !! In calculation of G',a2,' find determinant = '
     &	,g13.6,/)
	endif
C GET QAB IN TEMP. ITS SIZE IS (K1 X K2)
	call SUBMAT(QM,AB,TEMP,km,km,kdim,kdim)
C GET -QII*TEMP=GAB IN GAB
	call MATMUL(GAB,TEMP,GAB,K1,K1,K2,-ONE,
     & kr,kr,kdim,kdim,kr,kr)
c
	DEALLOCATE(temp)
	RETURN
	END

