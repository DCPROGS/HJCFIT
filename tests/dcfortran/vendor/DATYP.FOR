	SUBROUTINE DATYP(ARRAY,title,print1,KR,KC,KD1,KD2)
c
c 06/21/89 10:35am Lahey version; up to npr=6 cols before wrapping
c
C DATYPv is version of DATYP2 for dimension arrays, and 8 character title.
C VERSION THAT PRINTS EACH ROW ON 2 LINES IF MORE THAN 5 COLS
C**SIMPLIFIED VERSION OF DATYPE WITH NO FACILITY TO RENUMBER
C DOUBLE PRECISION VERSION OF ATYPE (CALL ALTERED SO KD1,KD2 SPECIFIED)
C SUBROUTINE TO TYPE SUBSET OF AN ARRAY.
C	TYPES ROW NUMBERS 1 TO KR.
C	TYPES COL NUMBERS 1 TO KC.
C	KD1,KD2=DIMENSIONS OF ARRAY DECLARED IN CALLING PROG
c	Title (4 char) also typed first
c	Prints too if PRINT=.true.
C
	real*8 ARRAY(KD1,KD2)
	character*8 title
	logical print1
	logical discprt
	common/dp/discprt
C
	IF(KR.LE.KD1.AND.KC.LE.KD2) GOTO 12
	print 11
11	FORMAT( ' ***ARRAY TOO BIG TO TYPE*** ',/)
	RETURN
12	CONTINUE
c
c Define number of values per row (was 5, now 6) (but still needs to
c be set numerically in format statements at present)
	npr=6
c
c Type title
	if(title.ne.'        ') then
	   print 20, title
	   if(print1) write(7,20)title
         if(discprt) write(8,20)title
20	   format(/,1x,A8,' =')
	endif
C TYPE COL NUMBERS
	IF(KC.LE.npr) print 100,(J,J=1,KC)
	IF(KC.GT.npr) print 101,(J,J=1,KC)
	if(.not.print1) goto 21
	IF(KC.LE.npr) write(7,100)(J,J=1,KC)
	IF(KC.GT.npr) write(7,101)(J,J=1,KC)
100	FORMAT(' J=',2X,6(I2,11X))
101	FORMAT(' J=',2X,6(I2,11X),/,5X,4(I2,11X))
21	continue
	if(.not.discprt) goto 211
	IF(KC.LE.npr) write(8,100)(J,J=1,KC)
	IF(KC.GT.npr) write(8,101)(J,J=1,KC)
211	continue
c
	DO 1 I=1,KR
	print 110,I
	IF(KC.LE.npr) print 5,(ARRAY(I,J),J=1,KC)
	IF(KC.GT.npr) print 6,(ARRAY(I,J),J=1,KC)
	if(.not.print1) goto 3
	write(7,110)I
	IF(KC.LE.npr) write(7,5)(ARRAY(I,J),J=1,KC)
	IF(KC.GT.npr) write(7,6)(ARRAY(I,J),J=1,KC)
3	if(.not.discprt) goto 1
	write(8,110)I
	IF(KC.LE.npr) write(8,5)(ARRAY(I,J),J=1,KC)
	IF(KC.GT.npr) write(8,6)(ARRAY(I,J),J=1,KC)
1	CONTINUE		!end of loop
C
110	FORMAT(I3)
5	FORMAT(1X,6G13.6)
6	FORMAT(1X,6G13.6,/,3X,4G13.6)
C
	RETURN
	END

