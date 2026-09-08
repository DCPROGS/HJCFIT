	subroutine SUBMAT3d(AMAT,m,AB,QAB,ka1,ka2,ka3,ks1,ks2)

c SUBMAT3d is version of SUBMAT3 allows all array dimensions to be specified
c SUBMAT3 is version of SUBMAT that produces the AB section of the 3D matrix
c AMAT(m) in top left corner of the 2D output matrix, QAB.
c NB Either F (=5) or T (=7) should define all shut states identically
c if kD=0 (burst progs), though only latter is OK if kD not zero (clusters).
c *** Sep-89: subset (C+D) added (#8) for use  for runs of single bursts
c in the 2 channel paper [for calcs, in BCHAN2.FOR, the sets have
c different letters from those given below so integers FG etc redefined
c appropriately in the calling prog: E=(B+C), F=(C+D), G=(B+C+D)].
C DOUBLE PRECISION SUBROUTINE TO EXTRACT A SUBMATRIX,QAB,
C FROM THE MATRIX QM.
C	KQ1,KQ2=DECLARED DIM OF QM IN MAIN PROG
C	KS1,KS2=DECLARED DIM OF QAB IN MAIN PROG
C	AB=2 DIGIT INTEGER (11,12 ETC) SUCH THAT 1ST
C		DIGIT=CODE FOR ROWS, AND 2ND=CODE FOR COLS THAT ARE
C		TO BE INCLUDED IN THE SUBMATRIX
C 1,2,3,4= CODES FOR A,B,C,D RESP. 5=F=(B+C). 6=E=(A+B). 7=T=(B+C+D).
c 8=H=(C+D), 9=G=(A+B+C)
C
c	real*8 amat(10,10,10),QAB(10,10),zero
	real*8 amat(ka1,ka2,ka3),QAB(ks1,ks2),zero
	integer AB
	COMMON/KBLK/KA,KB,KC,KD
C
C FIRST ZERO QAB SO ALL BUT REQ ELEMENTS ARE ZERO
	ZERO=0.0D0
	DO i=1,ks1
	   DO j=1,ks2
		QAB(i,j)=ZERO
	   enddo
	enddo
C FIRST SEPARATE THE DIGITS OF AB
	NR=AB/10		!LEFT HAND DIGIT=ROW
	NC=MOD(AB,10)		!RIGHT HAND DIGIT=COL
c	print 41,nc,nr
c41	format(' nc,nr= ',2i4)
C
	kE=kA+kB
	kG=kA+kB+kC
	k=kA+kB+kC+kD
C NR1 TO NR2=ROW NUMBERS FOR SUBMATRIX
C NC1 TO NC2=COL NUMBERS FOR SUBMATRIX
C
	GOTO (1,2,3,4,5,6,7,8,9),NR		!ASSIGN ROWS
1	NR1=1		!A
	NR2=KA
	GOTO 100
2	NR1=KA+1	!B
	NR2=KE
	GOTO 100
3	NR1=KE+1	!C
	NR2=KG
	GOTO 100
4	NR1=KG+1	!D
	NR2=K
	GOTO 100
5	NR1=KA+1	!=F SECTION=(B+C) (=E for BCHAN2)
	NR2=KG
	GOTO 100
6	NR1=1		!=E SECTION=(A+B)
	NR2=KE
	GOTO 100
7	NR1=KA+1	!=T SECTION=(B+C+D) (=G for BCHAN2)
	NR2=K
	GOTO 100
8	nr1=ka+kb+1	!SECTION=H=(C+D) (=G for BCHAN2)
	nr2=ka+kb+kc+kd
	goto 100
9	nr1=1
	nr2=kG	!section G=A+B+C in C&H 1982
	goto 100
C
100	GOTO (101,102,103,104,105,106,107,108,109),NC	!ASSIGN COLS
101	NC1=1		!A
	NC2=KA
	GOTO 200
102	NC1=KA+1	!B
	NC2=KE
	GOTO 200
103	NC1=KE+1	!C
	NC2=KG
	GOTO 200
104	NC1=KG+1	!D
	NC2=K
	GOTO 200
105	NC1=KA+1	!=F SECTION=(B+C)  (=E for BCHAN2)
	NC2=KG
	GOTO 200
106	NC1=1		!=E SECTION=(A+B)
	NC2=KE
	GOTO 200
107	NC1=KA+1	!=T SECTION=(B+C+D) (=G for BCHAN2)
	NC2=K
	GOTO 200
108	nc1=ka+kb+1	!=H SECTION=(C+D) (=F for BCHAN2)
	nc2=ka+kb+kc+kd
	goto 200
109	nc1=1
	nc2=kG	!section G=A+B+C in C&H 1982
c	goto 200
C
200	do 300 i=nr1,nr2
	do 300 j=nc1,nc2
300	QAB(i-nr1+1,j-nc1+1)=amat(i,j,m)
	RETURN
	END


