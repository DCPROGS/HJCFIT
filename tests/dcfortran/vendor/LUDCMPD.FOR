      SUBROUTINE LUDCMPd(A,N,NP,INDX,D,nerr)
c Double precision version of Numerical Recipes routine
c
c Modif 05/17/04 08:55pm for better precautions against singular matrix
c  (as of this date, only matinv2.for calls this subroutine)
c Extra argument added,
c  nerr=0	OK
c  nerr=1   singular
c Modified 02/09/06 02:36pm to set TINY=1.d-154 (NB in single precision version
c from book, tiny=1.e-20, NOT 1.e-37)
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      PARAMETER (NMAX=100,TINY=1.0E-20)
c      PARAMETER (NMAX=200,TINY=1.0d-307)
      PARAMETER (NMAX=200,TINY=1.0d-154)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
c
	nerr=0
      D=1.0d0
      DO 12 I=1,N
        AAMAX=0.0d0
        DO 11 J=1,N
          IF (dABS(A(I,J)).GT.AAMAX) AAMAX=dABS(A(I,J))
11      CONTINUE
c        IF (AAMAX.EQ.0.d0) PAUSE 'Singular matrix.'
        IF (dabs(AAMAX).le.tiny) then
c		PAUSE 'Singular matrix.'
		nerr=1
		RETURN
	  endif
        VV(I)=1.0d0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.d0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*dABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1.d0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.d0) A(N,N)=TINY
      RETURN
      END

