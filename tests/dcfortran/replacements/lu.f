c LUDCMPd / lubksbd -- LU decomposition with partial pivoting, and the
c corresponding solve.
c
c DCFORTRAN's LUDCMPD.FOR and LUBKSBD.FOR say in their own comments that they
c are real*8 versions of the Numerical Recipes routines, whose licence does
c not allow redistribution, so they are not vendored here
c (tests/dcfortran/README.md). These are written from the definition of
c Gaussian elimination with partial pivoting -- the right-looking form LAPACK
c calls dgetf2 -- and keep the interface MATINV2.FOR and DETERM2 use.
c
c One deliberate difference. The Numerical Recipes routine chooses its pivot
c by *implicit scaling*: it divides each candidate by the largest element of
c its row before comparing. This one compares the candidates themselves, which
c is ordinary partial pivoting. Both are backward stable and both give an
c exact LU of a permutation of the matrix; on a badly scaled matrix they can
c choose different pivots and so round differently in the last bits. That is
c the one place where this engine can disagree with the original, and it is
c measured rather than assumed -- see test_engines_agree.py.

c ----------------------------------------------------------------------
c LU decomposition, in place.
c
c   A(NP,NP)  on entry the N x N matrix; on exit U on and above the diagonal
c             and the multipliers of L below it, L having a unit diagonal
c   INDX(N)   INDX(j) is the row swapped into position j at step j
c   D         +1 or -1, the sign of the permutation, so that the determinant
c             is D times the product of the diagonal of U
c   nerr      0, or 1 if the matrix is singular

      SUBROUTINE LUDCMPd(A,N,NP,INDX,D,nerr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N)
c As in the original: the smallest pivot that is treated as nonzero. Well
c below anything a Q matrix produces, so this is a guard, not a threshold.
      PARAMETER (TINY=1.0d-154)

      nerr=0
      D=1.0d0
      do j=1,N
c the pivot: the largest remaining element of column j
         imax=j
         big=dabs(A(j,j))
         do i=j+1,N
            if(dabs(A(i,j)).gt.big) then
               big=dabs(A(i,j))
               imax=i
            endif
         enddo
         if(big.le.TINY) then
            nerr=1
            RETURN
         endif
         INDX(j)=imax
         if(imax.ne.j) then
            do k=1,N
               dum=A(imax,k)
               A(imax,k)=A(j,k)
               A(j,k)=dum
            enddo
            D=-D
         endif
c eliminate below it
         do i=j+1,N
            A(i,j)=A(i,j)/A(j,j)
            do k=j+1,N
               A(i,k)=A(i,k)-A(i,j)*A(j,k)
            enddo
         enddo
      enddo
      RETURN
      END

c ----------------------------------------------------------------------
c Solve A x = b, given the factors LUDCMPd left in A. b is overwritten with
c the solution. MATINV2 calls this once per column of the identity to invert
c a matrix.

      SUBROUTINE lubksbd(A,N,NP,INDX,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N),B(N)

c forward substitution, undoing the row swaps as it goes
      do i=1,N
         ll=INDX(i)
         s=B(ll)
         B(ll)=B(i)
         do j=1,i-1
            s=s-A(i,j)*B(j)
         enddo
         B(i)=s
      enddo

c back substitution
      do i=N,1,-1
         s=B(i)
         do j=i+1,N
            s=s-A(i,j)*B(j)
         enddo
         B(i)=s/A(i,i)
      enddo
      RETURN
      END
