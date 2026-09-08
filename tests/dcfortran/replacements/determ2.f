c DETERM2 -- determinant of a matrix by LU decomposition, with the result
c scaled by powers of 10^10 so that it cannot overflow or underflow.
c
c DCFORTRAN's determ2.for holds two routines: DETERM2 itself, which is
c DCPROGS' own, and an appended copy of the Numerical Recipes LUDCMP, whose
c licence does not allow redistribution. Because the two share a file, the
c file is not vendored (tests/dcfortran/README.md) and DETERM2 is written out
c again here, calling the LUDCMPd in lu.f. Its behaviour is the original's:
c the same decade scaling, the same ndscale convention, the same treatment of
c the 1 x 1 case.
c
c   A(ndim,ndim)  the n x n matrix; not destroyed
c   det, ndscale  determinant = det * 10^(10*ndscale)
c
c DETWA.FOR and DETWF.FOR use it for the determinant of W(s), whose roots are
c the asymptotic time constants, so the scaling matters: those determinants
c run over many orders of magnitude.

      subroutine DETERM2(A,n,ndim,det,ndscale)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 A(ndim,ndim)
      allocatable::A1,indx
      real*8 A1(:,:)
      integer indx(:)

      if(n.eq.1) then
         det=A(1,1)
         ndscale=0
         RETURN
      endif

c a copy, so the caller's matrix survives
      allocate(A1(n,n),indx(n))
      do i=1,n
         do j=1,n
            A1(i,j)=A(i,j)
         enddo
      enddo

      ndscale=0
      call LUDCMPd(A1,n,n,indx,D,nerr)
      if(nerr.ne.0) then
         det=0.0d0
         deallocate(A1,indx)
         RETURN
      endif

c The determinant is D times the product of the diagonal of U. Multiply it up
c one factor at a time, taking out a decade whenever the next multiplication
c would leave the range of a real*8, and counting the decades in ndscale.
      do i=1,n
1        di=dlog(dabs(D))
         ai=dlog(dabs(A1(i,i)))
         if(di+ai.gt.308.d0) then
            ndscale=ndscale+1
            D=D*1.d-10
            goto 1
         endif
2        di=dlog(dabs(D))
         ai=dlog(dabs(A1(i,i)))
         if(di+ai.lt.-308.d0) then
            ndscale=ndscale-1
            D=D*1.d10
            goto 2
         endif
         D=D*A1(i,i)
      enddo

      det=D
      deallocate(A1,indx)
      return
      end
