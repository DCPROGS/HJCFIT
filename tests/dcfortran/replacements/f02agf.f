c F02AGF -- eigenvalues and eigenvectors of a real unsymmetric matrix,
c with NAG's Mark 2 interface, computed by EISPACK.
c
c DCFORTRAN's own F02AGF is NAG's, and carries NAG's copyright, so it is not
c vendored here (tests/dcfortran/README.md). This provides the same entry
c point over the public-domain EISPACK routines in eispack/, which compute
c the same thing by the same algorithm: NAG's Mark 2 chain and EISPACK are
c both Fortran translations of the same ALGOL procedures from Wilkinson &
c Reinsch, Handbook for Automatic Computation vol II (1971) -- NAG's own
c comment lines name them, DIRHES, DIRTRANS, HQR2 and CDIV, which are
c elmhes, eltran, hqr2 and cdiv here.
c
c Interface, unchanged from NAG's, since QMAT5.FOR calls it:
c
c   A(IA,N)     the matrix; destroyed (QMAT5 passes a copy)
c   RR,RI       real and imaginary parts of the eigenvalues
c   VR(IVR,N)   real part of the eigenvectors, one per column
c   VI(IVI,N)   imaginary part
c   INT(N)      workspace
c   IFAIL       0 on success
c
c On scaling: QMAT5.FOR rescales every column so its top element is 1, and
c then forms the spectral matrices as A(m) = EM(:,m) EN(m,:) with EN =
c inv(EM). Scaling a column of EM scales the corresponding row of its inverse
c by the reciprocal, so the spectral matrices -- and therefore the likelihood
c -- do not depend on how the eigenvectors are normalised here at all. The
c normalisation below is for conditioning, not for agreement.

      SUBROUTINE F02AGF(A,IA,N,RR,RI,VR,IVR,VI,IVI,INT,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,N),VR(IVR,N),VI(IVI,N),RR(N),RI(N),INT(N)
      allocatable::H,Z
      real*8 H(:,:),Z(:,:)

c A local copy at leading dimension N, because EISPACK takes one leading
c dimension for both the matrix and the eigenvectors where NAG takes IA, IVR
c and IVI separately.
      allocate(H(N,N),Z(N,N))
      do j=1,N
         do i=1,N
            H(i,j)=A(i,j)
         enddo
      enddo

c NAG's F02AGF calls its chain with low=1 and igh=N, i.e. without balancing,
c so this does too.
      call elmhes(N,N,1,N,H,INT)
      call eltran(N,N,1,N,H,INT,Z)
      call hqr2(N,N,1,N,H,RR,RI,Z,ierr)

      IFAIL=ierr
      if(ierr.ne.0) then
         deallocate(H,Z)
         RETURN
      endif

c hqr2 returns the eigenvectors packed: a real eigenvalue's vector is its own
c column of Z; a complex pair, which hqr2 returns with the positive imaginary
c part first, has its real part in column i and its imaginary part in column
c i+1. Unpack that into the separate VR and VI the NAG interface promises,
c the second member of each pair being the conjugate of the first.
      i=1
      do while(i.le.N)
         if(RI(i).eq.0.0d0) then
            do j=1,N
               VR(j,i)=Z(j,i)
               VI(j,i)=0.0d0
            enddo
            i=i+1
         else
            do j=1,N
               VR(j,i)=Z(j,i)
               VI(j,i)=Z(j,i+1)
               VR(j,i+1)=Z(j,i)
               VI(j,i+1)=-Z(j,i+1)
            enddo
            i=i+2
         endif
      enddo

c Scale each eigenvector by its largest modulus. See the note above: this
c cancels out of everything QMAT5 goes on to compute.
      do i=1,N
         big=0.0d0
         do j=1,N
            t=dsqrt(VR(j,i)**2+VI(j,i)**2)
            if(t.gt.big) big=t
         enddo
         if(big.gt.0.0d0) then
            do j=1,N
               VR(j,i)=VR(j,i)/big
               VI(j,i)=VI(j,i)/big
            enddo
         endif
      enddo

      deallocate(H,Z)
      RETURN
      END
