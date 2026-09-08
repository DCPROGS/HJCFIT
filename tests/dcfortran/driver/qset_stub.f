c QSET_HJC replacement: hand the likelihood a Q matrix directly.
c
c The real QSET_HJC builds Q from theta through HJCFIT's whole topology and
c constraint machinery -- irate/jrate, the EBLK equality constraints, the
c micro-reversibility cycles, the EC50 constraint, the concentration
c insertion. None of that is what we are testing: the parameterisation and
c the constraints have already been eliminated as explanations, by paired
c refits that gave identical answers. What is left to test is the
c **likelihood given Q**.
c
c So this stub ignores theta entirely and copies a Q matrix supplied through
c common/qdgiven/, with the concentrations already in it, exactly as SCALCS
c hands the same matrix to the C++ Log10Likelihood. It then sets the
c diagonals the way SETDIAG does in the original (QSET_HJC.FOR, last lines),
c so the two implementations see the same generator.
c
c Call HJCLIK with kfit=0 and the theta loop above this point does nothing.

      subroutine QSET_HJC(jset,theta,QT,QD,kfit,k)
      implicit none
      integer jset,kfit,k
      real*8 theta(200),QT(100,100),QD(100,100)

      integer kg
      real*8 QGIVEN(100,100)
      common/qdgiven/QGIVEN,kg

      integer i,j
      real*8 s

      do i=1,kg
         do j=1,kg
            QD(i,j)=QGIVEN(i,j)
            QT(i,j)=QGIVEN(i,j)
         enddo
      enddo
c diagonals, as SETDIAG does
      do i=1,kg
         s=0.0d0
         do j=1,kg
            if(j.ne.i) s=s-QD(i,j)
         enddo
         QD(i,i)=s
         QT(i,i)=0.0d0
      enddo
      k=kg
      return
      end
