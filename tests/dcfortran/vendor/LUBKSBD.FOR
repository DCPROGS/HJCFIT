      SUBROUTINE lubksbd(a,n,np,indx,b)
c real*8 version of Numerical Recipes LUBKSB.for
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER n,np,indx(n)
      real*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      real*8 sum
c
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #1,15(15..
