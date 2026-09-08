	subroutine ROOT_FB(sa,sb,GfuncX,kX,s1,s2,nroot,kXm,detwx,nerr)
c To use Frank Ball's method to find suitable starting guesses for each root,
c returned in s1(m),s2(m), the upper and lower limits for BISECTion,
c which is done after this call.  Each pair, s1(i),s2(i) should contain exactly
c one root.
c
c This routine replaces DC's original ROOTSCH, and is based on
c Ball's matlab function resolve.m
c
c The input is an interval sa, sb and uses repeated calls to SPLIT to find
c nroot=kX disjoint subintervals each of which contains exactly one HJC root
c (used as initial guesses for bisection of roots)
c
c Modified 08/04/03 06:54am to cope with errors
c (1) Check that initial interval wide enough, and if not widen it
c (2) Need to cope with roots that are multiple (identical within rounding
c	error of calculations) -presumably this will mean that can't get down
c 	to one root ber interval, however narrow the interval is made -need a
c	minimum interval, such that anything smaller runs into numerical problems
c
c ng1 and ng2 arrays not used or returned so remove them
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 s1(100),s2(100)	!lower and upper guesses for subsequent bisection
	logical debprt,dprt
	common/detw2/tres,km,ndum			!for HJCLIK,DETWA,DETWF
	common/deb1/debprt
	logical discprt
	common/dp/discprt
c internal
c	integer ng1(100),ng2(100)	!corresponding g(s), internal array
	real*8 sv1(100),sv2(100)
	integer nv1(100),nv2(100)
c
	External GfuncX
c For check only
	logical check
	real*8 detwX
	common/det4/ndscale,idsign
	EXTERNAL detwx
c
	common/timer/idebugt,itick,itlast,ndeb
c
	check=.true.
	dprt=discprt.and.debprt	 !put problems in root-finding on disc as well as screen
	nerr=0
	if(idebugt.eq.1) then
	   call TIMER(itickr)
	endif

c Ball's resolve.m had nga and ngb as input arguments, but easier here to
c calculate them here
c NB smin=1.d-20 can give error 4 in GFUNCF/A which have
c	s1=-s*tres
c	if(dabs(s1).gt.709.d0) then
c	   nerr=4
c	   RETURN
c	endif
c so set smin=709/tres to prevent this happening
c (tres was not available in ROOT_FB so added (04/23/04 10:54am)
c	common/detw2/tres,km,nerr			!for HJCLIK,DETWA,DETWF
c This limits fastest asymptotic root to 0.14 microseconds at tres=100 microsec
c not an important limitation
c
	smin=-500.0d0/tres	!fastest = 1 microsecond (the default ratemax)
c	smin=-100.0d0/tres	!fastest = 1 microsecond (the default ratemax)
c	smin=-709.0d0/tres
c	smin=-1.d20
	smax=-1.d-20
	nsmax=1000		!max number of splits (nsplit=15 in test run)
10	call GfuncX(sa,nga)	!calculate g(sa)
	if(nga.gt.0) then		!should be no roots below sa
	   sa=sa*4.d0
	   if(sa.ge.smin) goto 10
	   goto 99			!error
	endif
11	call GfuncX(sb,ngb)	!calculate g(sb)
	if(ngb.lt.kX) then	!all roots should be below sb
	   sb=sb/4.d0
	   if(sb.le.smax) goto 11
	   goto 99			!error
	endif
c
	nroot=0	!number of roots found
c
c No need for split if there is only one root (in fact split does not work
c in this case)
	if(kX.eq.1) then
	   nroot=1
	   nsplit=0
	   if(ngb-nga.eq.1) then
		nroot=1
		s1(nroot)=sa	!output
		s2(nroot)=sb	!output
	   endif
	   goto 99
	endif
c
	nroot=0	!number of roots found
c	nlive=1
	ntodo=1
	ndone=kX-ntodo
c Ball's intlive(1,:) = int -= input sa, sb (so for output, s1=first col of
c intlive and  s2() = second col)
c Ball's glive(1,:) = input nga, ngb
c This notation was based on a definition of an interval as 'live' if it was
c unresolved (ie contained two or more roots) amd 'dead' when resolved
c (contains one root). Change this notation here -rather than 'nlive' use
c ntodo=kX-ndone.  Probably can't use output arrays for everything -need to
c define equivalent of intlive(:,:) and glive(:,:) for intermediates
c These are both (kX x 2) arrays -
c replace intlive(kX,2) with sv1(kX),sv2(kX)
c replace glive(kX,2) nv1(kX), nv2(kX)
	sv1(1)=sa
	sv2(1)=sb
	nv1(1)=nga
	nv2(1)=ngb
	nsplit=0	!counts iterations
c
c	do while(ntodo.gt.0)
	do while(ndone.lt.kX)
	   sa=sv1(ntodo)
	   sb=sv2(ntodo)
	   nga=nv1(ntodo)
	   ngb=nv2(ntodo)
	   call SPLIT(sa,sb,nga,ngb,GfuncX,
     &	sa1,sb1,sa2,sb2,nga1,ngb1,nga2,ngb2,nerrs)
	   if(nerrs.ne.0) then
      	print 4,sa
      	if(dprt) write(8,4) sa
4		format(/,
     &  ' ERROR: unable to split interval in root_fb: root = ',g13.6)
c		pause
		nerr=2
		goto 99
	   endif
	   nsplit=nsplit+	1
	   if(nsplit.gt.nsmax) goto 99
	   ntodo=ntodo-1
	   ndone=ndone+1
c 	Does either or both of the two subintervals output from SPLIT
c     contain only one root?
	   if(ngb1-nga1.eq.1) then
		nroot=nroot+1
		s1(nroot)=sa1	!output
		s2(nroot)=sb1	!output
c		ng1(nroot)=nga1
c		ng2(nroot)=ngb1
	   else
		ntodo=ntodo+1
		ndone=ndone-1
c==		intlive(ntodo,:)=int1
c==		glive(ntodo,:)=g1
		sv1(ntodo)=sa1
		sv2(ntodo)=sb1
		nv1(ntodo)=nga1
		nv2(ntodo)=ngb1
	   endif
c
	   if(ngb2-nga2.eq.1) then
		nroot=nroot+1
		s1(nroot)=sa2	!output
		s2(nroot)=sb2	!output
c		ng1(nroot)=nga2
c		ng2(nroot)=ngb2
	   else
		ntodo=ntodo+1
		ndone=ndone-1
		sv1(ntodo)=sa2
		sv2(ntodo)=sb2
		nv1(ntodo)=nga2
		nv2(ntodo)=ngb2
	   endif
	enddo
c Sort s1, s2 into ascending order
	call SORT2d(s1,s2,kX,100,.true.)
c
c Check that W(s) has different sign on each side of root
	if(check) then
	   do m=1,nroot
		d1=detwX(s1(m))
		id1=idsign	!sign of determinant from common
		d2=detwX(s2(m))
		id2=idsign	!sign of determinant from common
		if(id1*id2.eq.1) then
		   print 9,m,s1(m),d1,s2(m),d2
		   if(discprt) write(8,9) m,s1(m),d1,s2(m),d2
9		   format(' Problem with root #',i3,' in ROOT_FB',/,
     &	   ' det(W(',g13.6,')) = ',g13.6,/,
     &	   ' det(W(',g13.6,')) = ',g13.6,/,
     &	   ' have same sign')
		   nroot=nroot-1
		   nerr=1
c		   pause
		endif
	   enddo
	endif
c
99	if(nroot.lt.kX) then
         print 2,nroot,kX,nsplit
         if(dprt) write(8,2) nroot,kX,nsplit
2	   format(
     &  ' Only ',i2,' out of ',i2,' roots located in ROOT_FB (',i4,
     & ' steps)')
	else
         print 3,kX,nsplit
         if(dprt) write(8,3) kX,nsplit
3	   format(
     & ' Guesses for all ',i2,' roots located in ROOT_FB (',i4,
     & ' steps)')
	endif
	if(idebugt.eq.1) then
c=	   itlast=itick
	   ndeb=ndeb+1
	   call TIMER(itick)
	   print 401,ndeb,10*(itick-itickr)
	   if(discprt) write(8,3) ndeb,10*(itick-itickr)
401	   format(
     &   1x,i5,': Time in root_fb (ms) = ',i9)
	endif
	RETURN
	end

	subroutine SPLIT(sa,sb,nga,ngb,GfuncX,
     &  sa1,sb1,sa2,sb2,nga1,ngb1,nga2,ngb2,nerrs)
c Based on Frank Ball's matlab function, split.m
c Takes an interval sa, sb that contains at least two HJC roots and uses the
c bisection method to find two disjoint subintervals, each of which
c contains at least one root.
c Repeated calls to split from ROOTS_FB() locate all initial guesses in form of
c pairs of s values each of which enloses one root exactly.
c
c Modif 12/09/03 03:31pm by adding nerrs, because case found where
c sa,sb,sc become identical when a rate became very fast
c (all = -3136276.05446948 ie approx 0.3 microsec), so sticks in the do-while
c loop. It isn't clear how to define 'identical' so just set upper limit
c (say 1000?) to number of times the loop is done

c Input
c sa, sb =  the original interval (must contain at least two roots)
c nga = number of eigenvalues of H(sa) that are less than sa
c ngb = number of eigenvalues of H(sb) that are less than sb
c GfuncX = subroutine name -in call substitute GfuncA or GfuncF
c
c Output
c sa1,sb1 =one subinterval (includes at least one root)
c sa2,sb2 =other subinterval (includes at least one root)
c nga1,ngb1 = ng for sa1, sb1
c nga2,ngb2 = ng for sa2, sb2
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	logical end
	External GfuncX
c
	ntrymax=1000
	ntry=0
	nerrs=0
c
	end=.false.
	do while(.not.end)
	   sc=(sa+sb)*0.5d0
	   call GfuncX(sc,ngc)	!calculate g(sc)
	   if(ngc.eq.nga) then
		sa=sc
	   else if(ngc.eq.ngb) then
		sb=sc
	   else
		end=.true.
		sa1=sa	!interval 1
		sb1=sc	!interval 1
		sa2=sc	!interval 2
		sb2=sb	!interval 2
		nga1=nga	!surely sa has changed since input, but nga hasn't!
		ngb1=ngc
		nga2=ngc
		ngb2=ngb
	   endif
	   ntry=ntry+1
c	   if(ntry.gt.10) then
c		print 1,ntry
c1		format(' split problem in root_fb: attempt # ',i5)
c	   endif
	   if(ntry.gt.ntrymax) then
		nerrs=1	!return error
		end=.true.
	   endif
	enddo
	RETURN
	end

