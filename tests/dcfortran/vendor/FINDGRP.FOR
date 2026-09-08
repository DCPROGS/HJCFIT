	subroutine FINDGRP(detWX,kX,s1,d1,s2,d2,isignf,
     &	sout1,sout2,sout,detout,nerr)
c Used in HJCASYMP in place of bisection when we have a group of very
c closely spaced roots that are, pro tem, being treated as identical
c When this ocuurs the number in the group has to be even so the DETW(s)
c has the same sign on each side of the group, and bisection cannot
c therefore be used.
c Strategy is to use course search between guesses to look for minimum
c (if detw() is positive) or maximum (if it is negative).  Then interpolate
c to get min/max more closely?
c
c  When looking for minimum, stop search if a maximum is reached (and
c vice versa), otherwise it is possible that a coarse search will find
c minimum close to next root (depending on mesh size and how close s1, s2
c are to adjacent roots. In fact can stop loop as soon as min found (i.e
c previous slope negative and next postive.
c
c If search encounters a change in sign of detw(s) then can immediately
c do bisection to get (one of) the group of roots
c If no change in sign detected then have two options
c   (a) do finer search min/max untl sign change IS found, then bisect
c   (b) Get closer estimate of where min max is and use position of this
c      for all roots in the group
c  (only the latter will work if there really are identical roots)
c
c  Input:
c 	s1, s2 = values that encompass the group, such that detw(s1) and detw(s2)
c	have the same sign
c	d1, d2 =detW(s1), detW(s2) calculated before entry
c  Output:
c     isignf=1 change of sign found between detW(sout1) and detW(sout2)
c     isignf=-1 no change of sign found -use minimum
c     isignf=0 detw=0 found so sout=exact solution
c	sout= (approx to) value of s at the stationary point (if no sign
c			change detected in detw() (use if isignf=-1), or exact
c			solution if isignf=0
c	detout= detW(sout) if sout is used (isignf=-1)
c     sout1, sout2= use when sign change IS found (isignf=1) -these are suitable
c	     		starting values for bisection
c
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	logical first
c from detwf, detwa
	common/det4/ndscale,idsign
	EXTERNAL detwx
c
	zero=0.0d0
	slo=s1
	shi=s2
c	d1=DETWX(slo)
c	d2=DETWX(shi)
	if(d1.gt.zero.and.d2.gt.zero) then		!expect min
	   goto 1			!look for min
	else if(d1.lt.zero.and.d2.lt.zero) then	!expect max
	   goto 2			!look for max
	elseif(d1.gt.zero.and.d2.lt.zero) then
	   print 71,d1,d2
71	   format('  Error in FINDGRP: d1, d2 = ',2g13.6)
	elseif(d1.gt.zero.and.d2.lt.zero) then
	   print 71,d1,d2
	else if(d1.eq.zero) then
	   isignf=0
	   sout=s1		!exact solution
	   detout=d1
	else if(d2.eq.zero) then
	   isignf=0
	   sout=s2		!exact solution
	   detout=d2
	endif
c
c Look for minimum
1	continue
	np=100
	dlogs=(dlog10(-slo)-dlog10(-shi))/dfloat(np - 1)
	s=slo
	det=d1
	detmin=1.d308
	first=.true.
	istep=0
c start loop
11	continue
	istep=istep+1
	slast=s		!record value of previous s value
c
	s=slo/(10**(dfloat(istep)*dlogs))
	detlast=det		!record value of previous determinant
	det=DETWX(s)
	if(istep.gt.1) sloplast=slope
	slope=(det-detlast)/(s-slast)
	if(nerr.eq.10) then	!determinant overflowed
	   goto 9		!deallocate and RETURN
	endif
	if(det.lt.zero) then
	   isignf=1
	   sout1=slast
	   sout2=s
	   goto 9		!return for bisection between sout1,sout2
	else if(det.eq.zero) then
	   isignf=0
	   sout=s		!exact solution
	   detout=det
	   goto 9		!return with exact solution
	endif
	if(det.lt.detmin) then
	   detmin=det
	   smlast=slast
	   smin=s
	   snext=slo/(10**(dfloat(istep+1)*dlogs))
	endif
c leave loop if minimum passed or np reached
	if(.not.(slope.gt.zero.and.sloplast.lt.zero).AND.
     &	istep.lt.np) goto 11
c If reach here then all np steps found without sign change being found so
c best estimate so far is smin. Do finer search between points each side of
c smin.  If arrive here again after 2nd search then still no sign change
c found so use smin value in sout
	if(.not.first) then
	   isignf=-1	!use minimum in sout
	   sout=smin
	   detout=detmin
	   goto 9	!leave with best yet
	else
	   first=.false.		!doing 2nd search
	   slo=smlast
	   shi=snext
	   np=50
	   dlogs=(dlog10(-slo)-dlog10(-shi))/dfloat(np - 1)
	   s=slo
	   det=DETWX(s)
	   istep=0
	   isignf=-1	!use minimum in sout
         goto 11 	!do fine search
	endif
c
c Look for MAXIMUM
2	continue
	np=100
	dlogs=(dlog10(-slo)-dlog10(-shi))/dfloat(np - 1)
	s=slo
	det=d1
	detmax=-1.d308
	first=.true.
	istep=0
c start loop
21	continue
	istep=istep+1
	slast=s		!record value of previous s value
c
	s=slo/(10**(dfloat(istep)*dlogs))
	detlast=det		!record value of previous determinant
	det=DETWX(s)
	if(istep.gt.1) sloplast=slope
	slope=(det-detlast)/(s-slast)
	if(nerr.eq.10) then	!determinant overflowed
	   goto 9		!deallocate and RETURN
	endif
	if(det.gt.zero) then
	   isignf=1
	   sout1=slast
	   sout2=s
	   goto 9		!return for bisection between sout1,sout2
	else if(det.eq.zero) then
	   isignf=0
	   sout=s		!exact solution
	   detout=det
	   goto 9		!return with exact solution
	endif
	if(det.gt.detmax) then
	   detmax=det
	   smlast=slast
	   smax=s
	   snext=slo/(10**(dfloat(istep+1)*dlogs))
	endif
c leave loop if maximum passed or np reached
	if(.not.(slope.lt.zero.and.sloplast.gt.zero).AND.
     &	istep.lt.np) goto 21
c If reach here then all np steps found without sign change being found so
c best estimate so far is smax. Do finer search between points each side of
c smin.  If arrive here again after 2nd search then still no sign change
c found so use smax value in sout
	if(.not.first) then
	   isignf=-1	!use minimum in sout
	   sout=smax
	   detout=detmax
	   goto 9	!leave with best yet
	else
	   first=.false.		!doing 2nd search
	   slo=smlast
	   shi=snext
	   np=50
	   dlogs=(dlog10(-slo)-dlog10(-shi))/dfloat(np - 1)
	   s=slo
	   det=DETWX(s)
	   istep=0
	   isignf=-1	!use minimum in sout
         goto 21 	!do fine search
	endif

c Leave
9	continue
	RETURN
	end

