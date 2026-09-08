	subroutine PDFOUTd(ptitle,LI,LJ,w1,eigen,ncomp,am1,sd1,
     & kdim,pon,fprint,discprt)
c PDFOUTd is same as PDFOUTs but
c  (1) inputs real*8 values eigen(1/sec) and amplitudes w1(), rather
c	 than area,tau (need to do this to cope with case where eigenvalue=0)
c	 Now all real*8 internally, but outputs mean and SD as real*4 am1,sd1.
c  (2) Copes with cases with one or more zero eigenvalues. Get these eg in
c	 case of pdf of shut times ater step to zero conc, where there may be no
c	 more openings after the step so shut time is infinite, i.e. get delta
c	 function at t=infinity (LT of this is zero so part of the area seems
c	 to disappear, so area of this bit must be found by difference (unless
c	 there is some subtle way of taking the limit that preserves the area
c	 of the delta function?))
c  (3) does not print PTITLE if blank
c
c PDFOUTs is version of PDFOUT1 that takes SINGLE precision tau (in ms)
c and areas, as output by PDFOPEN() etc. Mean (AM) and SD also in msec.
C PDFOUT1 is combination of PDOUT2 which takes directly the coeffs (w1) and
c eigenvalues (r1) in double precision to print PDF, and PDFOUT which
c prints f(2,3), f(2) or just f(t) for condional dists. Also has overall title
c (up to 70 char) in ptitle.
c To print pdfs, mean and SD in TCHAN2 and return mean and SD (real*8)
C	L1=i LJ=j gives f(t;i,j)
C	L1=i LJ=-1 gives f(t;i)
C	L1=-1 LJ=-1 gives f(t)
C	L1=-2 LJ=-1 gives f(inf)
C
c FPRINT=true  -types/prints everything
c FPRINT=false  -types everything, but prints only mean and sd (if pon())
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real*8 eigen(kdim),w1(kdim),x2
	real*8 tau(100),area(100)
	real*4 am1,sd1			!for output
	logical inf(100),infin		!tau is 'infinite'
	logical pon,fprint,discprt
c	character*70 ptitle
	character*(*) ptitle
c
c Calc tau,area
c Here, atot=total area for components with non-zero eigevalues
	nbig=0
	zero=0.0d0
	atot=zero
	do 20 m=1,ncomp
	x2=dabs(eigen(m))
	if(x2.ge.1.d-20) then
	   inf(m)=.false.
	   tau(m)=1.0d3/x2
	   area(m)=w1(m)/x2
	   atot=atot+area(m)
	else
	   inf(m)=.true.		!tau is 'infinite'
	   nbig=nbig+1
	endif
20	continue
C Type LHS
	nl=NBLANK(ptitle)
	if(nl.gt.0) print *,CHARNB(ptitle)
	IF(LI.EQ.-1.AND.LJ.EQ.-1) print 1
	IF(LI.GE.0.AND.LJ.LT.0) print 2 ,LI
	IF(LI.GE.0.AND.LJ.GE.0) print 3 ,LI,LJ
	IF(LI.EQ.-2.AND.LJ.EQ.-1) print 4
	if(pon) then
	   if(nl.gt.0) write(7,*) ptitle
	   IF(LI.EQ.-1.AND.LJ.EQ.-1) write(7,1)
	   IF(LI.GE.0.AND.LJ.LT.0) write(7,2)LI
	   IF(LI.GE.0.AND.LJ.GE.0) write(7,3)LI,LJ
	   IF(LI.EQ.-2.AND.LJ.EQ.-1) write(7,4)
	endif
	if(discprt) then
	   if(nl.gt.0) write(8,*) ptitle
         if(LI.EQ.-2.AND.LJ.EQ.-1) write(8,4)
         if(LI.EQ.-1.AND.LJ.EQ.-1) write(8,1)
         if(LI.GE.0.AND.LJ.LT.0) write(8,2)LI
         if(LI.GE.0.AND.LJ.GE.0) write(8,3)LI,LJ
	endif
1	FORMAT( 1X,'f(t) = ')
2	FORMAT( 1X,'f(t;',I2,') =    ')
3	FORMAT( 1X,'f(t;',I2,',',I2,') = ')
4	FORMAT( 1X,'f(inf) = ')
50	format(//,1x,a70)
c
C TEST INDETERMINACY- in old version w(m) were set to zero in this case;
c here test for area=zero (but setting a=0.0 in indeterminate cases not
c yet done in any calling progs)
	do 10 m=1,ncomp
10	if(area(m).gt.1.d-20) goto 11	!out of loop
	print 12
	if(pon) write(7,12)
	if(discprt) write(8,12)
12	FORMAT(' INDETERMINATE: all areas are zero in PDFOUTd')
	RETURN
11	continue
c
c
	am=zero
	sum2=zero
	print 13
	if(pon) write(7,13)
	if(discprt) write(8,13)
13	format(
     & ' term   coeff (W)   rate const (1/sec)     area        tau (ms)'
     & )
	do 41 m=1,ncomp
	  infin=.false.		!pdf has no delta function
	  if(.not.inf(m)) then
		am=am + area(m)/dabs(eigen(m))
		sum2=sum2 + area(m)/(eigen(m)*eigen(m))
		print 14,m,w1(m),dabs(eigen(m)),area(m),tau(m)
		if(pon.and.fprint) write(7,14) m,w1(m),dabs(eigen(m)),
     &		area(m),tau(m)
	      if(discprt) write(8,14) m,w1(m),dabs(eigen(m)),
     &		area(m),tau(m)
14		format(i3,3x,g13.6,1x,g13.6,7x,g13.6,1x,g13.6)
	  else
		infin=.true.		!if any rate=0
		if(nbig.eq.1) then	!if only one zero eig print area here
		   ar=1.0d0-atot		!total area for all non-zero eigenvalues
		   area(m)=ar
		   print 140,m,w1(m),dabs(eigen(m)),ar
		   if(pon.and.fprint) write(7,140) m,w1(m),dabs(eigen(m)),ar
	         if(discprt) write(8,140) m,w1(m),dabs(eigen(m)),ar
140		   format(i3,3x,g13.6,1x,g13.6,7x,g13.6,'  >3.e12 years')
		else
		   print 141,m,w1(m),dabs(eigen(m))
		   if(pon.and.fprint) write(7,141) m,w1(m),dabs(eigen(m))
	         if(discprt) write(8,141) m,w1(m),dabs(eigen(m))
141		   format(i3,3x,g13.6,1x,g13.6,21x,' >3.e12 years')
		endif
	  endif
41	continue		!end of m=1,ncomp loop
c
c
	if(.not.infin) then		!usual case (no zero eigenvalues)
	   print 15,atot
	   if(pon.and.fprint) write(7,15) atot
         if(discprt) write(8,15) atot
15	   format(
     & '       Total area (for non-zero eigenvalues) = ',g13.6)
c
	   sd=2.0d0*sum2 - (am*am)
	   if(sd.gt.0.0d0) then
		sd=dsqrt(sd)
		am=am*1.d3			!in msec
		sd=sd*1.d3			!in msec
		print 16,am,sd,sd/am
		if(pon) write(7,16) am,sd,sd/am
	      if(discprt) write(8,16) am,sd,sd/am
16		format(
     & ' Mean (ms) = ',g13.6,'   SD= ',g13.6,'   SD/mean = ',g13.6)
	   else
		am=am*1.d3			!in msec
		print 17,am
		if(pon) write(7,17) am
	      if(discprt) write(8,17) am
17		format(
     & ' Mean (ms) = ',g13.6,'but variance is negative')
	   endif
	else if(nbig.ge.1) then
	   am=zero 		!for mean of non-zero eig components
	   s=zero
	   do 161 m=1,ncomp
		if(.not.inf(m)) then
		   am=am+area(m)*tau(m)
		   s=s+area(m)
		endif
161	   continue
c	   if(nbig.ge.2) then
	      ar=1.0d0-atot		!total area for all non-zero eigenvalues
            print 142,ar
            if(pon) write(7,142) ar
            if(discprt) write(8,142) ar
142	      format(
     & ' Total area for all zero eigenvalues (inf lifetimes) = ',g13.6)
c	   endif
	   if(s.gt.1.d-10) then
	      print 162,am/s
	      if(pon) write(7,162) am/s
	      if(discprt) write(8,162) am/s
162		format(/,
     &	' Mean for components with non-zero eigenvalues = ',g13.6)
	   endif
	   if(atot.gt.1.d-10) then
		print 18,atot
		if(pon.and.fprint) write(7,18) atot
	      if(discprt) write(8,18) atot
18		format(
     &	 ' Total area for non-zero eigenvalues = ',g13.6,/,
     &       '   Fractions of this area for each component:')
		do 163 m=1,ncomp
		 if(.not.inf(m).and.dabs(area(m)).gt.1.d-20) then
      	   print 164,m,area(m)/atot
	         if(pon) write(7,164) m,area(m)/atot
	         if(discprt) write(8,164) m,area(m)/atot
164		   format(
     &	  '     relative area for component ',i2,' = ',g13.6)
		 endif
163		continue
	   endif
	endif
c For output make real*4 values
	am1=sngl(am)
	sd1=sngl(sd)
c
	RETURN
	END



