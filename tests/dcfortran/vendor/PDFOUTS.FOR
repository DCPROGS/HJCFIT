	subroutine PDFOUTs(ptitle,LI,LJ,area,tau,ncomp,am,sd,
     & kdim,pon,fprint,discprt)
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
c
c Modified 01/17/03 10:18am to make SD calculation real*8 internally to
c prevent underflows
c
	real tau(kdim),area(kdim)
	logical pon,fprint,discprt
	real*8 atot,sd2,sum2,a,am2,w,x
c	character*70 ptitle
	character*(*) ptitle
c
C Type LHS
c need to rmove leading blanks form ptitle?
c	do 51 i=1,70
c	i1=i
c	ip=ichar(ptitle(i:i))
c51	if(ip.ne.32.and.ip.ne.0) goto 52
c52	i1=i1-1
c	ns=nblank(ptitle)
c	print 50,ptitle(i1:ns)
c	print 103
c103	format(/)
	print *,CHARNB(ptitle)
	IF(LI.EQ.-1.AND.LJ.EQ.-1) print 1
	IF(LI.GE.0.AND.LJ.LT.0) print 2 ,LI
	IF(LI.GE.0.AND.LJ.GE.0) print 3 ,LI,LJ
	IF(LI.EQ.-2.AND.LJ.EQ.-1) print 4
	if(pon) then
c	write(7,50) ptitle(i1:ns)
c	write(7,103)
	write(7,*) ptitle
	IF(LI.EQ.-1.AND.LJ.EQ.-1) write(7,1)
	IF(LI.GE.0.AND.LJ.LT.0) write(7,2)LI
	IF(LI.GE.0.AND.LJ.GE.0) write(7,3)LI,LJ
	IF(LI.EQ.-2.AND.LJ.EQ.-1) write(7,4)
	endif
	if(discprt) then
c	write(8,50) ptitle(i1:ns)
c	write(8,103)
	write(8,*) ptitle
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
	DO 10 M=1,ncomp
10	if(area(m).gt.1.e-20) goto 11	!out of loop
	print 12
	if(pon) write(7,12)
	write(8,12)
12	FORMAT(' INDETERMINATE: all areas are zero in PDFOUTs')
	RETURN
11	continue
c
c
	atot=0.0d0
	am2=0.0d0
	sum2=0.0d0
	print 13
	if(pon) write(7,13)
	if(discprt) write(8,13)
13	format(
     &' term   coeff (W)   rate const (1/sec)    area (%)     tau (ms)'
     & )
	DO 41 m=1,ncomp
	   x=1000.d0/dble(tau(m))	!lambda (1/sec)
	   a=dble(area(m))
	   w=x*a			!coeff
	   atot=atot+a
	   am2=am2 + a/x
	   sum2=sum2 + a/(x*x)
	   print 14,m,w,x,100.*a,tau(m)
	   if(pon.and.fprint) write(7,14) m,w,x,100.*a,tau(m)
         if(discprt) write(8,14) m,w,x,100.*a,tau(m)
14	   format(i3,3x,g13.6,1x,g13.6,7x,f11.6,2x,g13.6)
41	continue
c
	print 15,atot
	if(pon.and.fprint) write(7,15) atot
      if(discprt) write(8,15) atot
15	format(
     & '       Total area= ',g13.6)
c
	sd2=sqrt(2.0d0*sum2 - (am2*am2))
	am2=am2*1000.d0		!in msec
	sd2=sd2*1000.d0		!in msec
	am=sngl(am2)
	sd=sngl(sd2)
	print 16,am,sd,sd/am
	if(pon) write(7,16) am,sd,sd/am
      if(discprt) write(8,16) am,sd,sd/am
16	format(
     & ' Mean (ms) = ',g13.6,'   SD= ',g13.6,'   SD/mean = ',g13.6)
c
	RETURN
	END

