	subroutine QSET_HJC(jset,theta,QT,QD,kfit,k)
c
c Make a Q matrix for data set #jset from the current set of parameters:
c for new version of HJCFIT that takes multiple data sets.
c QD and QT are made from theta() (leaving fixed values in QT unchanged)
c QT has no concentrations etc, QD() has rates multiplied by appropriate conc
c
c NB micro rev parameters (and diagonals) are NOT set in QT -FIXED 12/20/01 05:37pm
c so that micro rev (but not diagonals) ARE set in QT as well as QD
c
c QT and QD are in common/qdblk
c conc(i,j) in common/cblk
c
c Modif 06/14/06 02:38pm for fast block correction (see SETFBLK)
c NB the blocker correction appears ONLY in QD(), which has, eg,
c alpha2/(1.0+B/KB) as actual transition rate.  As with aother connc dep
c transition rates, the conc is not present in theta0() or QT()
c
c Modif 09/02/04 03:27pm so that ec50 constraints calculate the ec50 using the
c right conc of other ligands when nlig>1.  Safest way will to duplicate the
c whole conc(i,j) matrix in conc1(i,j) and alter only the latter by inserting
c appropriate conc of other ligands, from conc_ec1 and conc_ec2
c
c Modif 10/09/03 11:01am for constraining a 2nd ec50
c
c Modif 11/27/01 02:59pm so that if called with k negative, this is signal
c to do everything EXCEPT put concentrations into output QD (and sign of k
c corrected)
c And QT,QD now made arguments (common/qdblk now split into qdblk1 and qdblk2)

c Modif 08/15/01 11:45am so that fitted values are inserted in QT as well
c as in QD (non-fitted, fixed values should be in QT before the fit,
c and never change)
c
c Modif 08/15/01 09:11am for common/indblk, so when nsetq>0 (for indep
c  models only) the rest of Q(i,j) are set at the end (in common with getqd,
c  qdefine and prtmod)
c
c Modif 07/12/00 12:17pm to add constraint by specified EC50
c
c Modif 03/20/00 05:42am Convert theta, simplex hjclik etc to real*8
c
c Modif 03/13/00 04:16pm
c Now calculate constraints BEFORE conc incorporated, wish is more
c convenient for fixing a binding constant K, by constraining k(-1)=K*k(+1)
c
	real*8 theta(200)
	real*4 conc(10,10)
	real*8 QT(100,100),QD(100,100)	!QD has conc in
	integer IQf(100,100),irate1(200),jrate1(200)
	logical discprt,debprt,doconc
	common/dp/discprt
	common/deb1/debprt
	integer NSC(50),IM(50,100),JM(50,100),IX(100),JX(100)
	COMMON/EBLK/NEQ,IE(200),JE(200),IF(200),JF(200),EFAC(200)
c
	common/CBLK/nset,conc,jsetlast		!for HJCLIK,DISP etc
c=	common/CBLK/nset,concA,concB,jsetlast		!for HJCLIK,DISP etc
c=	common/QDBLK/QT,QD,npar,IQf,irate1,jrate1,nlig
c=	common/QDBLK1/QT,QD
	common/QDBLK2/npar,IQf,irate1,jrate1,nlig   !hjclik, qset_hjc only
	common/cpar/ncdep,IX,JX,x
	common/LIG/nligsav,IL(100)
	COMMON/MPAR/NCYC,NSC,IM,JM
	COMMON/indblk/nsetq,ieq(200),jeq(200),ifq(200),jfq(200),efacq(200)
c For ec50 constraint
	logical fixec50,prtec50,fixec51,samegues,dcmodel
	real*8 ec50,xqlo,xqhi
	real*4 curinfs,curmax,concmax
	common/ec/fixec50,imod0,ec501,i50,j50,m50,prtec50,
     & xqlo,xqhi,dcmodel
	real*8 qval,qval2,ec501,ec50out,ec5012,ec50out2,
     &	x1s2,x2s2,fx12,fx22
c define common/ec1/ for printing in simphjc (in DSPOUT)
c Modif 10/09/03 05:24pm for case where two ec50s are constrained
c Modif 09/05/04 08:29am so that the (first) fixed ec50 value is called
c ec501 in common/ec/ (previously it was called 'ec50' but danger
c that this will be overwritten, So now need copy with diff name (ec5011)
c for common/ec1/
	real*8 ec5011
	common/ec1/nfixec51,fixec51,qval,ec5011,ec50out,i501,j501,x1s,x2s,
     &fx1,fx2,qval2,ec5012,ec50out2,i5012,j5012,x1s2,x2s2,fx12,fx22  !for print in SIMPHJC
c
	real*8 x1,x2,epsx,xout,yout,fx1,fx2,qsav,gfac,pmax,temp
	real*8 x1s,x2s,one,fx
	integer Leq50(5)
c For penalty function
	logical penalty
	real*8 penfunc,penfac
	common/pen/penalty,penfunc,penfac 	!for hjclik, qset_hjc
c Commons added 01/27/03 12:48pm for prtec50
	logical monot
	real*8 cfacd,dgamma(100)
	character*20 ligname(10)
	common/ec2/monot,nlvar,vhold,vref,cfacd,dgamma !for ec50_hjc, qset_hjc
	common/ec3/ligname		!for qset_hjc
c Added 10/09/03 11:01am for constraining a 2nd ec50
	real*8 ec502,xqlo2,xqhi2
	real*4 conc_ec1(10),conc_ec2(10),concsav(10),conc1(10,10)
	common/ec4/nfixec50,nlvar2,ec502,i502,j502,m502,
     &     xqlo2,xqhi2,conc_ec1,conc_ec2	 !for 2nd fixed ec50
c
	logical debug,caplock,deb
c
	logical*4 fastblk
	real*4 aKB
	integer*4 ifb(20),jfb(20),mfb(20),iflig
	common/fblck/fastblk,aKB,nfblock,iflig,ifb,jfb ,mfb
c
	debug()=caplock()
c
c	if(debug()) prtec50=.not.prtec50     !toggle
	kdim=100			!size of QT, QD
	one=1.0d0
	penalty=.false.
	penfunc=0.0d0
	doconc=.true.
	if(k.lt.0) then
	   k=-k
	   doconc=.false.
	endif
c
c Copy conc(i,j) to conc1(i,j) so former stays unchanged as conc used
c for set #j
	do i=1,10
	   do j=1,10
		conc1(i,j)=conc(i,j)
	   enddo
	enddo
c
c
c  Fixed parameters are already in QT
c  First insert the kfit estimated param from theta()
c  Copy fixed values to QD (other values overwritten later)
	do i=1,k
	   do j=1,k
		QD(i,j)=QT(i,j)
	   enddo
	enddo
c
c Copy current parameter estimates from theta (still without conc) into QD
	do n=1,kfit
	   i=irate1(n)
	   j=jrate1(n)
	   QT(i,j)=theta(n)	!define QT too
	   QD(i,j)=theta(n)
	enddo
c
c  Next calculate  constrained values (before micro rev) (as in QSETD)
c ALSO now moved to calculate BEFORE conc incorporated, wish is more
c convenient for fixing a binding constant K, by constraining k(-1)=K*k(+1)
	if(neq.gt.0) then
	   do L=1,neq
		if(efac(L).gt.0.) then
		   temp=dble(efac(L))*QT(if(L),jf(L))
		   QT(ie(L),je(L))=temp
		   QD(ie(L),je(L))=temp
		else
		   temp=dble(-efac(L)) - QT(if(L),jf(L))
		   if(temp.lt.0.d0) then
			temp=1.d-5	!make sure not neg
			print 10 ,ie(L),je(L),temp
10			format(' Warning: Q(',i2,',',i2,') set negative by',
     &		' constraint -reset to ',g12.5)
		   endif
		   QT(ie(L),je(L))=temp
		   QD(ie(L),je(L))=temp
		endif
	   enddo
	endif
c
c Now fix one rate from EC50, if required
c NB Rates are in QD which has NOT yet got conc incorporated, as required
c by MODEC50.  Until conc section reached, QT and QD are same, and both
c are reset. Now work on QTthroughout and copy to QD ate end of ec50 section
c NB if any other parameter is constrained to be a multiple of Q(i50,j50)
c then this parameter must be altered too, before EC50 is calculated!
c NB EC50_HJC uses only conc(i,1), so temporarily overwrite conc(i,1)
c with conc_ec1/2 that contain conc of other ligands for ec50 constraints
	if(fixec50) then
	 do i=1,nlig
	    concsav(i)=conc1(i,1)
	 enddo
	 do ifix=1,nfixec50
	   if(ifix.eq.1) then
		nlv=nlvar
		i5=i50
		j5=j50
c		m5=m50
		ec5=ec501
		do i=1,nlig
		   conc1(i,1)=conc_ec1(i)	!conc of other ligands at which ec50 determined
		enddo
	   else if(ifix.eq.2) then
		nlv=nlvar2
		i5=i502
		j5=j502
c		m5=m502
		ec5=ec502
		x1s=x1s2
		x2s=x2s2
		fx1=fx12
		fx2=fx22
		do i=1,nlig
		   conc1(i,1)=conc_ec2(i)	!conc of other ligands at which ec50 determined
		enddo
	   endif
	   deb=debug()
	   dlog2=0.693147180559945d0
c	   gfac=2.d0	!factor for initial guesses
	   gfac=10.d0	!factor for initial guesses
	   qsav=QT(i5,j5)	!save original entry (?)
c       Before calculating any ec50, check whether any other rate
c       is constrained to be a multiple of QD(i50,j50)=theta0(m50)
c       (could be more than one, possibly, record them in Leq50()
c	  NB need a complete QD (except for conc) to calculate EC50 -so
c	  need to set MR also each time (don't need diagonals, because
c	  in MODEC50 just need the rate constants which are defined
c	  from the off-diag elements of QD
	   if(neq.gt.0) then
		n=0
		do L=1,neq
		   if(if(L).eq.i5.and.jf(L).eq.j5) then		!found one
			n=n+1
			Leq50(n)=L
		   endif
		enddo
		neq50=n
	   endif
c
c EC50 never printed here -now all in PRNTEC50, called from hjcfit.
c	   if(prtec50) then	!print ec50, pmax for input QD
cc          NB setting MR should not depend on conc, but use original
cc		conc(i,j) to do it, because may have zero conc in conc_ec1
c		call SETMR(QT,jset,conc,kdim)
c
c		if(dcmodel.and.nlig.eq.1) then
c	         call MODEC50(imod0,QT,npar,kdim,pmax,fx1)	!fx1=ec50 for x1
c		else
cc========NB jset is now a parameter, but conc_ec1 etc have been
cc========copied into j=1 so call with jset1=1
c		   jset1=1
c		   call EC50_HJC(fx1,curinfs,curmax,concmax,cur0,pop0,
c     &	    pmax,nlv,QT,jset1,conc1,k,kdim,nerr)
c		   if(nerr.ne.0) then
c			call BELL(1)
c			print 590,nerr
c			if(discprt) write(8,590) nerr
c590			format('  Error in bisection for EC50: nerr = ',i2)
c		   endif
c		endif
c		print 593,ligname(nlv),cur0,pop0
c		if(discprt) write(8,593) ligname(nlv),cur0,pop0
c593		format(
c     &      ' For ',a10,':',/,
c     &      ' at zero concentration, current (pA) = ',g12.5,
c     &	', Popen = ',g12.5)
c		print 6, 1.d6*fx1,pmax
c		if(discprt) write(8,6) 1.d6*fx1,pmax
cc6		format(' EC50 = ',g13.6,', Max current = ',f8.4)
c6		format(' EC50 = ',g13.6,', Max Popen = ',f8.4)
c		print 594,i5,j5,qt(i5,j5)
c594		format(' Q(',i2,',',i2,') = ',g13.6,/)
c	   endif	!end of if(prtec50)
c
	   x1=qsav/gfac
	   x2=qsav*gfac
	   ntrymax=10	!max number of attempts to find good x1, x2 guesses
	   ntry=0
	   nreset=0
2	   continue
	   ntry=ntry+1
	   if(x1.gt.x2) then		!ensure x2 > x1
		temp=x1
		x1=x2
		x2=temp
	   endif
	   if(x1.lt.xqlo) then
		x1=xqlo
		nreset=nreset+1
	   endif
	   if(x2.gt.xqhi) then
		x2=xqhi
		nreset=nreset+1
	   endif
c Next line aborts repeated tries with same initial guess, as can occur
c when one guess hits xqlo or xqhi. However does not seen to be very effective
c because, at least in example tested, hardly ever get beyond ntry=1. Repeated
c 'Error in qset_hjc' which show the same x1, x2 each time actually result from
c leaving qset_hjc, back via hjclik, to simphjc which then alters parameters a
c bit (mostly relections at line 344 in my example) -if this reduces q(i50,j50)
c after resetting to xqlo here get same guesses, x1, x2, here on next entry to
c hjclik/qset_hjc. Therefore define samegues here to prevent excessive printing
c below. To do this, need to remember the guess used from previous calls
c of qset_hjc from hjclik -do this by adding x1s,x2s to common/ec1/
	   x1t=dabs(x1-x1s)
	   x2t=dabs(x2-x2s)
		samegues=.false.
	   if(x1t.lt.1.d-2.and.x2t.lt.1.d-2) then
		samegues=.true.
		if(ntry.gt.1) ntry=ntrymax
	   endif
	   x1s=x1	!save to print in case of error
	   x2s=x2	!save to print in case of error
c Test initial guesses. If they are bad, now use improved strategy to get good
c guesses (a) when func has neg slope: if x1,x2 both give fx>ec50 then guesses
c to low so make new x1= old x2 and new x2=gfac*new x1, and so on -in this case
c can keep gfac constant and x1, x2 should both creep up until they contain
c the answer
c PROBLEMS found here
c (a) if a parameter goes negative it is set in HJCLIK to 1.d-15 (at present)
c and this results in silly values for fx1 and fx2 -no sensible solution
c possible
c (b) in this case may have fx1=fx2 essentially, which original code
c did not allow for
	   QT(i5,j5)=x1
	   if(neq50.gt.0) then	!reset others that are constrained
		do n=1,neq50
		   L=Leq50(n)
		   QT(ie(L),je(L))=dble(efac(L))*QT(i5,j5)
		enddo
	   endif
	   call SETMR(QT,jset,conc,kdim)
	   if(dcmodel.and.nlig.eq.1) then
		call MODEC50(imod0,QT,npar,kdim,pmax,fx1)	!fx1=ec50 for x1
	   else
c NB jset is now a parameter of EC50_HJC, but conc_ec1 etc have been
c copied into j=1 so call with jset1=1 here
		jset1=1
		call EC50_HJC(fx1,curinfs,curmax,concmax,cur0,pop0,
     &	 pmax,nlv,QT,jset1,conc1,k,kdim,nerr)
		if(nerr.ne.0) then
		   call BELL(1)
		   print 590,nerr
		   if(discprt) write(8,590) nerr
590		   format('  Error in bisection for EC50: nerr = ',i2)
		endif
	   endif
c
	   QT(i5,j5)=x2
	   if(neq50.gt.0) then	!reset others that are constrained
		do n=1,neq50
		   L=Leq50(n)
		   QT(ie(L),je(L))=dble(efac(L))*QT(i5,j5)
		enddo
	   endif
	   call SETMR(QT,jset,conc,kdim)
	   if(dcmodel.and.nlig.eq.1) then
		call MODEC50(imod0,QT,npar,kdim,pmax,fx2)	!fx2=ec50 for x2
	   else
		jset1=1
		call EC50_HJC(fx2,curinfs,curmax,concmax,cur0,pop0,
     &	 pmax,nlv,QT,jset1,conc1,k,kdim,nerr)
		if(nerr.ne.0) then
		   call BELL(1)
		   print 590,nerr
		   if(discprt) write(8,590) nerr
c590		   format('  Error in bisection for EC50: nerr = ',i2)
		endif
	   endif
c
c Do guesses include required ec50?
c  Can get x1=x2=xqlo (or xqhi) -need to cope with this
	   if(dabs(x2/x1 - one).lt.5.d-3) then	!0.5 percent
		if(dabs(x1/xqlo-one) .lt.5.d-3) then	!0.5 percent
		   print 8,i5,j5,xqlo
8		   format(
     &' Guesses for q(',i2,',',i2,') have converged on lower limit = '
     &		,g13.6)
		else
		   print 9,i5,j5,xqhi
9		   format(
     &' Guesses for q(',i2,',',i2,') have converged on upper limit = '
     &		,g13.6)
		endif
		QT(i5,j5)=x1
		ec5011=1.d6*ec5
		ec50out=1.d6*fx1	!for print in simphjc
c==		goto 5
	   endif
	   if((fx1-ec5)*(fx2-ec5).gt.0.0) then		!same sign -guesses bad
		if(fx1.gt.ec5) then		!so fx2>ec50 too
		   if(fx1.gt.fx2) then		!neg slop (problem if not monotonic?)
			x1=x2
			x2=gfac*x2
		   else if(fx1.lt.fx2) then				!pos slope
			x2=x1
			x1=x1/gfac
		   else if(fx1.eq.fx2) then				!widen guesses
			x1=x1/gfac
			x2=x2*gfac
		   endif
		else			!fx1 < ec50 (and fx2 < ec50)
		   if(fx1.gt.fx2) then		!neg slop (problem if not monotonic?)
			x2=x1
			x1=x1/gfac
		   else if(fx1.lt.fx2) then				!pos slope
			x1=x2
			x2=x2*gfac
		   else if(fx1.eq.fx2) then				!widen guesses
			x1=x1/gfac
			x2=x2*gfac
		   endif
		endif
		if(x1.lt.xqlo) then
		   x1=xqlo
		   nreset=nreset+1
		endif
		if(x2.gt.xqhi) then
		   x2=xqhi
		   nreset=nreset+1
		endif
c need to go to error part if x1=x2 (e.g both=xqlo or both=xqhi)
		if(x1.ge.xqlo.and.x2.le.xqhi.and.
     &		dabs(x2-x1).gt.0.01d0) then
		   if(ntry.le.ntrymax.and.nreset.le.10) goto 2
		else
c		  If get many repeats of same guesses (x1 x2), e.g. if x1 hits xqlo,
c		  then next line prints for every function evaluation which slows
c		  program without conveying much
		   if(.not.samegues) then
		      print 3,ligname(nlv),1.d6*ec5,x1s,1.d6*fx1,
     &			x2s,1.d6*fx2
		      if(debprt) write(8,3) ligname(nlv),1.d6*ec5,
     &			x1s,1.d6*fx1,x2s,1.d6*fx2
3		      format(
     & ' ERROR in QSET_HJC: ec50 for ',a10,' should be ',g13.6,
     & ', but for',/,
     & '  rate = ',g13.6, ' ec50 = ',g13.6,' and for',/,
     & '  rate = ',g13.6, ' ec50 = ',g13.6)
		   endif
c		   STOP   !just carry on with wrong ec50 and hope it gets back later?
c		If initial guesses don't include result then use whichever of x1 and
c		x2 is the closer, and then skip bisection.  If fx2 nearer then QD
c		is already set, if not then reset using x1
		   if(dabs(fx1-ec5).le.dabs(fx2-ec5)) then	!x1 nearer
			QT(i5,j5)=x1s
			if(neq50.gt.0) then	!reset others that are constrained
			   do n=1,neq50
				L=Leq50(n)
				QT(ie(L),je(L))=dble(efac(L))*QT(i5,j5)
			   enddo
			endif
			call SETMR(QT,jset,conc,kdim)
			fx=fx1	!so fx is the closest to ec50
		   else			!x2 nearer
			QT(i5,j5)=x2s
			if(neq50.gt.0) then	!reset others that are constrained
			   do n=1,neq50
				L=Leq50(n)
				QT(ie(L),je(L))=dble(efac(L))*QT(i5,j5)
			   enddo
			endif
			call SETMR(QT,jset,conc,kdim)
			fx=fx2	!so fx is the closest to ec50
		   endif
c		Could also try setting
c		a penalty function here (dependent on difference between ec50 that
c		is required and the ec50 achieved (if this used, set the penalty to
c		0.00d0 exactly when solution IS found (even though solution has
c		 limited accuracy).  Say penalty is 1 log-likelihood unit for every
c		10% error in ec50.  Using penfunc=penfac*dabs(ec50-fx)/ec50 results
c	 	in bigger penalty for ec50 that is too big by factor of 2 than for
c		ec50 that is too small by factor of 2
		   penalty=.true.
c		   penfunc=penfac*dabs(ec5-fx)/ec5
		   if(fx.gt.ec5) then
			penfunc=penfac*dabs(ec5-fx)/ec5
		   else
			if(fx.gt.1.d-8) then
			   penfunc=penfac*dabs(ec5-fx)/fx
			else
			   penfunc=100.d0 	!if fx too small to divide
			endif
		   endif
		   if(deb) then
			print 7, penfunc
			if(debprt) write(8,7) penfunc
7			format(' Penalty function = ',f8.3)
		   endif
		   ec5011=1.d6*ec5
		   ec50out=1.d6*fx	!for print in simphjc
		   qval=-1.0d0	!to signal guesses bad
		   qval2=-1.0d0	!to signal guesses bad
		   x1s=x1
		   x2s=x2
		   yout=-1.d0
		   goto 5
		endif
	   endif	!end of bad guesses if((fx1-ec5)*(fx2-ec5).gt.0.0) (skip bisection bit via 'goto 5')
c
c Now bisection using x1, x2 as good guesses
	   if(fx1.ge.ec5) then
		temp=x1
		x1=x2
		x2=temp
	   endif
	   epsx=x1/10000.d0	!accuracy in qd()
	   if(deb) then
		print 31,1.d6*ec5,x1,1.d6*fx1,x2,1.d6*fx2,epsx
31		format(
     &	 ' Initial guesses in QSET_HJC: ec50 should be ',g13.6,/,
     &	 '  x1 = ',g13.6, ' ec50 = ',g13.6,' and for',/,
     &	 '  x2 = ',g13.6, ' ec50 = ',g13.6,/,
     &	 '  epsx = ',g13.6,/,
     &	 ' iter      x1       x2        xout        Yout      Yerr')
	   endif
	   nstep=int4(sngl(dlog(dabs(x1-x2)/epsx)/dlog2+0.5))
	   do ns=1,nstep
		xout=0.5d0*(x1+x2)
		QT(i5,j5)=xout
c       reset others that are constrained
		if(neq50.gt.0) then	!reset others that are constrained
		   do n=1,neq50
			L=Leq50(n)
			QT(ie(L),je(L))=dble(efac(L))*QT(i5,j5)
		   enddo
		endif
		call SETMR(QT,jset,conc,kdim)
		if(dcmodel.and.nlig.eq.1) then
		   call MODEC50(imod0,QT,npar,kdim,pmax,Yout)
		else
		   jset1=1
		   call EC50_HJC(Yout,curinfs,curmax,concmax,cur0,pop0,
     &		 pmax,nlv,QT,jset1,conc1,k,kdim,nerr)
		   if(nerr.ne.0) then
			call BELL(1)
			print 590,nerr
			if(discprt) write(8,590) nerr
c590			format('  Error in bisection for EC50: nerr = ',i2)
		   endif
		endif
		Yerr=Yout-ec5			!Yout=calc ec50 for the current theta
		if(deb) then
		   print 4,ns,x1,x2,xout,Yout,Yerr
4		   format(1x,i4,5g13.6)
		endif
		if(Yerr.eq.0.0d0) goto 98	!exact solution
		if(Yerr.lt.0.0d0) then
		   x1=xout
		else
		   x2=xout
		endif
	   enddo
98	   continue
	   if(deb) pause
c	   QD(i5,j5)=theta0(m5)
	   QT(i5,j5)=xout
c
c  An element of QT has been changed so re-calculate MR
	   call SETMR(QT,jset,conc,kdim)
c=	   QD(i5,j5)=xout
	   QD(i5,j5)=QT(i5,j5)	!conc inserted below
	   qsav=xout
c NB this is not one of the free parameters so it is not in theta(),
c only in theta0()
c Save values in common/ec1 for print in simphjc during iterations
c Jump to here? Or after settings below
5	   continue		!jump here if bisection skipped
	   fixec51=fixec50
	   nfixec51=nfixec50
	   if(ifix.eq.1) then
		qval=xout
		ec5011=1.d6*ec5
		ec50out=1.d6*yout
		i501=i5
		j501=j5
	   else if(ifix.eq.2) then
		qval2=xout
		ec5012=1.d6*ec5
		ec50out2=1.d6*yout
		i5012=i5
		j5012=j5
		x1s2=x1s
		x2s2=x2s
		fx12=fx1
		fx22=fx2
	   endif
c5	   continue		!jump here if bisection skipped
c and copy all values to QD (conc inserted below)
	   do i=1,k
		do j=1,k
		   QD(i,j)=QT(i,j)
		enddo
	   enddo
c===for debug
c	   print 32,ligname(nlv),1.d6*yout,1.d6*ec5
c32	   format(' In qset_hjc: EC50 for ',a10,' = ',g13.6,
c     &   ' (should be ',g13.6,')')
	   do i=1,nlig
		conc1(i,1)=concsav(i)	!restore conc(i,1)
	   enddo
	 enddo	!end of loop repeated if 2nd ec50 to be fixed (nfixec50=2)
	endif		!end of if(fixec50)
c
c  Could now set V-dep rates as in QSETD, but not in at present
c
c  Now set conc-dep rates AND define QT also, with all constraints
c etc incorporated but without conc
101	continue	!skip here to here if ec50 fails
	if(ncdep.gt.0) then
	   do L=1,ncdep
		i=ix(L)
		j=jx(L)
c Define NLIG=number of ligands and IL(i)=ligand type for ith conc-dep
c rate i=1,...,ncdep(needed only if NLIG>1)
		QT(i,j)=QD(i,j)
		if(doconc) then
		   x1=conc(IL(L),jset)
		   QD(i,j)=QD(i,j)*dble(x1)
		endif
	   enddo
	   jsetlast=jset	!record the conc last used to define QD
	endif
c
c Correct rates for fast channel block
	if(fastblk.and.doconc) then
	   do n=1,nfblock
		i=ifb(n)
		j=jfb(n)
		x1=dble(conc(iflig,jset)/aKB)
		if(mfb(n).gt.0) then
		   QD(i,j)=QD(i,j)/(1.0d0 + x1)
		else
		   QD(i,j)=QD(i,j)*x1/(1.0d0 + x1)
		endif
	   enddo
	endif
c
c  Correct rates for KMFAST case -NOT IN FOR NOW
c   Multiply rates out of state 1 by appropriate
c   occupancy to get true rate out of 1
c	if(kmfast) then
c	   do 7 J=2,K		!not 1,K
c	   if(abs(QD(1,J)).LT.1.d-20) GOTO 7	!NO ROUTE
c	   QT(1,J)=QT(1,J)*dble(PSTAR(KMCON(J-1)))
c	   QD(1,J)=QD(1,J)*dble(PSTAR(KMCON(J-1)))
c7	   continue
c	endif
c
C  Next micro rev- recalc those rates that are defined in terms of
c  the others by micro rev. NB already checked that no micro rev param
c  are already constrained, so no need to provide for this (cf QSETD)
c
c NB in case of fast block correction, MR should be obeyed for the effective
c rates, not the true ones, so don/t apply setmr to QT, but to QD
c then set QT by removing conc from QD.  QD with all conc set should
c contain the effective rates
c
	if(.not.fastblk) then
	   call SETMR(QT,jset,conc,kdim)
	   call SETMR(QD,jset,conc,kdim)
	else
	   call SETMR(QD,jset,conc,kdim)
	   call REMCONC(QD,QT,k,jset,conc,kdim)	!remove conc from QD
	endif
c NB micro rev parameters (and diagonals) are NOT set in QT
c
c12	continue
c
c For independent models, set all the other Q(i,j) that are
c the same basic npar rate constants
	if(nsetq.gt.0) then
	   do L=1,nsetq
		QT(ieq(L),jeq(L))=dble(efacq(L))*QT(ifq(L),jfq(L))
		QD(ieq(L),jeq(L))=dble(efacq(L))*QD(ifq(L),jfq(L))
	   enddo
	endif
c
c Lastly set diagonal elements in QD (not set in QT)
	call SETDIAG(QD,k,kdim)
	do i=1,k
	   QT(i,i)=0.d0
	enddo
c
	call CHECKQD(QD,k,jset,npar,nlig,nlvar,.true.,kdim)
c
	RETURN
	end

	subroutine SETDIAG(QD,k,kdim)
	real*8 QD(kdim,kdim)
	do i=1,k
	   QD(i,i)=0.0d0
	   do j=1,k
		if(j.ne.i) then
		   QD(i,i)=QD(i,i)-QD(i,j)
		endif
	   enddo
	enddo
	RETURN
	end

