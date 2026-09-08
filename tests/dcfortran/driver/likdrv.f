c LIKDRV -- evaluate the Fortran HJCLIK once, on inputs read from a file.
c
c Written for one purpose: to put the original Fortran likelihood and
c HJCFIT's C++ Log10Likelihood in front of the same Q matrix and the same
c record, and print both numbers. Everything else about HJCFIT -- the
c prompts, the graphics, the simplex, the topology and constraint machinery
c -- is left out, because all of it has already been eliminated as an
c explanation for the alpha2/beta2 bias.
c
c Input file (free format, one item per line or whitespace separated):
c
c   k kA kB kC kD
c   conc            [M]      concentration, only for the record
c   tres            [s]      dead time
c   tcrit           [s]      critical shut time; <=0 means no burst division
c   ichs                     1 to use CHS vectors, 0 for equilibrium
c   Q(1,1) ... Q(k,k)        row-major, concentrations already in,
c                            diagonals ignored (recomputed)
c   nint                     number of intervals
c   t(1) a(1) ... t(n) a(n)  [s] and [pA]: alternating record
c
c Output: one line, "HJCLIK <log10 likelihood> <nint> <k>".

      program LIKDRV
      implicit none

      integer kmx,ndmx
      parameter (kmx=100,ndmx=200000)

c ---- HJCLIK's own arguments
      integer kfit,nd1,nd2,kab,kAm,kFm,km
      real*8 theta(200)
      real*4 tint(ndmx,1),ampl(ndmx,1)
      integer*1 iprops(ndmx,1)
c NB km is the *declared* dimension of QD and QT (100), not the number of
c states -- Hjcfit1-09122003.for line 538, "km=100 !in COMMON -dimension of
c QD, QT etc". SUBMAT is told km as QD's leading dimension, so passing the
c state count instead makes it read QD with the wrong stride and QFF comes
c out zero. kAm and kFm, by contrast, are kA and kF exactly (lines 1984-5).
      real*8 Z00A(3,4,100),Z10A(3,4,100),Z11A(3,4,100)
      real*8 Z00F(4,3,100),Z10F(4,3,100),Z11F(4,3,100)
      real*8 XAF(3,4,3),XFA(4,3,4)
      real*8 QEXPQA(4,3),QEXPQF(3,4)
      real*8 alpha2(1),beta2(1),aloglik(1)
      real*8 HJCLIK,value
      external HJCLIK

c ---- the commons HJCLIK reads
      integer kA,kB,kC,kD
      common/KBLK/kA,kB,kC,kD

      integer Nint
      real*4 tcritc
      logical burst,chsvec,badend
      common/HJCBLK/Nint(10),tcritc(10),burst(10),chsvec(10),badend(10)

      integer nset,jsetlast
      real*4 conc
      common/CBLK/nset,conc(10,10),jsetlast

      real*8 tresd
      common/resblk/tresd(10)

      real*8 QT,QD
      common/QDBLK1/QT(100,100),QD(100,100)

      integer npar,IQf,irate1,jrate1,nlig
      common/QDBLK2/npar,IQf(100,100),irate1(200),jrate1(200),nlig

      integer ngp,nscal
      real*4 an
      logical first
      common/nblk/ngp(10),an(10),nscal(10),first

      logical oneset
      integer iset
      common/setblk/oneset,iset

      logical logfit
      common/logf/logfit

      logical discprt
      common/dp/discprt
      logical debprt
      common/deb1/debprt

      integer idebug,idebug1
      common/deb/idebug,idebug1

      logical abort
      common/abt/abort

      logical nodata
      common/ndata/nodata

      logical penalty
      real*8 penfunc,penfac
      common/pen/penalty,penfunc,penfac

      real*8 assmax,ratemax
      integer icdep
      common/amax/assmax,icdep(200),ratemax

      integer nab,jalpha,jbeta
      common/absave/nab,jalpha,jbeta

      integer ireset
      common/ires/ireset
      real*8 perfac
      common/pert/perfac

      real*8 absmin,thmin,step
      common/reset/absmin,thmin(200),step(200)

      integer idebugt,itick,itlast,ndeb
      common/timer/idebugt,itick,itlast,ndeb

      integer ncdep,IX,JX
      real*4 x
      common/cpar/ncdep,IX(100),JX(100),x

      integer nligsav,IL
      common/LIG/nligsav,IL(100)

      integer nfix,jfix
      common/FIXBLK/nfix,jfix(100)

      integer NEQ,IE,JE,IF,JF
      real*4 EFAC
      common/EBLK/NEQ,IE(200),JE(200),IF(200),JF(200),EFAC(200)

c ---- the Q matrix handed to the stubbed QSET_HJC
      integer kg
      real*8 QGIVEN
      common/qdgiven/QGIVEN(100,100),kg

c ---- locals
      integer k,i,j,n,ichs,iu
      real*8 cnc,tr,tc,qq,tt,aa
      character*256 fname

      iu=11
      call GETARG(1,fname)
      if(fname.eq.' ') fname='likin.txt'
      open(iu,file=fname,status='old')

      read(iu,*) k,kA,kB,kC,kD
      read(iu,*) cnc
      read(iu,*) tr
      read(iu,*) tc
      read(iu,*) ichs

      kg=k
      do i=1,k
         do j=1,k
            read(iu,*) qq
            QGIVEN(i,j)=qq
         enddo
      enddo

      read(iu,*) n
      if(n.gt.ndmx) stop 'too many intervals'
      do i=1,n
         read(iu,*) tt,aa
         tint(i,1)=sngl(tt*1.0d3)      ! HJCFIT holds intervals in ms
         ampl(i,1)=sngl(aa)
         iprops(i,1)=0
      enddo
      close(iu)

c ---- settings
      nset=1
      jsetlast=0
      conc(1,1)=sngl(cnc)
      nlig=1
      nligsav=1
      IL(1)=1

      Nint(1)=n
      tresd(1)=tr
c hjclik.for line 1009 tests "tint(in,jset).gt.tcrit(jset)" and is NOT
c guarded by burst(jset), so a tcrit of zero ends the group at every shut
c time and the record becomes one group per opening. HJCFIT means "do not
c divide" by setting tcrit enormous -- 3.1536e10 ms, one year
c (Hjcfit1-09122003.for line 1568) -- and the same convention is used here.
      if(tc.gt.0.0d0) then
         tcritc(1)=sngl(tc*1.0d3)      ! ms, as HJCFIT holds it
         burst(1)=.true.
      else
         tcritc(1)=3.1536e10           ! one year: no division
         burst(1)=.false.
      endif
      chsvec(1)=ichs.eq.1
      badend(1)=.false.

      first=.true.
      do i=1,10
         ngp(i)=0
         an(i)=0.0
         nscal(i)=0
      enddo

c oneset true makes HJCLIK print the asymptotic areas and the HJC mean open
c and shut times. Useful once -- it is how the missed-events machinery was
c checked against SCALCS, in scratch/asymptotic_areas_check.py -- but noise
c when the likelihood is being called in a loop.
      oneset=.false.
      iset=1
      logfit=.false.
      discprt=.false.
      debprt=.false.
      idebug=0
      idebug1=0
      abort=.false.
      nodata=.false.
      penalty=.false.
      penfunc=0.0d0
      penfac=0.0d0
      assmax=1.0d12
      ratemax=1.0d12
      do i=1,200
         icdep(i)=0
         thmin(i)=1.0d0
         step(i)=0.0d0
         theta(i)=0.0d0
      enddo
      absmin=0.0d0
c nab > kab, so the alpha2/beta2 recording branch is never taken
      nab=999
      jalpha=1
      jbeta=1
      ireset=0
      perfac=0.0d0
      idebugt=0
      itick=0
      itlast=0
      ndeb=0
      ncdep=0
      x=0.0
      nfix=0
      NEQ=0
      npar=0

c ---- call it
      kfit=0
      nd1=ndmx
      nd2=1
      kab=1
      kAm=kA
      kFm=kB+kC
      km=100

      value=HJCLIK(kfit,theta,
     & tint,ampl,iprops,nd1,nd2,
     & Z00A,Z10A,Z11A,Z00F,Z10F,Z11F,
     & XAF,XFA,QEXPQA,QEXPQF,
     & alpha2,beta2,aloglik,kab,
     & kAm,kFm,km)

      write(*,'(a,1x,g24.16,1x,i8,1x,i3)') 'HJCLIK',value,n,k
      end
