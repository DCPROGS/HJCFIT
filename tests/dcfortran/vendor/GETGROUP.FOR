	subroutine GETGROUP(xobs,yobs,w,jset,nj,niobs,njset,nxval,
     &  xgrp,ybar,sdy,sdm,ngrp,SSwg,ndf,SStot,SScol,ymean,n0)
c
C===check that OK if ALL values in a group are omitted
c Subroutine for CVFIT to
c (a) Check for repeated X values in Xval(i,jset) i.e.
c	for replicate Y measurements
c (b) If found, calculate the 'within group' error sum of squares and
c	return it as SSwg, dfwg=degrees of freedom -now all values are
c	calculated with weights in w(i,j)
c (c) Usually, if there are replicate Y values in data at given x, all
c 	the w(i,j) will be 1.0. But not necessarily -eg could specify arbitrary
c 	weights, or might even have legitimate weights. But hard to imagine
c	a case in which different w(i,j) would be legitimate at the same x value
c     NB In fact the SD of the group cannot be defined unless SD is the same
c     for all of the group (all Y at a given x)
c	This subroutine provides an SD for the replicates at each x value
c	and this could be used as weight (same for all Y at given x; NB since
c	all replicates are in data we want s(y) as weight for each y, not s(ybar)).
c	This is then the same as taking the mean before entering data, and
c	using the SDM as weight for the mean, or SDY (defined here) for the
c	weight attached to a single Y value
c
c Indexes
c  i=1,..,nj(jset): for all obs in the current set (as in Yobs(i,jset), w(i,jset)
c  k=1,..,nxval: for sum over groups (nxval=# of different X values)
c
c This assumes that x,y,w arrays have been sorted into ascending order
c of X values BEFORE this subroutine is called (as is done in
c latest CVFIT)
c NB call with nxval(jset), or with scalar e.g. nxv

	real*4 xobs(niobs,njset)
	real*4 yobs(niobs,njset)
	real*4 w(niobs,njset)
	integer nj(njset)
	real*4 xgrp(niobs,njset),ybar(niobs,njset),sdy(niobs,njset)
	real*4 sdm(niobs,njset)
	integer ngrp(niobs,njset)
	allocatable:: igrp,ngrp0
	integer igrp(:),ngrp0(:)
c
	allocatable wk		!for weights for kth group, i=1,nj(j)
	allocatable wkj		!for weights for kth group, k=1,nxval
	real*4 wk(:),wkj(:)
	allocatable omit		!to mark obs with zero weight for omission
	logical omit(:)
c

	n=nj(jset)
	n0=0		!counts excluded zero-weight values
c
	ALLOCATE(igrp(n),ngrp0(n))
	ALLOCATE(wk(n),wkj(n),omit(n))
	nxval=1
c
c Locate the groups of Y values at same X
c Define igrp(j)=index of 1st point in each group
c Note that this is defined regardless of whether some values have zero weight
	igrp(1)=1	!index for 1st member of group 1
	xgrp(1,jset)=Xobs(1,jset)	!x value for each group
	do i=2,n
	   if(Xobs(i,jset).ne.Xobs(i-1,jset)) then
		nxval=nxval+1
		igrp(nxval)=i	!1st grp=1,..igrp(1)
		xgrp(nxval,jset)=Xobs(i,jset)
	   endif
	enddo
c Calc number in each group
	do k=1,nxval-1
	   ngrp(k,jset)=igrp(k+1)-igrp(k)
	enddo
	ngrp(nxval,jset)=n - igrp(nxval) + 1
c
c Calculate unweighted mean of all Y, and set omit(i)
	s=0.0		!grand total
	do i=1,n
	   wt=w(i,jset)
	   if(wt.gt.1.e-37) then
		omit(i)=.false.
		y=Yobs(i,jset)
		s=s + y
	   else
		omit(i)=.true.
		n0=n0+1
	   endif
	enddo
	n1=n-n0
	if(n1.ge.1) then
	   ymean=s/float(n1)		!overall unweighted mean
	else
	   ymean=0.0
	endif
c
c Find mean of each group
c (all members of group have same SD so no weights, except weight=0 used to
c exclude values)
	do k=1,nxval
	   s=0.0
	   i1=igrp(k)
	   i2=i1+ngrp(k,jset)-1
	   nj0=0		!count zero weights in current group
	   do i=i1,i2
		wt=w(i,jset)
		if(omit(i)) then
		   nj0=nj0+1
		else
		   s=s + Yobs(i,jset)
		endif
	   enddo
	   ngrp0(k)=nj0	!record number of zero weight values in jth group
	   nj1=ngrp(k,jset)-nj0
	   if(nj1.ge.1) then
		enj=float(nj1)
		ybar(k,jset)=s/float(nj1)	!get weighted mean
	   else
		ybar(k,jset)=0.
		sdy(k,jset)=1.e36		!i.e.zero weight (see cvin2)
	   endif
	enddo
c
c Find SD of each group, and pooled SS within groups
c (all members of group have same SD so no weights, except weight=0 used to
c exclude values)
	SSwg=0.0
	do k=1,nxval
	   s=0.0
	   i1=igrp(k)
	   i2=i1+ngrp(k,jset)-1
	   nj0=0		!count zero weights in current group
	   do i=i1,i2
		if(omit(i)) then
		   nj0=nj0+1
		else
		   s=s + (Yobs(i,jset) - ybar(k,jset))**2
		endif
	   enddo
	   SSwg=SSwg + s
	   if(ngrp(k,jset)-nj0.gt.1) then
		en=float(ngrp(k,jset)-nj0)
		s=s/(en-1.)				!variance of y
		if(s.ge.0.0) then
		   sdy(k,jset)=sqrt(s)
		   sdm(k,jset)=sqrt(s/en)	!variance of ybar
		   wk1=1.0/s
		else
		   sdy(k,jset)=-2.0
		   sdm(k,jset)=-2.0
		   wk1=0.0
		endif
	   else
		sdy(k,jset)=-1.0
		sdm(k,jset)=-1.0
		wk1=0.0
	   endif
	   wkj(k)=wk1
	   do i=i1,i2
		if(omit(i)) then
		   wk(i)=0.0
		else
		   wk(i)=wk1	!weight same for each Y in the group
		endif
	   enddo
	enddo

c Now have SD for each group
c Can next calculate various SSD with or without these SD as weights
c
c Calculate weighted overall mean
	s=0.0		!grand total
	sw=0.0
	n0=0
	do i=1,n
	   if(omit(i)) then
		n0=n0+1
	   else
		s=s + wk(i)*Yobs(i,jset)
		sw=sw + wk(i)
	   endif
	enddo
	if(sw.gt.0.) then
	   ymeanw=s/sw
	else
	   ymeanw=0.0
	endif
c SStot (unweighted)
	SStot=0.0
	do i=1,n
	   if(.not.omit(i)) then
		y=Yobs(i,jset)
		SStot=SStot+(y-ymean)**2
	   endif
	enddo
c
c SScol unweighted
c (could have a whole column missing if weights in main prog were set
c to 0 for all the Y at one X value -in this case wkj(k) should be 0
	SScol=0.0
	do k=1,nxval
	   enj=float(ngrp(k,jset)-ngrp0(k))
	   if(wkj(k).gt.0.) then
		SScol=SScol + enj*((ybar(k,jset)-ymean)**2)	!ssd between cols
	   endif
	enddo
c
c NOW WEIGHTED VALUES =TO BE REVISED
cc Calculate weighted overall mean
c
cc Now calculate weighted total SSD
c	SStotw=0.0
c	do k=1,nxval
c	   sw=0.0
c	   sd=sdy(k,jset)
c	   if((sd.le.0.).or.(sd.gt.3.e18)) then
c		===no values for this group
c	   else
c		wk=1.0/(sd*sd)	!same for all x in the group=weight fot kth set
c		i1=igrp(k)
c		i2=i1+ngrp(k,jset)-1
c		nj0=0		!count zero weights in current group
c		do i=i1,i2
c		   wt=w(i,jset)
c		   if(wt.gt.1.e-37) then
c			y=Yobs(i,jset)
c===need ymeanw, not ymean here?
c			SStotw=SStotw+wk*(y-ymean)**2
c			sw=sw+wk
c		   else
c			nj0=nj0+1
c		  endif
c		enddo
c		ngrp0(k)=nj0	!record number of zero weight values in jth group
c	   endif
c	enddo
c==start of old
c Calc mean for each group
c
c	do j=1,nxval
c	   ybar(j,jset)=0.0
c	   sw=0.0
c	   i1=igrp(j)
c	   i2=i1+ngrp(j,jset)-1
c	   nj0=0		!count zero weights in current group
c	   do i=i1,i2
c		wt=w(i,jset)
c		if(wt.gt.1.e-37) then
c		   sw=sw + wt
c		   ybar(j,jset)=ybar(j,jset) + wt*Yobs(i,jset)
c		else
c		   nj0=nj0+1
c		endif
c	   enddo
c	   ngrp0(j)=nj0	!record number of zero weight values in jth group
c	   nj1=ngrp(j,jset)-nj0
c	   if(nj1.ge.1) then
c		enj=float(nj1)
c		ybar(j,jset)=ybar(j,jset)/sw	!get weighted mean
c	   else
c		ybar(j,jset)=0.
c		sdy(j,jset)=1.e36		!i.e.zero weight (see cvin2)0
c	   endif
cc
cc   NB w(j,jset) is based on SD for each group (of replicate y at same x)
cc   so have one value for each x, but it is not defined yet -done below!!!
cc==	   SScol=SScol + w(j,jset)*enj*((ybar(j,jset)-ymean)**2)	!ssd between cols
c	enddo
c
c
c If there are replicate Y values at same x, use then to calculate
c (a) the sdy for each group
c	 sdy(j)=SUM(j=1,nxval)[Sum(i=1,n(j)){w(j)[(y(i,j)-ybar(j)]**2]}/df
c (b) the 'within X' estimate of error.
c 	 SSwg=[Sum(i=1,n(j))[y(i,j)-ybar(j)]**2]/(n(j)*(n(j)-1)]
c
c	SSwg=0.0
c	do j=1,nxval
c	   s=0.0
c	   i1=igrp(j)
c	   i2=i1+ngrp(j,jset)-1
c	   sw=0.0
c	   nj0=0		!count zero weights in current group
c	   do i=i1,i2
c		wt=w(i,jset)
c		if(wt.gt.1.e-37) then
c		   sw=sw + wt
c		   s2=wt*((Yobs(i,jset) -ybar(j,jset))**2)
c		   s=s + s2
c		   SSwg=SSwg + s2
c		else
c		   nj0=nj0+1
c		endif
c	   enddo
c	   if(ngrp(j,jset)-nj0.gt.1) then
c		en=float(ngrp(j,jset)-nj0)
c		s=s/(en-1)		!variance of y
c		if(s.ge.0.0) then
c		   sdy(j,jset)=sqrt(s)
c		   sdm(j,jset)=sqrt(s/en)
c		else
c		   sdy(j,jset)=-2.0
c		   sdm(j,jset)=-2.0
c		endif
c	   else
c		sdy(j,jset)=-1.0
c		sdm(j,jset)=-1.0
c	   endif
c	enddo
c
c
c  Finally get SScol (between x values) using the sdy() as weights, as
c long as they are defined
c	SScol=0.0
c	do j=1,nxval
c	   sj=sdy(j,jset)
c	   if(sj.gt.1.e-30.and.sj.lt.3.e18) then
c		wj=1.0/(sj*sj)
c		SScol=SScol + wj*enj*((ybar(j,jset)-ymean)**2)	!ssd between cols
c	   endif
cc If any zero weight values are omitted, adjust number in group accordingly
c	   ngrp(j,jset)=ngrp(j,jset)-ngrp0(j)
c	enddo
cc
cc===========at the end can calulate weighted mean using the SD for each set
c===========of Y values at each x
c=====
c
c	sw=0.0
c	do i=1,n
c	   wt=w(i,jset)
c	   if(wt.gt.1.e-37) then
c		y=Yobs(i,jset)
c		s=s + wt*y
c		sw=sw + wt
c	   else
c		n0=n0+1
c	   endif
c	enddo
c	ymean=s/sw		!overall weighted mean
c====
c
c Define df for within group error
	ndf=0
	do j=1,nxval
	   nd1=ngrp(j,jset) - ngrp0(j)	!don't count zero weight values
	   ndf=ndf + nd1-1	     		!n-1 df for each group
	enddo
c
c Before returning, reset ngrp(k) to reflect omitted values
c If any zero weight values are omitted, adjust number in group accordingly
	do k=1,nxval
	   ngrp(k,jset)=ngrp(k,jset)-ngrp0(k)
	enddo
	DEALLOCATE(igrp,ngrp0,wk,wkj,omit)
c
	RETURN
	end

