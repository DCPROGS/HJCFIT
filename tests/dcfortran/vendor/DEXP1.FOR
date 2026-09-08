	real*8 function dexp1(x)
c Double precision version of EXP1, to prevent floating underflow
c 05/13/02 08:24pm added common/dexp/ to allow easy test for when set
c to 0.0d0
c
	real*8 x
	logical zero
	common/dex/zero
c
	zero=.false.
	if(x.lt.-708.d0) then
		dexp1=0.0d0
		zero=.true.
	else
		dexp1=dexp(x)
	endif
	RETURN
	end


