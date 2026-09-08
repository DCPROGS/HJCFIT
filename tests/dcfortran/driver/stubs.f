c Stubs for the DOS/Lahey library routines HJCLIK and its callees reach for.
c
c None of these affects a likelihood: they are a clock, a bell, a keyboard
c poll, an underflow mode switch and a random number generator used only on
c the parameter-resetting path that kfit=0 disables. They exist so the
c likelihood can be linked outside the DOS program.

      subroutine DCTIMER(itick)
      integer itick
      itick=0
      return
      end

      subroutine BELL(n)
      integer n
      return
      end

      subroutine UNDER0(mode)
c Lahey's underflow-to-zero switch. gfortran flushes denormals to zero under
c -ffpe-summary=none anyway, and underflow in these matrices only ever moves
c a quantity that is already negligible.
      character*(*) mode
      return
      end

      logical function CAPLOCK()
      CAPLOCK=.false.
      return
      end

      logical function KBHIT()
      KBHIT=.false.
      return
      end

      character function GETCH(ktype)
      integer*2 ktype
      ktype=0
      GETCH=char(0)
      return
      end

      real*8 function DRANDOM()
c Only reached on the reset-and-perturb path, which kfit=0 disables. Fixed
c value rather than a stream, so that any accidental use is obvious.
      DRANDOM=0.5d0
      return
      end

      subroutine RANPERT(theta,thmin,kfit,perfac)
c Argument order matches the original (RANPERT.FOR, and the call in
c hjclik.for line 533): theta, thmin, kfit, perfac. Reached only when the
c root finder reports nerr 7 or 8, which the driver treats as a failure
c rather than something to perturb away from.
      integer kfit
      real*8 theta(*),thmin(*),perfac
      return
      end

      subroutine CLS()
      return
      end

      subroutine LOCATE(i,j)
      integer i,j
      return
      end

      subroutine CLRKB()
      return
      end

c ---------------------------------------------------------------------
c Three routines from the external "Spindrift" utility library, which is
c not in DCFORTRAN. All three are used only for formatting printed
c titles, on branches (PDFOUTD, PDFOUTS, and the idebug printouts) that
c the likelihood path never takes.

      integer function NBLANK(s)
c Position of the last non-blank character; 0 for an all-blank string.
c Implemented properly rather than stubbed, since it is cheap and a
c wrong value could make a caller index outside a string.
      character*(*) s
      integer i
      NBLANK=0
      do i=len(s),1,-1
         if(s(i:i).ne.' ') then
            NBLANK=i
            return
         endif
      enddo
      return
      end

      function CHARNB(s)
c The original returns the string trimmed. Its callers (PDFOUTD.FOR line
c 58, PDFOUTS.FOR line 38) do not declare it, so by implicit typing it is
c REAL there and the result is only ever handed to a list-directed
c print. Returning a real keeps the types consistent with those call
c sites; the routine is unreachable on the likelihood path.
      character*(*) s
      real CHARNB
      CHARNB=0.0
      return
      end
