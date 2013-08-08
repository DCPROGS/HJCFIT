/* Copyright (c) 2001, 2002 Enthought, Inc.
 * All rights reserved.
 *
 * Copyright (c) 2003-2012 SciPy Developers.
 * All rights reserved.
**/

/* Written by Charles Harris charles.harris@sdl.usu.edu 
 * Taken from Scipy and adapted for use in dcprogs.
 * */

#include <DCProgsConfig.h>
#include "errors.h"
#include "brentq.h"


/*

  At the top of the loop the situation is the following:

    1. the root is bracketed between _xstart and _xend
    2. _xstart is the most recent estimate
    3. xp is the previous estimate
    4. |fp| < |fb|

  The order of _xstart and xp doesn't matter, but assume xp < _xend. Then _xstart lies to
  the right of xp and the assumption is that _xstart is increasing towards the root.
  In this situation we will attempt quadratic extrapolation as long as the
  condition

  *  |fa| < |fp| < |fb|

  is satisfied. That is, the function value is decreasing as we go along.
  Note the 4 above implies that the right inequlity already holds.

  The first check is that _xstart is still to the left of the root. If not, _xend is
  replaced by xp and the interval reverses, with _xend < _xstart. In this situation
  we will try linear interpolation. That this has happened is signaled by the
  equality _xend == xp;

  The second check is that |fa| < |fb|. If this is not the case, we swap
  _xstart and _xend and resort to bisection.

*/

namespace DCProgs {

  MSWINDOBE std::tuple<t_real, t_uint, t_uint>
    brentq( std::function<t_real(t_real)> const &_function, 
            t_real _xstart, t_real _xend,
            t_real _xtol, t_real _rtol, t_uint _itermax ) {

    t_real xpre = _xstart, xcur = _xend;
    t_real xblk = 0.0, fpre, fcur, fblk = 0.0, spre = 0.0, scur = 0.0, sbis, tol;
    t_real stry, dpre, dblk;
    t_uint i;

    t_uint function_calls = 0, iterations = 0;

    fpre = _function(xpre);
    fcur = _function(xcur);
    function_calls = 2;
    if (fpre*fcur > 0) throw errors::Mass("Interval does not bracket a root.");
    if (fpre == 0) return std::make_tuple(xpre, iterations, function_calls);
    if (fcur == 0) return std::make_tuple(xcur, iterations, function_calls);
    iterations = 0;
    for(i = 0; i < _itermax; i++) {
        ++iterations;
        if (fpre*fcur < 0) {
          xblk = xpre;
          fblk = fpre;
          spre = scur = xcur - xpre;
        }
        if (std::abs(fblk) < std::abs(fcur)) {
          xpre = xcur; xcur = xblk; xblk = xpre;
          fpre = fcur; fcur = fblk; fblk = fpre;
        }

        tol = _xtol + _rtol*std::abs(xcur);
        sbis = (xblk - xcur)/2;
        if (fcur == 0 || std::abs(sbis) < tol)
          return std::make_tuple(xcur, iterations, function_calls);

        if (std::abs(spre) > tol && std::abs(fcur) < std::abs(fpre)) {
          if (xpre == xblk) {
              /* interpolate */
              stry = -fcur*(xcur - xpre)/(fcur - fpre);
          }
          else {
              /* extrapolate */
              dpre = (fpre - fcur)/(xpre - xcur);
              dblk = (fblk - fcur)/(xblk - xcur);
              stry = -fcur*(fblk*dblk - fpre*dpre)
                     /(dblk*dpre*(fblk - fpre));
          }
          if (2*std::abs(stry) < std::min(std::abs(spre), 3*std::abs(sbis) - tol)) {
              /* good short step */
              spre = scur; scur = stry;
          } else {
              /* bisect */
              spre = sbis; scur = sbis;
          }
        }
        else {
          /* bisect */
          spre = sbis; scur = sbis;
        }

        xpre = xcur; fpre = fcur;
        if (std::abs(scur) > tol) xcur += scur;
        else xcur += (sbis > 0 ? tol : -tol);

        fcur = _function(xcur);
        ++function_calls;
    }
    throw errors::MaxIterations("Could not converge search for root in brentq.");
  }
}
