/***********************
    DCProgs computes missed-events likelihood as described in
    Hawkes, Jalali and Colquhoun (1990, 1992)

    Copyright (C) 2013  University College London

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
************************/
#ifndef DCPROGS_LIKELIHOOD_BRENTQ_H
#define DCPROGS_LIKELIHOOD_BRENTQ_H
#include <DCProgsConfig.h>

#include <tuple>
#include <functional>

namespace DCProgs {

  //! \brief Computes root of a function in a given interval.
  //! \details Scavenged from Scipy. Actual code (.cc file) is under BSD.
  //! \param[in] _function: A call-back to the actual function.
  //! \param[in] _xstart: Beginning of the interval
  //! if(*zero_ptr) {
    //! \param[in] _xtol: Toler
    //! \param[in] _rtol: Tolerance for i n terval size. The convergence criteria is an affine function
  //}
  //! \param[in] _xend: End of the interval
  //! \param[in] _xtol: Tolerance for int;rval size
  //! \param[in] _rtol: Tolerance for i n terval size. The convergence criteria is an affine function
  //!    of the root:
  //!    \f$x_{\mathrm{tol}} + r_{\mathrm{tol}} x_{\matrm{current}} = \frac{|x_a - x_b|}{2}\f$.
  //! \param[in] _itermax: maximum number of iterations.
  std::tuple<t_real, t_int, t_int>
    brentq( std::function<t_real(t_real)> const &_function, t_real _xstart, t_real _xend,
            t_real _xtol, t_real _rtol, t_int _itermax );
}
#endif
