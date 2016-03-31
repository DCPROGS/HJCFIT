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

#ifndef DCPROGS_LIKELIHOOD_APPROX_SURVIVOR_H
#define DCPROGS_LIKELIHOOD_APPROX_SURVIVOR_H

#include <DCProgsConfig.h>

#include <functional>
#include <memory>

#include "asymptotes.h"
#include "root_finder.h"

namespace DCProgs {

  //! \brief Implementation of asymptotic missed events.
  //! \details This object merely puts together two Asymptotes objects.
  class MSWINDOBE ApproxSurvivor {
    public:
      //! Type of the function used to minimize roots.
      typedef std::function<std::vector<Root>(DeterminantEq const &)> t_RootFinder;
      //! Initializes approximate survivor functor.
      //! \param[in] _af: Determinant equation for open->shut transitions
      //! \param[in] _roots_af: Roots of _af equation
      //! \param[in] _fa: Determinant equation for shut->open transitions
      //! \param[in] _roots_fa: Roots of _fa equation
      ApproxSurvivor(DeterminantEq const &_af, std::vector<Root> const &_roots_af, 
                     DeterminantEq const &_fa, std::vector<Root> const &_roots_fa );
      //! Initializes approx survivor functor.
      //! \param[in] _matrix: Transition matrix
      //! \param[in] _tau: resolution/max length missed events
      //! \param[in] _findroots: A functor with which to find all roots.
      //!                        This function should take a DeterminantEq as its sole argument and
      //!                        return a std::vector<DCProgs::RootInterval>
      ApproxSurvivor(QMatrix const &_matrix, t_real _tau, t_RootFinder const &_findroots);
      //! Move constructor
      ApproxSurvivor   (ApproxSurvivor &&_c) 
                     : asymptotes_af_(std::move(_c.asymptotes_af_)),
                       asymptotes_fa_(std::move(_c.asymptotes_fa_)) {}
      //! \brief Initializes missed-events functor.
      //! \param[in] _qmatrix Transition matrix
      //! \param[in] _tau resolution/max length missed events
      //! \param[in] _xtol Tolerance for interval size
      //! \param[in] _rtol Tolerance for interval size. The convergence criteria is an affine
      //!            function of the root:
      //!   \f$x_{\mathrm{tol}} + r_{\mathrm{tol}} x_{\mathrm{current}} = \frac{1}{2}|x_a - x_b|\f$.
      //! \param[in] _itermax maximum number of iterations for any of the three steps.
      //! \param[in] _lowerbound Lower bound of the interval bracketing all roots. If None, the
      //!            lower bound is obtained from find_lower_bound_for_roots().
      //! \param[in] _upperbound Upper bound of the interval bracketing all roots. If None, the
      //!            upper bound is obtained from find_upper_bound_for_roots().
      ApproxSurvivor( QMatrix const &_qmatrix, t_real _tau, t_real _xtol=1e-12, t_real _rtol=1e-12,
                      t_uint _itermax=100, t_real _lowerbound=quiet_nan, 
                      t_real _upperbound=quiet_nan );

      //! Open to close transitions 
      t_srmatrix af(t_real _t) const { return asymptotes_af_->operator()(_t); }
      //! Close to open transitions
      t_srmatrix fa(t_real _t) const { return asymptotes_fa_->operator()(_t); }
      //! Number of exponential components for af
      t_uint nb_af_components() const { return static_cast<t_uint>(asymptotes_af_->size()); }
      //! Number of exponential components for fa
      t_uint nb_fa_components() const { return static_cast<t_uint>(asymptotes_fa_->size()); }
      //! AF exponential components
      Asymptotes::t_MatrixAndRoot const & get_af_components(t_int i) const {
        if(i < 0) i += nb_af_components();
        if(i < 0 or static_cast<t_uint>(i) >= nb_af_components())
          throw errors::Index("AF component index out of range.");
        return (*asymptotes_af_)[i]; 
      }
      //! FA exponential components
      Asymptotes::t_MatrixAndRoot const & get_fa_components(t_int i) const {
        if(i < 0) i += nb_fa_components();
        if(i < 0 or static_cast<t_uint>(i) >= nb_fa_components())
          throw errors::Index("FA component index out of range.");
        return (*asymptotes_fa_)[i]; 
      }
    protected:
#     ifndef HAS_CXX11_UNIQUE_PTR
        //! Type of the pointers holding recursion interfaces.
        typedef std::auto_ptr<Asymptotes> t_AsymptotesPtr;
#     else
        //! Type of the pointers holding recursion interfaces.
        typedef std::unique_ptr<Asymptotes> t_AsymptotesPtr;
#     endif
      //! Pointer to AF recursion interface
      t_AsymptotesPtr asymptotes_af_;
      //! Pointer to FA recursion interface
      t_AsymptotesPtr asymptotes_fa_;
  };

  //! Dumps Survivor function equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, ApproxSurvivor const &_self);
}
#endif 
