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

#ifndef DCPROGS_ROOT_FINDER_H
#define DCPROGS_ROOT_FINDER_H

#include "DCProgsConfig.h"

#include <vector>

#include "determinant_equation.h"

namespace DCProgs {


  //! Defines an interval which one or more roots.
  struct MSWINDOBE RootInterval {
    //! Start of the interval (lower value).
    t_real start;
    //! End of the interval (max value).
    t_real end;
    //! Number of roots in interval.
    //! This is likely just a guess.
    t_uint multiplicity;

    //! Constructor.
    RootInterval   (t_real _start, t_real _end, t_uint _mult) noexcept
                 : start(_start), end(_end), multiplicity(_mult) {}
    //! Default Constructor.
    RootInterval() noexcept : start(0), end(0), multiplicity(0) {};
  };

  //! Defines a root, including mutliplicity.
  struct MSWINDOBE Root {
    //! Start of the interval (lower value).
    t_real root;
    //! Number of roots in interval.
    t_uint multiplicity;

    //! Constructor.
    Root(t_real _root, t_uint _mult) noexcept : root(_root), multiplicity(_mult) {}
    //! Default Constructor.
    Root() noexcept : root(0), multiplicity(0) {};
  };

  //! \brief Figures out interval where roots can be found.
  //! \param[in] _det The determinantal equation
  //! \param[in] _mins A valid lower bound, or DCProgs::quiet_nan. In the latter case, the lower
  //!            bound is determined using DCProgs::find_lower_bound_for_roots().
  //! \param[in] _maxs  A valid upper bound, or DCProgs::quiet_nan. In the latter case, the upper
  //!            bound is determined using DCProgs::find_upper_bound_for_roots().
  //! \param[in] _tolerance Minimum size of intervals. Below that, roots are expected to be
  //!            multiples.
  std::vector<RootInterval> MSWINDOBE
    find_root_intervals( DeterminantEq const &_det, t_real _mins=quiet_nan, t_real _maxs=quiet_nan,
                         t_real _tolerance=1e-8 );

  //! \brief Figures out an lower bound for root finding.
  //! \details Proceeds by computing the eigenvalues, then setting lower bound to somewhat lower than
  //!          the lowest eigenvalue. It then checks that the eigenvalues of the matrix computed at
  //!          that value, and so on and so forth. The algorithm stops when the lowest eigenvalue is
  //!          higher than the current bound.
  //! \param[in] _det The determinantal equation
  //! \param[in] _start Value where to start looking for lower bound.
  //! \param[in] _alpha factor by which to set new lower bound:
  //!            \f$s_{n+1} = min(\epsilon_i) + \alpha (s_N - min(\epsilon_i))\f$.
  //! \param[in] _itermax Maximum number of iterations.
  t_real MSWINDOBE find_lower_bound_for_roots(DeterminantEq const &_det, t_real _start=0e0,
                                              t_real _alpha=5e0, t_uint _itermax=100);

  //! \brief Figures out an upper bound for root finding.
  //! \param[in] _det The determinantal equation
  //! \param[in] _start Value where to start looking for lower bound.
  //! \param[in] _alpha factor by which to set new lower bound:
  //!            \f$s_{n+1} = min(\epsilon_i) + \alpha (s_N - min(\epsilon_i))\f$.
  //! \param[in] _itermax Maximum number of iterations.
  t_real MSWINDOBE find_upper_bound_for_roots(DeterminantEq const &_det, t_real _start=0e0,
                                              t_real _alpha=5e0, t_uint _itermax=100);
  //! \brief Finds roots via brute force search
  //! \details Computes all values between mins and maxs, for a given resolution.
  //!          If determinant changes sign between two values, or if it comes to within tolerance of
  //!          zero, then computes eigenvalues of H to determine possible multiplicity.
  //! \param[in] _det The determinantal equation
  //! \param[in] _resolution resolution at which computes values in interval.
  //! \param[in] _mins A valid lower bound, or DCProgs::quiet_nan. In the latter case, the lower
  //!            bound is determined using DCProgs::find_lower_bound_for_roots().
  //! \param[in] _maxs  A valid upper bound, or DCProgs::quiet_nan. In the latter case, the upper
  //!            bound is determined using DCProgs::find_upper_bound_for_roots().
  //! \param[in] _tolerance Tolerance below which the value of the determinant is considered
  //!            "close to zero".
  std::vector<RootInterval> MSWINDOBE 
    find_root_intervals_brute_force(DeterminantEq const &_det, 
                                    t_real _resolution = 1e-1,
                                    t_real _mins = quiet_nan,
                                    t_real _maxs = quiet_nan,
                                    t_real _tolerance = 1e-1);

  //! \brief Finds root using brentq and find_root_intervals.
  //! \details Tries and computes the roots of an input determinantal equation.
  //! This is a three fold process:
  //!
  //! 1. Find an interval that contains all roots/eigenvalues 
  //! 1. Split it into smaller intervals until each contains only one root
  //! 1. Optimize over this interval to find an accurate position for the root
  //!
  //| If _upperbound is larger than _lowerbound, then those bounds are used for the first
  //| step. Otherwise, DCProgs will try and determine the lower bound by iteratively checking
  //| the lowest eigenvalue :math:`\\\\epsilon_i^s` of :math:`H(s_i)`, where :math:`s_i` is the
  //| guess at iteration :math:`i`. If the lower eigenvalue is lower than :math:`s_i`, than
  //| :math:`s_{i+1} = \\\\epsilon_i^s + \\\\alpha(\\epsilon_i^s - s_i)` is created. The upper
  //| bound is obtained in a similar fashion, with an additional step where it is iteratively
  //| increased until the determinant is positive
  //
  //| The second step bisects the interval from the first step until intervals backeting a
  //| single eigenvalue is obtained, or until a resilution limit is achieved. In the latter
  //| case, the multiplicity is set to the number of eigenvalues in the interval.
  //|
  //! The last step is carried out by brentq().
  //!
  //! \param[in] _det The determinantal equation for which to compute the roots.
  //! \param[in] _resolution Size of sieve with which to look for roots.
  //! \param[in] _mins A valid lower bound, or DCProgs::quiet_nan. In the latter case, the lower
  //!            bound is determined using DCProgs::find_lower_bound_for_roots().
  //! \param[in] _maxs  A valid upper bound, or DCProgs::quiet_nan. In the latter case, the upper
  //!            bound is determined using DCProgs::find_upper_bound_for_roots().
  //! \param[in] _tolerance Tolerance to use for figuring out a root is in a sieve hole.
  std::vector<Root> MSWINDOBE find_roots( DeterminantEq const &_det, 
                                          t_real _xtol = 1e-8,
                                          t_real _rtol = 1e-8,
                                          t_uint _itermax = 100,
                                          t_real _lowerbound = quiet_nan,
                                          t_real _upperbound = quiet_nan );
}

#endif
