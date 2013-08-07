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
  //! \param[in] _det: The determinantal equation
  //! \param[in] _mins: A valid lower bound. All roots should be above that lower bound. 
  //!                   If _mins > _maxs, then tries to determine the lower bound using
  //!                   find_lower_bound_for_root.
  //! \param[in] _maxs: A valid upper bound. All roots should be below that upper bound. 
  //! \param[in] _tolerance: Minimum size of intervals. Below that, roots are expected to be
  //!                        multiples.
  std::vector<RootInterval> MSWINDOBE
    find_root_intervals( DeterminantEq const &_det, t_real _mins=1e8, t_real _maxs=1e1,
                         t_real _tolerance=1e-8 );

  //! \brief Figures out an lower bound for root finding.
  //! \details Proceeds by computing the eigenvalues, then setting lower bound to somewhat lower than
  //!          the lowest eigenvalue. It then checks that the eigenvalues of the matrix computed at
  //!          that value, and so on and so forth. The algorithm stops when the lowest eigenvalue is
  //!          higher than the current bound.
  //! \param[in] _det: The determinantal equation
  //! \param[in] _start: Value where to start looking for lower bound.
  //! \param[in] _alpha: factor by which to set new lower bound:
  //!                    \f$s_{n+1} = min(\epsilon_i) + \alpha (s_N - min(\epsilon_i))\f$.
  //! \param[in] _itermax: Maximum number of iterations.
  t_real MSWINDOBE find_lower_bound_for_roots(DeterminantEq const &_det, t_real _start=0e0,
                                              t_real _alpha=5e0, t_uint _itermax=100);

  //! \brief Figures out an upper bound for root finding.
  //! \param[in] _det: The determinantal equation
  //! \param[in] _start: Value where to start looking for lower bound.
  //! \param[in] _alpha: factor by which to set new lower bound:
  //!                    \f$s_{n+1} = min(\epsilon_i) + \alpha (s_N - min(\epsilon_i))\f$.
  //! \param[in] _itermax: Maximum number of iterations.
  t_real MSWINDOBE find_upper_bound_for_roots(DeterminantEq const &_det, t_real _start=0e0,
                                              t_real _alpha=5e0, t_uint _itermax=100);
  //! \brief Finds roots via brute force search
  //! \details Computes all values between mins and maxs, for a given resolution.
  //!          If determinant changes sign between two values, or if it comes to within tolerance of
  //!          zero, then computes eigenvalues of H to determine possible multiplicity.
  //! \param[in] _det: The determinantal equation
  //! \param[in] _resolution: resolution at which computes values in interval.
  //! \param[in] _mins: A valid lower bound. All roots should be above that lower bound. 
  //!                   If _mins > _maxs, then tries to determine the lower bound using
  //!                   find_lower_bound_for_root.
  //! \param[in] _maxs: A valid upper bound. All roots should be below that upper bound. 
  //! \param[in] _tolerance: Tolerance below which the value of the determinant is considered
  //!                        "close to zero".
  std::vector<RootInterval> MSWINDOBE 
    find_root_intervals_brute_force(DeterminantEq const &_det, 
                                    t_real _resolution = 1e-1,
                                    t_real _mins = 1e8,
                                    t_real _maxs   = 0e0,
                                    t_real _tolerance = 1e-1);

  //! Finds root using brentq and find_root_intervals.
  std::vector<Root> MSWINDOBE find_roots( DeterminantEq const &_det, 
                                          t_real _xtol = 1e-8,
                                          t_real _rtol = 1e-8,
                                          t_uint _itermax = 100 );
}

#endif
