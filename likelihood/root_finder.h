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
    t_int multiplicity;

    //! Constructor.
    RootInterval   (t_real _start, t_real _end, t_int _mult) noexcept
                 : start(_start), end(_end), multiplicity(_mult) {}
    //! Default Constructor.
    RootInterval() noexcept : start(0), end(0), multiplicity(0) {};
  };

  //! Defines a root, including mutliplicity.
  struct MSWINDOBE Root {
    //! Start of the interval (lower value).
    t_real root;
    //! Number of roots in interval.
    t_int multiplicity;

    //! Constructor.
    Root(t_real _root, t_int _mult) noexcept : root(_root), multiplicity(_mult) {}
    //! Default Constructor.
    Root() noexcept : root(0), multiplicity(0) {};
  };

  //! \brief Figures out interval where roots can be found.
  //! \details It starts from reasonable interval.
  //!          If that interval does not contain all eigenvalues of H, then it tries and increase
  //!          it until an overlapping interval is found.
  //!          It then proceeds by bisecting the interval until all roots are isolated.
  //!          If a bisection falls below a given convergence criteria, the root is deemed
  //!          degenerate.
  //! \param[in] _det: The determinantal equation
  //! \param[in] _mins: A valid lower bound. All roots should be above that lower bound. 
  //!                   If _mins > _maxs, then tries to determine the lower bound using
  //!                   find_lower_bound_for_root.
  //! \param[in] _maxs: A valid upper bound. All roots should be below that upper bound. 
  //! \param[in] _tolerance: Minimum size of intervals. Below that, roots are expected to be
  //!                        multiples.
  std::vector<RootInterval> MSWINDOBE
    find_root_intervals(DeterminantEq const &_det, t_real _mins = 1e8, t_real _maxs = 1e1,
                        t_real _tolerance = 1e-8);

  //! \brief Figures out an lower bound for root finding.
  //! \details Proceeds by computing the eigenvalues, then setting lower bound to somewhat lower than
  //!          the lowest eigenvalue. It then checks that the eigenvalues of the matrix computed at
  //!          that value, and so on and so forth. The algorithm stops when the lowest eigenvalue is
  //!          higher than the current bound. Then it checks that the determinant  is of the correct
  //!          sign.
  //! \param[in] _det: The determinantal equation
  //! \param[in] _start: Value where to start looking for lower bound.
  //! \param[in] _alpha: factor by which to set new lower bound:
  //!                    \f$s_{n+1} = min(\epsilon_i) + \alpha (s_N - min(\epsilon_i))\f$.
  //! \param[in] _itermax: Maximum number of iterations.
  t_real MSWINDOBE find_lower_bound_for_roots(DeterminantEq const &_det, t_real _start=0e0,
                                              t_real _alpha=5e0, t_int _itermax=100);

  //! \brief Figures out an upper bound for root finding.
  //! \param[in] _det: The determinantal equation
  //! \param[in] _start: Value where to start looking for lower bound.
  //! \param[in] _alpha: factor by which to set new lower bound:
  //!                    \f$s_{n+1} = min(\epsilon_i) + \alpha (s_N - min(\epsilon_i))\f$.
  //! \param[in] _itermax: Maximum number of iterations.
  t_real MSWINDOBE find_upper_bound_for_roots(DeterminantEq const &_det, t_real _start=0e0,
                                              t_real _alpha=5e0, t_int _itermax=100);
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
}

#endif
