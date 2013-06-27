#include <DCProgsConfig.h>

#include <iostream>

#include "root_finder.h"

namespace DCProgs {


  // Some functions we need here and nowhere else.
  namespace {

     // Compute number of roots above  s.
     t_int getUpper(DeterminantEq const &_det, t_real _s) {
       // Compute matrix
       t_rmatrix const H(_det.H(_s));

       // checks that it is valid.
       // If condition below is true, then means at least one element is NaN.
       if((H.array() != H.array()).any())
         throw errors::Mass("Encountered NaN in H matrix when trying to figure out roots.");

       // Computes eigenvalues of midpoint.
       Eigen::EigenSolver<t_rmatrix> eigsolver(H);
       if(eigsolver.info() != Eigen::Success) 
         throw errors::Mass("Could not solve eigenvalue problem.");

       // Checks we have no complex eigenvalues.
       if((eigsolver.eigenvalues().array().imag().abs() > 1e-8).any())
         throw errors::ComplexEigenvalues("when computing interval for roots.");

       t_rvector const eigs_real = eigsolver.eigenvalues().real();

       // compute number of roots
       return (eigs_real.array() >= _s).count();
     }

     // Actual bisecting algorithm.
     // Checks number of roots above and below midpoint. 
     // Based on that information, figures out whether to bisect some more or whether an interval
     // was found.
     void step_(DeterminantEq const &_det, 
                t_int const _nroots, t_real const _mins, 
                t_real const _maxs, t_real const _tolerance,
                std::vector<RootInterval> &_intervals) {
       
       t_real const mids = (_maxs + _mins) * 0.5;
       // compute number of roots
       t_int const upper = getUpper(_det, mids);
       t_int const lower = _nroots - upper;

       // Checks whether we have an interval, or whether we should bisect some more.
       auto check_and_set = [&](t_int _iroots, t_real _min, t_real _max) {

         if(_iroots == 1) 
           _intervals.push_back({_min, _max, _iroots});
         else if (_iroots > 1) {
           if ((_max - _min) < _tolerance)  
             _intervals.push_back({_min, _max, _iroots});
           else step_(_det, _iroots, _min, _max, _tolerance, _intervals);
         }
       };

       check_and_set(upper, mids, _maxs);
       check_and_set(lower, _mins, mids);
     }
  }

  std::vector<RootInterval> MSWINDOBE
    find_root_intervals(DeterminantEq const &_det, t_real _mins, t_real _maxs, t_real _tolerance) {

    
    t_real mins = _mins > _maxs ? find_lower_bound_for_roots(_det, _maxs): _mins;

    // Increases size of intervals if number of roots above boundary is unexpected.
    // This works both for min and max, since we can change the sign of alpha.
    auto find_boundary = [&](t_real _s, t_int _nbroots, t_real _alpha) {
      t_int i = 0;
      while( getUpper(_det, _s) != _nbroots ) {
        _s += _alpha;
        if( (++i) > 100 ) 
          throw errors::Runtime("Could not find interval with all roots.");
      }
      return _s;
    };

    // Sets up the size of the intervals to look for.
    t_int const nroots = _det.get_nbroots();
    t_real const maxs = find_boundary(_maxs, 0,  0.1 * (_maxs - _mins));
    mins = find_boundary(mins, nroots, -0.1 * (_maxs - mins));

    // Now calls a recurrent function to bisect intervals until all roots are accounted for.   
    std::vector<RootInterval> intervals;
    intervals.reserve(nroots);
    step_(_det, nroots, mins, maxs, _tolerance, intervals);
    return intervals;
  }

  t_real find_lower_bound_for_roots(DeterminantEq const &_det, t_real _start,
                                    t_real _alpha, t_int _itermax) {

    t_real minroot = _start;
    for(t_int i(0); i < _itermax; ++i) {

       t_rmatrix const H(_det.H(minroot));

       // checks that it is valid.
       // If condition below is true, then means at least one element is NaN.
       if((H.array() != H.array()).any()) {
         minroot = 0.9 * minroot;
         continue;
       }

       // Computes eigenvalues of midpoint.
       Eigen::EigenSolver<t_rmatrix> eigsolver(H);
       if(eigsolver.info() != Eigen::Success) 
         throw errors::Mass("Could not solve eigenvalue problem.");

       // Checks we have no complex eigenvalues.
       if((eigsolver.eigenvalues().array().imag().abs() > 1e-8).any())
         throw errors::ComplexEigenvalues("when computing interval for roots.");

       t_real const minimum(eigsolver.eigenvalues().real().minCoeff());
       if(minimum > minroot) return minroot;
       std::cout  << "minimum: " << minimum << " -- previous: " << minroot << std::endl;
       minroot = minimum + _alpha * (minimum - minroot);
    }
    throw errors::Runtime("Reached maximum number of iterations.");
  }
}
