#include <DCProgsConfig.h>

#include <iostream>

#include "root_finder.h"

namespace DCProgs {


  // Some functions we need here and nowhere else.
  namespace {

     // Compute number of roots in interval.
     t_rvector getEigenvalues(DeterminantEq const &_det, t_real _s) {
       // Compute matrix
       t_rmatrix const H(_det.H(_s));

       // checks that it is valid.
       // If condition below is true, then means at least one element is NaN.
       if(not (H.array() == H.array()).all()) throw errors::NaN("when computing matrix H.");

       // Computes eigenvalues of midpoint.
       Eigen::EigenSolver<t_rmatrix> eigsolver(H);
       if(eigsolver.info() != Eigen::Success) 
         throw errors::Mass("Could not solve eigenvalue problem.");

       // Checks we have no complex eigenvalues.
       if((eigsolver.eigenvalues().array().imag().abs() > 1e-8).any())
         throw errors::ComplexEigenvalues("when computing interval for roots.");

       t_rvector const eigs = eigsolver.eigenvalues().real();
       return eigsolver.eigenvalues().real();
     }
     
     // Compute number of roots in interval.
     t_int getInterval(DeterminantEq const &_det, t_real _min, t_real _max) {
       // Compute matrix
       t_rvector const eigs(getEigenvalues(_det, _min));
       // compute number of roots in interval
       return ( (eigs.array() >= _min) and (eigs.array() < _max)).count();
     }
     // Compute number of roots above input
     t_int getUpper(DeterminantEq const &_det, t_real _s) {
       // Compute matrix
       t_rvector const eigs(getEigenvalues(_det, _s));
       // compute number of roots in interval
       return (eigs.array() >= _s).count();
     }

     // Actual bisecting algorithm.
     // Checks number of roots above and below midpoint. 
     // Based on that information, figures out whether to bisect some more or whether an interval
     // was found.
     void step_(DeterminantEq const &_det, 
                t_real _mins, t_real _maxs, t_real const _tolerance,
                t_int _higher_than_min, t_int _higher_than_max,
                std::vector<RootInterval> &_intervals) {
       
       t_real const mids = (_maxs + _mins) * 0.5;
       t_rvector const eigenvalues = getEigenvalues(_det, mids);

       // compute number of roots
       t_int const higher_than_mid = (eigenvalues.array() >  mids).count();

       // This functor checks whether to bisect some more or whether an 
       auto check_and_set = [&](t_real _min, t_real _max, t_int _imin, t_int _imax) {

         t_int const nroots = _imin - _imax;
         if(nroots == 1) _intervals.push_back({_min, _max, 1});
         else if(nroots != 0) {
           if(_max - _min < _tolerance) {
             t_int const s = (_min + _max) * 0.5;
             t_rvector const eigenvalues = getEigenvalues(_det, s);
             t_int const multiplicity = ((eigenvalues.array() - s) < _tolerance).count();
             if(multiplicity == 0) throw errors::Runtime("Found interval with zero roots.");
             _intervals.push_back({_min, _max, multiplicity});
           } else step_(_det, _min, _max, _tolerance, _imin, _imax, _intervals);
         } 
       };

       // Checks whether we have an interval, or whether we should bisect some more.
       check_and_set(_mins, mids, _higher_than_min, higher_than_mid);
       check_and_set(mids, _maxs, higher_than_mid, _higher_than_max);
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
    step_(_det, mins, maxs, _tolerance, _det.get_nbroots(), 0, intervals);
    return intervals;
  }

  t_real find_lower_bound_for_roots(DeterminantEq const &_det, t_real _start,
                                    t_real _alpha, t_int _itermax) {

    t_real minroot = _start;
    for(t_int i(0); i < _itermax; ++i) {

       t_rmatrix const H(_det.H(minroot));

       // checks that it is valid.
       // If condition below is true, then means at least one element is NaN.
       if(not (H.array() == H.array()).all()) {
         if(std::abs(minroot) < 1e-1) {
           std::ostringstream sstr;
           sstr << "when computing matrix H(" << minroot << "):\n"
                << numpy_io(H) << std::endl;
           throw errors::NaN(sstr.str());
         }
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
       minroot = minimum + _alpha * (minimum - minroot);
    }
    throw errors::Runtime("Reached maximum number of iterations.");
  }
}
