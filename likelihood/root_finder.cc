#include <DCProgsConfig.h>

#include <iostream>

#include "root_finder.h"

namespace DCProgs {


  // Some functions we need here and nowhere else.
  namespace {

     // Compute number of roots in interval.
     t_cvector getEigenvalues(DeterminantEq const &_det, t_real _s) {
       // Compute matrix
       t_rmatrix const H(_det.H(_s));

       // checks that it is valid.
       // If condition below is true, then means at least one element is NaN.
       if(not (H.array() == H.array()).all()) throw errors::NaN("when computing matrix H.");

       // Computes eigenvalues of midpoint.
       Eigen::EigenSolver<t_rmatrix> eigsolver(H);
       if(eigsolver.info() != Eigen::Success) {
         std::ostringstream sstr("Could not solve eigenvalue problem.\n");
         sstr << _det << "\n"; 
         throw errors::Mass(sstr.str());
       }

       return eigsolver.eigenvalues();
     }
     
     //! Count the number of roots at a given point.
     t_int getMultiplicity(DeterminantEq const &_det, t_real _s, t_real _tolerance) {
       t_cvector const eigenvalues = getEigenvalues(_det, _s);
       // NOTE: Newer versions of eigen have and "and" operator, so the
       // following could be simplified in the future.
       t_int result(0);
       t_cvector::Scalar const * i_data = &eigenvalues(0);
       t_cvector::Scalar const * const i_data_end = i_data + eigenvalues.size();
       for(; i_data != i_data_end; ++i_data) 
         if( std::abs(i_data->imag()) < _tolerance
             and std::abs(i_data->real() - _s) < _tolerance ) ++result;
       return result;
     }
     
     // Compute number of roots above input
     t_int getUpper(DeterminantEq const &_det, t_real _s) {
       // Compute matrix
       t_rvector const eigs(getEigenvalues(_det, _s).real());
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
       t_cvector const eigenvalues = getEigenvalues(_det, mids);

       // compute number of roots
       t_int const higher_than_mid = (eigenvalues.array().real() >  mids).count();

       // This functor checks whether to bisect some more or whether an 
       auto check_and_set = [&](t_real _min, t_real _max, t_int _imin, t_int _imax) {

         t_int const nroots = _imin - _imax;
         if(nroots == 1) {
           if(_det(_min) * _det(_max) <= 0e0) _intervals.emplace_back(_min, _max, 1);
         } else if(nroots != 0) {
           if(_max - _min < _tolerance) {
             t_real const s = (_min + _max) * 0.5;
             // count number of eigenvalues that are equal to s *and* real.
             t_int const multiplicity = getMultiplicity(_det, s, _tolerance);
             // Multiplicity == 0 corresponds to complex poles? 
             if(multiplicity != 0) _intervals.emplace_back(_min, _max, multiplicity);
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

    // Now calls a recurrent function to bisect intervals until all roots are accounted for.   
    std::vector<RootInterval> intervals;
    step_(_det, mins, _maxs, _tolerance, _det.get_nbroots(), 0, intervals);
    return intervals;
  }

  t_real MSWINDOBE find_lower_bound_for_roots(DeterminantEq const &_det, t_real _start,
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

     // // Checks we have no complex eigenvalues.
     // if((eigsolver.eigenvalues().array().imag().abs() > 1e-8).any())
     //   throw errors::ComplexEigenvalues("when computing interval for roots.");

       t_real const minimum(eigsolver.eigenvalues().real().minCoeff());
       if(minimum > minroot) return minroot;
       minroot = minimum + _alpha * (minimum - minroot);
    }
    throw errors::Runtime("Reached maximum number of iterations.");
  }

  std::vector<RootInterval> find_root_intervals_brute_force(DeterminantEq const &_det, 
                                                            t_real _resolution,
                                                            t_real _mins,
                                                            t_real _maxs,
                                                            t_real _root_tolerance) {
    if(_mins > _maxs) _mins = find_lower_bound_for_roots(_det, _maxs);
    if(_resolution < 1e-12) throw errors::Domain("Resolution cannot be negative or null");
    if(_mins + 2e0 * _resolution > _maxs) _resolution = (_maxs - _mins) / 10e0;
 
    std::vector<RootInterval> intervals;
    t_real const half_step = 0.5 * _resolution;
    t_real previous = _det(_mins);
    t_real s(_mins + _resolution);
    do {
      t_real current = _det(s);
 
      // Checks we have sensible values. 
      if(not (DCPROGS_ISNAN(current) or DCPROGS_ISNAN(previous)))
      {
        // Sign changed. There should be at least one root.
        if(current * previous < 0e0) {
          // Tries and figures out the multiplicity.
          // It should be at least one and odd. Hence, falls back to one if result is even.
          t_cvector const eigenvalues = getEigenvalues(_det, s - half_step); 
          t_int const multiplicity = getMultiplicity(_det, s - half_step, _resolution);
          intervals.emplace_back(s - _resolution, s, multiplicity % 2 == 1 ? multiplicity: 1);
 
        } else if( std::abs(current) < _root_tolerance) {
          // Tries and figures whether we are skimming the x axis, or whether we will cross it
          // later.
          t_real s_next = s;
          t_real value_next(current);
          bool skimmed_out = false, crossed = false;
          t_real minimum_value(current), minimum_s(s); 
          do {
            s_next += _resolution;
            value_next = _det(s_next);
            crossed = value_next * current < 0;
            skimmed_out = (not crossed) and std::abs(value_next) > _root_tolerance;
            if(std::abs(value_next) < std::abs(minimum_value)) {
              minimum_value = value_next;
              minimum_s = s_next;
            }
          } while( (not skimmed_out) and (not crossed) );
          
          // If we crossed we want to add it as a potential root.
          if(crossed) {
            t_int const multiplicity = getMultiplicity(_det, s_next - half_step, _resolution);
            intervals.emplace_back( s_next - half_step, s_next + half_step, 
                                    multiplicity % 2 == 0 ? multiplicity: 1);
          } else if(skimmed_out) {
            t_int const multiplicity = getMultiplicity(_det, minimum_s, 2*_resolution);
            intervals.emplace_back( minimum_s - _resolution, minimum_s + _resolution, 
                                    multiplicity % 2 == 1 ? multiplicity: 1 );
          }
          // Adjust current value and s position.
          current = value_next;
          s = s_next;
        }
      }
      // Move to next point.
      previous = current;
      s += _resolution;
    } while (s < _maxs);
 
    // return result.
    return intervals;
  }
}
