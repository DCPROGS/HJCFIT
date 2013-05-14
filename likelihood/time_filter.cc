#include <DCProgsConfig.h>

#include <iostream>
#include <vector>
#include "errors.h"
#include "time_filter.h"

namespace DCProgs {

  namespace {
    template<class T> t_rvector time_filter_impl(Eigen::DenseBase<T> const & _series, t_real _tau) {
      
      long const n(_series.size());
      auto const durations = _series.tail(n-1) - _series.head(n-1);
      long i(0);
      long const nbIntervals = durations.size();

      // Finds starting interval (first detectable interval).
      for(; i < nbIntervals && durations(i) < _tau; i += 2);
      if(i >= nbIntervals) return t_rvector::Zero(0);
      if(i == nbIntervals - 1) return _series.tail(2);

      // Now construct filtered time series
      std::vector<t_real> result;
      result.reserve(durations.size());
      result.push_back(_series(i));
      result.push_back(durations(i));
      ++i;
      while (i < nbIntervals - 1) {
        t_real const delta_t(durations(i));
        // Normal case: interval is detected.
        if(delta_t >= _tau) result.push_back(delta_t); 
        else {
          // interval is not detected. We must add the length of this interval and the next to the
          // current interval.
          ++i;
          result.back() += delta_t + durations(i);
        }
        ++i;
      } // for loop
      // Case for last interval is specialized to avoid lots of in-loop if statements.
      if(i == nbIntervals - 1) {
        if(durations(i) >= _tau) result.push_back(durations(i)); 
        else result.back() += durations(i);
      }
  
      // Now we can figure out result:
      t_rvector vector(result.size());
      vector(0) = result.front();
      for(size_t j(1); j < vector.size(); ++j)
        vector(j) = vector(j-1) + result[j];
      return vector;
    }
  }
  t_rvector time_filter(t_rvector const & _series, t_real _tau) {
     return time_filter_impl(_series, _tau);
  }

}

  extern "C" int time_filter_wrapper(int _N, DCProgs::t_real const _series[], DCProgs::t_real _tau, DCProgs::t_real _out[]) {
    try {
      Eigen::Map<DCProgs::t_rvector const> const map(_series, _N);
      
      DCProgs::t_rvector const result(DCProgs::time_filter_impl(map, _tau));
      std::copy(&result(0), (&result(0)) + result.size(), _out);
      return int(result.size());
    } catch (...) { return 0; }
  }
