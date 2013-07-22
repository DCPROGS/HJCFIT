#include <DCProgsConfig.h>

#include <iostream>
#include <vector>
#include "errors.h"
#include "time_filter.h"

namespace DCProgs {

  namespace {
    template<class T> t_rvector interval_filter_impl(Eigen::DenseBase<T> const & _intervals, t_real _tau) {
      
      typename Eigen::DenseBase<T>::Index i(0);
      typename Eigen::DenseBase<T>::Index const nbIntervals = _intervals.size();

      // Finds starting interval (first detectable interval).
      for(; i < nbIntervals && _intervals(i) < _tau; i += 2);
      if(i >= nbIntervals) return t_rvector::Zero(0);
      if(i == nbIntervals - 1) return t_rvector::Zero(0);

      // Now construct filtered time series
      std::vector<t_real> result;
      result.reserve(_intervals.size()-i);
      result.push_back(_intervals(i));
      ++i;
      while (i < nbIntervals - 1) {
        t_real const delta_t(_intervals(i));
        // Normal case: interval is detected.
        if(delta_t >= _tau) result.push_back(delta_t); 
        else {
          // interval is not detected. We must add the length of this interval and the next to the
          // current interval.
          ++i;
          result.back() += delta_t + _intervals(i);
        }
        ++i;
      } // for loop
      // Case for last interval is specialized to avoid lots of in-loop if statements.
      if(i == nbIntervals - 1) {
        if(_intervals(i) >= _tau) result.push_back(_intervals(i)); 
        else result.back() += _intervals(i);
      }
  
      // Now we can figure out result:
      t_rvector vector(result.size());
      std::copy(result.begin(), result.end(), &vector(0));
      return vector;
    }

    template<class T> t_rvector time_filter_impl(Eigen::DenseBase<T> const & _series, t_real _tau) {
      
      typename Eigen::DenseBase<T>::Index const n(_series.size());
      auto const durations = _series.tail(n-1) - _series.head(n-1);
      typename Eigen::DenseBase<T>::Index i(0);
      typename Eigen::DenseBase<T>::Index const nbIntervals = durations.size();

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
      for(decltype(vector.size()) j(1); j < vector.size(); ++j)
        vector(j) = vector(j-1) + result[j];
      return vector;
    }
  }
  t_rvector MSWINDOBE time_filter(t_rvector const & _series, t_real _tau) {
     return time_filter_impl(_series, _tau);
  }
}
