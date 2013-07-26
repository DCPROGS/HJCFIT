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

#include <DCProgsConfig.h>

#include <iostream>
#include <vector>
#include <tuple>
#include "errors.h"
#include "time_filter.h"

namespace DCProgs {

  namespace {
    template<class T> std::tuple<std::vector<t_real>, t_int>
      interval_filter_impl(Eigen::DenseBase<T> const & _intervals, t_real _tau) {
      
        typename Eigen::DenseBase<T>::Index i(0);
        typename Eigen::DenseBase<T>::Index const nbIntervals = _intervals.size();
  
        // Finds starting interval (first *open* detectable interval).
        for(; i < nbIntervals && _intervals(i) < _tau; i += 2);
        if(i >= nbIntervals) return std::make_tuple(std::vector<t_real>(), nbIntervals); 
        if(i == nbIntervals - 1) return std::make_tuple(std::vector<t_real>(), nbIntervals); 
  
        // Save i index
        t_int const initial = i;
        // Now construct filtered time series
        std::vector<t_real> result;
        result.reserve(_intervals.size()-i);
        result.push_back(_intervals(i));
        while (++i < nbIntervals - 1) {
          t_real const delta_t(_intervals(i));
          // Normal case: interval is detected.
          if(delta_t >= _tau) result.push_back(delta_t); 
          else
            // interval is not detected. We must add the length of this interval and the next to the
            // current interval.
            result.back() += delta_t + _intervals(++i);
        }
        // Case for last interval is specialized to avoid lots of in-loop if statements.
        if(i == nbIntervals - 1) {
          if(_intervals(i) >= _tau) result.push_back(_intervals(i)); 
          else result.back() += _intervals(i);
        }
        return std::make_tuple(std::move(result), initial);
      }
  } // anonymous namesace


  t_rvector MSWINDOBE time_filter(t_rvector const & _series, t_real _tau) {
    t_rvector::Index const n(_series.size());
    std::tuple<std::vector<t_real>, t_int> const intervals(
        interval_filter_impl(_series.tail(n-1) - _series.head(n-1), _tau)
    );
    if(std::get<0>(intervals).size() == 0) return t_rvector::Zero(0);
    t_rvector vector(std::get<0>(intervals).size()+1);
    vector(0) = _series(std::get<1>(intervals));
    std::vector<t_real> :: const_iterator i_first = std::get<0>(intervals).begin();
    for(t_rvector::Index j(1); j < vector.size(); ++j, ++i_first)
      vector(j) = vector(j-1) + *i_first;
    return vector;
  }

  // Filters an incoming list of intervals.
  t_rvector MSWINDOBE interval_filter(t_rvector const & _intervals, t_real _tau) {
    t_rvector::Index const n(_intervals.size());
    std::tuple<std::vector<t_real>, t_int> const intervals
       = interval_filter_impl(_intervals, _tau);

    t_rvector result;
    result = t_rvector::Map(&std::get<0>(intervals)[0], std::get<0>(intervals).size());
    return result;
  }
}
