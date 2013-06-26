#ifndef DCPROGS_TIME_FILTER_H
#define DCPROGS_TIME_FILTER_H

#include <DCProgsConfig.h>

namespace DCProgs {

  //! \brief Filters an incoming time serie.
  //! \param[in] _series: Array of times.
  //! \param[in] _tau: Critical resolution below which to remove events.
  //! \return: A  vector with the intervals smaller than _tau filtered out, and the adjacent
  //!          intervals joined.
  t_rvector MSWINDOBE time_filter(t_rvector const & _series, t_real _tau);
}

//! \brief Wrapper around the time filter function.
//! \param[in] _N: Number of times in the time series.
//! \param[in] _series: Array of times.
//! \param[in] _tau: Critical resolution below which to remove events.
//! \param[inout] _out: A pre-allocated array of times containing the filtered series. It must be at
//!                     least of length _N.
//! \return: The number of significant values in _out.
extern "C" int time_filter_wrapper(int _N, 
                                   DCProgs::t_real const _series[], 
                                   DCProgs::t_real _tau,
                                   DCProgs::t_real _out[]);
#endif 
