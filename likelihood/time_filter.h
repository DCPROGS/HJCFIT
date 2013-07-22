#ifndef DCPROGS_TIME_FILTER_H
#define DCPROGS_TIME_FILTER_H

#include <DCProgsConfig.h>

namespace DCProgs {

  //! \brief Filters an incoming time series.
  //! \param[in] _series: Array of times.
  //! \param[in] _tau: Critical resolution below which to remove events.
  //! \return: A time series such that the intervals smaller than _tau filtered out, and the
  //!          adjacent intervals joined.
  t_rvector MSWINDOBE time_filter(t_rvector const & _series, t_real _tau);
  //! \brief Filters an incoming list of intervals.
  //! \param[in] _intervals: Array of time intervals.
  //! \param[in] _tau: Critical resolution below which to remove events.
  //! \return: A  vector with the intervals smaller than _tau filtered out, and the adjacent
  //!          intervals joined.
  t_rvector MSWINDOBE interval_filter(t_rvector const & _intervals, t_real _tau);
}

#endif 
