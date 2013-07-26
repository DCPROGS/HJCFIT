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
