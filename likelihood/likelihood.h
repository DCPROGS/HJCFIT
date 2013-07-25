#ifndef DCPROGS_LIKELIHOOD_H
#define DCPROGS_LIKELIHOOD_H

#include <DCProgsConfig.h>

#include <unsupported/Eigen/MatrixFunctions>

namespace DCProgs {

  //! Computes likelihood of a time series.
  //! \param[in] _begin: First interval in the time series. This must be an "open" interval.
  //! \param[in] _end: One past last interval.
  //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
  //!                member function, where the argument is the length of an open or shut interval.
  //! \param[in] _initial: initial occupancies.
  //! \param[in] _final: final occupancies.
  template<class T_INTERVAL_ITERATOR, class T_G>
    t_real chained_likelihood(T_G const & _g, T_INTERVAL_ITERATOR _begin, T_INTERVAL_ITERATOR _end, 
                              t_initvec const &_initial, t_rvector const &_final ) {
      if( (_end - _begin) % 2 != 1 )
        throw errors::Domain("Expected a burst with odd number of intervals");
      t_initvec current = _initial * _g.af(static_cast<t_real>(*_begin));
      while(++_begin != _end) {
        current = current * _g.fa(static_cast<t_real>(*_begin));
        current = current * _g.af(static_cast<t_real>(*(++_begin)));
      }
      return current * _final;
    }
  //! Computes likelihood of a time series.
  //! \param[in] _intervals: Time intervals, starting and ending with an "open" interval.
  //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
  //!                member function, where the argument is the length of an open or shut interval.
  //! \param[in] _initial: initial occupancies.
  //! \param[in] _final: final occupancies.
  template<class T_G>
    t_real chained_likelihood(T_G const & _g, t_rvector const &_intervals, 
                              t_initvec const &_initial, t_rvector const &_final ) {
      return chained_likelihood( _g, &_intervals(0), &_intervals(_intervals.size()-1) + 1, 
                                 _initial, _final );
    }
}

#endif 

