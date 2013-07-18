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
    t_real chained_likelihood(T_INTERVAL_ITERATOR _begin, T_INTERVAL_ITERATOR _end, 
                              T_G const & _g, t_initvec const &_initial,
                              t_rvector const &_final ) {
      t_rvector current = _final;
      bool is_AF = true;
      for(; _begin != _end; ++_begin, is_AF = not is_AF)
        current = is_AF ? _g.af(static_cast<t_real>(*_begin)): _g.fa(static_cast<t_real>(*_begin));
      return _initial * current;
    }
  //! Computes likelihood of a time series.
  //! \param[in] _intervals: Time intervals, starting and ending with an "open" interval.
  //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
  //!                member function, where the argument is the length of an open or shut interval.
  //! \param[in] _initial: initial occupancies.
  //! \param[in] _final: final occupancies.
  template<class T_G>
    t_real chained_likelihood(t_rvector const &_intervals, T_G const & _g,
                              t_initvec const &_initial, t_rvector const &_final ) {
      return chained_likelihood( &_intervals(0), &_intervals(_intervals.size()-1) + 1, 
                                 _g, _initial, _final );
    }
}

#endif 

