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

#ifndef DCPROGS_LIKELIHOOD_H
#define DCPROGS_LIKELIHOOD_H

#include <DCProgsConfig.h>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
#if defined(_OPENMP)
#include <omp.h>
#endif

#include "qmatrix.h"
#include "errors.h"

namespace DCProgs {

  //! Computes likelihood of a time series.
  //! \param[in] _begin First interval in the time series. This must be an "open" interval.
  //! \param[in] _end One past last interval.
  //! \param[in] _g The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
  //!                member function, where the argument is the length of an open or shut interval.
  //! \param[in] _initial initial occupancies.
  //! \param[in] _final final occupancies.
  template<class T_INTERVAL_ITERATOR, class T_G>
    t_real chained_likelihood( T_G const & _g, T_INTERVAL_ITERATOR _begin, T_INTERVAL_ITERATOR _end, 
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

  //! \brief Computes log10-likelihood of a time series.
  //! \details Adds a bit of trickery to take care of exponent. May make this a bit more stable.
  //! \param[in] _begin First interval in the time series. This must be an "open" interval.
  //! \param[in] _end One past last interval.
  //! \param[in] _g The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
  //!                member function, where the argument is the length of an open or shut interval.
  //! \param[in] _initial initial occupancies.
  //! \param[in] _final final occupancies.
  template<class T_G>
    t_real chained_log10_likelihood( T_G const & _g, const t_Burst burst,
                                     t_initvec const &_initial, t_rvector const &_final,
                                     t_int const threads, bool const openmphighlevel) {
      auto _begin = burst.begin();
      auto _end = burst.end();
      t_int const intervals = _end - _begin;
      if( (intervals) % 2 != 1 )
        throw errors::Domain("Expected a burst with odd number of intervals");
      t_initvec current = _initial * _g.af(static_cast<t_real>(*_begin));
      t_int exponent(0);
      const t_int cols = current.cols();
      const auto identity = t_rmatrix::Identity(cols, cols);
      std::vector<t_rmatrix> current_vec(threads, identity);
      std::vector<t_int> exponents(threads, 0);
      bool openmplowlevel = (intervals>100) && (!openmphighlevel);
      #pragma omp parallel default(none), shared(_g, current_vec, exponents), if(openmplowlevel)
      {
        t_int thread;
        #if defined(_OPENMP)
          thread = omp_get_thread_num();
        #else
          thread = 0;
        #endif
        exponents[thread] = 0;
        current_vec[thread] = identity;
        #pragma omp for schedule(static)
        for(t_int j=1; j<intervals-1; j=j+2) {
          current_vec[thread] = current_vec[thread] * _g.fa(static_cast<t_real>(burst[j]));
          current_vec[thread] = current_vec[thread] * _g.af(static_cast<t_real>(burst[j+1]));
          t_real const max_coeff = current_vec[thread].array().abs().maxCoeff();
          if(max_coeff > 1e20) {
            current_vec[thread]  *= 1e-20;
            exponents[thread] += 20;
          } else if(max_coeff < 1e-20) {
            current_vec[thread]  *= 1e+20;
            exponents[thread] -= 20;
          }
        }
      }
      for (t_int tmpexp : exponents)
        exponent += tmpexp;
      for (auto tmpcurrent : current_vec)
        current = current * tmpcurrent;
      return std::log10(current * _final) + exponent;
    }


  //! \brief Likelihood of a set of bursts
  //! \details This functor takes as input a Q-matrix and returns the likelihood for a given set of
  //!          bursts. It is, in practice, a convenience object with which to perform likelihood
  //!          optimization.
  //!
  //!          At each call, it creates a likelihood object MissedEventsG using the input parameters
  //!          given during initialization, and the Q-matrix given on input. It then returns the
  //!          log10 likelihood for the set of bursts given on input.
  class MSWINDOBE Log10Likelihood {
    public:
      //! Set of bursts for which to compute the likelihood.
      t_Bursts bursts;
      //! Number of open states.
      t_uint nopen;
      //! Max length of missed events
      t_real tau;
      //! \brief \f$t_{\mathrm{crit}}\f$. 
      //! \details If negative or null, will use equilibrium occupancies rather than CHS
      //!          occupancies.
      t_real tcritical;
      //! Number of intervals for which to compute exact result.
      t_uint nmax;
      //! Tolerance for root finding.
      t_real xtol;
      //! Tolerance for root finding.
      t_real rtol;
      //! Maximum number of iterations for root finding.
      t_uint itermax;
      //! Lower bound bracketing all roots.
      t_real lower_bound;
      //! Upper bound bracketing all roots.
      t_real upper_bound;

      t_int omp_num_threads;
      //! Constructor
      //! \param[in] _bursts A vector of bursts. Each burst is a vector of intervals, starting with
      //!            an open interval. The intervals should be prefiltered for the maximum
      //!            missed-event length.
      //! \param[in] _nopen Number of open states. The open states must be in the top left corner of
      //!            the Q-matrix.
      //! \param[in] _tau Maximum length of the missed events
      //! \param[in] _tcritical Parameter for CHS vectors (see \cite colquhoun:1996) if positive. If
      //!            negative, then equilibrium occupancies will be used as initial and final
      //!            states (as in \cite colquhoun:1982)
      //! \param[in] _nmax The exact missed-event likelihood will be computed for 
      //!            \f$ t < n_{\mathrm{max}} \tau\f$
      //! \param[in] _xtol Tolerance criteria for brentq().
      //! \param[in] _rtol Tolerance criteria for brentq().
      //! \param[in] _itermax Maximum number of iteration when calling brentq().
      //! \param[in] _lowerbound Lower bound of the interval bracketing all roots. If None, the
      //!            lower bound is obtained from find_lower_bound_for_roots().
      //! \param[in] _upperbound Upper bound of the interval bracketing all roots. If None, the
      //!            upper bound is obtained from find_upper_bound_for_roots().
      Log10Likelihood   ( t_Bursts const &_bursts, t_uint _nopen, t_real _tau,
                          t_real _tcritical=-1e0, t_uint _nmax=2, t_real _xtol=1e-10,
                          t_real _rtol=1e-10, t_uint _itermax=100,
                          t_real _lowerbound=quiet_nan,
                          t_real _upperbound=quiet_nan ) 
                      : bursts(_bursts), nopen(_nopen), tau(_tau), tcritical(_tcritical),
                        nmax(_nmax), xtol(_xtol), rtol(_rtol), itermax(_itermax),
                        lower_bound(_lowerbound), upper_bound(_upperbound) {
                          #pragma omp parallel default(none)
                          {
                            #pragma omp single
                            {
                              #if defined(_OPENMP)
                              omp_num_threads = omp_get_num_threads();
                              #else
                              omp_num_threads = 1;
                              #endif
                            }
                          }
                        }
     
      //! \brief Computes likelihood for each burst
      //! \return a DCProgs::t_rvector 
      t_rvector vector(t_rmatrix const &_Q) const { return vector(QMatrix(_Q, nopen)); }
      //! \brief Computes likelihood for each burst
      //! \return a DCProgs::t_rvector 
      t_rvector vector(QMatrix const &_Q) const;
      //! Log-likelihood 
      t_real operator()(t_rmatrix const &_Q) const { return operator()(QMatrix(_Q, nopen)); }
      //! Log-likelihood 
      t_real operator()(QMatrix const &_Q) const;
  };
  //! Dumps likelihood to stream.
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, Log10Likelihood const & _self);
  //! Dumps burst to stream.
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, t_Burst const & _self);
  //! Dumps bursts to stream.
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, t_Bursts const & _self);
}

#endif 
