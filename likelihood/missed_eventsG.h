#ifndef DCPROGS_LIKELIHOOD_MISSED_EVENT_H
#define DCPROGS_LIKELIHOOD_MISSED_EVENT_H

#include <DCProgsConfig.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "exact_survivor.h"
#include "approx_survivor.h"

namespace DCProgs {

  //! \brief Implementation of recursion for exact missed-event G function
  //! \details Implements the exact-missed event probability calculations, as detailed in Hawkes,
  //! Jalali, and Colquhoun (1990). Specifically, this is equation 3.2.
  class MSWINDOBE MissedEventsG : protected ExactSurvivor, protected ApproxSurvivor {
    public:
      //! Initializes missed events G functor.
      //! \param[in] _af: Determinant equation for af.
      //! \param[in] _fa: Determinant equation for af.
      //! \param[in] _roots_af: Roots of determinant equation for af
      //! \param[in] _roots_fa: Roots of determinant equation for fa
      //! \param[in] _nmax: Switches to asymptotic values after \f$t\geq n_{\mathrm{max}}\tau\f$.
      MissedEventsG  ( DeterminantEq const &_af,
                       std::vector<Root> const &_roots_af, 
                       DeterminantEq const &_fa,
                       std::vector<Root> const &_roots_fa,
                       t_int _nmax=2 )
                    : ExactSurvivor(_af.get_qmatrix(), _af.get_tau()),
                      ApproxSurvivor(_af, _roots_af, _fa, _roots_fa),
                      nmax_(_nmax), tmax_(_af.get_tau()*_nmax),
                      af_factor_( _af.get_qmatrix().af()
                                  * (_af.get_tau() * _af.get_qmatrix().ff()).exp() ),
                      fa_factor_( _fa.get_qmatrix().af()
                                  * (_af.get_tau() * _fa.get_qmatrix().ff()).exp() ) {}
      //! Initializes missed events functor.
      //! \param[in] _qmatrix: Transition matrix
      //! \param[in] _tau: resolution/max length missed events
      //! \param[in] _findroots: A functor with which to find all roots.
      //!                        This function should take a DeterminantEq as its sole argument and
      //!                        return a std::vector<RootIntervals>
      MissedEventsG   ( QMatrix const &_qmatrix, t_real _tau, 
                        t_RootFinder const &_findroots, t_int _nmax=2 )
                    : ExactSurvivor(_qmatrix, _tau),
                      ApproxSurvivor(_qmatrix, _tau, _findroots), 
                      nmax_(_nmax), tmax_(_tau*_nmax),
                      af_factor_(_qmatrix.af() * (_tau * _qmatrix.ff()).exp()),
                      fa_factor_(_qmatrix.fa() * (_tau * _qmatrix.aa()).exp()) {}

      //! Open to close transitions 
      t_rmatrix af(t_real _t) const {
        return survivor_af(_t - ExactSurvivor::get_tau()) * af_factor_; 
      }
      //! Close to open transitions
      t_rmatrix fa(t_real _t) const {
        return survivor_af(_t - ExactSurvivor::get_tau()) * fa_factor_; 
      }
      //! Probability of no shut times detected between 0 and t.
      t_rmatrix survivor_af(t_real _t) const {
        return _t > tmax_ ? ApproxSurvivor::af(_t): ExactSurvivor::af(_t);
      }
      //! Probability of no open times detected between 0 and t.
      t_rmatrix survivor_fa(t_real _t) const {
        return _t > tmax_ ? ApproxSurvivor::af(_t): ExactSurvivor::af(_t);
      }

      //! Sets \f$t\geq n_{\mathrm{max}}\tau\f$
      void  set_nmax(t_int _n) { nmax_ = _n; tmax_ = _n * ExactSurvivor::get_tau(); }
      //! When to switch to asymptotic values
      t_int  get_nmax() const { return nmax_; }
      //! Gets the value of missed event resolution;
      t_real get_tau() const { return ExactSurvivor::get_tau(); }

      //! \f$Q_{AF}e^{-Q_{FF}\tau} \f$
      t_rmatrix const & get_af_factor() const { return af_factor_; }
      //! \f$Q_{FA}e^{-Q_{AA}\tau} \f$
      t_rmatrix const & get_fa_factor() const { return fa_factor_; }


    protected:
      //! Switches to asymptotic values for \f$t\geq n_{\mathrm{max}}\tau\f$.
      t_int nmax_;
      //! Max length of missed events.
      t_int tmax_;
      //! \f$Q_{AF}e^{-Q_{FF}\tau} \f$
      t_rmatrix af_factor_;
      //! \f$Q_{FA}e^{-Q_{AA}\tau} \f$
      t_rmatrix fa_factor_;
  };
}

#endif 

