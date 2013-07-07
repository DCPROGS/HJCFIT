#ifndef DCPROGS_LIKELIHOOD_MISSED_EVENT_H
#define DCPROGS_LIKELIHOOD_MISSED_EVENT_H

#include <DCProgsConfig.h>

#include <tuple>
#include <map>
#include <vector>
#include <memory>


#include "exactG.h"
#include "approxG.h"

namespace DCProgs {

  //! \brief Implementation of recursion for exact missed-event G function
  //! \details Implements the exact-missed event probability calculations, as detailed in Hawkes,
  //! Jalali, and Colquhoun (1990). Specifically, this is equation 3.2.
  class MSWINDOBE MissedEventG : protected ExactSurvivor, protected ApproxSurvivor {
    public:
      //! Initializes exact G functor.
      //! \param[in] _af: Determinant equation for af.
      //! \param[in] _fa: Determinant equation for af.
      //! \param[in] _roots_af: Roots of determinant equation for af
      //! \param[in] _roots_fa: Roots of determinant equation for fa
      //! \param[in] _nmax: Switches to asymptotic values after \f$t\geq n_{\mathrm{max}}\tau\f$.
      ApproxG( DeterminantEq const &_af,
               DeterminantEq const &_fa,
               std::vector<Root> const &_roots_af, 
               std::vector<Root> const &_roots_fa,
               t_int _nmax=2 );

      //! Open to close transitions 
      t_rmatrix af(t_real _t) const { return survivor_af(_t) * factor_af_; }
      //! Close to open transitions
      t_rmatrix fa(t_real _t) const { return survivor_af(_t) * factor_af_; }
      //! Probability of no shut times detected between 0 and t.
      t_rmatrix survivor_af(t_real _t) const {
        return _t > tmax_ ? ApproxSurvivor::af(_t): ExactG::af(_t);
      }
      //! Probability of no open times detected between 0 and t.
      t_rmatrix survivor_fa(t_real _t) const {
        return _t > tmax_ ? ApproxSurvivor::af(_t): ExactG::af(_t);
      }

      //! Sets \f$t\geq n_{\mathrm{max}}\tau\f$
      void  set_nmax(t_int _n) {
        if(_n < 0)
          throw errors::DomainError("Switch over between exact and asymptotic"
                                    " cannot be smaller than 0.");
        else if(_n > 10) 
          throw errors::DomainError("Switch over between exact and asymptotic"
                                    " cannot be larger than 10.");
        nmax_ = _n; tmax_ = _n * ExactG::get_tau(); 
      }
      //! When to switch to asymptotic values
      t_int  get_nmax() const { return _nmax; }
      //! Gets the value of missed event resolution;
      t_real get_tau() const { return ExactG::get_tau(); }


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

