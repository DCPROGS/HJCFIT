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
  class MSWINDOBE MissedEventG : protected ExactG {
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
      t_rmatrix af(t_real _t) const {
        return _t > tmax_ ? asymptotic_af_->operator()(_t) * factor_af_: ExactG::af(_t);
      }
      //! Close to open transitions
      t_rmatrix fa(t_real _t) const {
        return _t > tmax_ ? asymptotic_fa_->operator()(_t) * factor_fa_: ExactG::af(_t);
      }
      //! Probability of no shut times detected between 0 and t.
      t_rmatrix R_af(t_real _t) const {
        return _t > tmax_ ? asymptotic_af_->operator()(_t): ExactG::af(_t);
      }
      //! Probability of no open times detected between 0 and t.
      t_rmatrix R_fa(t_real _t) const {
        return _t > tmax_ ? asymptotic_fa_->operator()(_t): ExactG::af(_t);
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
      //! \brief Implementation of recursion for exact missed-event G function
      //! \details This is an interface to the function recursion_formula.  In practice, this object
      //!          needs not be called directly. Rather the public interface (which is about
      //!          computing the likelihood for an event of duration t) is in the containing class
      //!          ExactG.
      class MSWINDOBE RecursionInterface;

#     ifndef HAS_CXX11_UNIQUE_PTR
        //! Type of the pointers holding recursion interfaces.
        typedef std::auto_ptr<Asymptotes> t_AsymptotesPtr;
#     else
        //! Type of the pointers holding recursion interfaces.
        typedef std::unique_ptr<Asymptotes> t_AsymptotesPtr;
#     endif
      //! Switches to asymptotic values for \f$t\geq n_{\mathrm{max}}\tau\f$.
      t_int nmax_;
      //! Max length of missed events.
      t_int tmax_;
      //! Pointer to AF recursion interface
      t_AsymptotesPtr asymptotes_af_;
      //! Pointer to FA recursion interface
      t_AsymptotesPtr asymptotes_fa_;
  };
}

#endif 

