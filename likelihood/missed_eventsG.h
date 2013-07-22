#ifndef DCPROGS_LIKELIHOOD_MISSED_EVENT_H
#define DCPROGS_LIKELIHOOD_MISSED_EVENT_H

#include <DCProgsConfig.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "exact_survivor.h"
#include "approx_survivor.h"

namespace DCProgs {

  class MissedEventsG;

  //! CHS matrices \f$H_{FA}\f$
  t_rmatrix CHS_matrix_Hfa(MissedEventsG const &, t_real);
  //! CHS matrices \f$H_{AF}\f$
  t_rmatrix CHS_matrix_Haf(MissedEventsG const &, t_real);


  //! \brief Implementation of recursion for exact missed-event G function
  //! \details Implements the exact-missed event probability calculations, as detailed in Hawkes,
  //! Jalali, and Colquhoun (1990). Specifically, this is equation 3.2.
  class MSWINDOBE MissedEventsG : protected ExactSurvivor, protected ApproxSurvivor {
     
      friend t_rmatrix CHS_matrix_Hfa(MissedEventsG const &, t_real);
      friend t_rmatrix CHS_matrix_Haf(MissedEventsG const &, t_real);
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
                      laplace_a_(new LaplaceSurvivor(_af.get_qmatrix())),
                      laplace_f_(new LaplaceSurvivor(_fa.get_qmatrix())),
                      nmax_(_nmax), tmax_(nmax_ * get_tau()),
                      af_factor_( _af.get_qmatrix().af()
                                  * (_af.get_tau() * _af.get_qmatrix().ff()).exp() ),
                      // _fa is already transpose of _af, so it is indeed _fa.matrix.af * e^...
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
                      laplace_a_(new LaplaceSurvivor(_qmatrix)),
                      laplace_f_(new LaplaceSurvivor(_qmatrix.transpose())),
                      nmax_(_nmax), tmax_(_tau*t_real(_nmax)),
                      af_factor_(_qmatrix.af() * (_tau * _qmatrix.ff()).exp()),
                      fa_factor_(_qmatrix.fa() * (_tau * _qmatrix.aa()).exp()) {}

      //! Open to close transitions 
      t_rmatrix af(t_real _t) const {
        return survivor_af(_t - ExactSurvivor::get_tau()) * af_factor_; 
      }
      //! Close to open transitions
      t_rmatrix fa(t_real _t) const {
        return survivor_fa(_t - ExactSurvivor::get_tau()) * fa_factor_; 
      }
      //! Probability of no shut times detected between 0 and t.
      t_rmatrix survivor_af(t_real _t) const {
        return _t > tmax_ ? ApproxSurvivor::af(_t): ExactSurvivor::af(_t);
      }
      //! Probability of no open times detected between 0 and t.
      t_rmatrix survivor_fa(t_real _t) const {
        return _t > tmax_ ? ApproxSurvivor::fa(_t): ExactSurvivor::fa(_t);
      }

      //! Sets \f$t\geq n_{\mathrm{max}}\tau\f$
      void  set_nmax(t_int _n) { 
        if(_n <= 0) throw errors::Domain("n should be strictly positive.");
        nmax_ = _n; tmax_ = t_real(_n) * ExactSurvivor::get_tau(); 
      }
      //! When to switch to asymptotic values
      t_int  get_nmax() const { return nmax_; }
      //! Gets the value of missed event resolution;
      t_real get_tau() const { return ExactSurvivor::get_tau(); }
      //! Tmax is the maximum time after which to switch to approximate calculations.
      t_real get_tmax() const { return tmax_; }

      //! \f$Q_{AF}e^{-Q_{FF}\tau} \f$
      t_rmatrix const & get_af_factor() const { return af_factor_; }
      //! \f$Q_{FA}e^{-Q_{AA}\tau} \f$
      t_rmatrix const & get_fa_factor() const { return fa_factor_; }

      //! Exact laplace of AF
      t_rmatrix laplace_af(t_real _s) const {
        return laplace_a_->operator()(_s, get_tau()) * std::exp(-_s*get_tau()) * af_factor_;
      }
      //! Exact laplace of FA
      t_rmatrix laplace_fa(t_real _s) const {
        return laplace_f_->operator()(_s, get_tau()) * std::exp(-_s*get_tau()) * fa_factor_;
      }
      //! Returns current QMatrix
      QMatrix const & get_qmatrix() const { return laplace_a_->get_qmatrix(); }

    protected:
#     ifndef HAS_CXX11_UNIQUE_PTR
        //! Type of the pointers holding laplace object.
        typedef std::auto_ptr<LaplaceSurvivor> t_LaplacePtr;
#     else
        //! Type of the pointers holding laplace object.
        typedef std::unique_ptr<LaplaceSurvivor> t_LaplacePtr;
#     endif
      //! Laplace Survivor function \f$^{A}R(s)\f$.
      t_LaplacePtr laplace_a_;
      //! Laplace Survivor function \f$^{F}R(s)\f$.
      t_LaplacePtr laplace_f_;
      //! Switches to asymptotic values for \f$t\geq n_{\mathrm{max}}\tau\f$.
      t_int nmax_;
      //! Max length of missed events.
      t_real tmax_;
      //! \f$Q_{AF}e^{-Q_{FF}\tau} \f$
      t_rmatrix af_factor_;
      //! \f$Q_{FA}e^{-Q_{AA}\tau} \f$
      t_rmatrix fa_factor_;
  };

  //! Dumps Missed-Events likelihood to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, MissedEventsG const &_self);
}

#endif 
