#ifndef DCPROGS_LIKELIHOOD_DETERMINANT_EQUATION_H
#define DCPROGS_LIKELIHOOD_DETERMINANT_EQUATION_H

#include <DCProgsConfig.h>
#include <ostream>
#include "laplace_survivor.h"

namespace DCProgs {

  //! A functor to compute the W matrix, so as to find its roots.
  //! \detail The whole implementation is done w.r.t. to AF transitions. 
  //!         However, in practice, this is sufficient to compute FA transitions as well, by messing
  //!         with the input matrix.
  class MSWINDOBE DeterminantEq : public LaplaceSurvivor {
    
    public:
      //! Constructor.
      DeterminantEq(QMatrix const & _qmatrix, t_real _tau) : LaplaceSurvivor(_qmatrix), tau_(_tau) {};
      //! Constructor.
      DeterminantEq(DeterminantEq const & _c) : LaplaceSurvivor(_c), tau_(_c.tau_) {};

      //! Computes the determinant \f$\mathrm{det}(sI - H(s))\f$
      //! \param[in] _s: Value of the laplacian scale.
      t_real operator()(t_real _s) const { return DeterminantEq::operator()(_s, tau_); }
      //! Computes the determinant \f$\mathrm{det}(sI - H(s, \tau))\f$
      //! \param[in] _s: Value of the laplacian scale.
      t_real operator()(t_real _s, t_real _tau) const {
        return (_s * this->id_() - LaplaceSurvivor::H(_s, _tau)).determinant();
      }
      DeterminantEq transpose() const { return DeterminantEq(get_qmatrix().transpose(), tau_); }

      //! Computes \f$sI - Q_{AA} - Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}\f$
      //! \param[in] _s: Value of the laplacian scale.
      t_rmatrix H(t_real _s) const { return LaplaceSurvivor::H(_s, tau_); }
      using LaplaceSurvivor::H;

      //! Derivative along _s
      t_rmatrix s_derivative(t_real _s) const { 
        return LaplaceSurvivor::s_derivative(_s, tau_); 
      }
      using LaplaceSurvivor::s_derivative;
      //! Max length of missed events.
      t_real get_tau() const { return tau_; }
      //! Max length of missed events.
      void set_tau(t_real _tau) { tau_ = _tau; }
    
    protected:
      //! Max length of missed events
      t_real tau_;
  };

  //! Dumps Determinantal equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, DeterminantEq const &_self);
}
#endif 
