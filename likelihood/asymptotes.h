#ifndef DCPROGS_LIKELIHOOD_ASYMPTOTES_H
#define DCPROGS_LIKELIHOOD_ASYMPTOTES_H

#include <DCProgsConfig.h>
#include "state_matrix.h"

namespace DCProgs {


  //! A functor to compute asymptotic missed event G.
  //! \detail The whole implementation is done w.r.t. to AF transitions. 
  //!         However, in practice, this is sufficient to compute FA transitions as well, by messing
  //!         with the input matrix.
  class DeterminantEq {
    
    public:
      //! Constructor. 
      //! \param[in] _matrix: The transition state matrix for which to compute
      //!                     \f$^eG_{AF}(t\righarrow\infty)\f$
      //! \param[in] _tau: Missed event resolution.
      //! \param[in] _doopen: Whether to do AF or FA.
      DeterminantEq(StateMatrix & _matrix, t_real _tau, bool _doopen=true);

      //! Computes \f$Q_{AA} + Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}\f$
      inline t_rmatrix H(t_real _s) const {
        return matrix_.aa() + matrix_.af() * this->integral_(_s) * matrix_.fa();
      }

      //! Computes the determinant \f$\mathrm{det}(sI - H(s))\f$
      inline t_real operator()(t_real _s) const { 
        return (_s * this->id_() - H(_s)).determinant();
      }
      //! Derivative along _s
      inline t_rmatrix s_derivative(t_real _s) const;

    protected:
      //! Computes integral \f$\int_0^\tau\partial\,t\ e^{(Q_{FF} - sI)t}\f$
      t_rmatrix integral_(t_real _s) const;
      //! Just the identity, just to write shorter code.
      inline auto id_() const ->decltype(t_rmatrix::Identity(1, 1)) 
        { return t_rmatrix::Identity(matrix_.nopen, matrix_.nopen); }

    protected:
      //! Time below which events are missed
      t_real tau_;
      //! The transition state matrix on which to act.
      StateMatrix matrix_;
      //! The eigenvalues of the ff matrix. Computed once.
      t_rvector ff_eigenvalues_;
      //! The eigenvectors of the ff matrix. Computed once.
      t_rmatrix ff_eigenvectors_;
      //! Hard coded static constant zero.
      constexpr static t_real ZERO = 1e-12;
  };

}
#endif 
