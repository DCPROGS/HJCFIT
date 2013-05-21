#ifndef DCPROGS_LIKELIHOOD_ASYMPTOTES_H
#define DCPROGS_LIKELIHOOD_ASYMPTOTES_H

#include <DCProgsConfig.h>

namespace DCProgs {


  //! A functor to compute asymptotic missed event G.
  //! \detail The whole implementation is done w.r.t. to AF transitions. 
  //!         However, in practice, this is sufficient to compute FA transitions as well, by messing
  //!         with the input matrix.
  class Asymptotes {
    
    public:
      //! Constructor. 
      //! \param[in] _matrix: The transition state matrix for which to compute
      //!                     \f$^eG_{AF}(t\righarrow\infty)\f$
      //! \param[in] _tau: Missed event resolution.
      //! \param[in] _doopen: Whether to do AF or FA.
      Asymptotes(StateMatrix & _matrix, t_real _tau, bool _doopen=true);

      inline t_rmatrix H(t_real _s) const {
        return _Qaa + _Qaf * this->integral_(s) * _Qfa;
      }

      inline t_rmatrix operator(t_real _s) const { 
        return (_s * t_rmatrix::Identity(matrix_.nopen, matrix_.nopen) - H(_s)).determinant();
      }


    protected:
      //! Computes integral \f$\int_0^\tau\partial\,t\ e^{(Q_{FF} - sI)t}\f$
      t_rmatrix integral_(t_real _s);
      //! Computes \f$S_{FF}^{*}=I - e^{(Q_{FF} - sI)\tau}\f$
      t_rmatrix Sff_star(t_real _s);

    protected:
      //! The transition state matrix on which to act.
      StateMatrix matrix_;
      //! The eigenvalues of the ff matrix. Computed once.
      t_rvector ff_eigenvalues;
      //! The eigenvectors of the ff matrix. Computed once.
      t_rmatrix ff_eigenvectors;
  };


  //! \brief Computes H(s) with no frills
  //! \details This is meant to go in some interface later on. Just not sure what it should look
  //!          like just yet. The equation is 
  //!          \f\[Q_{AA}+ Q_{AF} \left[\int_0^\tau\partial\,t\ e^{(Q_{FF} - sI)t}\right]Q_{FA}\f\].
  //! \param[in] _Qaa: The "aa" partition of the transition matrix
  //! \param[in] _Qaf: The "af" partition of the transition matrix
  //! \param[in] _eigenvalues_ff: The eigenvalues of the "ff" partition of the transition matrix
  //! \param[in] _trans_ff: The transform to diagonalise the "ff" partition of the transition matrix
  //! \param[in] _invtrans_ff: The inverse transform to diagonalise the "ff" partition of the
  //!                          transition matrix. We should have
  //!                          \f$Q_{FF} = T \cdot \mathrm{diagm}(\epsilon) \cdot T^{-1}\f$.
  //! \param[in] _Qfa: The "fa" partition of the transition matrix
  //! \param[in] _s: the value for which to compute H(s).
  //! \param[in] _tau: critical time.
  //! \param[in] _zero: Value below which a value is considered == 0. Used when evaluating integral.
  t_rmatrix compute_H(t_rmatrix const &_Qaa, 
                      t_rmatrix const &_Qaf,
                      t_rvector const &_eigenvalues_ff,
                      t_rmatrix const &_trans_ff,
                      t_rmatrix const &_invtrans_ff,
                      t_rmatrix const &_Qfa,
                      t_real const _s,
                      t_real const _tau, 
                      t_real const _zero = 1e-12);

}
#endif 
