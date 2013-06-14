#ifndef DCPROGS_LIKELIHOOD_ASYMPTOTES_H
#define DCPROGS_LIKELIHOOD_ASYMPTOTES_H

#include <DCProgsConfig.h>
#include <ostream>
#include "state_matrix.h"

namespace DCProgs {

  class DeterminantEq;

  //! Dumps Determinantal equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream&, DeterminantEq const &);

  //! A functor to compute asymptotic missed event G.
  //! \detail The whole implementation is done w.r.t. to AF transitions. 
  //!         However, in practice, this is sufficient to compute FA transitions as well, by messing
  //!         with the input matrix.
  class MSWINDOBE DeterminantEq {
    
    friend std::ostream& operator<<(std::ostream&, DeterminantEq const &);

    public:
      //! Constructor. 
      //! \param[in] _matrix: The transition state matrix for which to compute
      //!                     \f$^eG_{AF}(t\rightarrow\infty)\f$
      //! \param[in] _tau: Missed event resolution.
      //! \param[in] _doopen: Whether to do AF or FA.
      DeterminantEq(StateMatrix const & _matrix, t_real _tau, bool _doopen=true);
      //! Copy constructor
      DeterminantEq   (DeterminantEq const & _c)
                    : tau_(_c.tau_), matrix_(_c.matrix_), ff_eigenvalues_(_c.ff_eigenvalues_),
                      ff_eigenvectors_(_c.ff_eigenvectors_), 
                      ff_eigenvectors_inv_(_c.ff_eigenvectors_inv_) {}

      //! Computes \f$Q_{AA} + Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}\f$
      //! \param[in] _s: Value of the laplacian scale.
      inline t_rmatrix H(t_real _s) const {
        return matrix_.aa() + matrix_.af() * this->integral_(_s) * matrix_.fa();
      }
      //! Computes \f$Q_{AA} + Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}\f$
      //! \param[in] _s: Value of the laplacian scale.
      //! \param[in] _tau: Value of tau for duration of call.
      inline t_rmatrix H(t_real _s, t_real _tau) const {
        return DeterminantEq(*this, _tau).H(_s);
      }

      //! Computes the determinant \f$\mathrm{det}(sI - H(s))\f$
      //! \param[in] _s: Value of the laplacian scale.
      inline t_real operator()(t_real _s) const { 
        return (_s * this->id_() - H(_s)).determinant();
      }
      //! Computes the determinant \f$\mathrm{det}(sI - H(s))\f$
      //! \param[in] _s: Value of the laplacian scale.
      //! \param[in] _tau: Value of tau for duration of call.
      inline t_real operator()(t_real _s, t_real _tau) const { 
        return (_s * this->id_() - H(_s, _tau)).determinant();
      }
      //! Derivative along _s
      t_rmatrix s_derivative(t_real _s) const;
      //! Derivative along _s
      inline t_rmatrix s_derivative(t_real _s, t_real _tau) const {
        return DeterminantEq(*this, _tau).s_derivative(_s);
      }

      //! Get resolution
      t_real get_tau() const { return tau_; }
      //! Set resolution
      void set_tau(t_real const &_tau) { tau_ = _tau; }

    protected:
      //! Computes integral \f$\int_0^\tau\partial\,t\ e^{(Q_{FF} - sI)t}\f$
      t_rmatrix integral_(t_real _s) const;
      //! Just the identity, just to write shorter code.
      inline auto id_() const ->decltype(t_rmatrix::Identity(1, 1)) 
        { return t_rmatrix::Identity(matrix_.nopen, matrix_.nopen); }


    private:
      //! Copy constructor for changing tau in constant functions.
      DeterminantEq   (DeterminantEq const & _c, t_real _tau)
                    : tau_(_tau), matrix_(_c.matrix_), ff_eigenvalues_(_c.ff_eigenvalues_),
                      ff_eigenvectors_(_c.ff_eigenvectors_), 
                      ff_eigenvectors_inv_(_c.ff_eigenvectors_inv_) {}

    protected:
      //! Time below which events are missed
      t_real tau_;
      //! The transition state matrix on which to act.
      StateMatrix matrix_;
      //! The eigenvalues of the ff matrix. Computed once.
      t_rvector ff_eigenvalues_;
      //! The eigenvectors of the ff matrix. Computed once.
      t_rmatrix ff_eigenvectors_;
      //! The inverse eigenvectors of the ff matrix. Computed once.
      t_rmatrix ff_eigenvectors_inv_;
#     ifdef HAS_CXX11_CONSTEXPR
        //! Hard coded static constant zero.
        constexpr static t_real ZERO = 1e-12;
#     else
        //! Hard coded static constant zero.
        const static t_real ZERO;
#     endif
  };

}
extern "C" void * create_determinant_eq(int _n0, int _n1, double *_matrix,  int _nopen, double _tau, bool _doopen);
extern "C" void delete_determinant_eq(void *_self);
extern "C" double call_determinant_eq(void *_self, double _s);
extern "C" char const * str_determinant_eq(void *_self);
#endif 
