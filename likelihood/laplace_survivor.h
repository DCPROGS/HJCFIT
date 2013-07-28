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

#ifndef DCPROGS_LIKELIHOOD_LAPLACE_SURVIVOR_H
#define DCPROGS_LIKELIHOOD_LAPLACE_SURVIVOR_H

#include <DCProgsConfig.h>
#include <ostream>
#include "qmatrix.h"

namespace DCProgs {


  //! \brief Survivor functions $^{A}R(s)$ in Laplace space. 
  class MSWINDOBE LaplaceSurvivor {
    
    public:
      //! Constructor. 
      //! \param[in] _qmatrix: The transition state matrix for which to compute
      //!                     \f$^eG_{AF}(t\rightarrow\infty)\f$
      //! \param[in] _tau: Missed event resolution.
      LaplaceSurvivor(QMatrix const & _qmatrix);
      //! Copy constructor
      LaplaceSurvivor   (LaplaceSurvivor const & _c)
                      : qmatrix_(_c.qmatrix_), ff_eigenvalues_(_c.ff_eigenvalues_),
                        ff_eigenvectors_(_c.ff_eigenvectors_), 
                        ff_eigenvectors_inv_(_c.ff_eigenvectors_inv_) {}

      //! Computes \f$sI - Q_{AA} - Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}\f$
      //! \param[in] _s: Value of the laplacian scale.
      t_rmatrix H(t_real _s, t_real _tau) const {
        return qmatrix_.aa() + qmatrix_.af() * this->integral_(_s, _tau) * qmatrix_.fa();
      }
      //! Computes the determinant \f$\mathrm{det}(sI - H(s))\f$
      //! \param[in] _s: Value of the laplacian scale.
      t_rmatrix operator()(t_real _s, t_real _tau) const {
        return (_s * id_() - H(_s, _tau)).inverse(); 
      }
      //! Computes the matrix \f$W=\mathrm{det}(sI - H(s, \tau))\f$
      //! \param[in] _s: Value of the laplacian scale.
      t_rmatrix W(t_real _s, t_real _tau) const {
        return _s * this->id_() - LaplaceSurvivor::H(_s, _tau);
      }
      //! Derivative along of W along s
      t_rmatrix s_derivative(t_real _s, t_real _tau) const;

      //! \brief Returns the Q matrix.
      //! \details This is strictly a read-only function since changing the matrix has fairly far
      //! ranging implications.
      QMatrix const &get_qmatrix() const { return qmatrix_; }
      //! Get expected number of open-states.
      t_uint get_nopen() const { return qmatrix_.nopen; }


      //! Returns \f$^{F}R(s)\f$.
      LaplaceSurvivor transpose() const { return LaplaceSurvivor(qmatrix_.transpose()); }

      //! Returns eigenvalues of ff block.
      t_cvector get_ff_eigenvalues() const { return ff_eigenvalues_; }

    protected:
      //! Computes integral \f$\int_0^\tau\partial\,t\ e^{(Q_{FF} - sI)t}\f$
      t_rmatrix integral_(t_real _s, t_real _tau) const;
      //! Just the identity, just to write shorter code.
      inline auto id_() const ->decltype(t_rmatrix::Identity(1, 1)) 
        { return t_rmatrix::Identity(qmatrix_.nopen, qmatrix_.nopen); }

    protected:
      //! The transition state matrix on which to act.
      QMatrix qmatrix_;
      //! The eigenvalues of the ff matrix. Computed once.
      t_cvector ff_eigenvalues_;
      //! The eigenvectors of the ff matrix. Computed once.
      t_cmatrix ff_eigenvectors_;
      //! The inverse eigenvectors of the ff matrix. Computed once.
      t_cmatrix ff_eigenvectors_inv_;
#     ifdef HAS_CXX11_CONSTEXPR
        //! Hard coded static constant zero.
        constexpr static t_real ZERO = 1e-12;
#     else
        //! Hard coded static constant zero.
        const static t_real ZERO;
#     endif
  };

  
  //! Dumps survivor function equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, LaplaceSurvivor const &_self);
}
#endif 
