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

#ifndef DCPROGS_LIKELIHOOD_DETERMINANT_EQUATION_H
#define DCPROGS_LIKELIHOOD_DETERMINANT_EQUATION_H

#include <DCProgsConfig.h>
#include <ostream>
#include "laplace_survivor.h"

namespace DCProgs {

  //! \brief A functor to compute the W matrix, so as to find its roots.
  //! \details The whole implementation is done w.r.t. to AF transitions. However, in practice,
  //!          this is sufficient to compute FA transitions as well, by messing with the input
  //!          matrix.
  class MSWINDOBE DeterminantEq : public LaplaceSurvivor {
    
    public:
      //! Constructor.
      DeterminantEq(QMatrix const & _qmatrix, t_real _tau) : LaplaceSurvivor(_qmatrix), tau_(_tau) {};
      //! Constructor.
      DeterminantEq(DeterminantEq const & _c) : LaplaceSurvivor(_c), tau_(_c.tau_) {};
      //! Constructor
      template<class T> 
        DeterminantEq   (Eigen::DenseBase<T> const &_matrix, t_uint _nopen, t_real _tau)
                      : LaplaceSurvivor(QMatrix(_matrix, _nopen)), tau_(_tau) {};

      //! Computes the determinant \f$\mathrm{det}(sI - H(s))\f$
      //! \param[in] _s Value of the laplacian scale.
      t_real operator()(t_real _s) const { return DeterminantEq::operator()(_s, tau_); }
      //! Computes the determinant \f$\mathrm{det}(sI - H(s, \tau))\f$
      //! \param[in] _s Value of the laplacian scale
      //! \param[in] _tau Maximum length of missed events
      t_real operator()(t_real _s, t_real _tau) const {
        return (_s * this->id_() - LaplaceSurvivor::H(_s, _tau)).determinant();
      }
      //! \brief Determinant equation for transpose matrix
      //! \details In other words, if looking at open states, then returns equation for shut states,
      //!          and vice-versa.
      DeterminantEq transpose() const { return DeterminantEq(get_qmatrix().transpose(), tau_); }

      //! Computes \f$sI - Q_{AA} - Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}\f$
      //! \param[in] _s Value of the laplacian scale.
      t_rmatrix H(t_real _s) const { return LaplaceSurvivor::H(_s, tau_); }
      using LaplaceSurvivor::H;

      //! Derivative of W along s
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

  //! Dumps Determinant equation to stream
  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, DeterminantEq const &_self);
}
#endif 
