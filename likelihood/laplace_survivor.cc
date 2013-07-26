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

#include <DCProgsConfig.h>

#include <sstream>
#include <iostream>

#include "laplace_survivor.h"

namespace DCProgs {

# ifdef HAS_CXX11_CONSTEXPR
    // Only the God of linkers knows why we need this declaration twice.
    constexpr t_real LaplaceSurvivor :: ZERO;
# else
    // Only the God of linkers knows why we need this declaration twice.
    const t_real LaplaceSurvivor :: ZERO = 1e-12;
# endif

  LaplaceSurvivor :: LaplaceSurvivor   (QMatrix const & _qmatrix)
                                     : qmatrix_(_qmatrix), ff_eigenvalues_(),
                                       ff_eigenvectors_() {
    Eigen::EigenSolver<t_rmatrix> eigsolver(qmatrix_.ff());
    if(eigsolver.info() != Eigen::Success)  {
      std::ostringstream sstr("Could not solve eigenvalue problem.");
      sstr << numpy_io(qmatrix_.ff()) << "\n";
      throw errors::Mass(sstr.str());
    }
    ff_eigenvalues_ = eigsolver.eigenvalues();
    ff_eigenvectors_ = eigsolver.eigenvectors();
    ff_eigenvectors_inv_ = ff_eigenvectors_.inverse();
  }

  t_rmatrix LaplaceSurvivor :: integral_(t_real _s, t_real _tau) const {
 
    t_cvector const alpha = ff_eigenvalues_.array() - _s;
    t_cmatrix diagonal = t_cmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) 
      diagonal(i, i) = std::abs(alpha(i)) > ZERO ? 
                         (std::exp(alpha(i) * _tau) - t_real(1e0)) / alpha(i):
                         _tau; 
    t_cmatrix const result = ff_eigenvectors_ * diagonal * ff_eigenvectors_inv_;
    if((result.imag().array().abs() > t_real(1e-8)).any())
      throw errors::ComplexEigenvalues("Integral calculation yielded complex values.\n");
    return result.real();
  }

  t_rmatrix LaplaceSurvivor :: s_derivative(t_real _s, t_real _tau) const { 

    t_cvector const alpha = ff_eigenvalues_.array() - _s;
    t_cmatrix diagonal = t_cmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) {
      if(std::abs(alpha(i)) > ZERO) {
        t_complex const invalpha = t_real(1e0) / alpha(i);
        diagonal(i, i) = invalpha * ((invalpha - _tau) * std::exp(alpha(i)*_tau) - invalpha);
      } else diagonal(i, i) = -_tau * _tau * t_real(0.5);
    }

    t_cmatrix const integral = ff_eigenvectors_ * diagonal * ff_eigenvectors_inv_;
    if((integral.imag().array().abs() > t_real(1e-8)).any())
      throw errors::ComplexEigenvalues("Integral calculation yielded complex values.\n");
    return this->id_() - qmatrix_.af() * integral.real() * qmatrix_.fa(); 
  }

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, LaplaceSurvivor const & _self) {
    return _stream << "Survivor function in Laplace space:\n"
                   << "===================================\n\n" 
                   << "  * Transition Rate matrix:\n" << numpy_io(_self.get_qmatrix().matrix) << "\n"
                   << "  * Number of 'A' states: " << _self.get_qmatrix().nopen << "\n"
                   << "  * FF eigenvalues: " << _self.get_ff_eigenvalues().transpose() << "\n";
  }
}
