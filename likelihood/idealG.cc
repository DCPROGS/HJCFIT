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

#include "DCProgsConfig.h"

#include <iostream>
#include "idealG.h"

namespace DCProgs {
  void IdealG::set(t_rmatrix const &_Q, t_uint const &_nopen) {
  
    if(_Q.rows() != _Q.cols()) throw errors::Domain("Transition matrix is not square.");
    if(_nopen + 1 > static_cast<t_uint>(_Q.rows()))
      throw errors::Domain("There should be at least one shut state.");
    if(_nopen < 1) throw errors::Domain("There should be at least one open state.");
    this->matrix = _Q;
    // Enforces row constraints.
    for(t_rmatrix::Index i(0); i < _Q.rows(); ++i) 
    {
      this->matrix(i, i) = 0;
      this->matrix(i, i) = -this->matrix.row(i).sum();
    }
    this->nopen = _nopen;
  }

  t_rmatrix IdealG::laplace_af(t_real s) const {
    t_rmatrix const Qaa( QMatrix::aa() );
    auto const stuff = s * t_rmatrix::Identity(Qaa.rows(), Qaa.rows()) - Qaa;
    Eigen::FullPivLU<t_rmatrix> pivotLU(stuff);
    if(not pivotLU.isInvertible()) {
      std::ostringstream sstr; 
      sstr << *this << "\n ***** " << s << " is an eigenvalue of Qaa.";
      throw errors::NotInvertible(sstr.str());
    }
    return pivotLU.inverse() * QMatrix::af();
  }
  t_rmatrix IdealG::laplace_fa(t_real s) const {
    t_rmatrix const Qff( QMatrix::ff() );
    auto const stuff = s * t_rmatrix::Identity(Qff.rows(), Qff.rows()) - Qff;
    Eigen::FullPivLU<t_rmatrix> pivotLU(stuff);
    if(not pivotLU.isInvertible()) {
      std::ostringstream sstr; 
      sstr << *this << "\n  ***** " << s << " is an eigenvalue of Qff.";
      throw errors::NotInvertible(sstr.str());
    }
    return pivotLU.inverse() * QMatrix::fa();
  }
 
  MSWINDOBE std::ostream & operator<< (std::ostream &_stream, IdealG const &_mat) {
    return _stream << "Ideal Likelihood:\n" 
                   << "=================\n\n" 
                   << "  * nopen: "  << _mat.get_nopen() << "\n"
                   << "  * matrix: " << DCProgs::numpy_io(_mat.get_matrix()) << "\n";
  }
}
