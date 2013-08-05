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

#include <iostream>
#include <sstream>
#include "qmatrix.h"

namespace DCProgs {
  QMatrix QMatrix::transpose() const {
  
    QMatrix result;
    t_uint const nclosed(nshut());
    result.nopen = nclosed;
    result.matrix.resize(matrix.rows(), matrix.cols());
    result.matrix.topLeftCorner(nclosed, nclosed) = matrix.bottomRightCorner(nclosed, nclosed);
    result.matrix.topRightCorner(nclosed, nopen) = matrix.bottomLeftCorner(nclosed, nopen);
    result.matrix.bottomLeftCorner(nopen, nclosed) = matrix.topRightCorner(nopen, nclosed);
    result.matrix.bottomRightCorner(nopen, nopen) = matrix.topLeftCorner(nopen, nopen);
    return result;
  }

  std::tuple<t_cvector, t_cmatrix> QMatrix::eigenstuff() const {
     Eigen::EigenSolver<t_rmatrix> eigsolver(matrix.transpose());
     if(eigsolver.info() != Eigen::Success) {
       std::ostringstream sstr("Could not solve eigenvalue problem.\n");
       sstr << *this << "\n";
       throw errors::Mass(sstr.str());
     }
     t_cvector const eigs = eigsolver.eigenvalues();
     t_cmatrix const vecs = eigsolver.eigenvectors();
     return std::make_tuple(eigs, vecs.transpose());
  }

  MSWINDOBE std::ostream & operator<< (std::ostream &_stream, QMatrix const &_mat) {
    std::ostringstream sstr;
    sstr << "Transition matrix with " << _mat.nopen << " open states:";
    _stream << sstr.str() << "\n";
    for(std::string::size_type i(0); i < sstr.str().size(); ++i) _stream << '-';
    return _stream << "\n" << DCProgs::numpy_io(_mat.matrix) << "\n";
  }
}
