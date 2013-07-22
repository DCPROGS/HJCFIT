#include <iostream>
#include "qmatrix.h"

namespace DCProgs {
  QMatrix QMatrix::transpose() const {
  
    QMatrix result;
    t_int const nclosed(nshut());
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
    return _stream << "Transition matrix with " << _mat.nopen
                   << " open states:\n" << DCProgs::numpy_io(_mat.matrix) << "\n";
  }
}