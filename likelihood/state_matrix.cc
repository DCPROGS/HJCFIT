#include <iostream>
#include "state_matrix.h"

namespace DCProgs {
  StateMatrix StateMatrix::transpose() const {
  
    StateMatrix result;
    t_int const nclose(matrix.rows() - nopen);
    result.nopen = nclose;
    result.matrix.resize(matrix.rows(), matrix.cols());
    result.matrix.topLeftCorner(nclose, nclose) = matrix.bottomRightCorner(nclose, nclose);
    result.matrix.topRightCorner(nclose, nopen) = matrix.bottomLeftCorner(nclose, nopen);
    result.matrix.bottomLeftCorner(nopen, nclose) = matrix.topRightCorner(nopen, nclose);
    result.matrix.bottomRightCorner(nopen, nopen) = matrix.topLeftCorner(nopen, nopen);
    return result;
  }

  std::tuple<t_cvector, t_cmatrix> StateMatrix::eigenstuff() const {
     Eigen::EigenSolver<t_rmatrix> eigsolver(matrix.transpose());
     if(eigsolver.info() != Eigen::Success) 
       throw errors::Mass("Could not solve eigenvalue problem.");
     t_cvector const eigs = eigsolver.eigenvalues();
     t_cmatrix const vecs = eigsolver.eigenvectors();
     return std::make_tuple(eigs, vecs.transpose());
  }
}
