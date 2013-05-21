#include <DCProgsConfig.h>

#include <iostream>
#include "asymptotes.h"

namespace DCProgs {

  // Only the God of linkers knows why we need this declaration twice.
  constexpr t_real DeterminantEq :: ZERO;

  DeterminantEq :: DeterminantEq   (StateMatrix & _matrix, t_real _tau, bool _doopen)
                                 : tau_(_tau), matrix_(_matrix), ff_eigenvalues_(),
                                   ff_eigenvectors_() {
    if(not _doopen) {
      t_int const nopen = matrix_.nopen;
      t_int const nclose = matrix_.matrix.rows() - nopen;
      t_rmatrix const aa = matrix_.aa();
      t_rmatrix const af = matrix_.af();
      t_rmatrix const fa = matrix_.fa();
      t_rmatrix const ff = matrix_.ff();
      matrix_.nopen = nclose;
      matrix_.aa() = aa;
      matrix_.af() = af;
      matrix_.ff() = ff;
      matrix_.fa() = fa;
    }
    Eigen::EigenSolver<t_rmatrix> eigsolver(matrix_.ff());
    if(eigsolver.info() != Eigen::Success) 
        throw errors::Mass("Could not solve eigenvalue problem.");
    ff_eigenvalues_ = eigsolver.eigenvalues().real();
    ff_eigenvectors_ = eigsolver.eigenvectors().real();
  }

  t_rmatrix DeterminantEq :: integral_(t_real _s) const {
 
    t_rvector const alpha = ff_eigenvalues_.array() - _s;
    auto const is_not_eig = alpha.array().abs() > ZERO;
    t_rmatrix diagonal = t_rmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) 
      diagonal(i, i) = is_not_eig(i) ?  (std::exp(alpha(i) * tau_) - 1e0) / alpha(i): tau_; 
    return ff_eigenvectors_ * diagonal * ff_eigenvectors_.transpose();
  }

  t_rmatrix DeterminantEq::s_derivative(t_real _s) const { 

    t_rvector const alpha = ff_eigenvalues_.array() - _s;
    auto const is_not_eig = alpha.array().abs() > ZERO;
    t_rmatrix diagonal = t_rmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) {
      if(is_not_eig(i)) {
        t_real const invalpha = 1e0 / alpha(i);
        diagonal(i, i) = invalpha * ((invalpha - tau_) * std::exp(alpha(i)*tau_) - invalpha);
      } else diagonal(i, i) = -tau_ * tau_ * 0.5;
    }
    return this->id_() +
           matrix_.af() * ff_eigenvectors_ * diagonal * ff_eigenvectors_.transpose() * matrix_.fa(); 
  }
}
