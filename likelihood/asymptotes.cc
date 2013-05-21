#include <DCProgsConfig.h>

#include <iostream>
#include "asymptotes.h"

namespace DCProgs {

  Asymptotes :: Asymptotes   (StateMatrix & _matrix, t_real _tau, bool _doopen=true)
                           : matrix_(_matrix), tau_(_tau) {
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
    std::tie(ff_eigenvalues, ff_eigenvectors) = matrix_.eigenstuff();
  }

  t_rmatrix Asymptotes :: integral_(t_real _s) {
 
    t_rvector const alpha = ff_eigenvalues - _s;
    t_rmatrix diagonal = t_rmatrix::Zero(ff_eigenvalues.size(), ff_eigenvalues.size());
    for(t_int i(0); i < ff_eigenvalues.size(); ++i) {
      if(std::abs(alpha(i)) > _zero)
        diagonal(i, i) =  (std::exp(alpha(i) * tau_) - 1e0) / alpha(i);
    }
    return ff_eigenvectors_ * diagonal * ff_eigenvectors_.transpose();
  }

  t_rmatrix Asymptotes :: Sff_star(t_real _s) {

    t_rmatrix result = t_rmatrix::Zero(ff_eigenvalues.size(), ff_eigenvalues.size());
    result.diagonal() = 1e0 - ((ff_eigenvalues - _s) * tau_).exp();
    return ff_eigenvectors_ * result * ff_eigenvectors_.transpose();
  }

}
