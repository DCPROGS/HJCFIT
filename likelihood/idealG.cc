#include <iostream>
#include "idealG.h"

namespace DCProgs {
  void IdealG::set(t_rmatrix const &_Q, t_int const &_nopen) {
  
    if(_Q.rows() != _Q.cols()) throw errors::Domain("Transition matrix is not square.");
    if(_nopen > _Q.rows() - 1)
      throw errors::Domain("There should be at least one shut state.");
    if(_nopen < 1) throw errors::Domain("There should be at least one open state.");
    this->matrix = _Q;
    // Enforces row constraints.
    for(size_t i(0); i < _Q.rows(); ++i) 
    {
      this->matrix(i, i) = 0;
      this->matrix(i, i) = -this->matrix.row(i).sum();
    }
    this->nopen = _nopen;
  }

  t_rmatrix IdealG::ff(t_real t) const {
    long const N{this->matrix.rows() - this->nopen};
    return t_rmatrix::Zero(N, N);
  }
  t_rmatrix IdealG::laplace_af(t_real s) const {
    t_rmatrix const Qaa( StateMatrix::aa() );
    t_rmatrix const Qaf( StateMatrix::af() );
    auto const stuff = s * t_rmatrix::Identity(Qaa.rows(), Qaa.rows()) - Qaa;
    Eigen::FullPivLU<t_rmatrix> pivotLU(stuff);
    if(not pivotLU.isInvertible()) throw errors::NotInvertible("Found pole of laplacian.");
    return pivotLU.inverse() * Qaf;
  }
  t_rmatrix IdealG::laplace_fa(t_real s) const {
    t_rmatrix const Qff( StateMatrix::ff() );
    t_rmatrix const Qfa( StateMatrix::fa() );
    Eigen::FullPivLU<t_rmatrix> pivotLU(s * t_rmatrix::Identity(Qff.rows(), Qff.rows()) - Qff);
    if(not pivotLU.isInvertible()) throw errors::NotInvertible("Found pole of laplacian.");
    return pivotLU.inverse() * Qfa;
  }
}
