#include "idealG.h"

namespace DCProgs {
  void IdealG::set(t_rmatrix const &_Q, t_int const &_nopen) {
  
    if(_Q.rows() != _Q.cols()) throw errors::Domain("Transition matrix is not square.");
    if(_Q.rows() < _nopen)
      throw errors::Domain("Number of open-states greater than size of Q matrix.");
    if(_nopen < 0) throw errors::Domain("Negative number of open states.");
    this->matrix = _Q;
    // Enforces row constraints.
    for(size_t i(0); i < _Q.rows(); ++i) 
    {
      this->matrix(i, i) = 0;
      this->matrix(i, i) = -this->matrix.row(i).sum();
    }
    this->nopen = _nopen;
  }
}
