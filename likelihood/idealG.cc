#include "idealG.h"

namespace DCProgs {
  void IdealG::set(t_rmatrix const &_Q, t_int const &_nopen) {
  
    if(_Q.rows() != _Q.cols()) throw errors::Domain("Transition matrix is not square.");
    if(_Q.rows() < _nopen)
      throw errors::Domain("Number of open-states greater than size of Q matrix.");
    if(_nopen < 0) throw errors::Domain("Negative number of open states.");
    Q_ = _Q;
    // Enforces row constraints.
    for(size_t i(0); i < _Q.rows(); ++i) 
    {
      Q_(i, i) = 0;
      Q_(i, i) = -Q_.row(i).sum();
    }
    nopen_ = _nopen;
  }
}
