#include "idealG.h"

namespace DCProgs {
  void IdealG::set(t_rmatrix const &_Q, t_partition const &_open) {
  
    if(_Q.rows() != _Q.cols()) throw errors::Domain("Transition matrix is not square.");
    if(_Q.rows() != _open.rows())
      throw errors::Domain("Size of partition and Q matrix do not match.");
    Q_ = _Q;
    // Enforces row constraints.
    for(size_t i(0); i < _Q.rows(); ++i) 
    {
      Q_(i, i) = 0;
      Q_(i, i) = -Q_.row(i).sum();
    }
    open_ = _open;
  }
}
