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
    for(t_int i(0); i < _Q.rows(); ++i) 
    {
      this->matrix(i, i) = 0;
      this->matrix(i, i) = -this->matrix.row(i).sum();
    }
    this->nopen = _nopen;
  }

  t_rmatrix IdealG::laplace_af(t_real s) const {
    t_rmatrix const Qaa( StateMatrix::aa() );
    auto const stuff = s * t_rmatrix::Identity(Qaa.rows(), Qaa.rows()) - Qaa;
    Eigen::FullPivLU<t_rmatrix> pivotLU(stuff);
    if(not pivotLU.isInvertible()) {
      std::ostringstream sstr; 
      sstr << *this << "\nQaa has eigenvalue of 1.0";
      throw errors::NotInvertible(sstr.str());
    }
    return pivotLU.inverse() * StateMatrix::af();
  }
  t_rmatrix IdealG::laplace_fa(t_real s) const {
    t_rmatrix const Qff( StateMatrix::ff() );
    auto const stuff = s * t_rmatrix::Identity(Qff.rows(), Qff.rows()) - Qff;
    Eigen::FullPivLU<t_rmatrix> pivotLU(stuff);
    if(not pivotLU.isInvertible()) {
      std::ostringstream sstr; 
      sstr << *this << "\nQff has eigenvalue of 1.0";
      throw errors::NotInvertible(sstr.str());
    }
    return pivotLU.inverse() * StateMatrix::fa();
  }
 
  MSWINDOBE std::ostream & operator<< (std::ostream &_stream, IdealG const &_mat) {
    return _stream << "Ideal Likelihood:\n" 
                   << "-----------------" 
                   << "  - nopen: "  << _mat.get_nopen()
                   << "  - matrix: " << DCProgs::numpy_io(_mat.get_matrix()) << "\n";
  }
}
