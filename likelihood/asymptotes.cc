#include <DCProgsConfig.h>

#include <sstream>
#include <iostream>

#include "asymptotes.h"

namespace DCProgs {

# ifdef HAS_CXX11_CONSTEXPR
    // Only the God of linkers knows why we need this declaration twice.
    constexpr t_real DeterminantEq :: ZERO;
# else
    // Only the God of linkers knows why we need this declaration twice.
    const t_real DeterminantEq :: ZERO = 1e-12;
# endif

  DeterminantEq :: DeterminantEq   (StateMatrix const & _matrix, t_real _tau, bool _doopen)
                                 : tau_(_tau), matrix_(_matrix), ff_eigenvalues_(),
                                   ff_eigenvectors_() {
    if(not _doopen) {
      t_int const nopen = matrix_.nopen;
      t_int const nclose = matrix_.matrix.rows() - nopen;
      t_rmatrix const aa = matrix_.ff();
      t_rmatrix const af = matrix_.fa();
      t_rmatrix const fa = matrix_.af();
      t_rmatrix const ff = matrix_.aa();
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
    ff_eigenvectors_inv_ = ff_eigenvectors_.inverse();
  }

  t_rmatrix DeterminantEq :: integral_(t_real _s) const {
 
    t_rvector const alpha = ff_eigenvalues_.array() - _s;
    t_rmatrix diagonal = t_rmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) 
      diagonal(i, i) = std::abs(alpha(i)) > ZERO ?  (std::exp(alpha(i) * tau_) - 1e0) / alpha(i): tau_; 
    return ff_eigenvectors_ * diagonal * ff_eigenvectors_inv_;
  }

  t_rmatrix DeterminantEq::s_derivative(t_real _s) const { 

    t_rvector const alpha = ff_eigenvalues_.array() - _s;
    t_rmatrix diagonal = t_rmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) {
      if(std::abs(alpha(i)) > ZERO) {
        t_real const invalpha = 1e0 / alpha(i);
        diagonal(i, i) = invalpha * ((invalpha - tau_) * std::exp(alpha(i)*tau_) - invalpha);
      } else diagonal(i, i) = -tau_ * tau_ * 0.5;
    }
    return this->id_() +
           matrix_.af() * ff_eigenvectors_ * diagonal * ff_eigenvectors_inv_ * matrix_.fa(); 
  }

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, DeterminantEq const & _self) {
    
    return _stream << "Determinant equation:\n"
                   << "=====================\n\n" 
                   << "  * Transition Rate matrix\n" << _self.matrix_.matrix << "\n"
                   << "  * Number of open states " << _self.matrix_.nopen << "\n"
                   << "  * eigenvalues: " << _self.ff_eigenvalues_.transpose() << "\n";
  }

}
extern "C" void * create_determinant_eq(int _n0, int _n1, double *_matrix,  int _nopen, double _tau,
                                        bool _doopen) {
  using namespace DCProgs;
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> map(_matrix, _n0, _n1);
  StateMatrix states(map.cast<t_rmatrix::Scalar>(), _nopen);
  void *result = (void*) new DeterminantEq(states, _tau, _doopen);
  return result;
}
extern "C" void delete_determinant_eq(void *_self) {
  using namespace DCProgs;
  DeterminantEq *self = static_cast<DeterminantEq*>(_self);
  delete self;
}
extern "C" double call_determinant_eq(void *_self, double _s) {
  return (*static_cast<DCProgs::DeterminantEq*>(_self))(_s);
}

extern "C" char const * str_determinant_eq(void *_self) {
  std::ostringstream sstr;
  sstr << (*static_cast<DCProgs::DeterminantEq*>(_self));
  return sstr.str().c_str();
}
