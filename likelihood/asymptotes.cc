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
    if(eigsolver.info() != Eigen::Success)  {
      std::ostringstream sstr("Could not solve eigenvalue problem.");
      sstr << numpy_io(matrix_.ff()) << "\n";
      throw errors::Mass(sstr.str());
    }
    ff_eigenvalues_ = eigsolver.eigenvalues();
    ff_eigenvectors_ = eigsolver.eigenvectors();
    ff_eigenvectors_inv_ = ff_eigenvectors_.inverse();
  }

  t_rmatrix DeterminantEq :: integral_(t_real _s) const {
 
    t_cvector const alpha = ff_eigenvalues_.array() - _s;
    t_cmatrix diagonal = t_cmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) 
      diagonal(i, i) = std::abs(alpha(i)) > ZERO ?  (std::exp(alpha(i) * tau_) - 1e0) / alpha(i): tau_; 
    t_cmatrix const result = ff_eigenvectors_ * diagonal * ff_eigenvectors_inv_;
    if((result.imag().array().abs() > 1e-8).any())
      throw errors::ComplexEigenvalues("Integral calculation yielded complex values.\n");
    return result.real();
  }

  t_rmatrix DeterminantEq::s_derivative(t_real _s) const { 

    t_cvector const alpha = ff_eigenvalues_.array() - _s;
    t_cmatrix diagonal = t_cmatrix::Zero(ff_eigenvalues_.size(), ff_eigenvalues_.size());
    for(t_int i(0); i < ff_eigenvalues_.size(); ++i) {
      if(std::abs(alpha(i)) > ZERO) {
        t_complex const invalpha = 1e0 / alpha(i);
        diagonal(i, i) = invalpha * ((invalpha - tau_) * std::exp(alpha(i)*tau_) - invalpha);
      } else diagonal(i, i) = -tau_ * tau_ * 0.5;
    }
    t_cmatrix const integral = ff_eigenvectors_ * diagonal * ff_eigenvectors_inv_;
    if((integral.imag().array().abs() > 1e-8).any())
      throw errors::ComplexEigenvalues("Integral calculation yielded complex values.\n");
    return this->id_() + matrix_.af() * integral.real() * matrix_.fa(); 
  }

  MSWINDOBE std::ostream& operator<<(std::ostream& _stream, DeterminantEq const & _self) {
    
    return _stream << "Determinant equation:\n"
                   << "=====================\n\n" 
                   << "  * Transition Rate matrix:\n" << numpy_io(_self.matrix_.matrix) << "\n"
                   << "  * Number of 'A' states: " << _self.matrix_.nopen << "\n"
                   << "  * Tau: " << _self.tau_ << "\n"
                   << "  * FF eigenvalues: " << _self.ff_eigenvalues_.transpose() << "\n";
  }
}
