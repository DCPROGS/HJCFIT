#include <DCProgsConfig.h>

#include <sstream>
#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>

#include "errors.h"
#include "exactG.h"

namespace DCProgs {

  void ExactG :: set(StateMatrix const &_matrix, t_real _tau) {
    // Two step process. Otherwise, reset would catch any exception thrown. 
    RecursionInterface afinterface(_matrix, _tau, true);
    RecursionInterface fainterface(_matrix, _tau, false);
    recursion_af_.reset(new RecursionInterface(std::move(afinterface)));
    if(not recursion_af_.get()) throw errors::Runtime("Could not initialize unique_ptr");
    recursion_fa_.reset(new RecursionInterface(std::move(fainterface)));
    if(not recursion_fa_.get()) throw errors::Runtime("Could not initialize unique_ptr");

    tau_ = _tau;
    af_factor_ = _matrix.af() * (_tau * _matrix.ff()).exp();
    fa_factor_ = _matrix.fa() * (_tau * _matrix.aa()).exp();
  }


  ExactG :: RecursionInterface::RecursionInterface( StateMatrix const & _matrix,
                                                    t_real _tau, bool _doAF ) {
                   
    // Sets matrix depending on whether this is AF or FA stuff.
    StateMatrix const transitions = _doAF ? _matrix: _matrix.transpose();

    // Solves eigenvalue problem
    Eigen::EigenSolver<t_rmatrix> eigsolver(transitions.matrix);
    if(eigsolver.info() != Eigen::Success) 
        throw errors::Mass("Could not solve eigenvalue problem.");

    // Initializes eigenvalues
    eigenvalues_ = -eigsolver.eigenvalues().real();

    // Initializes recursion formula for m == l == 0
    t_rmatrix const eigenvectors = eigsolver.eigenvectors().real();
    t_rmatrix const eigenvectors_inv = eigenvectors.inverse();
    for(t_int i(0); i < nbeigvals(); ++i) {
      auto left = eigenvectors.col(i).head(transitions.nopen);
      auto right = eigenvectors_inv.row(i).head(transitions.nopen);
      coeff_map_[ std::make_tuple(i, 0, 0) ] = left * right;
    }

    // Computes all Di values
    t_rmatrix const exponential_factor = (_tau * transitions.ff()).exp() * transitions.fa();
    for(t_int i(0); i < nbeigvals(); ++i) {
      auto left = eigenvectors.col(i).head(transitions.nopen);
      auto right = eigenvectors_inv.row(i).tail(transitions.nclose());
      dvalues_.push_back((left * right) * exponential_factor);
    }
    // set number of open states
    nopen = transitions.nopen;
  }

  // Recursion element i, m, l.
  ExactG::RecursionInterface::t_element
    ExactG::RecursionInterface::operator()(t_int _i, t_int _m, t_int _l) {

      assert(_i >= 0 and _m >= 0 and _l >= 0);
      assert(_i < nbeigvals());

      // Checks for existence in cache.
      t_key const key(_i, _m, _l);
      std::map<t_key, t_element>::const_iterator const i_found = coeff_map_.find(key);
      if(i_found != coeff_map_.end()) return i_found->second;

      // Otherwise compute it from recursion
      t_int const nopen_ = this->nopen;
      auto zero = [&nopen_]() { return t_element::Zero(nopen_, nopen_); };
      t_element const result = recursion_formula(*this, _i, _m, _l, zero);
      return result;
    }

  namespace {
    
    //! \brief Computes \f$B_{im}(t) = \sum_{r=0}^m C_{imr}t^r\f$
    //! \details See Theorem below equation 3.12
    template<class T> 
      typename T::t_element B_im_of_t(T &_C, t_int _i, t_int _m, t_real _t) {
        t_rmatrix result = _C(_i, _m, 0);
        t_real t(_t);
        for(t_int r(1); r <= _m; ++r, t *= _t) result += _C(_i, _m, r) * t;
        return result;
      }
    //! \brief Computes \f$M_m(t) = \sum_{i=1}^k B_{im}(t) e^{-\lambda_i t}\f$
    //! \details See Theorem below equation 3.12
    template<class T> 
      typename T::t_element M_m_of_t(T &_C, t_int _m, t_real _t) {

        // Eigenvalues == 0 are special (eq. solutions) and should be avoided.
        t_int i(0);
        for(; i < _C.nbeigvals(); ++i) 
          if( std::abs(_C.get_eigvals(i)) > 1e-10 ) break;

        t_rmatrix result = B_im_of_t(_C, i, _m, _t) * std::exp(-_C.get_eigvals(i)*_t);
        for(++i; i < _C.nbeigvals(); ++i)
          if( std::abs(_C.get_eigvals(i)) > 1e-10 ) 
            result += B_im_of_t(_C, i, _m, _t) * std::exp(-_C.get_eigvals(i)*_t);
        return result;
      }
    
    //! \brief Computes \f$R(t) = \sum_{m=0}^{t-m\tau>0} M_{m}(t-m\tau)\f$
    //! \details Equation 3.12 of HJC (1990)
    template<class T> 
      typename T::t_element R_of_t(T &_C, t_real _t, t_real _tau) {

        assert(_t > 0);

        t_real current_t(_t);
        t_rmatrix result = M_m_of_t(_C, 0, _t);
        t_int m=1;
        
        for(current_t -= _tau; current_t > 0; current_t -= _tau, ++m)  {
          if(m % 2 == 0) result += M_m_of_t(_C, m, current_t);
          else           result -= M_m_of_t(_C, m, current_t);
        }
        return result;
      }

  }

  t_rmatrix ExactG :: af(t_real _t) const {
    assert(_t > tau_);
    return R_of_t(*recursion_af_, _t - tau_, tau_) * af_factor_; 
  }
  t_rmatrix ExactG :: fa(t_real _t) const {
    assert(_t > tau_);
    return R_of_t(*recursion_fa_, _t - tau_, tau_) * fa_factor_; 
  }
}
