/***********************
    DCProgs computes missed-events likelihood as described in
    Hawkes, Jalali and Colquhoun (1990, 1992)

    Copyright (C) 2013  University College London

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
************************/

#include <DCProgsConfig.h>

#include <sstream>
#include <iostream>
#include <algorithm>
#include <vector>

#include <unsupported/Eigen/MatrixFunctions>

#include "errors.h"
#include "exact_survivor.h"

namespace DCProgs {

 
# if defined(_DEBUG) || defined(DEBUG)
    // anonymous namespace -- functions inside are not expected in the shared library.
    namespace {
      // Computes lowest two eigenvalues,
      // Checks they are real
      // Checks they are positive
      // Checks they are not equal,
      // Checks the lowest is below tolerance,
      void sanity_check(t_cvector const &_eigs, t_real _tolerance=1e-8) {
    
        // NOTE: neg_eig will be -_eigs.real(). Watch sign. 
        std::vector<t_real> neg_eigs(_eigs.size());
        std::vector<t_real> :: iterator i_first = neg_eigs.begin();
        t_cvector::Scalar const * ptr_eig = &_eigs[0];
        t_cvector::Scalar const * const ptr_eig_end = ptr_eig + _eigs.size();
        for(; ptr_eig != ptr_eig_end; ++ptr_eig, ++i_first) {
          if(std::abs(ptr_eig->imag()) > _tolerance) 
            throw errors::ComplexEigenvalues("Transition matrix has complex eigenvalues.");
          *i_first = -ptr_eig->real();
        }
    
        std::partial_sort(neg_eigs.begin(), neg_eigs.begin() + 2, neg_eigs.end());
        if( neg_eigs.front() < -_tolerance )
          throw errors::Mass("Expected eigenvalues of transition matrix to be negative.");
        if( neg_eigs[0] > _tolerance)
          throw errors::Mass("Kernel of transition matrix is of dimension 0.");
        if( neg_eigs[1] - neg_eigs[0] < _tolerance)
          throw errors::Mass("Kernel of transition matrix is larger than 1.");
      }
    }
#  endif

  void ExactSurvivor :: set(QMatrix const &_qmatrix, t_real _tau) {
    if(_tau <= 0e0) throw errors::Domain("The resolution time tau cannot be zero or negative.");

    // Two step process. Otherwise, reset would catch any exception thrown. 
    RecursionInterface afinterface(_qmatrix, _tau, true);
    RecursionInterface fainterface(_qmatrix, _tau, false);
    recursion_af_.reset(new RecursionInterface(std::move(afinterface)));
    if(not recursion_af_.get()) throw errors::Runtime("Could not initialize unique_ptr");
    recursion_fa_.reset(new RecursionInterface(std::move(fainterface)));
    if(not recursion_fa_.get()) throw errors::Runtime("Could not initialize unique_ptr");

    tau_ = _tau;
  }


  ExactSurvivor :: RecursionInterface::RecursionInterface( QMatrix const & _qmatrix,
                                                           t_real _tau, bool _doAF ) {
                   
    // Sets matrix depending on whether this is AF or FA stuff.
    QMatrix const transitions = _doAF ? _qmatrix: _qmatrix.transpose();

    // Solves eigenvalue problem
    Eigen::EigenSolver<t_rmatrix> eigsolver(transitions.matrix);
    if(eigsolver.info() != Eigen::Success) 
        throw errors::Mass("Could not solve eigenvalue problem.");

    // Initializes eigenvalues
#   if defined(_DEBUG) || defined(DEBUG)
      sanity_check(eigsolver.eigenvalues());
#   endif
    eigenvalues_ = -eigsolver.eigenvalues().real();

    // Initializes recursion formula for m == l == 0
    t_rmatrix const eigenvectors = eigsolver.eigenvectors().real();
    t_rmatrix const eigenvectors_inv = eigsolver.eigenvectors().inverse().real();
    for(t_rvector::Index i(0); i < eigenvalues_.size(); ++i) {
      auto left = eigenvectors.col(i).head(transitions.nopen);
      auto right = eigenvectors_inv.row(i).head(transitions.nopen);
      coeff_map_[ std::make_tuple(i, 0, 0) ] = left * right;
    }

    // Computes all Di values
    t_rmatrix const exponential_factor = (_tau * transitions.ff()).exp() * transitions.fa();
    for(t_rvector::Index i(0); i < eigenvalues_.size(); ++i) {
      auto left = eigenvectors.col(i).head(transitions.nopen);
      auto right = eigenvectors_inv.row(i).tail(transitions.nshut());
      dvalues_.push_back((left * right) * exponential_factor);
    }
    // set number of open states
    nopen = transitions.nopen;
  }

  // Recursion element i, m, l.
  ExactSurvivor::RecursionInterface::t_element
    ExactSurvivor::RecursionInterface::operator()(t_uint _i, t_uint _m, t_uint _l) {

      assert(_i >= 0 and _m >= 0 and _l >= 0);
      assert(_i < nbeigvals());

      // Checks for existence in cache.
      t_key const key(_i, _m, _l);
      std::map<t_key, t_element>::const_iterator const i_found = coeff_map_.find(key);
      if(i_found != coeff_map_.end()) return i_found->second;

      // Otherwise compute it from recursion
      return recursion_formula(*this, _i, _m, _l);
    }

  namespace {
    
    //! \brief Computes \f$B_{im}(t) = \sum_{r=0}^m C_{imr}t^r\f$
    //! \details See Theorem below equation 3.12
    template<class T> 
      typename T::t_element B_im_of_t(T &_C, t_uint _i, t_uint _m, t_real _t) {
        t_rmatrix result = _C(_i, _m, 0);
        t_real t(_t);
        for(t_uint r(1); r <= _m; ++r, t *= _t) result += _C(_i, _m, r) * t;
        return result;
      }
    //! \brief Computes \f$M_m(t) = \sum_{i=1}^k B_{im}(t) e^{-\lambda_i t}\f$
    //! \details See Theorem below equation 3.12
    template<class T> 
      typename T::t_element M_m_of_t(T &_C, t_uint _m, t_real _t) {

        t_rmatrix result = B_im_of_t(_C, 0, _m, _t) * std::exp(-_C.get_eigvals(0)*_t);
        for(t_uint i(1); i < _C.nbeigvals(); ++i)
          result += B_im_of_t(_C, i, _m, _t) * std::exp(-_C.get_eigvals(i)*_t);
        return result;
      }
    
    //! \brief Computes \f$R(t) = \sum_{m=0}^{t-m\tau>0} M_{m}(t-m\tau)\f$
    //! \details Equation 3.12 of HJC (1990)
    template<class T> 
      typename T::t_element R_of_t(T &_C, t_real _t, t_real _tau) {

        t_real current_t(_t);
        t_rmatrix result = M_m_of_t(_C, 0, _t);
        t_uint m=1;
        
        for(current_t -= _tau; current_t > 0; current_t -= _tau, ++m)  {
          if(m % 2 == 0) result += M_m_of_t(_C, m, current_t);
          else           result -= M_m_of_t(_C, m, current_t);
        }
        return result;
      }

  }

  t_rmatrix ExactSurvivor :: af(t_real _t) const {
    if(_t < 0e0) return recursion_af_->zero();
    return R_of_t(*recursion_af_, _t, tau_);
  }
  t_rmatrix ExactSurvivor :: fa(t_real _t) const {
    if(_t < 0e0) return recursion_fa_->zero();
    return R_of_t(*recursion_fa_, _t, tau_);
  }

  t_rmatrix ExactSurvivor :: recursion_af(t_uint _i, t_uint _m, t_uint _l) const {
    if(_i >= recursion_af_->nbeigvals())  
      throw errors::Index("i index should be smaller than the number of eigenvalues.");
    if(_l > _m) throw errors::Index("l index should be smaller than m index.");
    return  recursion_af_->operator()(_i, _m, _l);
  }
  t_rmatrix ExactSurvivor :: recursion_fa(t_uint _i, t_uint _m, t_uint _l) const {
    if(_i >= recursion_fa_->nbeigvals())  
      throw errors::Index("i index should be smaller than the number of eigenvalues.");
    if(_l > _m) throw errors::Index("l index should be smaller than m index.");
    return  recursion_fa_->operator()(_i, _m, _l);
  }
  t_rmatrix ExactSurvivor :: D_af(t_uint _i) const {
    if(_i >= recursion_af_->nbeigvals())  
      throw errors::Index("i index should be smaller than the number of eigenvalues.");
    return  recursion_af_->getD(_i);
  }
  t_rmatrix ExactSurvivor :: D_fa(t_uint _i) const {
    if(_i >= recursion_fa_->nbeigvals())  
      throw errors::Index("i index should be smaller than the number of eigenvalues.");
    return  recursion_fa_->getD(_i);
  }
  t_rvector ExactSurvivor::eigenvalues_af() const { return recursion_af_->eigenvalues(); }
  t_rvector ExactSurvivor::eigenvalues_fa() const { return recursion_af_->eigenvalues(); }

  ExactSurvivor& ExactSurvivor :: operator=(ExactSurvivor &&_c) {
    recursion_af_ = std::move(_c.recursion_af_);
    recursion_fa_ = std::move(_c.recursion_fa_);
    tau_ = _c.tau_;
    return *this;
  }
}
