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

#ifndef DCPROGS_LIKELIHOOD_RECURSION_FORMULA_H
#define DCPROGS_LIKELIHOOD_RECURSION_FORMULA_H

#include <DCProgsConfig.h>

#include "errors.h"

namespace DCProgs {

  //| \cond
  namespace details {
    template<class T> 
      typename T::t_element general(T & _C, t_uint _i, t_uint _m, t_uint _l);
    template<class T> 
      typename T::t_element lzero(T & _C, t_uint _i, t_uint _m);
  }
  //| \endcond

  //! \brief Obtains _C[_i, _m, _l] if prior terms are known.
  //! \details This function implements the recursion with as few requirements as possible. The
  //!          objective is to make the recursion clearer and testing easier. It also means we
  //!          separate concerns, such as caching prior results, or the possibility of using
  //!          arbitrary precision types.
  //!
  //!          Implements eq 3.18 from Hawkes, Jalali, Colquhoun (1990).
  //! \tparam T A type with the following form:
  //!    \code{.cpp}
  //!      class T {
  //!        // Type of the element
  //!        typedef t_element;
  //!
  //!        // Returns (prior) element in recursion
  //!        t_element operator()(t_uint _i, t_uint _j, t_uint _m);
  //!        // Returns D objects, e.g. $A_{iAF}e^{Q_{FF}\tau}Q_{FA}$.
  //!        t_element getD(t_uint _i) const;
  //!        // Returns specific eigenvalue of $Q$.
  //!        t_real get_eigval(t_uint _i) const;
  //!        // Returns number of eigenvalues.
  //!        t_uint nbeigval(t_uint _i) const;
  //!      };
  //!    \endcode
  //! \param _C The object over which the recursion is perfomed.
  //! \param _i Index to the eigenvalues
  //! \param _m Index to the time interval \f$(m-1)\tau < t < m\tau\f$
  //! \param _l An integer \f$0 \leq l \leq m\f$
  template<class T> 
    typename T::t_element recursion_formula(T & _C, t_uint _i, t_uint _m, t_uint _l) {
      
      assert(_m >= 0);
      assert(_i >= 0);
      assert(_l >= 0);
      assert(_l <= _m);
      assert(_i < _C.nbeigvals());
      // first, deal with _m == 0 and _l == 0 case.
      if(_m == 0 and _l == 0) return _C(_i, 0, 0);
      // then deals with two _l = 0 case
      if(_l == 0) return details::lzero(_C, _i, _m);
      // then deals with _l == _m
      if(_l == _m) return _C.getD(_i) * _C(_i, _m-1, _m-1) / t_real(_m);

      return details::general(_C, _i, _m, _l);
    }



  //| \cond
  namespace details {

    // This is basically an iterator that loops simultaneously over r and eigenvalues terms in
    // C(i, m, 0). 
    // The advantage of this kind of implementation is that we do not need to initialise return
    // argument. Other implementations would. See history of this file.
    template<class T> 
      class LZero {
        public:
          LZero   (T &_C, t_uint _i, t_uint _m, t_real _tolerance=1e-8) 
                : C_(_C), i_(_i), m_(_m), j_(-1), r_(_m), tolerance_(_tolerance) {}
          bool next() {
            // This is outer loop over eigenvalues.
            if(r_ >= m_) {
         
              ++j_;
              if(j_ >= C_.nbeigvals()) return false;
              else if(i_ == j_) return next();
         
              t_real const diff_lambda(C_.get_eigvals(j_) - C_.get_eigvals(i_));
              if(std::abs(diff_lambda) < tolerance_) return next();
         
              // set some values before moving on to inner loop over r_.
              inv_diff_lambda_ = 1e0 / diff_lambda;
              factor_ = inv_diff_lambda_;
              r_ = 0;
            } else {  // increment r_ loop.
              ++r_;
              if(r_ >= m_) return next(); // reached end of r loop.
              factor_ *= inv_diff_lambda_ * r_;
            }
            return true;
          };
          // Shouldn't be called if next returns false.
          typename T::t_element current() const {
            assert(j_ >= 0);
            assert(i_ >= 0);
            assert(m_ - 1 >= 0);
            assert(r_ >= 0);
            assert(i_ < C_.nbeigvals());
            assert(j_ < C_.nbeigvals());
            if(r_ % 2 == 0) 
              return (C_.getD(i_) * C_(j_, m_-1, r_) + C_.getD(j_) * C_(i_, m_-1, r_)) * factor_;
            return (C_.getD(i_) * C_(j_, m_-1, r_) - C_.getD(j_) * C_(i_, m_-1, r_)) * factor_;
          };
        protected:
          T & C_;
          t_uint i_, m_, j_, r_;
          t_real tolerance_;
          t_real inv_diff_lambda_;
          t_real factor_;
      };
    //! Recursion formula for l = 0
    template<class T> 
      typename T::t_element lzero(T & _C, t_uint _i, t_uint _m) {

        LZero<T> functor(_C, _i, _m); 
        if(not functor.next()) throw errors::Runtime("Expected to have something to do.");
        typename T::t_element result = functor.current();
        while(functor.next()) result += functor.current();
        return result;
      }

    //! Recursion formula for l != 0 and l != m
    template<class T> 
      typename T::t_element general(T & _C, t_uint _i, t_uint _m, t_uint _l) {

        typename T::t_element result(_C.getD(_i) * _C(_i, _m-1, _l-1) / t_real(_l));
        for(t_uint j(0); j < _C.nbeigvals(); ++j) {
          if(_i == j) continue;
 
          t_real const diff_lambda(_C.get_eigvals(_i)-_C.get_eigvals(j));
          if(std::abs(diff_lambda) < 1e-12) continue;
 
          t_real const diff_lambda_inv(1e0/diff_lambda); 
          t_real factor(diff_lambda_inv); 
          typename T::t_element intermediate = _C(_i, _m-1, _l) * factor;
          for(t_uint r(_l+1); r < _m; ++r) {
            factor *= diff_lambda_inv * r;
            intermediate += _C(_i, _m-1, r) * factor;
          }
          result -= _C.getD(j) * intermediate;
        } // loop over j
  
        return result;
      }
  }
  //| \endcond
}
#endif 
