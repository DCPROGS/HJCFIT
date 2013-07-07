#ifndef DCPROGS_LIKELIHOOD_RECURSION_FORMULA_H
#define DCPROGS_LIKELIHOOD_RECURSION_FORMULA_H

#include <DCProgsConfig.h>

namespace DCProgs {

  namespace details {
    template<class T, class T_ZERO> 
      typename T::t_element general(T & _C, t_int _i, t_int _m, t_int _l, T_ZERO const &_zero);
    template<class T, class T_ZERO> 
      typename T::t_element lzero(T & _C, t_int _i, t_int _m, T_ZERO const &_zero);
  }

  //! \brief Obtains _C[_i, _m, _l] if prior terms are known.
  //! \details This function implements the recursion with as few requirements as possible. The
  //!          objective is to make the recursion clearer and testing easier. It also means we
  //!          separate concerns, such as caching prior results, or the possibility of using
  //!          arbitrary precision types.
  //!
  //!          Implements eq 3.18 from Hawkes, Jalali, Colquhoun (1990).
  //! \tparam T: A type with the following form:
  //!    \code{.cpp}
  //!      class T {
  //!        //! Type of the element
  //!        typedef t_element;
  //!
  //!        //! Returns (prior) element in recursion
  //!        t_element operator()(t_int _i, t_int _j, t_int _m);
  //!        //! \brief Returns D objects, e.g. \f$A_{iAF}e^{Q_{FF}\tau}Q_{FA}\f$.
  //!        auto getD(t_int _i) const;
  //!        //! Returns specific eigenvalue of \f$Q\f$.
  //!        t_real get_eigval(t_int _i) const;
  //!        //! Returns number of eigenvalues.
  //!        t_int nbeigval(t_int _i) const;
  //!      };
  //!    \endcode
  //! \tparam T_ZERO: Type of a functor.
  //! \param _C: The object over which the recursion is perfomed.
  //! \param _i: An integer
  //! \param _m: An integer
  //! \param _l: An integer
  //! \parma _zero: A functor used to initialise intermediate objects.
  template<class T, class T_ZERO> 
    typename T::t_element recursion_formula( T & _C, t_int _i, t_int _m, t_int _l,
                                             T_ZERO const &_zero ) {
      
      // first, deal with _m == 0 and _l == 0 case.
      if(_m == 0 and _l == 0) return _C(_i, 0, 0);
      // then deals with two _l = 0 case
      if(_l == 0) return details::lzero(_C, _i, _m, _zero);
      // then deals with _l == _m
      if(_l == _m) return _C.getD(_i) * _C(_i, _m-1, _m-1) / t_real(_m);

      return details::general(_C, _i, _m, _l, _zero);
    }



  namespace details {

    //! Recursion formula for l = 0
    template<class T, class T_ZERO> 
      typename T::t_element lzero(T & _C, t_int _i, t_int _m, T_ZERO const &_zero) {

        typename T::t_element result = _zero();
        auto Di = _C.getD(_i);
        t_real const lambda_i = _C.get_eigvals(_i); 

        for(t_int j(0); j < _C.nbeigvals(); ++j) {
          if(_i == j) continue;

          t_real const lambda_j = _C.get_eigvals(j);
          t_real const diff_lambda(lambda_j - lambda_i);
          if(std::abs(diff_lambda) < 1e-12) continue;


          t_real const lambda_invdiff(1e0/diff_lambda); 
          auto Dj = _C.getD(j);
          t_real factor(lambda_invdiff); 

          t_real sign(1);
          for(t_int r(0); r < _m; ++r, sign = -sign) {
            result += (Di * _C(j, _m-1, r) + sign * Dj * _C(_i, _m-1, r)) * factor;
            factor *= lambda_invdiff * (r+1);
          }
        } // loop over j
  
        return result;
      }

    //! Recursion formula for l != 0 and l != m
    template<class T, class T_ZERO> 
      typename T::t_element general(T & _C, t_int _i, t_int _m, t_int _l, T_ZERO const &_zero) {

        typename T::t_element result(_C.getD(_i) * _C(_i, _m-1, _l-1) / t_real(_l));
        for(t_int j(0); j < _C.nbeigvals(); ++j) {
          if(_i == j) continue;
 
          t_real const diff_lambda(_C.get_eigvals(_i)-_C.get_eigvals(j));
          if(std::abs(diff_lambda) < 1e-12) continue;
 
          t_real const diff_lambda_inv(1e0/diff_lambda); 
          t_real factor(diff_lambda_inv); 
          typename T::t_element intermediate = _C(_i, _m-1, _l) * factor;
          for(t_int r(_l+1); r < _m; ++r) {
            factor *= diff_lambda_inv * r;
            intermediate += _C(_i, _m-1, r) * factor;
          }
          result -= _C.getD(j) * intermediate;
        } // loop over j
  
        return result;
      }
  }
}
#endif 
