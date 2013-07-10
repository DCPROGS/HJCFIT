#ifndef DCPROGS_LIKELIHOOD_ASYMPTOTES_H
#define DCPROGS_LIKELIHOOD_ASYMPTOTES_H

#include <DCProgsConfig.h>
#include <ostream>
#include "root_finder.h"

namespace DCProgs {


  //! A functor to compute asymptotic missed event G.
  //! \detail From the knowledged of roots and eigenvectors, figures out how to compute the
  //!         asymptotic values of missed events R.
  class MSWINDOBE Asymptotes {
    
    public:
      //! Holds a pair defining each exponential function.
      typedef std::pair<t_rmatrix, t_real> t_MatrixAndRoot;
      //| Holds all data relating to this functor.
      typedef std::vector<t_MatrixAndRoot> t_MatricesAndRoots;

      //! Constructor. 
      Asymptotes(t_MatricesAndRoots const &_values ) : matrices_and_roots_(_values) {} 
      //! Creates functor from equation and roots.
      Asymptotes(DeterminantEq const &_equation, std::vector<Root> const &_roots);

      //! Computes \f$\sum_i R_i e^{-\frac{t}{s_i}}\f$.
      //! \details The \f$s_i\f$ are the roots. \f$R_i\f$ matrices are weighted
      //! projections of the eigenvectors corresponding to the roots:
      //! \f$R_i = \frac{c_i\cross r_i}{r_i \cdot W'(s_i) \cdot c_i}\f$.
      t_rmatrix operator()(t_real _t) const;

      //! Access to matrices and roots
      //! The matrices are \f$^AR_i = \frac{c_i\cdot r_i}{r_i \cdot W'(s_i) \cdot c_i}\f$, where
      //! \f$c_i\f$ and \f$r_i\f$ are the left and right eigenvectors of \f$H(s_i)\f$ with
      //! eigenvalue $s_i$. \f$W'(s) = \left.\frac{d W(s)}{d s}\right|_{s=s_i}\f$.
      t_MatrixAndRoot const & operator[](t_int _i) const {
        if(_i < 0) _i += matrices_and_roots_.size();
        if(_i < 0 or _i >= matrices_and_roots_.size())
          throw errors::Index("Index to matrices and roots out-of-range.");
        return matrices_and_roots_[_i]; 
      }
      //! Number of matrices and roots.
      t_int size() { return matrices_and_roots_.size(); }
    protected:
      //! Holds the weight and the exponent of the exponential functions.
      t_MatricesAndRoots matrices_and_roots_;
  };

  //! \brief Matrix with which to compute \f$H_{FA}\f$ for  the CHS vectors.
  //! \cite hawkes:1996.
  t_rmatrix tcrit_H_FA_matrix( Asymptotes const &_asymptotes,
                               QMatrix const &_qmatrix,
                               t_real _eta, t_real _tcrit);
  
}
#endif 
