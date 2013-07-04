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

    protected:
      //! Holds the weight and the exponent of the exponential functions.
      t_MatricesAndRoots matrices_and_roots_;
  };

  
}
#endif 
