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

#ifndef DCPROGS_LIKELIHOOD_ASYMPTOTES_H
#define DCPROGS_LIKELIHOOD_ASYMPTOTES_H

#include <DCProgsConfig.h>
#include <ostream>
#include "root_finder.h"

namespace DCProgs {


  //! \brief A functor to compute asymptotic missed event G.
  //! \details From the knowledged of roots and eigenvectors, figures out how to compute the
  //!          asymptotic values of missed events R.
  class MSWINDOBE Asymptotes {
    
    public:
      //! \brief Holds a pair defining each exponential function.
      //! \details - The first item is the weigh (as an matrix) of the exponential.
      //! - The second item is the exponent of the exponential.
      typedef std::pair<t_rmatrix, t_real> t_MatrixAndRoot;

      //| \typedef  std::vector<t_MatrixAndRoot> t_MatricesAndRoots
      //! \brief Container holding the parameters for each exponential
      typedef std::vector<t_MatrixAndRoot> t_MatricesAndRoots;

      //! Constructor. 
      Asymptotes(t_MatricesAndRoots const &_values ) : matrices_and_roots_(_values) {} 
      //! Creates functor from equation and roots.
      Asymptotes(DeterminantEq const &_equation, std::vector<Root> const &_roots);

      //! Computes \f$\sum_i R_i e^{-\frac{t}{s_i}}\f$.
      //! \details The \f$s_i\f$ are the roots. \f$R_i\f$ matrices are weighted
      //! projections of the eigenvectors corresponding to the roots:
      //! \f$R_i = \frac{c_i\times r_i}{r_i \cdot W'(s_i) \cdot c_i}\f$.
      t_rmatrix operator()(t_real _t) const;

      //! Access to matrices and roots
      //! The matrices are \f$^AR_i = \frac{c_i\cdot r_i}{r_i \cdot W'(s_i) \cdot c_i}\f$, where
      //! \f$c_i\f$ and \f$r_i\f$ are the left and right eigenvectors of \f$H(s_i)\f$ with
      //! eigenvalue $s_i$. \f$W'(s) = \left.\frac{d W(s)}{d s}\right|_{s=s_i}\f$.
      t_MatrixAndRoot const & operator[](t_int _i) const {
        if(_i < 0) _i += static_cast<t_int>(matrices_and_roots_.size());
        if(_i < 0 or _i >= static_cast<t_int>(matrices_and_roots_.size())) {
          std::ostringstream sstr; 
          sstr << "Index to matrices and roots out-of-range: " << _i
               << " vs. " << matrices_and_roots_.size() << ".";
          throw errors::Index(sstr.str());
        }
        return matrices_and_roots_[_i]; 
      }
      //! Number of matrices and roots.
      t_int size() const { return static_cast<t_int>(matrices_and_roots_.size()); }
    protected:
      //! Holds the weight and the exponent of the exponential functions.
      t_MatricesAndRoots matrices_and_roots_;
  };

  //! \brief Partial computation of \f$H_{FA}\f$ for CHS vectors.
  //! \details The object is to implement Eq. 5.10 from \cite colquhoun:1996. In order to do this,
  //!          we abstract the part that comes from the likelihood and compute in this function 
  //!          \f$-\sum_{i=1}^{k_F}{}^{A}R_i\frac{1}{s_i} e^{(t_{\mathrm{crit}}-\tau)s_i}\f$. The
  //!          additional factor \f$Q_{AF}e^{Q_{FF}\tau}\f$ is computed in MissedEventsG.
  //! \note This is somewhat outside the remit of Asymptotes, although the calculations are similar
  //!       for good reasons. In any case, we keep this function outside the class itself, so as to
  //!       not confuse the purpose of the class itself.
  t_rmatrix MSWINDOBE partial_CHS_matrix( Asymptotes const &_asymptotes,
                                          t_real _tau, t_real _tcrit );
}
#endif 
