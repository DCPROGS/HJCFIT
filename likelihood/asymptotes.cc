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

#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>

#include "asymptotes.h"

namespace DCProgs {

  t_srmatrix Asymptotes :: operator()(t_real _t) const {

    t_MatricesAndRoots :: const_iterator i_first = matrices_and_roots_.begin();
    t_MatricesAndRoots :: const_iterator const i_end = matrices_and_roots_.end();

    auto function = [_t](t_MatrixAndRoot const &_mat_and_root) -> t_srmatrix {
      return std::get<0>(_mat_and_root) * std::exp(_t * std::get<1>(_mat_and_root));
    };

    t_srmatrix result = function(*i_first);
    for(++i_first; i_first != i_end; ++i_first) result += function(*i_first);
    return result;
  }


  Asymptotes :: Asymptotes   (DeterminantEq const &_equation, std::vector<Root> const &_roots)
                           : matrices_and_roots_() {

    verify_qmatrix(_equation.get_qmatrix());
    matrices_and_roots_.reserve(_roots.size());
    for(Root const & root: _roots) {

      t_srmatrix const H(_equation.H(root.root));
      t_srmatrix const derivative(_equation.s_derivative(root.root));
      Eigen::JacobiSVD<t_srmatrix> svd(H - root.root * t_srmatrix::Identity(H.rows(), H.cols()),
                                      Eigen::ComputeThinU|Eigen::ComputeThinV);

      t_rvector const abs_singval( svd.singularValues().array().abs() );

      if(static_cast<t_uint>(abs_singval.size()) < root.multiplicity)
        throw errors::Mass("Requesting more roots than there are singular values.");

      // Figures out lowest singular values
      // To do this, creates a vector of indices that is partially sorted.
      t_int i(0);
      std::vector<t_int> indices( abs_singval.size() );
      std::generate(indices.begin(), indices.end(), [&i]() { return i++; });
      auto comparison = [abs_singval](t_int a, t_int b) {
        return abs_singval(a) < abs_singval(b);
      };
      std::partial_sort(indices.begin(), indices.begin() + root.multiplicity,
                        indices.end(), comparison);

      // Following is 2.29 from Colquhounm Hawkes, Srodzinski (1996)
      auto single_root_function = [&svd, &derivative, &H, &root](t_int _index) -> t_rmatrix { 
         auto c_i = svd.matrixV().col(_index) / svd.matrixV().col(_index).sum() ;
         auto r_i = svd.matrixU().col(_index).transpose() / svd.matrixU().col(_index).sum();
         return c_i * r_i / (r_i * derivative * c_i);
      };
      // Now loop over all degenerate roots.
      std::vector<t_int> :: const_iterator i_first = indices.begin();
      std::vector<t_int> :: const_iterator const i_end = i_first + root.multiplicity;
      t_rmatrix Ri = single_root_function(*i_first);
      for(++i_first; i_first != i_end; ++i_first) Ri += single_root_function(*i_first);

      // Finally, add to vector of matrices and roots
      matrices_and_roots_.emplace_back(std::move(Ri), root.root);
    }
  }

  // Matrix with which to compute \f$H_{FA}\f$ for  the CHS vectors.
  t_srmatrix MSWINDOBE partial_CHS_matrix( Asymptotes const &_asymptotes,
                                          t_real _tau, t_real _tcrit ) {

    auto function = [&_asymptotes, &_tcrit, &_tau](t_int i) -> t_rmatrix {
      Asymptotes::t_MatrixAndRoot const & mat_and_root = _asymptotes[i];
      t_rmatrix const & Ri = std::get<0>(mat_and_root);
      t_real const root = std::get<1>(mat_and_root);
      return -Ri * std::exp(root * (_tcrit - _tau)) / root;
    };
    t_srmatrix result = function(0);
    for(t_int i(1); i < _asymptotes.size(); ++i) result += function(i);
    return result;
  }

}
