#include <DCProgsConfig.h>

#include <iostream>
#include "asymptotes.h"

namespace DCProgs {

  t_rmatrix Asymptotes :: operator()(t_real _t) const {

    t_MatricesAndRoots :: const_iterator i_first = matrices_and_roots_.begin();
    t_MatricesAndRoots :: const_iterator const i_end = matrices_and_roots_.end();

    auto function = [_t](t_MatrixAndRoot const &_mat_and_root) -> t_rmatrix {
      return std::get<0>(_mat_and_root) * std::exp(_t * std::get<1>(_mat_and_root));
    };

    t_rmatrix result = function(*i_first);
    for(++i_first; i_first != i_end; ++i_first) result += function(*i_first); 
    return result;
  }


  Asymptotes :: Asymptotes   (DeterminantEq const &_equation, std::vector<Root> const &_roots) 
                           : matrices_and_roots_() {

    matrices_and_roots_.reserve(_roots.size());

    for(Root const & root: _roots) {

      t_rmatrix const H(_equation.H(root.root));
      t_rmatrix const derivative(_equation.s_derivative(root.root)); 
      Eigen::JacobiSVD<t_rmatrix> svd(H - root.root * t_rmatrix::Identity(H.rows(), H.cols()), 
                                      Eigen::ComputeThinU|Eigen::ComputeThinV);

      t_rvector const abs_singval( svd.singularValues().array().abs() );

      if(abs_singval.size() < root.multiplicity) 
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
         auto c_i = svd.matrixV().col(_index);
         auto r_i = svd.matrixU().col(_index).transpose();
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

}
