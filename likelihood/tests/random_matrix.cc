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

#include "DCProgsConfig.h"

#include <math.h>
#include <cmath>
#include <ctime>

#include <iostream>

#include "random_matrix.h"

namespace DCProgs {
  
  std::mt19937 & global_mersenne() {
    std::mt19937 static mersenne;
    bool static is_first_entry = true;
    if(is_first_entry) {
  #   ifdef HAS_CXX11_RANDOM_DEVICE
        std::random_device rd;
        mersenne.seed(rd()); 
  #   else 
        mersenne.seed(static_cast<unsigned int>(std::time(nullptr))); 
  #   endif
      is_first_entry = false;
    }
    return mersenne;
  }


  //! Creates a fake Qmatrix
  t_rmatrix rate_matrix(t_uint nmin, t_uint nmax, t_real large, t_real zeroprob) {
    
    std::mt19937& mersenne = global_mersenne();
    typedef std::uniform_int_distribution<t_uint> t_idist;
    typedef std::uniform_real_distribution<t_real> t_rdist;
    t_rdist rdist(0, 1);
    t_uint const N = t_idist(nmin, nmax)(mersenne);
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_uint i(0); i < N; ++i) 
      for(t_uint j(i+1); j < N; ++j) {
        if(rdist(mersenne) > zeroprob)  continue;
        result(i, j) = rdist(mersenne);
        result(j, i) = rdist(mersenne); 
        if(rdist(mersenne) > large) {
          if(rdist(mersenne) > 0.5) result(i, j) *= 3e3;
          else result(j, i) *=  3e3;
        } 
    }
    for(t_uint i(0); i < N; ++i) result(i, i) -= result.row(i).sum();
    return result;
  }

  t_rmatrix nonnan_rate_matrix(t_uint nmin, t_uint nmax, t_real large, t_real zeroprob) {
    t_rmatrix result; 
    do { result = rate_matrix(nmin, nmax, large, zeroprob); }
    while(DCPROGS_ISNAN(result.determinant()));
    return result;
  }
  t_rmatrix nonsingular_rate_matrix(t_uint nmin, t_uint nmax, t_real large, t_real zeroprob) {
    t_rmatrix result;
    do { result = nonnan_rate_matrix(nmin, nmax, large, zeroprob); }
    while(std::abs(result.determinant()) < 1e-4);
    return result;
  }

  QMatrix qmatrix(t_uint nmin, t_uint nmax, t_real large, t_real zeroprob) {
     typedef std::uniform_int_distribution<t_uint> t_idist;
     t_rmatrix matrix = rate_matrix(nmin, nmax, large, zeroprob);
     return QMatrix(matrix, t_idist(2, matrix.rows()-2)(global_mersenne()));
  }
  QMatrix nonnan_qmatrix(t_uint nmin, t_uint nmax, t_real large, t_real zeroprob) {
     typedef std::uniform_int_distribution<t_uint> t_idist;
     t_rmatrix matrix = nonnan_rate_matrix(nmin, nmax, large, zeroprob);
     return QMatrix(matrix, t_idist(2, matrix.rows()-2)(global_mersenne()));
  }
  QMatrix nonsingular_qmatrix(t_uint nmin, t_uint nmax, t_real large, t_real zeroprob) {
     typedef std::uniform_int_distribution<t_uint> t_idist;
     t_rmatrix matrix = nonsingular_rate_matrix(nmin, nmax, large, zeroprob);
     return QMatrix(matrix, t_idist(2, matrix.rows()-2)(global_mersenne()));
  }
  t_rvector random_vector(t_uint nmin, t_uint nmax) {
    std::mt19937& mersenne = global_mersenne();
    typedef std::uniform_int_distribution<t_uint> t_idist;
    typedef std::uniform_real_distribution<t_real> t_rdist;
    t_rdist rdist(0, 1);
    t_uint const N = t_idist(nmin, nmax)(mersenne);
    t_rvector result = t_rvector::Zero(N);
    for(t_uint i(0); i < N; ++i) result(i) = rdist(mersenne);
    return result;
  }
}
