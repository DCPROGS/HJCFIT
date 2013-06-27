#include "DCProgsConfig.h"

#include <math.h>

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
  t_rmatrix qmatrix(int nmin, int nmax, t_real large, t_real zeroprob) {
    
    std::mt19937& mersenne = global_mersenne();
    typedef std::uniform_int_distribution<t_int> t_idist;
    typedef std::uniform_real_distribution<t_real> t_rdist;
    t_rdist rdist(0, 1);
    t_int const N = t_idist(nmin, nmax)(mersenne);
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int i(0); i < N; ++i) 
      for(t_int j(i+1); j < N; ++j) {
        if(rdist(mersenne) > zeroprob)  continue;
        result(i, j) = rdist(mersenne);
        result(j, i) = rdist(mersenne); 
        if(rdist(mersenne) > large) {
          if(rdist(mersenne) > 0.5) result(i, j) *= 3e3;
          else result(j, i) *=  3e3;
        } 
    }
    for(t_int i(0); i < N; ++i) result(i, i) -= result.row(i).sum();
    return result;
  }

  t_rmatrix nonnan_qmatrix(int nmin, int nmax, t_real large, t_real zeroprob) {
    t_rmatrix result; 
    do { result = qmatrix(nmin, nmax, large, zeroprob); }
    while(std::isnan(result.determinant()));
    return result;
  }
  t_rmatrix nonsingular_qmatrix(int nmin, int nmax, t_real large, t_real zeroprob) {
    t_rmatrix result;
    do { result = nonnan_qmatrix(nmin, nmax, large, zeroprob); }
    while(std::abs(result.determinant()) < 1e-4);
    return result;
  }

  StateMatrix state_matrix(int nmin, int nmax, t_real large, t_real zeroprob) {
     typedef std::uniform_int_distribution<t_int> t_idist;
     t_rmatrix matrix = qmatrix(nmin, nmax, large, zeroprob);
     return StateMatrix(matrix, t_idist(2, matrix.rows()-2)(global_mersenne()));
  }
  StateMatrix nonnan_state_matrix(int nmin, int nmax, t_real large, t_real zeroprob) {
     typedef std::uniform_int_distribution<t_int> t_idist;
     t_rmatrix matrix = nonnan_qmatrix(nmin, nmax, large, zeroprob);
     return StateMatrix(matrix, t_idist(2, matrix.rows()-2)(global_mersenne()));
  }
  StateMatrix nonsingular_state_matrix(int nmin, int nmax, t_real large, t_real zeroprob) {
     typedef std::uniform_int_distribution<t_int> t_idist;
     t_rmatrix matrix = nonsingular_qmatrix(nmin, nmax, large, zeroprob);
     return StateMatrix(matrix, t_idist(2, matrix.rows()-2)(global_mersenne()));
  }
}
