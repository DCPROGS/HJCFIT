#ifndef DCPROGS_TEST_RANDOM_MATRIX_H
#define DCPROGS_TEST_RANDOM_MATRIX_H

#include "DCProgsConfig.h"
#include "../state_matrix.h"

#include <random>

namespace DCProgs {

  std::mt19937 & global_mersenne();
  //! Computes a fake qmatrix. 
  //! Tries to give it some structure, with a fair number of zeros, as well as large and small
  //! numbers.
  t_rmatrix qmatrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a non-nan qmatrix.
  //! Checks the determinant is not NaN.
  t_rmatrix nonnan_qmatrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a fake_Qmatrix that is not singular.
  t_rmatrix nonsingular_qmatrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);

  //! Computes a fake state matrix. 
  //! Tries to give it some structure, with a fair number of zeros, as well as large and small
  //! numbers.
  StateMatrix state_matrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a non-nan state matrix.
  //! Checks the determinant is not NaN.
  StateMatrix nonnan_state_matrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a fake_Qmatrix that is not singular.
  StateMatrix nonsingular_state_matrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
}

#endif
