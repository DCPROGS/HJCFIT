#ifndef DCPROGS_TEST_RANDOM_MATRIX_H
#define DCPROGS_TEST_RANDOM_MATRIX_H

#include "DCProgsConfig.h"
#include "../qmatrix.h"

#include <random>

namespace DCProgs {

  std::mt19937 & global_mersenne();
  //! Computes a fake rate_matrix. 
  //! Tries to give it some structure, with a fair number of zeros, as well as large and small
  //! numbers.
  t_rmatrix rate_matrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a non-nan rate_matrix.
  //! Checks the determinant is not NaN.
  t_rmatrix nonnan_rate_matrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a fake_Qmatrix that is not singular.
  t_rmatrix nonsingular_rate_matrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);

  //! Computes a fake state matrix. 
  //! Tries to give it some structure, with a fair number of zeros, as well as large and small
  //! numbers.
  QMatrix qmatrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a non-nan state matrix.
  //! Checks the determinant is not NaN.
  QMatrix nonnan_qmatrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a fake_Qmatrix that is not singular.
  QMatrix nonsingular_qmatrix(int nmin = 5, int nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
}

#endif
