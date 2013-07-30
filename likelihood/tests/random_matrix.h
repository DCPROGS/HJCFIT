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
  t_rmatrix rate_matrix(t_uint nmin = 5, t_uint nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a non-nan rate_matrix.
  //! Checks the determinant is not NaN.
  t_rmatrix nonnan_rate_matrix(t_uint nmin = 5, t_uint nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a fake_Qmatrix that is not singular.
  t_rmatrix nonsingular_rate_matrix(t_uint nmin = 5, t_uint nmax = 20, t_real large=0.5,
                                    t_real zeroprob=0.3);

  //! Computes a fake state matrix. 
  //! Tries to give it some structure, with a fair number of zeros, as well as large and small
  //! numbers.
  QMatrix qmatrix(t_uint nmin = 5, t_uint nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a non-nan state matrix.
  //! Checks the determinant is not NaN.
  QMatrix nonnan_qmatrix(t_uint nmin = 5, t_uint nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Computes a fake_Qmatrix that is not singular.
  QMatrix nonsingular_qmatrix(t_uint nmin = 5, t_uint nmax = 20, t_real large=0.5, t_real zeroprob=0.3);
  //! Vector of random numbers 
  t_rvector random_vector(t_uint nmin=5, t_uint nmax=20);
}

#endif
