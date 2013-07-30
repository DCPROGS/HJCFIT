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
#include <random>
#include <memory>
#include <iostream>
#include <gtest/gtest.h>
#include "../idealG.h"
using namespace DCProgs;

// Max exponential time is 5*t (t=1)
size_t const nexponents = 5;
// Random matrix have zero if n in [zeros[0], zeros[1][ is > zeros[2]
t_uint const zeros[3] = {0, 5, 2}; 
// Min max size of random matrices.
size_t const matsizes[2] = {2, 10};
// Min max numbers in random matrices.
t_real const randreal[2] = {-1e0, 1e0};
// Number of random matrices to tests
size_t const nrands = 100;


// Try and make sure that blocks are either zero, or the exponential-like form they have.
class Exponentiation
   : public ::testing::TestWithParam<QMatrix> { 
   protected:
     IdealG idealg;
};

TEST_P(Exponentiation, af){
  idealg.set(GetParam());
  t_rmatrix exponential = GetParam().aa().exp();
  t_rmatrix current = GetParam().af();
  for(size_t i(0); i < nexponents; ++i, current = exponential * current) {
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = (idealg.af(t_real(i)) - current).array().abs();
    EXPECT_TRUE((diff < 1e-8).all()); 
  }
}
TEST_P(Exponentiation, fa){
  idealg.set(GetParam());
  t_rmatrix exponential = GetParam().ff().exp();
  t_rmatrix current = GetParam().fa();
  for(size_t i(0); i < nexponents; ++i, current = exponential * current) {
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = (idealg.fa(t_real(i)) - current).array().abs();
    EXPECT_TRUE((diff < 1e-8).all()); 
  }
}

void add_data(std::vector<QMatrix> &_container, t_rmatrix const &_matrix) {
  for(int i(1); i < _matrix.rows()-1; ++i)
  { _container.emplace_back(_matrix, i); break; }
}
template<class T> t_rmatrix random_matrix(T && _rng) {
  typedef std::uniform_real_distribution<t_real> t_rdist;
  typedef std::uniform_int_distribution<t_uint> t_idist;
  t_rdist __rnumbers(randreal[0], randreal[1]);
  t_idist __matsize(matsizes[0], matsizes[1]);
  t_idist __isnotzero(zeros[0], zeros[1]);

  auto matsize = [&] { return __matsize(_rng); };
  auto isnotzero = [&] { return __isnotzero(_rng) < zeros[2]; };
  auto rnumbers = [&] { return isnotzero() ? __rnumbers(_rng): 0; };

  t_uint N = matsize();
  t_rmatrix Q(N, N);
  for(t_uint i(0); i < N; ++i) {
    for(t_uint j(0); j < N; ++j) 
      Q(i, j) = rnumbers();
    Q(i, i) = 0e0;
    Q(i, i) = -Q.row(i).sum();
  }
  return Q;
}

std::shared_ptr<std::vector<QMatrix>> create_container() {

  std::shared_ptr<std::vector<QMatrix>> result(new std::vector<QMatrix>());

  t_rmatrix Q(5, 5);
  Q <<  -3050,        50,  3000,      0,    0,
        2./3., -1502./3.,     0,    500,    0, 
           15,         0, -2065,     50, 2000, 
            0,     15000,  4000, -19000,    0, 
            0,         0,    10,      0,  -10;
  add_data(*result, Q);

  std::mt19937 mersenne;
  for(size_t i(0); i < nrands; ++i) add_data(*result, random_matrix(mersenne));
  return result;
}

std::shared_ptr<std::vector<QMatrix>> testcases = create_container();
INSTANTIATE_TEST_CASE_P(IdealG, Exponentiation, ::testing::ValuesIn(*testcases));


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

