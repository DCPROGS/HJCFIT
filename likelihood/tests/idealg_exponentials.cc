#include <random>
#include <memory>
#include <iostream>
#include <gtest/gtest.h>
#include "../idealG.h"
using namespace DCProgs;

// Max exponential time is 5*t (t=1)
size_t const nexponents = 5;
// Random matrix have zero if n in [zeros[0], zeros[1][ is > zeros[2]
size_t const zeros[3] = {0, 5, 2}; 
// Min max size of random matrices.
size_t const matsizes[2] = {2, 10};
// Min max numbers in random matrices.
t_real const randreal[2] = {-1e0, 1e0};
// Number of random matrices to tests
size_t const nrands = 100;


// Try and make sure that blocks are either zero, or the exponential-like form they have.
class Exponentiation
   : public ::testing::TestWithParam<StateMatrix> { 
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

TEST_P(Exponentiation, aa){
  idealg.set(GetParam());
  EXPECT_EQ(idealg.laplace_aa(0).rows(), idealg.get_nopen());
  EXPECT_EQ(idealg.laplace_aa(0).cols(), idealg.get_nopen());
  EXPECT_EQ(idealg.laplace_aa(1).rows(), idealg.get_nopen());
  EXPECT_EQ(idealg.laplace_aa(1).cols(), idealg.get_nopen());
  EXPECT_TRUE((idealg.aa(0).array().abs() < 1e-8).all());
  EXPECT_TRUE((idealg.aa(1).array().abs() < 1e-8).all());
}
TEST_P(Exponentiation, ff){
  idealg.set(GetParam());
  EXPECT_EQ(idealg.laplace_ff(0).rows(), idealg.get_Q().rows() - idealg.get_nopen());
  EXPECT_EQ(idealg.laplace_ff(0).cols(), idealg.get_Q().rows() - idealg.get_nopen());
  EXPECT_EQ(idealg.laplace_ff(1).rows(), idealg.get_Q().rows() - idealg.get_nopen());
  EXPECT_EQ(idealg.laplace_ff(1).cols(), idealg.get_Q().rows() - idealg.get_nopen());
  EXPECT_TRUE((idealg.ff(0).array().abs() < 1e-8).all());
  EXPECT_TRUE((idealg.ff(1).array().abs() < 1e-8).all());
}


void add_data(std::vector<StateMatrix> &_container, t_rmatrix const &_matrix) {
  for(int i(1); i < _matrix.rows()-1; ++i)
    _container.push_back(StateMatrix(_matrix, i));
}
template<class T> t_rmatrix random_matrix(T && _rng) {
  typedef std::uniform_real_distribution<t_real> t_rdist;
  typedef std::uniform_int_distribution<t_int> t_idist;
  t_rdist __rnumbers{randreal[0], randreal[1]};
  t_idist __matsize{matsizes[0], matsizes[1]};
  t_idist __isnotzero{zeros[0], zeros[1]};

  auto matsize = [&] { return __matsize(_rng); };
  auto isnotzero = [&] { return __isnotzero(_rng) < zeros[2]; };
  auto rnumbers = [&] { return isnotzero() ? __rnumbers(_rng): 0; };

  t_int N = matsize();
  t_rmatrix Q(N, N);
  for(size_t i(0); i < N; ++i) {
    for(size_t j(0); j < N; ++j) 
      Q(i, j) = rnumbers();
    Q(i, i) = 0e0;
    Q(i, i) = -Q.row(i).sum();
  }
  return Q;
}

std::shared_ptr<std::vector<StateMatrix>> create_container() {

  std::shared_ptr<std::vector<StateMatrix>> result(new std::vector<StateMatrix>);

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

std::shared_ptr<std::vector<StateMatrix>> testcases = create_container();
INSTANTIATE_TEST_CASE_P(IdealG, Exponentiation, ::testing::ValuesIn(*testcases));
