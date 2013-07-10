#include <random>
#include <memory>
#include <iostream>
#include <gtest/gtest.h>
#include "../idealG.h"
using namespace DCProgs;

std::mt19937& rng() {
  std::mt19937 static mersenne;
  return mersenne; 
}
// Max exponential time is 5*t (t=1)
size_t const nexponents = 5;
// Random matrix have zero if n in [zeros[0], zeros[1][ is > zeros[2]
t_int const zeros[3] = {0, 5, 2}; 
// Min max size of random matrices.
size_t const matsizes[2] = {2, 10};
// Min max numbers in random matrices.
t_real const randreal[2] = {0e0, 1e4};
// Number of random matrices to tests
size_t const nrands = 10;
// Number of random scales to tests per matrix
size_t const nscales = 10;



// Try and make sure that blocks are either zero, or the exponential-like form they have.
class Laplacian
   : public ::testing::TestWithParam<QMatrix> { 
   protected:
     IdealG idealg;
};

TEST_P(Laplacian, af) {
  idealg.set(GetParam());
  typedef std::uniform_real_distribution<t_real> t_rdist;
  t_rdist __rnumbers(randreal[0], randreal[1]);
  auto rnumbers = [&] { return __rnumbers(rng()); };

  t_rmatrix const aa = GetParam().aa();
  t_rmatrix const af = GetParam().af();
  t_int const nrows = aa.rows();
  t_rmatrix const id = t_rmatrix::Identity(nrows, nrows); 
  for(size_t i(0); i < nscales; ++i) {
    t_real const s(rnumbers()); 
    t_rmatrix const factor = s * id - aa;
    t_rmatrix laplace;
    try { laplace = idealg.laplace_af(s); }
    catch(errors::NotInvertible) { 
      Eigen::FullPivLU<t_rmatrix> pivotLU(factor);
      EXPECT_FALSE(pivotLU.isInvertible());
      continue;
    }
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = ((s*id - aa) * laplace - af).array().abs();
    EXPECT_TRUE((diff < 1e-8).all()); 
  }
}

TEST_P(Laplacian, fa) {
  idealg.set(GetParam());
  typedef std::uniform_real_distribution<t_real> t_rdist;
  t_rdist __rnumbers(randreal[0], randreal[1]);
  auto rnumbers = [&] { return __rnumbers(rng()); };

  t_rmatrix const ff = GetParam().ff();
  t_rmatrix const fa = GetParam().fa();
  t_int const nrows = ff.rows();
  t_rmatrix const id = t_rmatrix::Identity(nrows, nrows); 
  for(size_t i(0); i < nscales; ++i) {
    t_real const s(rnumbers()); 
    t_rmatrix const factor = s * id - ff;
    t_rmatrix laplace;
    try { laplace = idealg.laplace_fa(s); }
    catch(errors::NotInvertible) { 
      Eigen::FullPivLU<t_rmatrix> pivotLU(factor);
      EXPECT_FALSE(pivotLU.isInvertible());
      continue;
    }
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = ((s*id - ff) * laplace - fa).array().abs();
    EXPECT_TRUE((diff < 1e-8).all()); 
  }
}

void add_data(std::vector<QMatrix> &_container, t_rmatrix const &_matrix) {
  for(int i(1); i < _matrix.rows()-1; ++i)
    _container.push_back(QMatrix(_matrix, i));
}
t_rmatrix random_matrix() {
  typedef std::uniform_real_distribution<t_real> t_rdist;
  typedef std::uniform_int_distribution<t_int> t_idist;
  t_rdist __rnumbers(randreal[0], randreal[1]);
  t_idist __matsize(matsizes[0], matsizes[1]);
  t_idist __isnotzero(zeros[0], zeros[1]);

  auto matsize = [&] { return __matsize(rng()); };
  auto isnotzero = [&] { return __isnotzero(rng()) < zeros[2]; };
  auto rnumbers = [&] { return isnotzero() ? __rnumbers(rng()): 0; };

  t_int N = matsize();
  t_rmatrix Q(N, N);
  for(t_int i(0); i < N; ++i) {
    for(t_int j(0); j < N; ++j) 
      Q(i, j) = rnumbers();
    Q(i, i) = 0e0;
    Q(i, i) = -Q.row(i).sum();
  }
  return Q;
}

std::shared_ptr<std::vector<QMatrix>> create_container() {

  std::shared_ptr<std::vector<QMatrix>> result(new std::vector<QMatrix>);

  t_rmatrix Q(5, 5);
  Q <<  -3050,        50,  3000,      0,    0,
        2./3., -1502./3.,     0,    500,    0, 
           15,         0, -2065,     50, 2000, 
            0,     15000,  4000, -19000,    0, 
            0,         0,    10,      0,  -10;
  add_data(*result, Q);

  for(size_t i(0); i < nrands; ++i) add_data(*result, random_matrix());
  return result;
}

std::shared_ptr<std::vector<QMatrix>> testcases = create_container();
INSTANTIATE_TEST_CASE_P(IdealG, Laplacian, ::testing::ValuesIn(*testcases));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

