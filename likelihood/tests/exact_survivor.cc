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
#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "../exact_survivor.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class ExactSurvivorTest : public ::testing::Test {
  
  public:
  virtual void SetUp() {
    Q.resize(5, 5);
    Q <<  -3050,        50,  3000,      0,    0,
          2./3., -1502./3.,     0,    500,    0, 
             15,         0, -2065,     50, 2000, 
              0,     15000,  4000, -19000,    0, 
              0,         0,    10,      0,  -10;
    
    nopen = 2; nclose = 3;

    Eigen::EigenSolver<t_rmatrix> eigsolver(Q);
    if(eigsolver.info() != Eigen::Success) 
        throw errors::Mass("Could not solve eigenvalue problem.");
    eigenvalues = -eigsolver.eigenvalues().real();
    eigenvectors = eigsolver.eigenvectors().real();
    eigenvectors_inv = eigenvectors.inverse();
  }

  t_rmatrix N0(t_real _t, t_real _tau, bool _a = true)  {
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int i(0); i < eigenvalues.size(); ++i) 
      result += get_ci00(i, _a) * std::exp(-eigenvalues(i)*_t);
    return result;
  }
  t_rmatrix N1(t_real _t, t_real _tau, bool _a = true)  {
    t_real const t(_t - _tau);
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int i(0); i < eigenvalues.size(); ++i) 
      result += (get_ci10(i, _tau, _a)  + get_ci11(i, _tau, _a) * _t) 
                * std::exp(-eigenvalues(i)*_t);
    return result;
  }

  t_rmatrix get_di(t_int _i, t_real _tau, bool _a = true) {
    t_rvector const left = _a ? eigenvectors.col(_i).head(nopen): eigenvectors.col(_i).tail(nclose);
    t_rvector const right = _a ? eigenvectors_inv.row(_i).tail(nclose):
                                 eigenvectors_inv.row(_i).head(nopen);
    t_rmatrix const exponent = _a ? Q.bottomRightCorner(nclose, nclose):
                                    Q.topLeftCorner(nopen, nopen);
    t_rmatrix const endvec = _a ? Q.bottomLeftCorner(nclose, nopen):  
                                  Q.topRightCorner(nopen, nclose);
    return (left * right.transpose()) * (_tau * exponent).exp() * endvec;
  }
  t_rmatrix get_ci00(t_int _i, bool _a = true) {
    t_rvector const left = _a ? eigenvectors.col(_i).head(nopen):
                                eigenvectors.col(_i).tail(nclose);
    t_rvector const right = _a ? eigenvectors_inv.row(_i).head(nopen):
                                 eigenvectors_inv.row(_i).tail(nclose);
    return left * right.transpose();
  }
  t_rmatrix get_ci10(t_int _i, t_real _tau, bool _a = true) {
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int j(0); j < eigenvalues.size(); ++j) {
      t_real const delta(eigenvalues(j) - eigenvalues(_i));
      if(std::abs(delta) > 1e-8) {
        result += (get_di(_i, _tau, _a) * get_ci00(j, _a) + get_di(j, _tau, _a) * get_ci00(_i, _a)) 
                  / delta;
      } 
    }
    return result;
  }
  t_rmatrix get_ci11(t_int _i, t_real _tau, bool _a = true) {
    return get_di(_i, _tau, _a) * get_ci00(_i, _a);
  }
  protected:
    t_int nopen, nclose;
    t_rmatrix Q;
    t_rvector eigenvalues;
    t_rmatrix eigenvectors;
    t_rmatrix eigenvectors_inv;
};

// Compares recursive implementation to the expanded one in this file.
TEST_F(ExactSurvivorTest, negative_times) {
  QMatrix qmatrix(Q, 2);
  ExactSurvivor survivor(qmatrix, 1e-4);

  EXPECT_TRUE( (survivor.af(-1e-5).array().abs() < 1e-8).all()  );
  EXPECT_EQ(survivor.af(-1e-5).rows(), 2);
  EXPECT_EQ(survivor.af(-1e-5).cols(), 2);
  EXPECT_TRUE( (survivor.fa(-1e-5).array().abs() < 1e-8).all()  );
  EXPECT_EQ(survivor.fa(-1e-5).rows(), 3);
  EXPECT_EQ(survivor.fa(-1e-5).cols(), 3);
}

// Compares recursive implementation to the expanded one in this file.
TEST_F(ExactSurvivorTest, first_interval) {
  std::cout.precision(15);
  QMatrix qmatrix(Q, 2);
  ExactSurvivor survivor(qmatrix, 1e-4);

  auto aR0 = [&](t_real _t, t_real _tau) { return N0(_t, _tau, true); };
  auto fR0 = [&](t_real _t, t_real _tau) { return N0(_t, _tau, false); };
  auto compare = [](t_rmatrix const &_a, t_rmatrix const & _b) {
    return ((_a - _b).array().abs() < 1e-10).all();
  };
  EXPECT_TRUE( compare(survivor.af(1e-5), aR0(1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-5), fR0(1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(3e-5), aR0(3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(3e-5), fR0(3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(5e-5), aR0(5e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(5e-5), fR0(5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.af(1e-5), aR0(5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.fa(1e-5), fR0(5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.af(1e-4 + 1e-5), aR0(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.fa(1e-4 + 1e-5), fR0(1e-4 + 5e-5, 1e-4)) );
}

TEST_F(ExactSurvivorTest, second_interval) {
  std::cout.precision(15);
  QMatrix qmatrix(Q, 2);
  ExactSurvivor survivor(qmatrix, 1e-4);

  auto aG0 = [&](t_real _t, t_real _tau) { return N0(_t, _tau, true); };
  auto aG1 = [&](t_real _t, t_real _tau) -> t_rmatrix {
    return N0(_t, _tau, true) - N1(_t - _tau, _tau, true); 
  };
  auto fG1 = [&](t_real _t, t_real _tau) -> t_rmatrix {
    return N0(_t, _tau, false) - N1(_t - _tau, _tau, false);
  };
  auto compare = [](t_rmatrix const &_a, t_rmatrix const & _b) {
    return ((_a - _b).array().abs() < 1e-10).all();
  };
  EXPECT_TRUE( compare(survivor.af(1e-4 + 1e-5), aG1(1e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-4 + 1e-5), fG1(1e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(1e-4 + 3e-5), aG1(1e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-4 + 3e-5), fG1(1e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(1e-4 + 5e-5), aG1(1e-4 + 5e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-4 + 5e-5), fG1(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.af(1e-4 + 1e-5), aG1(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.fa(1e-4 + 1e-5), fG1(1e-4 + 5e-5, 1e-4)) );
}

// Checks that Di matrices throw out-of-range
TEST_F(ExactSurvivorTest, out_of_range_Di) {

  QMatrix qmatrix(Q, 2);
  ExactSurvivor survivor(qmatrix, 1e-4);
  EXPECT_THROW(survivor.D_af(-1), errors::Index);
  EXPECT_THROW(survivor.D_af(5), errors::Index);
  EXPECT_THROW(survivor.D_fa(-1), errors::Index);
  EXPECT_THROW(survivor.D_fa(5), errors::Index);
}

TEST_F(ExactSurvivorTest, out_of_range_Ciml_matrices) {
  QMatrix qmatrix(Q, 2);
  ExactSurvivor survivor(qmatrix, 1e-4);
  
  EXPECT_THROW(survivor.recursion_af(-1, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_af(5, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(-1, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(5, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, 0, -1), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, 0, -1), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, -1, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, -1, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, 0, 1), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, 0, 1), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, 1, 2), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, 1, 2), errors::Index);
}

// Check that negative values are zero
TEST_F(ExactSurvivorTest, negative_t_is_zero) {
  QMatrix qmatrix(Q, 2);
  ExactSurvivor survivor(qmatrix, 1e-4);
  EXPECT_FALSE((survivor.fa(0).array().abs() < 1e-8).all());
  EXPECT_FALSE((survivor.af(0).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.af(-1e-8).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.fa(-1e-8).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.af(-1e-4).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.fa(-1e-4).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.af(-1e-2).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.fa(-1e-2).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.af(-1e-1).array().abs() < 1e-8).all());
  EXPECT_TRUE((survivor.fa(-1e-1).array().abs() < 1e-8).all());
}

// Check continuity at first interval.
TEST_F(ExactSurvivorTest, continuity_at_t_eq_tau) {

  QMatrix qmatrix(Q, 2);
  t_real const tau = 1e-4;
  ExactSurvivor survivor(qmatrix, tau);
  { 
    t_rmatrix const left = survivor.af(tau - 1e-8);
    t_rmatrix const right = survivor.af(tau + 1e-8);
    t_rmatrix const center = survivor.af(tau);
    EXPECT_FALSE((left.array().abs() < 1e-8).all());
    EXPECT_FALSE((right.array().abs() < 1e-8).all());
    EXPECT_TRUE(( ((left - right).array() / center.array()).abs() < 1e-3).all());
  }
  { 
    t_rmatrix const left = survivor.fa(tau - 1e-8);
    t_rmatrix const right = survivor.fa(tau + 1e-8);
    t_rmatrix const center = survivor.fa(tau);
    EXPECT_FALSE((left.array().abs() < 1e-8).all());
    EXPECT_FALSE((right.array().abs() < 1e-8).all());
    EXPECT_TRUE(( ((left - right).array() / center.array()).abs() < 1e-3).all());
  }
}

// Check continuity at second interval.
TEST_F(ExactSurvivorTest, continuity_at_t_eq_two_tau) {

  QMatrix qmatrix(Q, 2);
  t_real const tau = 1e-4;
  ExactSurvivor survivor(qmatrix, tau);
  { 
    t_rmatrix const left = survivor.af(2e0*tau - 1e-8);
    t_rmatrix const right = survivor.af(2e0*tau + 1e-8);
    t_rmatrix const center = survivor.af(2e0*tau);
    EXPECT_FALSE((left.array().abs() < 1e-8).all());
    EXPECT_FALSE((right.array().abs() < 1e-8).all());
    EXPECT_TRUE(( ((left - right).array() / center.array()).abs() < 1e-3).all());
  }
  { 
    t_rmatrix const left = survivor.fa(2e0*tau - 1e-8);
    t_rmatrix const right = survivor.fa(2e0*tau + 1e-8);
    t_rmatrix const center = survivor.fa(2e0*tau);
    EXPECT_FALSE((left.array().abs() < 1e-8).all());
    EXPECT_FALSE((right.array().abs() < 1e-8).all());
    EXPECT_TRUE(( ((left - right).array() / center.array()).abs() < 1e-3).all());
  }
}

// Check continuity at third interval.
TEST_F(ExactSurvivorTest, continuity_at_t_eq_three_tau) {

  QMatrix qmatrix(Q, 2);
  t_real const tau = 1e-4;
  ExactSurvivor survivor(qmatrix, tau);
  { 
    t_rmatrix const left = survivor.af(3e0*tau - 1e-8);
    t_rmatrix const right = survivor.af(3e0*tau + 1e-8);
    t_rmatrix const center = survivor.af(3e0*tau);
    EXPECT_FALSE((left.array().abs() < 1e-8).all());
    EXPECT_FALSE((right.array().abs() < 1e-8).all());
    EXPECT_TRUE(( ((left - right).array() / center.array()).abs() < 1e-3).all());
  }
  { 
    t_rmatrix const left = survivor.fa(3e0*tau - 1e-8);
    t_rmatrix const right = survivor.fa(3e0*tau + 1e-8);
    t_rmatrix const center = survivor.fa(3e0*tau);
    EXPECT_FALSE((left.array().abs() < 1e-8).all());
    EXPECT_FALSE((right.array().abs() < 1e-8).all());
    EXPECT_TRUE(( ((left - right).array() / center.array()).abs() < 1e-3).all());
  }
}

// Checks that likelihood for t=0 is identity (states are where they start at).
TEST_F(ExactSurvivorTest, at_t_equal_0_is_identity) {

  QMatrix qmatrix(Q, 2);
  t_real const tau = 1e-4;
  ExactSurvivor survivor(qmatrix, tau);
  EXPECT_TRUE(((survivor.af(0) - t_rmatrix::Identity(2, 2)).array().abs() < 1e-12).all());
  EXPECT_TRUE(((survivor.fa(0) - t_rmatrix::Identity(3, 3)).array().abs() < 1e-12).all());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

