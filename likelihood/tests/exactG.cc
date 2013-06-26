#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "../exactG.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class ExactGTest : public ::testing::Test {
  
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
    for(t_int i(0); i < eigenvalues.size(); ++i) {
      if(std::abs(eigenvalues(i) > 1e-12))  {
        result += get_ci00(i, _a) * std::exp(-eigenvalues(i)*_t);
      }
    }
    return result;
  }
  t_rmatrix N1(t_real _t, t_real _tau, bool _a = true)  {
    t_real const t(_t - _tau);
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int i(0); i < eigenvalues.size(); ++i) {
      if(std::abs(eigenvalues(i) > 1e-12)) 
        result += (get_ci10(i, _tau, _a)  + get_ci11(i, _tau, _a) * _t) * std::exp(-eigenvalues(i)*_t);
    }
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
TEST_F(ExactGTest, first_interval) {
  std::cout.precision(15);
  StateMatrix transitions(Q, 2);
  ExactG exactG(transitions, 1e-4);

  auto aG0 = [&](t_real _t, t_real _tau) {
    return (N0(_t - _tau, _tau, true) * transitions.af() * (_tau * transitions.ff()).exp()).eval();
  };
  auto fG0 = [&](t_real _t, t_real _tau) {
    return (N0(_t - _tau, _tau, false) * transitions.fa() * (_tau * transitions.aa()).exp()).eval();
  };
  auto compare = [](t_rmatrix const &_a, t_rmatrix const & _b) {
    return ((_a - _b).array().abs() < 1e-10).all();
  };
  EXPECT_TRUE( compare(exactG.af(1e-4 + 1e-5), aG0(1e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.fa(1e-4 + 1e-5), fG0(1e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.af(1e-4 + 3e-5), aG0(1e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.fa(1e-4 + 3e-5), fG0(1e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.af(1e-4 + 5e-5), aG0(1e-4 + 5e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.fa(1e-4 + 5e-5), fG0(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(exactG.af(1e-4 + 1e-5), aG0(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(exactG.fa(1e-4 + 1e-5), fG0(1e-4 + 5e-5, 1e-4)) );
}

TEST_F(ExactGTest, second_interval) {
  std::cout.precision(15);
  StateMatrix transitions(Q, 2);
  ExactG exactG(transitions, 1e-4);

  auto aG0 = [&](t_real _t, t_real _tau) {
    return (N0(_t - _tau, _tau, true) * transitions.af() * (_tau * transitions.ff()).exp()).eval();
  };
  auto aG1 = [&](t_real _t, t_real _tau) {
    return ( (N0(_t - _tau, _tau, true) - N1(_t - 2e0*_tau, _tau, true)) 
             * transitions.af() * (_tau * transitions.ff()).exp()).eval();
  };
  auto fG1 = [&](t_real _t, t_real _tau) {
    return ( (N0(_t - _tau, _tau, false) - N1(_t - 2e0*_tau, _tau, false))
             * transitions.fa() * (_tau * transitions.aa()).exp()).eval();
  };
  auto compare = [](t_rmatrix const &_a, t_rmatrix const & _b) {
    return ((_a - _b).array().abs() < 1e-10).all();
  };
  EXPECT_TRUE( compare(exactG.af(2e-4 + 1e-5), aG1(2e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.fa(2e-4 + 1e-5), fG1(2e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.af(2e-4 + 3e-5), aG1(2e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.fa(2e-4 + 3e-5), fG1(2e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.af(2e-4 + 5e-5), aG1(2e-4 + 5e-5, 1e-4)) );
  EXPECT_TRUE( compare(exactG.fa(2e-4 + 5e-5), fG1(2e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(exactG.af(2e-4 + 1e-5), aG1(2e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(exactG.fa(2e-4 + 1e-5), fG1(2e-4 + 5e-5, 1e-4)) );
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
