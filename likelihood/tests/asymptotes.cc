#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>
#include <unsupported/Eigen/MatrixFunctions>
#include "../asymptotes.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class DeterminantEqTest : public ::testing::Test {
  
  public:
  virtual void SetUp() {
    Q.resize(5, 5);
    Q <<  -3050,        50,  3000,      0,    0,
          2./3., -1502./3.,     0,    500,    0, 
             15,         0, -2065,     50, 2000, 
              0,     15000,  4000, -19000,    0, 
              0,         0,    10,      0,  -10;
  }
  protected:
    t_rmatrix Q;
};


// Missed event with t critical == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, FF_critical_resolution_is_zero) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 0, false); 
  EXPECT_TRUE( ((det.H(0).array() - states.ff().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(1).array() - states.ff().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(10).array() - states.ff().array()).abs() < 1e-8).all() );

  Eigen::EigenSolver<t_rmatrix> eigsolver(states.ff());
  t_rvector const eigs = eigsolver.eigenvalues().real();
  EXPECT_TRUE(std::abs(det(eigs(0)) / eigs(0)) < 1e-5);
  EXPECT_TRUE(std::abs(det(eigs(1)) / eigs(1)) < 1e-5);
  EXPECT_TRUE(std::abs(det(eigs(2)) / eigs(2)) < 1e-5);
}
// Missed event with t critical == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, FF_critical_resolution_is_zero_check_derivative) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 0, false); 
  t_rmatrix id = t_rmatrix::Identity(3, 3);
  EXPECT_TRUE( ((det.s_derivative(0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(1).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(10).array() - id.array()).abs() < 1e-8).all() );
}
// Missed event with t critical == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, AA_critical_resolution_is_zero) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 0, true); 
  EXPECT_TRUE( ((det.H(0).array() - states.aa().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(1).array() - states.aa().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(10).array() - states.aa().array()).abs() < 1e-8).all() );

  Eigen::EigenSolver<t_rmatrix> eigsolver(states.aa());
  t_rvector const eigs = eigsolver.eigenvalues().real();
  EXPECT_TRUE(std::abs(det(eigs(0)) / eigs(0)) < 1e-5);
  EXPECT_TRUE(std::abs(det(eigs(1)) / eigs(1)) < 1e-5);
}
// Missed event with t critical == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, AA_critical_resolution_is_zero_check_derivative) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 0, true); 
  t_rmatrix id = t_rmatrix::Identity(2, 2);
  EXPECT_TRUE( ((det.s_derivative(0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(1).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(10).array() - id.array()).abs() < 1e-8).all() );
}

// Test non-zero critical tau using numerical and enalytical derivatives.
TEST_F(DeterminantEqTest, from_tau_derivative) {
  StateMatrix const states(Q, 2);
  t_real dt = 1e-8;
  DeterminantEq detEpsilon(states, dt, false); 
  DeterminantEq det0(states, 0, false); 
  t_rmatrix approx = (detEpsilon.H(0) - det0.H(0)) / dt;
  t_rmatrix exact = states.fa() * states.af();
  t_real magnitude = approx.array().abs().maxCoeff(); 
  EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
               < 10 * std::sqrt(dt)).all() ) 
             << "Derivative for tau = 0 and s = 0\n";

  approx = (detEpsilon.H(1) - det0.H(1)) / dt;
  magnitude = approx.array().abs().maxCoeff(); 
  EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
               < 10 * std::sqrt(dt)).all() )
             << "Derivative for tau = 0 and s != 0\n";


  det0.set_tau(1e-2);
  detEpsilon.set_tau(det0.get_tau()+dt);

  exact = states.fa() * (det0.get_tau() * states.aa()).exp() * states.af();
  approx = (detEpsilon.H(0) - det0.H(0)) / dt;
  magnitude = approx.array().abs().maxCoeff(); 
  EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
               < 10 * std::sqrt(dt)).all() ) 
             << "Derivative for tau != 0 and s = 0\n";

  exact = states.fa() * std::exp(-det0.get_tau()) * (det0.get_tau() * states.aa()).exp() * states.af();
  approx = (detEpsilon.H(1) - det0.H(1)) / dt;
  magnitude = approx.array().abs().maxCoeff(); 
  EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
               < 10 * std::sqrt(dt)).all() ) 
             << "Derivative for tau != 0 and s != 0\n";
}

// Test s derivatives for non-zero critical tau using numerical and enalytical derivatives.
TEST_F(DeterminantEqTest, s_derivative_from_tau_derivative) {
  StateMatrix const states(Q, 2);
  t_real dt = 1e-4;
  DeterminantEq detEpsilon(states, dt, false); 
  DeterminantEq det0(states, 0, false); 
  t_rmatrix approx = (detEpsilon.s_derivative(0) - det0.s_derivative(0)) / dt;
  std::cout.precision(10);
  std::cout << approx << std::endl << std::endl;
  EXPECT_TRUE( (approx.array().abs() < 10 * std::sqrt(dt)).all() ) 
             << "s Derivative for tau = 0 and s = 0\n";

// approx = (detEpsilon.H(1) - det0.H(1)) / dt;
// magnitude = approx.array().abs().maxCoeff(); 
// EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
//              < 10 * std::sqrt(dt)).all() )
//            << "Derivative for tau = 0 and s != 0\n";
//
//
// det0.set_tau(1e-2);
// detEpsilon.set_tau(det0.get_tau()+dt);
//
// exact = states.fa() * (det0.get_tau() * states.aa()).exp() * states.af();
// approx = (detEpsilon.H(0) - det0.H(0)) / dt;
// magnitude = approx.array().abs().maxCoeff(); 
// EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
//              < 10 * std::sqrt(dt)).all() ) 
//            << "Derivative for tau != 0 and s = 0\n";
//
// exact = states.fa() * std::exp(-det0.get_tau()) * (det0.get_tau() * states.aa()).exp() * states.af();
// approx = (detEpsilon.H(1) - det0.H(1)) / dt;
// magnitude = approx.array().abs().maxCoeff(); 
// EXPECT_TRUE( ((approx - exact).array().abs() / magnitude
//              < 10 * std::sqrt(dt)).all() ) 
//            << "Derivative for tau != 0 and s != 0\n";
}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

