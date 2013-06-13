#include <iostream>
#include <type_traits>
#include <tuple>
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

//! Convergence criteria for taylor expansion.
template<class T_APPROX, class T_EXACT> 
   void taylor_convergence( T_APPROX &_approx, T_EXACT &_exact, 
                            t_real _s, t_real _tau, t_real _dtau,
                            std::string const &_message = "", t_real M = 10e0 ) {
     t_rmatrix const order2 = 0.5 * _dtau * (_exact(_s, _tau + _dtau) - _exact(_s, _tau));
     t_rmatrix const approx = _approx(_s, _tau + _dtau) - _approx(_s, _tau);
     t_rmatrix const exact  = _dtau * _exact(_s, _tau);
     t_rmatrix const diff = approx - exact;
     t_bmatrix condition =    (diff.array().abs() < M * order2.array().abs())
                           or (diff.array().abs() < 1e-12);
     EXPECT_TRUE( condition.all() ) << _message 
       << "Params:  tau=" << _tau << " -- s=" << _s << " -- dtau=" << _dtau << "\n"
       << "approx: \n" << approx << "\n"
       << "exact:  \n" << exact  << "\n"
       << "diff:  \n" << diff  << "\n"
       << "condition:  \n" << condition  << "\n"
       << "max (M=" << M << "):  \n" << (M * order2)  << "\n";
 }


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
  DeterminantEq determinant(states, 0, false); 
  auto approx = [&](t_real _s, t_real _tau) { return determinant.H(_s, _tau); };
  auto exact = [&](t_real _s, t_real _tau) -> t_rmatrix {
    return states.fa() * std::exp(-_s * _tau) * (_tau * states.aa()).exp() * states.af();
  };

  t_real svec[] = {0e0, 1e-1, 1e0};
  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dtaus[] = {1e-4, 1e-6, 1e-8};
  for(t_real s: svec) for(t_real tau: taus) for(t_real dtau: dtaus) {
    if(dtau >= tau) continue;
    taylor_convergence(approx, exact, s, tau, dtau, "Testing H\n");
  }
}

// Test s derivatives for non-zero critical tau using numerical and enalytical derivatives.
TEST_F(DeterminantEqTest, s_derivative_from_tau_derivative) {

  StateMatrix const states(Q, 2);
  DeterminantEq determinant(states, 0, false); 
  auto approx = [&determinant](t_real _s, t_real _tau) {
    return determinant.s_derivative(_s, _tau); 
  };
  auto exact = [&states](t_real _s, t_real _tau) -> t_rmatrix {
    return -_tau * std::exp(-_s * _tau) * states.fa() * (_tau * states.aa()).exp() * states.af();
  };

  t_real svec[] = {0e0, 1e-1, 1e0};
  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dtaus[] = {1e-4, 1e-6, 1e-8};
  for(t_real s: svec) for(t_real tau: taus) for(t_real dtau: dtaus) {
    if(dtau >= tau) continue;
    taylor_convergence(approx, exact, s, tau, dtau, "Testing s_derivative\n");
  }
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

