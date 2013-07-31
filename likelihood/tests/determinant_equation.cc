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
#include <tuple>
#include <gtest/gtest.h>
#include <unsupported/Eigen/MatrixFunctions>
#include "../determinant_equation.h"
using namespace DCProgs;

#ifdef HAS_CXX11_TYPETRAITS
  // Checks some assumption about eigen matrix types.
  static_assert( std::is_move_constructible<DeterminantEq>::value,
         	       "DeterminantEq is not move constructible." );  
  static_assert( std::is_move_assignable<DeterminantEq>::value, 
  	             "DeterminantEq is not move assignable." );  
#endif

#ifdef HAS_CXX11_TRIVIALTYPETRAITS
  static_assert( not std::is_trivially_move_constructible<DeterminantEq>::value,
  	             "DeterminantEq is trivially move constructible." );  
  static_assert( not std::is_trivially_move_assignable<DeterminantEq>::value, 
  	             "DeterminantEq is trivially move assignable." );  
#endif

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
                           t_real _x, t_real _dx,
                           std::string const &_message = "", t_real M = 10e0 ) {
    t_rmatrix const order2 = 0.5 * _dx * (_exact(_x + _dx) - _exact(_x));
    t_rmatrix const approx = _approx(_x + _dx) - _approx(_x);
    t_rmatrix const exact  = _dx * _exact(_x);
    t_rmatrix const diff = approx - exact;
    t_bmatrix condition =    (diff.array().abs() < M * order2.array().abs())
                          || (diff.array().abs() < 1e-11);
    EXPECT_TRUE( condition.all() ) << _message 
      << "Params:  x=" << _x << " -- dx=" << _dx << "\n"
      << "approx: \n" << approx << "\n"
      << "exact:  \n" << exact  << "\n"
      << "diff:  \n" << diff  << "\n"
      << "condition:  \n" << condition  << "\n"
      << "max (M=" << M << "):  \n" << (M * order2)  << "\n";
}


// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, FF_resolution_is_zero) {
  QMatrix const qmatrix(Q, 2);
  DeterminantEq det(qmatrix.transpose(), 0); 
  EXPECT_TRUE( ((det.H(0).array() - qmatrix.ff().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(1).array() - qmatrix.ff().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(10).array() - qmatrix.ff().array()).abs() < 1e-8).all() );

  Eigen::EigenSolver<t_rmatrix> eigsolver(qmatrix.ff());
  t_rvector const eigs = eigsolver.eigenvalues().real();
  EXPECT_TRUE(std::abs(det(eigs(0)) / eigs(0)) < 1e-5);
  EXPECT_TRUE(std::abs(det(eigs(1)) / eigs(1)) < 1e-5);
  EXPECT_TRUE(std::abs(det(eigs(2)) / eigs(2)) < 1e-5);
}
// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, FF_resolution_is_zero_check_derivative) {
  QMatrix const qmatrix(Q, 2);
  DeterminantEq det(qmatrix.transpose(), 0); 
  t_rmatrix id = t_rmatrix::Identity(3, 3);
  EXPECT_TRUE( ((det.s_derivative(0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(1).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(10).array() - id.array()).abs() < 1e-8).all() );
}
// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, AA_resolution_is_zero) {
  QMatrix const qmatrix(Q, 2);
  DeterminantEq det(qmatrix, 0); 
  EXPECT_TRUE( ((det.H(0).array() - qmatrix.aa().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(1).array() - qmatrix.aa().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.H(10).array() - qmatrix.aa().array()).abs() < 1e-8).all() );

  Eigen::EigenSolver<t_rmatrix> eigsolver(qmatrix.aa());
  t_rvector const eigs = eigsolver.eigenvalues().real();
  EXPECT_TRUE(std::abs(det(eigs(0)) / eigs(0)) < 1e-5);
  EXPECT_TRUE(std::abs(det(eigs(1)) / eigs(1)) < 1e-5);
}
// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(DeterminantEqTest, AA_resolution_is_zero_check_derivative) {
  QMatrix const qmatrix(Q, 2);
  DeterminantEq det(qmatrix, 0); 
  t_rmatrix id = t_rmatrix::Identity(2, 2);
  EXPECT_TRUE( ((det.s_derivative(0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(1).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((det.s_derivative(10).array() - id.array()).abs() < 1e-8).all() );
}

// Test non-zero resolution tau using numerical and analytical derivatives.
TEST_F(DeterminantEqTest, from_tau_derivative) {

  QMatrix const qmatrix(Q, 2);
  DeterminantEq determinant(qmatrix.transpose(), 0); 

  t_real svec[] = {0e0, 1e-1, 1e0};
  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dtaus[] = {1e-4, 1e-6, 1e-8};
  for(t_real s: svec) {
    auto approx = [&determinant, &s](t_real _tau) -> t_rmatrix { return determinant.H(s, _tau); };
    auto exact = [&qmatrix, &s](t_real _tau) -> t_rmatrix {
      return qmatrix.fa() * std::exp(-s * _tau) * (_tau * qmatrix.aa()).exp() * qmatrix.af();
    };
    std::ostringstream sstr;
    sstr << "Testing H with s=" << s << "\n";
    std::string const message = sstr.str();

    for(t_real tau: taus)
      for(t_real dtau: dtaus) 
        taylor_convergence(approx, exact, tau, dtau, sstr.str());
  }
}

// Test s derivatives for non-zero resolution tau using numerical and analytical derivatives.
TEST_F(DeterminantEqTest, s_derivative_from_tau_derivative) {

  QMatrix const qmatrix(Q, 2);
  DeterminantEq determinant(qmatrix.transpose(), 0); 

  t_real svec[] = {0e0, 1e-1, 1e0};
  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dtaus[] = {1e-4, 1e-6, 1e-8};
  for(t_real s: svec) {
    auto approx = [&determinant, &s](t_real _tau) { return determinant.s_derivative(s, _tau); };
    auto exact = [&qmatrix, &s](t_real _tau) -> t_rmatrix {
      return _tau * std::exp(-s * _tau) * qmatrix.fa() * (_tau * qmatrix.aa()).exp() * qmatrix.af();
    };
    
    std::ostringstream sstr;
    sstr << "Testing s_derivative with s=" << s << "\n";
    std::string const message = sstr.str();

    for(t_real tau: taus)
      for(t_real dtau: dtaus) 
        taylor_convergence(approx, exact, tau, dtau, sstr.str());
  }
}

// Test s derivative from numerical derivative of H.
TEST_F(DeterminantEqTest, s_derivative_from_H) {

  QMatrix const qmatrix(Q, 2);
  DeterminantEq determinant(qmatrix.transpose(), 0); 

  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real svec[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dsvec[] = {1e-4, 1e-6, 1e-8};
  for(t_real tau: taus) {

    auto approx = [&determinant, &tau](t_real _s) { return determinant.H(_s, tau); };
    auto exact = [&determinant, &tau](t_real _s) -> t_rmatrix { 
      t_rmatrix const result = determinant.s_derivative(_s, tau); 
      return t_rmatrix::Identity(result.rows(), result.cols()) - result;
    };
    
    std::ostringstream sstr;
    sstr << "Testing s_derivative with tau=" << tau << "\n";
    std::string const message = sstr.str();

    for(t_real s: svec)
      for(t_real ds: dsvec) 
        taylor_convergence(approx, exact, s, ds, sstr.str());
  }
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

