#include <iostream>
#include <type_traits>
#include <tuple>
#include <gtest/gtest.h>
#include <unsupported/Eigen/MatrixFunctions>
#include "../laplace_survivor.h"
using namespace DCProgs;

#ifdef HAS_CXX11_TYPETRAITS
  // Checks some assumption about eigen matrix types.
  static_assert( std::is_move_constructible<LaplaceSurvivor>::value,
  	             "LaplaceSurvivor is not move constructible." );  
  static_assert( std::is_move_assignable<LaplaceSurvivor>::value, 
  	             "LaplaceSurvivor is not move assignable." );  
#endif

#ifdef HAS_CXX11_TRIVIALTYPETRAITS
  static_assert( not std::is_trivially_move_constructible<LaplaceSurvivor>::value,
         	       "LaplaceSurvivor is trivially move constructible." );  
  static_assert( not std::is_trivially_move_assignable<LaplaceSurvivor>::value, 
  	             "LaplaceSurvivor is trivially move assignable." );  
#endif

// Sets up test with parameters from CH82, 1e-7 nM.
class LaplaceSurvivorTest : public ::testing::Test {
  
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
    t_bmatrix condition =   (diff.array().abs() < M * order2.array().abs()).eval()
                          + (diff.array().abs() < 1e-11).eval();
    EXPECT_TRUE( condition.all() ) << _message 
      << "Params:  x=" << _x << " -- dx=" << _dx << "\n"
      << "approx: \n" << approx << "\n"
      << "exact:  \n" << exact  << "\n"
      << "diff:  \n" << diff  << "\n"
      << "condition:  \n" << condition  << "\n"
      << "max (M=" << M << "):  \n" << (M * order2)  << "\n";
}


// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(LaplaceSurvivorTest, FF_resolution_is_zero) {
  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor R(qmatrix.transpose()); 
  EXPECT_TRUE( ((R.H(0, 0).array() - qmatrix.ff().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.H(1, 0).array() - qmatrix.ff().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.H(10, 0).array() - qmatrix.ff().array()).abs() < 1e-8).all() );

  auto value = [&qmatrix](t_real _s) -> t_rmatrix {
    return (_s * t_rmatrix::Identity(qmatrix.ff().rows(), qmatrix.ff().cols())
            - qmatrix.ff()).inverse(); 
  };
  EXPECT_TRUE( ((R(0, 0).array() - value(0).array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R(1, 0).array() - value(1).array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R(10, 0).array() - value(10).array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R(11, 0).array() - value(10).array()).abs() > 1e-8).any() );
}
// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(LaplaceSurvivorTest, FF_resolution_is_zero_check_derivative) {
  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor R(qmatrix.transpose()); 
  t_rmatrix id = t_rmatrix::Identity(3, 3);
  EXPECT_TRUE( ((R.s_derivative(0, 0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.s_derivative(1, 0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.s_derivative(10, 0).array() - id.array()).abs() < 1e-8).all() );
}
// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(LaplaceSurvivorTest, AA_resolution_is_zero) {
  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor R(qmatrix); 
  EXPECT_TRUE( ((R.H(0, 0).array() - qmatrix.aa().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.H(1, 0).array() - qmatrix.aa().array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.H(10, 0).array() - qmatrix.aa().array()).abs() < 1e-8).all() );

  auto value = [&qmatrix](t_real _s) -> t_rmatrix {
    return (_s * t_rmatrix::Identity(qmatrix.aa().rows(), qmatrix.aa().cols()) 
            - qmatrix.aa()).inverse(); 
  };
  EXPECT_TRUE( ((R(0, 0).array() - value(0).array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R(1, 0).array() - value(1).array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R(10, 0).array() - value(10).array()).abs() < 1e-8).all() );
}
// Missed event with t resolution == zero (e.g. no missed event)
TEST_F(LaplaceSurvivorTest, AA_resolution_is_zero_check_derivative) {
  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor R(qmatrix); 
  t_rmatrix id = t_rmatrix::Identity(2, 2);
  EXPECT_TRUE( ((R.s_derivative(0, 0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.s_derivative(1, 0).array() - id.array()).abs() < 1e-8).all() );
  EXPECT_TRUE( ((R.s_derivative(10, 0).array() - id.array()).abs() < 1e-8).all() );
}

// Test non-zero resolution tau using numerical and analytical derivatives.
TEST_F(LaplaceSurvivorTest, H_from_tau_derivative) {

  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor survivor(qmatrix.transpose()); 

  t_real svec[] = {0e0, 1e-1, 1e0};
  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dtaus[] = {1e-4, 1e-6, 1e-8};
  for(t_real s: svec) {
    auto approx = [&survivor, &s](t_real _tau) { return survivor.H(s, _tau); };
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
TEST_F(LaplaceSurvivorTest, s_derivative_from_tau_derivative) {

  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor survivor(qmatrix.transpose()); 

  t_real svec[] = {0e0, 1e-1, 1e0};
  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dtaus[] = {1e-4, 1e-6, 1e-8};
  for(t_real s: svec) {
    auto approx = [&survivor, &s](t_real _tau) { return survivor.s_derivative(s, _tau); };
    auto exact = [&qmatrix, &s](t_real _tau) -> t_rmatrix {
      return -_tau * std::exp(-s * _tau) * qmatrix.fa() * (_tau * qmatrix.aa()).exp() * qmatrix.af();
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
TEST_F(LaplaceSurvivorTest, s_derivative_from_H) {

  QMatrix const qmatrix(Q, 2);
  LaplaceSurvivor survivor(qmatrix.transpose()); 

  t_real taus[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real svec[] = {0e0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
  t_real dsvec[] = {1e-4, 1e-6, 1e-8};
  for(t_real tau: taus) {

    auto approx = [&survivor, &tau](t_real _s) { return survivor.H(_s, tau); };
    auto exact = [&survivor, &tau](t_real _s) -> t_rmatrix { 
      t_rmatrix const result = survivor.s_derivative(_s, tau); 
      return result - t_rmatrix::Identity(result.rows(), result.cols());
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

