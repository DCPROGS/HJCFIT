#include <iostream>
#include <utility>
#include <gtest/gtest.h>
#include <iostream>
#include "../qmatrix.h"
using namespace DCProgs;

#ifdef HAS_CXX11_TYPETRAITS
  // Checks some assumption about eigen matrix types.
  static_assert( std::is_move_constructible<t_rmatrix>::value,
        	       "t_rmatrix is not move constructible." );  
  static_assert( std::is_move_assignable<t_rmatrix>::value, 
  	             "t_rmatrix is not move assignable." );  
  
  // Checks some assumption about QMatrix type.
  static_assert( std::is_move_constructible<QMatrix>::value,
  	             "QMatrix is not move constructible." );  
  static_assert( std::is_move_assignable<QMatrix>::value, 
  	             "QMatrix is not move assignable." );  
#endif
#ifdef HAS_CXX11_TRIVIALTYPETRAITS
  static_assert( not std::is_trivially_move_constructible<t_rmatrix>::value,
  	             "t_rmatrix is trivially move constructible." );  
  static_assert( not std::is_trivially_move_assignable<t_rmatrix>::value, 
  	             "t_rmatrix is trivially move assignable." );  
  static_assert( not std::is_trivially_move_constructible<QMatrix>::value,
  	            "QMatrix is trivially move constructible." );  
  static_assert( not std::is_trivially_move_assignable<QMatrix>::value, 
  	             "QMatrix is trivially move assignable." );  
#endif

// Checks some assumption about return of aa, af, fa, ff functions.
static_assert( not std::is_same<decltype(std::declval<QMatrix>().aa()),  t_rmatrix>::value,
               "QMatrix's aa function does not return an expression." );
static_assert( not std::is_same<decltype(std::declval<QMatrix>().af()),  t_rmatrix>::value,
               "QMatrix's af function does not return an expression." );
static_assert( not std::is_same<decltype(std::declval<QMatrix>().ff()),  t_rmatrix>::value,
               "QMatrix's ff function does not return an expression." );
static_assert( not std::is_same<decltype(std::declval<QMatrix>().fa()),  t_rmatrix>::value,
               "QMatrix's fa function does not return an expression." );


// Sets up test with parameters from CH82, 1e-7 nM.
class QMatrixTest : public ::testing::Test {
  
  public:
  QMatrixTest() : Q(), matrix() {}

  virtual void SetUp() {
    Q.matrix.resize(5, 5);
    Q.matrix <<  -3050,        50,  3000,      0,    0,
                 2./3., -1502./3.,     0,    500,    0, 
                    15,         0, -2065,     50, 2000, 
                     0,     15000,  4000, -19000,    0, 
                     0,         0,    10,      0,  -10;
    Q.nopen = 2;
    matrix = Q.matrix;
  }
  protected:
    QMatrix Q;
    t_rmatrix matrix;
};

TEST_F(QMatrixTest, blocks){
  Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
    diff = (Q.matrix  - matrix).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
 
  { t_rmatrix aa(2, 2); aa << -3050, 50, 2./3., -1502./3.;
    t_rmatrix af(2, 3); af << 3000, 0, 0, 0, 500, 0;
    t_rmatrix fa(3, 2); fa << 15, 0, 0, 15000, 0, 0;
    t_rmatrix ff(3, 3); ff << -2065, 50, 2000, 4000, -19000, 0, 10, 0, -10;
    Q.nopen = 2;
    EXPECT_TRUE( ((Q.aa() - aa).array().abs() < 1e-8).all() );
    EXPECT_TRUE( ((Q.af() - af).array().abs() < 1e-8).all() );
    EXPECT_TRUE( ((Q.fa() - fa).array().abs() < 1e-8).all() );
    EXPECT_TRUE( ((Q.ff() - ff).array().abs() < 1e-8).all() ); }
 
  { t_rmatrix aa(3, 3); aa << -3050, 50, 3000, 2./3., -1502./3., 0, 15, 0, -2065;
    t_rmatrix af(3, 2); af << 0, 0, 500, 0, 50, 2000;
    t_rmatrix fa(2, 3); fa << 0, 15000, 4000, 0, 0, 10;
    t_rmatrix ff(2, 2); ff << -19000, 0, 0, -10;
    Q.nopen = 3;
    EXPECT_TRUE( ((Q.aa() - aa).array().abs() < 1e-8).all() );
    EXPECT_TRUE( ((Q.af() - af).array().abs() < 1e-8).all() );
    EXPECT_TRUE( ((Q.fa() - fa).array().abs() < 1e-8).all() );
    EXPECT_TRUE( ((Q.ff() - ff).array().abs() < 1e-8).all() ); }
}


TEST_F(QMatrixTest, eigenvalues){
  
  auto const eigenstuff = Q.eigenstuff();
  auto const & vectors = std::get<1>(eigenstuff);
  auto const & eigenvalues = std::get<0>(eigenstuff);
 
  EXPECT_TRUE( std::abs(std::imag(eigenvalues(0))) < 1e-8 );
  EXPECT_TRUE( std::abs(std::real(eigenvalues(0)) + 3.09352723698141e+03) < 1e-8 );
 
  EXPECT_TRUE( std::abs(std::imag(eigenvalues(1))) < 1e-8 );
  EXPECT_TRUE( std::abs(std::real(eigenvalues(1)) + 2.02211926949769e+03) < 1e-8 );
 
  EXPECT_TRUE( std::abs(std::imag(eigenvalues(2))) < 1e-8 );
  EXPECT_TRUE( std::abs(std::real(eigenvalues(2))) < 1e-8 );
 
  EXPECT_TRUE( std::abs(std::imag(eigenvalues(3))) < 1e-8 );
  EXPECT_TRUE( std::abs(std::real(eigenvalues(3)) + 1.94082022553873e+04) < 1e-8 );
 
  EXPECT_TRUE( std::abs(std::imag(eigenvalues(4))) < 1e-8 );
  EXPECT_TRUE( std::abs(std::real(eigenvalues(4)) + 1.01817904800281e+02) < 1e-8 );
 
  for(size_t i(0); i < 5; ++i) {
    auto vector = vectors.row(i);
    auto eigenvalue = eigenvalues(i);
    EXPECT_TRUE( ( (eigenvalue * vector - vector * Q.matrix).array().abs() < 1e-8).all() );
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

