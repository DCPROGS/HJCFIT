#include <iostream>
#include <utility>
#include <gtest/gtest.h>
#include "../idealG.h"
using namespace DCProgs;

// Checks some assumption about eigen matrix types.
static_assert( std::is_move_constructible<t_rmatrix>::value,
               "t_rmatrix is not move constructible." );  
static_assert( not std::is_trivially_move_constructible<t_rmatrix>::value,
               "t_rmatrix is trivially move constructible." );  
static_assert( std::is_move_assignable<t_rmatrix>::value, 
               "t_rmatrix is not move assignable." );  
static_assert( not std::is_trivially_move_assignable<t_rmatrix>::value, 
               "t_rmatrix is trivially move assignable." );  

// Checks some assumption about StateMatrix type.
static_assert( std::is_move_constructible<StateMatrix>::value,
               "StateMatrix is not move constructible." );  
static_assert( not std::is_trivially_move_constructible<StateMatrix>::value,
               "StateMatrix is trivially move constructible." );  
static_assert( std::is_move_assignable<StateMatrix>::value, 
               "StateMatrix is not move assignable." );  
static_assert( not std::is_trivially_move_assignable<StateMatrix>::value, 
               "StateMatrix is trivially move assignable." );  

// Checks some assumption about return of aa, af, fa, ff functions.
static_assert( not std::is_same<decltype(std::declval<StateMatrix>().aa()),  t_rmatrix>::value,
               "StateMatrix's aa function does not return an expression." );
static_assert( not std::is_same<decltype(std::declval<StateMatrix>().af()),  t_rmatrix>::value,
               "StateMatrix's af function does not return an expression." );
static_assert( not std::is_same<decltype(std::declval<StateMatrix>().ff()),  t_rmatrix>::value,
               "StateMatrix's ff function does not return an expression." );
static_assert( not std::is_same<decltype(std::declval<StateMatrix>().fa()),  t_rmatrix>::value,
               "StateMatrix's fa function does not return an expression." );


// Sets up test with parameters from CH82, 1e-7 nM.
class StateMatrixTest : public ::testing::Test {
  
  public:
  StateMatrixTest() {}

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
    StateMatrix Q;
    t_rmatrix matrix;
};

TEST_F(StateMatrixTest, blocks){
  Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
    diff = (Q.matrix  - matrix).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());

  Q.nopen = 2;
  diff = (matrix.topLeftCorner(2, 2) - Q.aa()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  diff = (matrix.topRightCorner(2, 2) - Q.af()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  diff = (matrix.bottomLeftCorner(2, 2) - Q.fa()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  diff = (matrix.bottomRightCorner(2, 2) - Q.ff()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());

  Q.nopen = 3;
  diff = (matrix.topLeftCorner(3, 3) - Q.aa()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  diff = (matrix.topRightCorner(3, 3) - Q.af()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  diff = (matrix.bottomLeftCorner(3, 3) - Q.fa()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  diff = (matrix.bottomRightCorner(3, 3) - Q.ff()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
}

