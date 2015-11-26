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
#include "DCProgsConfig.h"
#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "../asymptotes.h"
using namespace DCProgs;
#ifdef HAS_CXX11_TYPE_TRAITS
  // Checks some assumption about Asymptotes matrix types.
  static_assert( std::is_move_constructible<Asymptotes>::value,
        	       "Asymptotes is not move constructible." );  
  static_assert( std::is_move_assignable<Asymptotes>::value, 
        	       "Asymptotes is not move assignable." );  
#endif

//! Params are whether to do open or closed states, 
//! and the root of the determinant equation.
typedef std::tuple<bool, std::vector<Root>> t_Params;

class TestAsymptotes : public ::testing::TestWithParam<t_Params> {
  
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


// Checks the size of the matrices is correct.
TEST_P(TestAsymptotes, correct_size) {

  std::vector<Root> roots = std::get<1>(GetParam());
  QMatrix qmatrix = std::get<0>(GetParam()) ? QMatrix(Q, 2): QMatrix(Q, 2).transpose();
  DeterminantEq equation(qmatrix, 1e-4);
  Asymptotes asymptotes(equation, roots);

  t_rmatrix const result = asymptotes(0);
  EXPECT_EQ(result.rows(), std::get<0>(GetParam()) ? 2: 3);
  EXPECT_EQ(result.cols(), result.rows());
}

// Checks that left and right apply leave matrix untouched
TEST_P(TestAsymptotes, is_projection_matrix) {

  t_uint const nopen = std::get<0>(GetParam()) ? 2: 3;
  QMatrix qmatrix = std::get<0>(GetParam()) ? QMatrix(Q, 2): QMatrix(Q, 2).transpose();
  DeterminantEq equation(qmatrix, 1e-4);

  for(auto root: std::get<1>(GetParam())) {
    // Creating asymptote equation with single root.
    // Makes it possible to test that each R_i matrix is indeed a projection matrix onto eigenvalue
    // correspinding to root.
    std::vector<Root> roots(1, root);
    Asymptotes asymptotes(equation, roots);
  
    for(t_int i(0); i < 3; ++i) {
      // Following tests imply that asymptotes is a factor of the projection matrix of H for the
      // eigenvalue root.root. The loop over different times ensure this is the case for more than
      // one time, eg time independent result.
      t_rmatrix const result = asymptotes(t_real(i) * equation.get_tau());
      t_rmatrix const HmI = equation.H(root.root) 
                            - t_rmatrix::Identity(nopen, nopen) * root.root; 
      EXPECT_TRUE(((HmI * result).array().abs() < 1e-8).all()) 
                 << "right applied projection matrix does not yield zero \n" 
                 << result << std::endl << std::endl
                 << HmI * result << std::endl;
      EXPECT_TRUE(((result * HmI).array().abs() < 1e-8).all()) 
                 << "left applied projection matrix does not yield zero \n"
                 << result << std::endl << std::endl
                 << result * HmI << std::endl;
    }
  }
}

t_Params create_open_params() {
  std::vector<Root> roots;
  roots.emplace_back(-3045.285776037674, 1);
  roots.emplace_back(-162.92946543451328, 1);
  return std::make_tuple(true, roots);
}

t_Params create_closed_params() {
  std::vector<Root> roots;
  roots.emplace_back(-17090.192769236815, 1);
  roots.emplace_back(-2058.0812921673496, 1);
  roots.emplace_back(-0.24356535498785126, 1);
  return std::make_tuple(false, roots);
}

INSTANTIATE_TEST_CASE_P( ClassicMatrix, TestAsymptotes,
                         ::testing::Values(create_open_params(), create_closed_params()) );

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

