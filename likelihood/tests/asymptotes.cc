#include "DCProgsConfig.h"
#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "../asymptotes.h"
using namespace DCProgs;

//! Params are whether to do open or closed states, 
//! and the root of the determinantal equation.
typedef std::tuple<bool, std::vector<Root>> t_Params;
// Sets up test with parameters from CH82, 1e-7 nM.
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
  DeterminantEq equation(StateMatrix(Q, 2), 1e-4, std::get<0>(GetParam()));
  Asymptotes asymptotes(equation, roots);

  t_rmatrix const result = asymptotes(0);
  EXPECT_EQ(result.rows(), std::get<0>(GetParam()) ? 2: 3);
  EXPECT_EQ(result.cols(), result.rows());
}

// Checks that left and right apply leave matrix untouched
TEST_P(TestAsymptotes, is_projection_matrix) {

  t_int const nopen = std::get<0>(GetParam()) ? 2: 3;
  DeterminantEq equation(StateMatrix(Q, 2), 1e-4, std::get<0>(GetParam()));

  for(auto root: std::get<1>(GetParam())) {
    // Creating asymptote equation with single root.
    // Makes it possible to test that each R_i matrix is indeed a projection matrix onto eigenvalue
    // correspinding to root.
    std::vector<Root> roots(1, root);
    Asymptotes asymptotes(equation, roots);
  
    for(t_int i(0); i < 5; ++i) {
      // Following tests imply that asymptotes is a factor of the projection matrix of H for the
      // eigenvalue root.root. The loop over different times ensure this is the case for more than
      // one time, eg time independent result.
      t_rmatrix const result = asymptotes(t_real(i) * equation.get_tau());
      t_rmatrix const HmI = equation.H(root.root) 
                            - t_rmatrix::Identity(nopen, nopen) * root.root; 
      EXPECT_TRUE(((HmI * result).array().abs() < 1e-8).all()) 
                 << "right applied projection matrix does not yield zero \n" 
                 << HmI * result << std::endl;
      EXPECT_TRUE(((result * HmI).array().abs() < 1e-8).all()) 
                 << "left applied projection matrix does not yield zero \n"
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

