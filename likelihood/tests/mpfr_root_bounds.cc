// This is meant as an example where the default double precision rootfinding fails
// but the fallback to mpfr works. As such it may fail if rootfinding it tweeked to make 
// this work with doubles

#include <iostream>

#include "DCProgsConfig.h"
#include "../determinant_equation.h"
#include "../root_finder.h"
#include "../brentq.h"
#include <gtest/gtest.h>
using namespace DCProgs;

class MPFRROOTTest : public ::testing::Test {
  
  public:
  MPFRROOTTest() {}

  virtual void SetUp() {
    Q.resize(10, 10);
    Q << -2.63635924e+03, 0.00000000e+00, 0.00000000e+00, 2.63635924e+03, 0.00000000e+00, 0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,  -1.31139767e+03,   0.00000000e+00,   0.00000000e+00, 1.31139767e+03,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,  -6.82928534e+02,   0.00000000e+00, 0.00000000e+00,   6.82928534e+02,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00,
         9.44807825e+01,   0.00000000e+00,   0.00000000e+00,  -1.10439558e+06, 2.27110722e+04,   0.00000000e+00,   1.08159003e+06,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   1.27102751e+04,   0.00000000e+00,   1.99674768e+04, -5.81182261e+04,   1.13555361e+04,   0.00000000e+00,   1.40849381e+04,0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   2.36790649e+04,   0.00000000e+00, 2.99512152e+04,  -5.56479033e+04,   0.00000000e+00,   0.00000000e+00, 2.01762324e+03,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   9.77672433e+04, 0.00000000e+00,   0.00000000e+00,  -1.01997759e+05,   5.67571675e+02, 0.00000000e+00,   3.66294412e+03,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 1.86913065e+04,   0.00000000e+00,   7.32588824e+03,  -2.63009806e+04, 2.83785838e+02,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   3.93078069e+04,   0.00000000e+00,   1.09888324e+04, -5.02966392e+04,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, 0.00000000e+00,   0.00000000e+00,   8.51357513e+02,   0.00000000e+00, 0.00000000e+00,  -8.51357513e+02;
  }
  protected:
    t_rmatrix Q;
};

TEST_F(MPFRROOTTest, initialize){

  // // Define parameters.
  QMatrix const qmatrix(Q, 3);
  DeterminantEq det(qmatrix, 0.000030);
  // 
  // // Find upper and lower bound
  t_real upper_bound = find_upper_bound_for_roots(det.transpose(), 0.0, 5e0, 100);
  t_real lower_bound = find_lower_bound_for_roots(det.transpose(), 0.0, 5e0, 100);
  EXPECT_DOUBLE_EQ(upper_bound, 0.0);
  EXPECT_DOUBLE_EQ(lower_bound, -1801505.9230927574);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
