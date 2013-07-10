#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>
#include "../idealG.h"
#include "../occupancies.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class EquilibriumTest : public ::testing::Test {
  
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


TEST_F(EquilibriumTest, initial_states){
  IdealG idealg(Q, 2);
  t_initvec const phiA = occupancies(idealg);
  EXPECT_EQ(phiA.size(), 2);
  EXPECT_NEAR(phiA(0), 10./135., 1e-10);
  EXPECT_NEAR(phiA(1), 1 - 10./135., 1e-10);

  t_initvec const phiF = occupancies(idealg, false);
  EXPECT_EQ(phiF.size(), 3);
  EXPECT_NEAR(phiF(0), 10./135., 1e-10);
  EXPECT_NEAR(phiF(1), 1 - 10./135., 1e-10);
  EXPECT_NEAR(phiF(2), 0e0, 1e-10);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

