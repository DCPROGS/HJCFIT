#include <iostream>
#include <gtest/gtest.h>
#include "../idealG.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class IdealGTest : public ::testing::Test {
  
  public:
  IdealGTest() {}

  virtual void SetUp() {
    Q.resize(5, 5);
    Q <<  -3050,        50,  3000,      0,    0,
          2./3., -1502./3.,     0,    500,    0, 
             15,         0, -2065,     50, 2000, 
              0,     15000,  4000, -19000,    0, 
              0,         0,    10,      0,  -10;
  }
  protected:
    IdealG idealg;
    t_rmatrix Q;
};

TEST_F(IdealGTest, set){
  idealg.set(Q, 2);
  Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
    diff = (Q - idealg.get_Q()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  EXPECT_EQ(idealg.get_nopen(), 2);

  { t_rmatrix qq(5, 3);
    EXPECT_THROW(idealg.set(qq, 2), errors::Domain); }
  { t_rmatrix qq(3, 5);
    EXPECT_THROW(idealg.set(qq, 2), errors::Domain); }
  EXPECT_THROW(idealg.set(Q, 6), errors::Domain); 
  EXPECT_THROW(idealg.set(Q, -1), errors::Domain); 
 
  // Tests row constraints.
  for(size_t i(0); i < Q.rows(); ++i)
    EXPECT_DOUBLE_EQ(std::abs(idealg.get_Q().row(i).sum()), 0e0);
 
  // Test that row constraints always works.
  { t_rmatrix qq(Q);
    qq(1, 1) = 5e5;
    idealg.set(qq, 2);
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = (Q - idealg.get_Q()).array().abs();
    EXPECT_TRUE((diff < 1e-8).all());
  }
}
