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
    partition.resize(5);
    partition << true, true, false, false, false;
  }
  protected:
    IdealG idealg;
    t_rmatrix Q;
    t_partition partition;
};

TEST_F(IdealGTest, set){
  idealg.set(Q, partition);
  Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
    diff = (Q - idealg.get_Q()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  EXPECT_TRUE(partition == idealg.get_open_states());

  { t_rmatrix qq(5, 3);
    EXPECT_THROW(idealg.set(qq, partition), errors::Domain); }
  { t_rmatrix qq(3, 5);
    EXPECT_THROW(idealg.set(qq, partition), errors::Domain); }
  { t_partition pp(4);
    EXPECT_THROW(idealg.set(Q, pp), errors::Domain); }

  // Tests row constraints.
  for(size_t i(0); i < Q.rows(); ++i)
    EXPECT_DOUBLE_EQ(std::abs(idealg.get_Q().row(i).sum()), 0e0);

  // Test that row constraints always works.
  { t_rmatrix qq(Q);
    qq(1, 1) = 5e5;
    idealg.set(qq, partition);
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = (Q - idealg.get_Q()).array().abs();
    EXPECT_TRUE((diff < 1e-8).all());
  }
}
