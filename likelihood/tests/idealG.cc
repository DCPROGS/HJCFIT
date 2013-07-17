#include "DCProgsConfig.h"
#include <iostream>
#include <gtest/gtest.h>
#include "../idealG.h"
using namespace DCProgs;

#ifdef HAS_CXX11_TYPETRAITS
  // Checks some assumption about eigen matrix types.
  static_assert( std::is_move_constructible<IdealG>::value,
  	       "IdealG is not move constructible." );  
  static_assert( std::is_move_assignable<IdealG>::value, 
  	       "IdealG is not move assignable." );  
#endif

#ifdef HAS_CXX11_TRIVIALTYPETRAITS
  static_assert( not std::is_trivially_move_constructible<IdealG>::value,
  	       "IdealG is trivially move constructible." );  
  static_assert( not std::is_trivially_move_assignable<IdealG>::value, 
  	       "IdealG is trivially move assignable." );  
#endif

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

TEST_F(IdealGTest, initialize){
  idealg.set(Q, 2);
  Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic> diff = (Q - idealg.get_matrix()).array().abs();
  EXPECT_TRUE((diff < 1e-8).all());
  EXPECT_EQ(idealg.get_nopen(), 2);

  { t_rmatrix qq(5, 3);
    EXPECT_THROW(idealg.set(qq, 2), errors::Domain); }
// { t_rmatrix qq(3, 5);
//   EXPECT_THROW(idealg.set(qq, 2), errors::Domain); }
// EXPECT_THROW(idealg.set(Q, 6), errors::Domain); 
// EXPECT_THROW(idealg.set(Q, -1), errors::Domain); 
 
  // Tests row constraints.
  for(t_int i(0); i < Q.rows(); ++i)
    EXPECT_DOUBLE_EQ(std::abs(idealg.get_matrix().row(i).sum()), 0e0);
 
  // Test that row constraints always works.
  { t_rmatrix qq(Q);
    qq(1, 1) = 5e5;
    idealg.set(qq, 2);
    Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
      diff = (Q - idealg.get_matrix()).array().abs();
    EXPECT_TRUE((diff < 1e-8).all());
  }
}

TEST_F(IdealGTest, blocks){
  QMatrix qmatrix(Q, 2);
  idealg.set(Q, 2);

  // This test pretty much ensures that we are dealing with an exponential
  // At least over 10 integers. 
  { t_rmatrix exponential = qmatrix.aa().exp();
    t_rmatrix current = qmatrix.af();
    for(size_t i(0); i < 21; ++i, current = exponential * current) {
      Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
        diff = (idealg.af(t_real(i)) - current).array().abs();
      EXPECT_TRUE((diff < 1e-8).all()); 
    }
  }
  { t_rmatrix exponential = qmatrix.ff().exp();
    t_rmatrix current = qmatrix.fa();
    for(size_t i(0); i < 21; ++i, current = exponential * current) {
      Eigen::Array<t_real, Eigen::Dynamic, Eigen::Dynamic>
        diff = (idealg.fa(t_real(i)) - current).array().abs();
      EXPECT_TRUE((diff < 1e-8).all()); 
    }
  }
}

TEST_F(IdealGTest, laplacians) {

  idealg.set(Q, 2);

  t_rmatrix af0(2, 3);
  af0 << 0.98362802881467, 0.01637197118533,  0., 
         0.00130975769483, 0.99869024230517,  0.;
  EXPECT_TRUE(((idealg.laplace_af(0) - af0).array().abs() < 1e-8).all());
  t_rmatrix af1(2, 3);
  af1 << 0.98330558371655, 0.01633397979596, 0.,
         0.00130671838368, 0.99669944714923, 0.;
  EXPECT_TRUE(((idealg.laplace_af(1) - af1).array().abs() < 1e-8).all());

  t_rmatrix fa0(3, 2);
  fa0 << 0.27536231884058, 0.7246376811594,
         0.05797101449275, 0.94202898550724,
         0.27536231884058, 0.7246376811594;
  EXPECT_TRUE(((idealg.laplace_fa(0) - fa0).array().abs() < 1e-8).all());
  t_rmatrix fa1(3, 2);
  fa1 << 0.063213144351,  0.16634162505001,
         0.013307330004,  0.8244495816115,
         0.057466494865,  0.15121965913637;
  EXPECT_TRUE(((idealg.laplace_fa(1) - fa1).array().abs() < 1e-8).all());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

