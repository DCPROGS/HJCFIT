#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>
#include "../recursion_formula.h"
using namespace DCProgs;

t_real ZERO() { return 0e0; }
struct FullRecursion {
  typedef t_real t_element;

  FullRecursion() {
     D.resize(5);
     D << 1, 2, 3, 4, 5;
     eigvals.resize(5);
     eigvals << 1./1, 1./2., 1./3., 1./4., 1./5.;
     initials.resize(5);
     initials << 1, 1, 1, 1, 1;
  }

  t_element operator()(t_int _i, t_int _j, t_int _k) const {
    if(_j == 0 and _k == 0) return initials(_i);
    return recursion_formula(*this, _i, _j, _k, ZERO);
  }
  t_element getD(t_int _i) const { return D(_i); }
  t_element get_eigvals(t_int _i) const { return eigvals(_i); }
  t_int nbeigvals() const { return eigvals.size(); }

  t_rvector D, eigvals, initials;
};


TEST(TestRecursion, m_equal_l_InitialZero) {
  FullRecursion fullrecursion;
  fullrecursion.initials << 0, 0, 0, 0, 0;

  for(t_int i(0); i < 5; ++i) {
    for(t_int m(0); m < 10; ++m) 
      EXPECT_DOUBLE_EQ(recursion_formula(fullrecursion, i, m, m, ZERO), 0e0);
  }
}

t_int factorial(t_int _i) { return _i == 0 ? 1: _i * factorial(_i-1); }


TEST(TestRecursion, m_equal_l_InitialOnes) {
  FullRecursion fullrecursion;

  for(t_int i(0); i < 5; ++i) {
    EXPECT_DOUBLE_EQ(recursion_formula(fullrecursion, i, 0, 0, ZERO), fullrecursion.initials(i));
    for(t_int m(1); m < 10; ++m)  {
      EXPECT_NEAR( recursion_formula(fullrecursion, i, m, m, ZERO), 
                   std::pow(fullrecursion.D(i), m) / t_real(factorial(m)) * fullrecursion.initials(i),
                   1e-8 );
    }
  }
}

TEST(TestRecursion, m_equal_l_InitialRandom) {
  FullRecursion fullrecursion;
  for(t_int trial(0); trial < 500; ++trial) {
    fullrecursion.initials = t_rvector::Random(5);
  
    for(t_int i(0); i < 5; ++i) {
      EXPECT_DOUBLE_EQ(recursion_formula(fullrecursion, i, 0, 0, ZERO), fullrecursion.initials(i));
      for(t_int m(1); m < 10; ++m)  {
        EXPECT_NEAR( recursion_formula(fullrecursion, i, m, m, ZERO), 
                     std::pow(fullrecursion.D(i), m) / t_real(factorial(m)) * fullrecursion.initials(i),
                     1e-8 );
      }
    }
  }
}
 
TEST(TestRecursion, lzero_InitialZero) {
  FullRecursion fullrecursion;
  fullrecursion.initials << 0, 0, 0, 0, 0;

  for(t_int i(0); i < 5; ++i) {
    for(t_int m(0); m < 5; ++m) 
      EXPECT_DOUBLE_EQ(recursion_formula(fullrecursion, i, m, 0, ZERO), 0e0);
  }
}


class TestViaMatrix : public ::testing::TestWithParam<t_int> {
  public:
    struct Mocking {
      typedef t_real t_element;
      t_rvector D, eigvals;
      t_rmatrix valmat;
      t_int current_m;
  
      t_element operator()(t_int _i, t_int _j, t_int _k) const {
        EXPECT_EQ(_j, current_m) << "Unexpected m value in  recursion.";
        return valmat(_i, _k);
      }
      t_element getD(t_int _i) const { return D(_i); }
      t_element get_eigvals(t_int _i) const { return eigvals(_i); }
      t_int nbeigvals() const { return eigvals.size(); }
    };
  
    TestViaMatrix() {
       jay.D.resize(matsize);
       jay.eigvals.resize(matsize);
       for(t_int i(0); i < matsize; ++i) {
         jay.D(i) = i+1; 
         jay.eigvals(i) = 1e-1 * t_real(i+1);
       }
    
       for(t_int i(0); i < matsize; ++i) {
       
         t_rmatrix term0 = t_rmatrix::Zero(jay.nbeigvals(), mmax);
         t_rmatrix term1 = t_rmatrix::Zero(jay.nbeigvals(), mmax);
         for(t_int j(0); j < jay.nbeigvals(); ++j) {
           if(i == j) continue;
           t_real const lambda_diff = 1e0 / (jay.eigvals(j) - jay.eigvals(i));
           t_real factor(lambda_diff);
           for(t_int r(0), s(1); r < mmax; ++r, s = -s) {
             term0(j, r) = factor * jay.D(i);
             term1(j, r) = s * factor * jay.D(j);
             factor *= (r + 1) * lambda_diff;
           }
         }
         terms0.push_back(term0);
         terms1.push_back(term1);
       }
    }
    

  protected:
    Mocking jay;
    //! \f$\frac{D_i}{(\lambda_i-\lambda_j)^{r+1}}, with i std::vector first index
    std::vector<t_rmatrix> terms0;
    //! \f$\frac{D_j}{(\lambda_i-\lambda_j)^{r+1}}, with i std::vector first index
    std::vector<t_rmatrix> terms1;
    //! \f$\frac{D_j}{(\lambda_i-\lambda_j)^{r+1}}, with i std::vector first index
    t_int const mmax = 8;
    t_int const matsize = 5;
};

TEST_P(TestViaMatrix, lzero) {

  // Create the matrix by which to multiply.
  jay.valmat = t_rmatrix::Zero(jay.nbeigvals(), mmax);
  switch(GetParam()) {
    case 0: jay.valmat(0, 0) = 1; break;
    case 1: jay.valmat(1, 1) = 1; break;
    case 2: jay.valmat = t_rmatrix::Ones(jay.nbeigvals(), mmax); break;
    default: jay.valmat = t_rmatrix::Random(jay.nbeigvals(), mmax);
  }


  // Now do checks
  for(t_int i(0); i < matsize; ++i) {

    t_rmatrix const term0 = terms0[i].array() * jay.valmat.array();
    for(t_int m(1); m < mmax; ++m) {
      jay.current_m = m - 1; // Checks we don't make inner calls
      t_real const term1= (terms1[i].leftCols(m)
                           * jay.valmat.row(i).head(m).transpose()).array().sum();
      
      t_real const check = term0.topLeftCorner(term0.rows(), m).sum() - term1; 
      EXPECT_NEAR(recursion_formula(jay, i, m, 0, ZERO), check, std::abs(check) * 1e-8)
          << "i=" << i << " m=" << m;
                   
    }
  }
}

TEST_P(TestViaMatrix, general) {

  // Create the matrix by which to multiply.
  jay.valmat = t_rmatrix::Zero(jay.nbeigvals(), mmax);
  switch(GetParam()) {
    case 0: jay.valmat(0, 0) = 1; break;
    case 1: jay.valmat(1, 1) = 1; break;
    case 2: jay.valmat = t_rmatrix::Ones(jay.nbeigvals(), mmax); break;
    default: jay.valmat = t_rmatrix::Random(jay.nbeigvals(), mmax);
  }

  
  // Now do checks
  for(t_int i(0); i < matsize; ++i) {
    t_rmatrix const &term = terms1[i];
    auto valrow = jay.valmat.row(i);
    for(t_int m(2); m < mmax; ++m) {
      jay.current_m = m - 1; // Checks we don't make inner calls
      for(t_int l(1); l < m; ++l) {


        //! dcl: \f$D_i C_{i(m-1)(l-1)}/l\f$
        t_real dcl = jay.D(i) * jay.valmat(i, l-1) / t_real(l);
        //! Now add other term to dcl
        for(t_int j(0); j < matsize; ++j) {
          if(i == j) continue;
          t_real const factor = jay.D(j)  / term(j, l-1) / t_real(l);
          dcl -= (valrow.segment(l, m-l) * term.row(j).segment(l, m-l).transpose() * factor)(0, 0);
        }
          
        EXPECT_NEAR( recursion_formula(jay, i, m, l, ZERO), dcl, std::abs(dcl) * 1e-8 )
          << "i=" << i << " m=" << m << " l=" << l;
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P(SingleRecursion, TestViaMatrix, ::testing::Range(0, 300));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

