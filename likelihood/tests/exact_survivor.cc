#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "../exact_survivor.h"
using namespace DCProgs;

// Sets up test with parameters from CH82, 1e-7 nM.
class ExactSurvivorTest : public ::testing::Test {
  
  public:
  virtual void SetUp() {
    Q.resize(5, 5);
    Q <<  -3050,        50,  3000,      0,    0,
          2./3., -1502./3.,     0,    500,    0, 
             15,         0, -2065,     50, 2000, 
              0,     15000,  4000, -19000,    0, 
              0,         0,    10,      0,  -10;
    
    nopen = 2; nclose = 3;

    Eigen::EigenSolver<t_rmatrix> eigsolver(Q);
    if(eigsolver.info() != Eigen::Success) 
        throw errors::Mass("Could not solve eigenvalue problem.");
    eigenvalues = -eigsolver.eigenvalues().real();
    eigenvectors = eigsolver.eigenvectors().real();
    eigenvectors_inv = eigenvectors.inverse();
  }

  t_rmatrix N0(t_real _t, t_real _tau, bool _a = true)  {
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int i(0); i < eigenvalues.size(); ++i) {
      if(std::abs(eigenvalues(i) > 1e-12))  {
        result += get_ci00(i, _a) * std::exp(-eigenvalues(i)*_t);
      }
    }
    return result;
  }
  t_rmatrix N1(t_real _t, t_real _tau, bool _a = true)  {
    t_real const t(_t - _tau);
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int i(0); i < eigenvalues.size(); ++i) {
      if(std::abs(eigenvalues(i) > 1e-12)) 
        result += (get_ci10(i, _tau, _a)  + get_ci11(i, _tau, _a) * _t) * std::exp(-eigenvalues(i)*_t);
    }
    return result;
  }

  t_rmatrix get_di(t_int _i, t_real _tau, bool _a = true) {
    t_rvector const left = _a ? eigenvectors.col(_i).head(nopen): eigenvectors.col(_i).tail(nclose);
    t_rvector const right = _a ? eigenvectors_inv.row(_i).tail(nclose):
                                 eigenvectors_inv.row(_i).head(nopen);
    t_rmatrix const exponent = _a ? Q.bottomRightCorner(nclose, nclose):
                                    Q.topLeftCorner(nopen, nopen);
    t_rmatrix const endvec = _a ? Q.bottomLeftCorner(nclose, nopen):  
                                  Q.topRightCorner(nopen, nclose);
    return (left * right.transpose()) * (_tau * exponent).exp() * endvec;
  }
  t_rmatrix get_ci00(t_int _i, bool _a = true) {
    t_rvector const left = _a ? eigenvectors.col(_i).head(nopen):
                                eigenvectors.col(_i).tail(nclose);
    t_rvector const right = _a ? eigenvectors_inv.row(_i).head(nopen):
                                 eigenvectors_inv.row(_i).tail(nclose);
    return left * right.transpose();
  }
  t_rmatrix get_ci10(t_int _i, t_real _tau, bool _a = true) {
    t_int const N = _a ? nopen: nclose;
    t_rmatrix result = t_rmatrix::Zero(N, N);
    for(t_int j(0); j < eigenvalues.size(); ++j) {
      t_real const delta(eigenvalues(j) - eigenvalues(_i));
      if(std::abs(delta) > 1e-8) {
        result += (get_di(_i, _tau, _a) * get_ci00(j, _a) + get_di(j, _tau, _a) * get_ci00(_i, _a)) 
                  / delta;
      } 
    }
    return result;
  }
  t_rmatrix get_ci11(t_int _i, t_real _tau, bool _a = true) {
    return get_di(_i, _tau, _a) * get_ci00(_i, _a);
  }
  protected:
    t_int nopen, nclose;
    t_rmatrix Q;
    t_rvector eigenvalues;
    t_rmatrix eigenvectors;
    t_rmatrix eigenvectors_inv;
};

// Compares recursive implementation to the expanded one in this file.
TEST_F(ExactSurvivorTest, negative_times) {
  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  EXPECT_TRUE( (survivor.af(-1e-5).array().abs() < 1e-8).all()  );
  EXPECT_EQ(survivor.af(-1e-5).rows(), 2);
  EXPECT_EQ(survivor.af(-1e-5).cols(), 2);
  EXPECT_TRUE( (survivor.fa(-1e-5).array().abs() < 1e-8).all()  );
  EXPECT_EQ(survivor.fa(-1e-5).rows(), 2);
  EXPECT_EQ(survivor.fa(-1e-5).cols(), 2);
}

// Compares recursive implementation to the expanded one in this file.
TEST_F(ExactSurvivorTest, first_interval) {
  std::cout.precision(15);
  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  auto aR0 = [&](t_real _t, t_real _tau) { return N0(_t, _tau, true); };
  auto fR0 = [&](t_real _t, t_real _tau) { return N0(_t, _tau, false); };
  auto compare = [](t_rmatrix const &_a, t_rmatrix const & _b) {
    return ((_a - _b).array().abs() < 1e-10).all();
  };
  EXPECT_TRUE( compare(survivor.af(1e-5), aR0(1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-5), fR0(1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(3e-5), aR0(3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(3e-5), fR0(3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(5e-5), aR0(5e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(5e-5), fR0(5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.af(1e-5), aR0(5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.fa(1e-5), fR0(5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.af(1e-4 + 1e-5), aR0(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.fa(1e-4 + 1e-5), fR0(1e-4 + 5e-5, 1e-4)) );
}

TEST_F(ExactSurvivorTest, second_interval) {
  std::cout.precision(15);
  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  auto aG0 = [&](t_real _t, t_real _tau) { return N0(_t, _tau, true); };
  auto aG1 = [&](t_real _t, t_real _tau) -> t_rmatrix {
    return N0(_t, _tau, true) - N1(_t - _tau, _tau, true); 
  };
  auto fG1 = [&](t_real _t, t_real _tau) -> t_rmatrix {
    return N0(_t, _tau, false) - N1(_t - _tau, _tau, false);
  };
  auto compare = [](t_rmatrix const &_a, t_rmatrix const & _b) {
    return ((_a - _b).array().abs() < 1e-10).all();
  };
  EXPECT_TRUE( compare(survivor.af(1e-4 + 1e-5), aG1(1e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-4 + 1e-5), fG1(1e-4 + 1e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(1e-4 + 3e-5), aG1(1e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-4 + 3e-5), fG1(1e-4 + 3e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.af(1e-4 + 5e-5), aG1(1e-4 + 5e-5, 1e-4)) );
  EXPECT_TRUE( compare(survivor.fa(1e-4 + 5e-5), fG1(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.af(1e-4 + 1e-5), aG1(1e-4 + 5e-5, 1e-4)) );
  EXPECT_FALSE( compare(survivor.fa(1e-4 + 1e-5), fG1(1e-4 + 5e-5, 1e-4)) );
}

// Checks that Di matrices throw out-of-range
TEST_F(ExactSurvivorTest, out_of_range_Di) {

  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);
  EXPECT_THROW(survivor.D_af(-1), errors::Index);
  EXPECT_THROW(survivor.D_af(5), errors::Index);
  EXPECT_THROW(survivor.D_fa(-1), errors::Index);
  EXPECT_THROW(survivor.D_fa(5), errors::Index);
}

// Checks that Di matrices are what we expect them to be. 
TEST_F(ExactSurvivorTest, regression_Di_matrices) {

  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  t_rmatrix D0(2, 2); D0 << -2.766592241217363e-04, 1.337366986151818e+00, 1.199995862780172e-02,
                            -5.800763938001960e+01;
  EXPECT_TRUE( ((survivor.D_af(0) - D0).array().abs() < 1e-8).all() );

  t_rmatrix D1(2, 2); D1 << -3.400371754358680e+01, -9.964228726415260e+01, -1.285444415722078e-02,
                            -3.766782898643305e-02;
  EXPECT_TRUE( ((survivor.D_af(1) - D1).array().abs() < 1e-8).all() );

  t_rmatrix D2(2, 2); D2 << 33.906985979306185, 95.85937122179871, -0.709318227330929,
                            -2.005333039911058;
  EXPECT_TRUE( ((survivor.D_af(2) - D2).array().abs() < 1e-8).all() );

  t_rmatrix D3(2, 2); D3 << 2.277817593601981e-02, 2.139946009230033e+00, 6.359426652916392e-01,
                            5.974503720194507e+01;
  EXPECT_TRUE( ((survivor.D_af(3) - D3).array().abs() < 1e-8).all() );

  t_rmatrix D4(2, 2); D4 << 0.074230047568711, 0.305603046972027, 0.074230047568708,
                            0.305603046972018;
  EXPECT_TRUE( ((survivor.D_af(4) - D4).array().abs() < 1e-8).all() );
}

// Checks that Ci00 matrices are what we expect them to be. These are the initial matrices that
// yield the recursion.
TEST_F(ExactSurvivorTest, regression_Ci00_matrices) {

  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  t_rmatrix C000(2, 2); C000 << 1.455330874554811e-07, -4.734319903811114e-04, -6.312426538410201e-06,
                                2.053488119070023e-02;
  EXPECT_TRUE( ((survivor.recursion_af(0, 0, 0) - C000).array().abs() < 1e-8).all() );
  t_rmatrix C100(2, 2); C100 << 9.594482244561667e-01, 2.720255573646091e-02, 
                                3.627007431527878e-04, 1.028339719619780e-05;
  EXPECT_TRUE( ((survivor.recursion_af(1, 0, 0) - C100).array().abs() < 1e-8).all() );
  t_rmatrix C200(2, 2); C200 << 0.040510103210602, -0.063558925473074, -0.000847452339641,
                                0.001329622879932;
  EXPECT_TRUE( ((survivor.recursion_af(2, 0, 0) - C200).array().abs() < 1e-8).all() );
  t_rmatrix C300(2, 2); C300 << 1.669965911328917e-05, 3.496776614970140e-02, 4.662368819960321e-04,
                                9.762631769548786e-01;
  EXPECT_TRUE( ((survivor.recursion_af(3, 0, 0) - C300).array().abs() < 1e-8).all() );
  t_rmatrix C400(2, 2); C400 << 2.482714103057449e-05, 1.862035577292865e-03, 2.482714103057373e-05,
                                1.862035577292807e-03;
  EXPECT_TRUE( ((survivor.recursion_af(4, 0, 0) - C400).array().abs() < 1e-8).all() );
}

// Tests one bit of the regression, eg m=1 using m=0 
TEST_F(ExactSurvivorTest, regression_Ci10_matrices) {

  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  t_rmatrix C010(2, 2); C010 << 1.543216801935930e-08, -5.832161141244571e-05,
                                -7.776214854987692e-07, 2.881846635882015e-03;
  EXPECT_TRUE( ((survivor.recursion_af(0, 1, 0) - C010).array().abs() < 1e-8).all() );
  t_rmatrix C110(2, 2); C110 << -2.920314067555906e-02, 3.011280303325819e-02,
                                4.015040404434260e-04, 2.308011939325624e-05;
  EXPECT_TRUE( ((survivor.recursion_af(1, 1, 0) - C110).array().abs() < 1e-8).all() );
  t_rmatrix C210(2, 2); C210 << 0.029163203160271, -0.046601780057635, -0.000621357067435,
                                0.000992579868274;
  EXPECT_TRUE( ((survivor.recursion_af(2, 1, 0) - C210).array().abs() < 1e-8).all() );
  t_rmatrix C310(2, 2); C310 << 1.304381195661902e-05, 1.351401119241333e-02, 1.801868158988379e-04,
                               -7.948211173074680e-03;
  EXPECT_TRUE( ((survivor.recursion_af(3, 1, 0) - C310).array().abs() < 1e-8).all() );
  t_rmatrix C410(2, 2); C410 << 2.687827116302181e-05, 3.033287443375991e-03, 4.044383257834748e-05,
                                4.050704549525265e-03;
  EXPECT_TRUE( ((survivor.recursion_af(4, 1, 0) - C410).array().abs() < 1e-8).all() );
}


// Tests other bit of the regression
TEST_F(ExactSurvivorTest, regression_Ci11_matrices) {

  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);

  t_rmatrix C011(2, 2); C011 << -8.442071118049461e-06,  2.746280314831956e-02,
                                 3.661707086439929e-04, -1.191185663985985e+00;
  EXPECT_TRUE( ((survivor.recursion_af(0, 1, 1) - C011).array().abs() < 1e-8).all() );
  t_rmatrix C111(2, 2); C111 << -3.266094675374351e+01, -9.260126829437684e-01,
                                -1.234683577258275e-02, -3.500610868950084e-04;
  EXPECT_TRUE( ((survivor.recursion_af(1, 1, 1) - C111).array().abs() < 1e-8).all() );
  t_rmatrix C211(2, 2); C211 <<  1.292339253163696, -2.027634781642889,
                                -0.027035130421905,  0.04241716765587;
  EXPECT_TRUE( ((survivor.recursion_af(2, 1, 1) - C211).array().abs() < 1e-8).all() );
  t_rmatrix C311(2, 2); C311 <<  9.981021427566167e-04, 2.089946991412273e+00,
                                 2.786595988549777e-02, 5.834911732046284e+01;
  EXPECT_TRUE( ((survivor.recursion_af(3, 1, 1) - C311).array().abs() < 1e-8).all() );
  t_rmatrix C411(2, 2); C411 <<  9.430169806242189e-06, 7.072627354680796e-04,
                                 9.430169806241902e-06, 7.072627354680582e-04;
  EXPECT_TRUE( ((survivor.recursion_af(4, 1, 1) - C411).array().abs() < 1e-8).all() );
}

TEST_F(ExactSurvivorTest, out_of_range_Ciml_matrices) {
  StateMatrix transitions(Q, 2);
  ExactSurvivor survivor(transitions, 1e-4);
  
  EXPECT_THROW(survivor.recursion_af(-1, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_af(5, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(-1, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(5, 0, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, 0, -1), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, 0, -1), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, -1, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, -1, 0), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, 0, 1), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, 0, 1), errors::Index);
  EXPECT_THROW(survivor.recursion_af(0, 1, 2), errors::Index);
  EXPECT_THROW(survivor.recursion_fa(0, 1, 2), errors::Index);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

