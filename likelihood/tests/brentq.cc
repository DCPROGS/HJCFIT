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

#include <gtest/gtest.h>
#include "../brentq.h"
#include "random_matrix.h"

using namespace DCProgs;

class TestPolynomial : public ::testing::TestWithParam<int> {
  
  public:
  virtual void SetUp() {
    t_rvector randoms = DCProgs::random_vector();
    roots.resize(randoms.size());
    for(t_rvector::Index i(0); i < randoms.size(); ++i) roots(i) = randoms.head(i+1).sum();
    intervals.resize(roots.size()+1);
    intervals(0) = 0e0;
    for(t_rvector::Index i(0); i < randoms.size() - 1; ++i)
      intervals(i+1) = (roots(i) + roots(i+1)) * 0.5;
    intervals(roots.size()) = roots(roots.size()-1) + 1;
  }


  protected:
    t_rvector roots;
    t_rvector intervals;
};


// Checks the size of the matrices is correct.
TEST_P(TestPolynomial, random) {

  SetUp();
  t_rvector const roots = this->roots;
  auto lambda = [&roots](t_rvector::Scalar const &_t) -> t_rvector::Scalar {
    t_rvector::Scalar result(1);
    for(t_rvector::Index i(0); i < roots.size(); ++i) result *= (roots(i) - _t);
    return result;
  };
  EXPECT_TRUE(roots.size() > 0);
  for(t_rvector::Index i(0); i < roots.size(); ++i) {
    auto result = brentq(lambda, intervals(i), intervals(i+1), 1e-8, 1e-8);
    EXPECT_TRUE(std::abs(std::get<0>(result) - roots(i)) < 1e-7);
  }
}

INSTANTIATE_TEST_CASE_P(BrentQ, TestPolynomial, ::testing::Range(0, 100));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

