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

#include <memory>
#include <gtest/gtest.h>
#include "../likelihood.h"
#include "../missed_eventsG.h"
#include "../occupancies.h"

using namespace DCProgs;

class TestLikelihood : public ::testing::Test {
  
  public:
  virtual void SetUp() {
    qmatrix.matrix.resize(5, 5);
    qmatrix.matrix <<  -3050,        50,  3000,      0,    0,
                     2./3., -1502./3.,     0,    500,    0, 
                        15,         0, -2065,     50, 2000, 
                         0,     15000,  4000, -19000,    0, 
                         0,         0,    10,      0,  -10;
    qmatrix.nopen = 2;
    t_Bursts bursts;
    likelihood.reset(new Log10Likelihood(bursts, qmatrix.nopen, 1e-4));
  }


  protected:
    std::shared_ptr<Log10Likelihood> likelihood;
    QMatrix qmatrix;
};


// Checks that vector and operator() are coherent wrt one another.
TEST_F(TestLikelihood, vector_vs_real) {

  likelihood->bursts = t_Bursts(1, t_Burst(1, 2.5e-4) );
  t_real const result_real = (*likelihood)(qmatrix);
  t_rvector result_vector = likelihood->vector(qmatrix.matrix);

  EXPECT_EQ(result_vector.size(), 1);
  EXPECT_DOUBLE_EQ(result_vector(0), result_real);

  likelihood->bursts.push_back(likelihood->bursts.back());
  result_vector = likelihood->vector(qmatrix.matrix);
  EXPECT_EQ(result_vector.size(), 2);
  EXPECT_DOUBLE_EQ(result_vector(0), result_vector(1));

  likelihood->bursts.back().push_back(5.5e-4);
  likelihood->bursts.back().push_back(5.5e-4);
  result_vector = likelihood->vector(qmatrix.matrix);
  EXPECT_DOUBLE_EQ(result_vector(0), result_real);
  likelihood->bursts.erase(likelihood->bursts.begin());
  EXPECT_DOUBLE_EQ(result_vector(1), (*likelihood)(qmatrix.matrix));
}


// Checks that adding bursts is additive indeed.
TEST_F(TestLikelihood, additive_vs_bursts) {

  likelihood->bursts = t_Bursts(1, t_Burst(1, 1.5e-4) );
  t_real const result = (*likelihood)(qmatrix);

  likelihood->bursts.push_back(likelihood->bursts.back());
  EXPECT_TRUE( std::abs(2e0*result - (*likelihood)(qmatrix)) < 1e-8 * std::abs(result));
  EXPECT_FALSE( std::abs(result - (*likelihood)(qmatrix)) < 1e-8 * std::abs(result));

  likelihood->bursts.push_back(likelihood->bursts.back());
  EXPECT_TRUE( std::abs(3e0*result - (*likelihood)(qmatrix)) < 1e-8 * std::abs(result));
  EXPECT_FALSE( std::abs(result - (*likelihood)(qmatrix)) < 1e-8 * std::abs(result));

  likelihood->bursts.back().back() *= 1.5; 
  EXPECT_FALSE( std::abs(3e0*result - (*likelihood)(qmatrix)) < 1e-8 * std::abs(result));

  likelihood->bursts.back().back() = likelihood->bursts.front().front();
  EXPECT_TRUE( std::abs(3e0*result - (*likelihood)(qmatrix)) < 1e-8 * std::abs(result));
  likelihood->bursts.back().push_back(4.5e-4);
  likelihood->bursts.back().push_back(4.5e-4);
  t_real const allthree = (*likelihood)(qmatrix);
  EXPECT_FALSE( std::abs(3e0*result - allthree) < 1e-8 * std::abs(result));
  likelihood->bursts.erase(likelihood->bursts.begin());
  likelihood->bursts.erase(likelihood->bursts.begin());
  EXPECT_TRUE( std::abs(2e0*result + (*likelihood)(qmatrix) - allthree) < 1e-8 * std::abs(result));
}

// Checks that a manual test gives what we want.
TEST_F(TestLikelihood, manual_check) {

  likelihood->bursts = t_Bursts(1, t_Burst(1, 1.5e-4) );

  likelihood->tcritical = -1e0;
  MissedEventsG eG(qmatrix, 1e-4);
  t_initvec const initial = occupancies(eG);
  t_rvector const final = t_rvector::Ones(3,1);

  t_rmatrix matrix = eG.af( likelihood->bursts.back()[0] );

  t_real check = std::log10(initial * matrix * final);
  EXPECT_TRUE(std::abs((*likelihood)(qmatrix) - check) < 1e-8 * std::abs(check)) 
          << "check: " << check << "\n"
          << "likelihood: " << (*likelihood)(qmatrix) << "\n";

  likelihood->bursts.back().push_back(4.5e-4);
  likelihood->bursts.back().push_back(3.5e-4);
  matrix = eG.af( likelihood->bursts.back()[0] ) 
           * eG.fa( likelihood->bursts.back()[1] )
           * eG.af( likelihood->bursts.back()[2] );
  check = std::log10(initial * matrix * final);
  EXPECT_TRUE(std::abs((*likelihood)(qmatrix) - check) < 1e-8 * std::abs(check))
          << "check: " << check << "\n"
          << "likelihood: " << (*likelihood)(qmatrix) << "\n";
}

// Checks that a manual test gives what we want.
TEST_F(TestLikelihood, manual_check_CHS) {

  likelihood->bursts = t_Bursts(1, t_Burst(1, 1.5e-4) );

  likelihood->tcritical = 1e-3;
  MissedEventsG eG(qmatrix, 1e-4);
  t_initvec const initial_eq = occupancies(eG);
  t_rvector const final_eq = occupancies(eG, false).transpose();
  t_initvec const initial = CHS_occupancies(eG, likelihood->tcritical);
  t_rvector const final = CHS_occupancies(eG, likelihood->tcritical, false).transpose();

  t_rmatrix matrix = eG.af( likelihood->bursts.back()[0] );

  t_real check = std::log10(initial * matrix * final);
  EXPECT_TRUE(std::abs((*likelihood)(qmatrix) - check) < 1e-8 * std::abs(check)) 
          << "check: " << check << "\n"
          << "likelihood: " << (*likelihood)(qmatrix) << "\n";
  EXPECT_FALSE(std::abs(check - std::log10(initial_eq * matrix * final_eq))
                  < 1e-8 * std::abs(check));

  likelihood->bursts.back().push_back(4.5e-4);
  likelihood->bursts.back().push_back(3.5e-4);
  matrix = eG.af( likelihood->bursts.back()[0] ) 
           * eG.fa( likelihood->bursts.back()[1] )
           * eG.af( likelihood->bursts.back()[2] );
  check = std::log10(initial * matrix * final);
  EXPECT_TRUE(std::abs((*likelihood)(qmatrix) - check) < 1e-8 * std::abs(check))
          << "check: " << check << "\n"
          << "likelihood: " << (*likelihood)(qmatrix) << "\n";
  EXPECT_FALSE(std::abs(check - std::log10(initial_eq * matrix * final_eq))
                  < 1e-8 * std::abs(check));
}

// Checks throws on odd number of intervals
TEST_F(TestLikelihood, odd_intervals) {

  likelihood->bursts = t_Bursts(1, t_Burst(1, 1.5e-4) );
  EXPECT_NO_THROW((*likelihood)(qmatrix));
  likelihood->bursts.back().push_back(3.5e-4);
  EXPECT_THROW((*likelihood)(qmatrix), errors::Domain);
}

// Makes sure that code throws if qmatrix is too large to fit the stack
TEST_F(TestLikelihood, exceeds_stack_throws) {
  qmatrix.matrix.resize(dcprogs_stack_matrix+1,dcprogs_stack_matrix+1);
  likelihood->bursts = t_Bursts(1, t_Burst(1, 1.5e-4) );
  EXPECT_THROW((*likelihood)(qmatrix), errors::Domain);
}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
