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

