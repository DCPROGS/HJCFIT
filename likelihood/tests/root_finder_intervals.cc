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

#include "../root_finder.h"
#include "random_matrix.h"
using namespace DCProgs;

#ifdef HAS_CXX11_TYPE_TRAITS
  // Checks some assumption about RootInterval
  static_assert( std::is_move_constructible<RootInterval>::value,
  	             "RootInterval is not move constructible." );  
  static_assert( std::is_move_assignable<RootInterval>::value, 
        	       "RootInterval is not move assignable." );  
  static_assert( std::is_move_constructible<Root>::value,
  	             "Root is not move constructible." );  
  static_assert( std::is_move_assignable<Root>::value, 
        	       "Root is not move assignable." );  
#endif
static_assert(std::is_standard_layout<RootInterval>::value, "RootInterval is not a standard layout." );
static_assert(std::is_standard_layout<Root>::value, "Root is not a standard layout." );

#ifdef HAS_CXX11_TRIVIAL_TYPE_TRAITS
// Another MSWindows Fail(TM)
# ifndef MSVC
  static_assert( std::is_trivially_move_constructible<RootInterval>::value,
           "RootInterval is trivially move constructible." );  
  static_assert( std::is_trivially_move_assignable<RootInterval>::value, 
           "RootInterval is trivially move assignable." );  
  static_assert( std::is_trivially_move_constructible<Root>::value,
           "Root is trivially move constructible." );  
  static_assert( std::is_trivially_move_assignable<Root>::value, 
           "Root is trivially move assignable." );  
# endif
#endif

// Sets up test with parameters from CH82, 1e-7 nM.
class RootFinderIntervalsTest : public ::testing::Test {
  
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

std::ostream &operator<< (std::ostream &_stream, DCProgs::RootInterval const &_interval) {
  return _stream << "@(" << _interval.start << ", "
                 << _interval.end << ") = " << _interval.multiplicity;
}

// Checks can find zero for cannonical matrix, open states.
TEST_F(RootFinderIntervalsTest, open) {
  QMatrix const qmatrix(Q, 2);
  DeterminantEq det(qmatrix, 1e-4); 
  
  std::vector<RootInterval> intervals = find_root_intervals(det);
  EXPECT_EQ(std::vector<RootInterval>::size_type(2), intervals.size());

  std::sort(intervals.begin(), intervals.end(),
            [](RootInterval const &_a, RootInterval const &_b) { return _a.start < _b.end; });
  EXPECT_TRUE(intervals.front().end <= intervals.back().start 
              or std::abs(intervals.front().end - intervals.back().start) < 1e-12 );
  EXPECT_TRUE(intervals.front().start < -3045.285776);
  EXPECT_TRUE(intervals.front().end   > -3045.285776);
  EXPECT_TRUE(intervals.back().start < -162.929465);
  EXPECT_TRUE(intervals.back().end   > -162.929465);
}

// Checks can find zero for cannonical matrix, closed qmatrix.
TEST_F(RootFinderIntervalsTest, closed) {
  QMatrix const qmatrix(Q, 2);
  DeterminantEq det(qmatrix.transpose(), 1e-4); 
  
  std::vector<RootInterval> intervals = find_root_intervals(det);
  EXPECT_EQ(std::vector<RootInterval>::size_type(3), intervals.size());

  std::sort(intervals.begin(), intervals.end(),
            [](RootInterval const &_a, RootInterval const &_b) { return _a.start < _b.end; });
  EXPECT_TRUE(intervals[0].end <= intervals[1].start 
              or std::abs(intervals[0].end - intervals[1].start) < 1e-12 );
  EXPECT_TRUE(intervals[1].end <= intervals[2].start 
              or std::abs(intervals[1].end - intervals[2].start) < 1e-12 );
  EXPECT_TRUE(intervals[0].start < -17090.1927692368);
  EXPECT_TRUE(intervals[0].end   > -17090.1927692368);
  EXPECT_TRUE(intervals[1].start < -2058.08129216735);
  EXPECT_TRUE(intervals[1].end   > -2058.08129216735);
  EXPECT_TRUE(intervals[2].start < -0.243565355);
  EXPECT_TRUE(intervals[2].end   > -0.243565355);
}

class TestFindIntervals : public ::testing::TestWithParam<t_int> {
  public:
  TestFindIntervals() {};
};

TEST_P(TestFindIntervals, random_matrix) {

  typedef std::uniform_int_distribution<t_uint> t_idist;
  
  QMatrix qmatrix;
  try {

    // This is the meat of the test.
    // The rest if pfaff to catch complex eigenvalues and other errors.
    // First create an appropriate Q matrix.
    qmatrix.matrix = nonsingular_rate_matrix(5, 8);
    qmatrix.nopen = t_idist(2, qmatrix.matrix.rows()-2)(global_mersenne());

    // Then the determinant object
    DeterminantEq det(qmatrix, 1e-4);
    // Look for roots and sort them
    t_real const convergence = 1e-6;
    std::vector<RootInterval> intervals = find_root_intervals(det, quiet_nan, quiet_nan, convergence);
    std::sort(intervals.begin(), intervals.end(),
              [](RootInterval const &_a, RootInterval const &_b)
              { return _a.start < _b.end; });

    t_uint nroots = 0;
    // A few tests:
    //   - multiplicity should never be zero. What would be the point of an interval without
    //     roots ?
    //   - interval should have start < end, of course. 
    //   - if multiplicity is larger than 1, then interval should be smaller than convergence.
    //   - if multiplicity is _even_, then the sign of det(W(s)) _should not_ change between start and
    //   end of interval.
    //     If is it _odd_, then it _should_ change between start and end.
    //   - The multiplicity should add up to the size of aa matrix.
    for(RootInterval const &interval: intervals) {


      nroots += interval.multiplicity;
      EXPECT_TRUE(interval.multiplicity != 0);
      EXPECT_TRUE(interval.end > interval.start);
      if(interval.multiplicity > 1) {
        EXPECT_TRUE(interval.end - interval.start < convergence );
      } 

      t_int const start_sign = det(interval.start) > 0 ? 1: -1;
      t_int const end_sign = det(interval.end) > 0 ? 1: -1;
      EXPECT_TRUE(interval.multiplicity % 2 == 0?
                      start_sign == end_sign:
                      start_sign != end_sign );
    }
    EXPECT_EQ(static_cast<t_uint>(qmatrix.aa().rows()), det.get_nopen()) << qmatrix;
    EXPECT_TRUE(nroots > 0) << qmatrix;
  } catch(...) {
    std::cerr.precision(15);
    std::cerr << "Error for " << qmatrix << std::endl;
    throw;
  }
}

INSTANTIATE_TEST_CASE_P(random, TestFindIntervals, ::testing::Range(t_int(0), t_int(100)));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

