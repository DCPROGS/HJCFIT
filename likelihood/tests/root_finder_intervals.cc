#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>

#include "../root_finder.h"
using namespace DCProgs;

#ifdef HAS_CXX11_TYPETRAITS
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

#ifdef HAS_CXX11_TRIVIALTYPETRAITS
  static_assert( std::is_trivially_move_constructible<RootInterval>::value,
  	       "RootInterval is trivially move constructible." );  
  static_assert( std::is_trivially_move_assignable<RootInterval>::value, 
  	       "RootInterval is trivially move assignable." );  
  static_assert( std::is_trivially_move_constructible<Root>::value,
  	       "Root is trivially move constructible." );  
  static_assert( std::is_trivially_move_assignable<Root>::value, 
  	       "Root is trivially move assignable." );  
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

std::ostream &operator<< (std::ostream &_stream, DCProgs::RootInterval &_interval) {
  return _stream << "@(" << _interval.start << ", "
                 << _interval.end << ") = " << _interval.multiplicity;
}

// Checks can find zero for cannonical matrix, open states.
TEST_F(RootFinderIntervalsTest, open) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 1e-4, true); 
  
  std::vector<RootInterval> intervals = find_root_intervals(det, -1e6);
  EXPECT_EQ(intervals.size(), 2);

  std::sort(intervals.begin(), intervals.end(),
            [](RootInterval const &_a, RootInterval const &_b) { return _a.start < _b.end; });
  EXPECT_TRUE(intervals.front().end <= intervals.back().start 
              or std::abs(intervals.front().end - intervals.back().start) < 1e-12 );
  EXPECT_TRUE(intervals.front().start < -3045.285776);
  EXPECT_TRUE(intervals.front().end   > -3045.285776);
  EXPECT_TRUE(intervals.back().start < -162.929465);
  EXPECT_TRUE(intervals.back().end   > -162.929465);
}

// Checks can find zero for cannonical matrix, closed states.
TEST_F(RootFinderIntervalsTest, closed) {
  StateMatrix const states(Q, 2);
  DeterminantEq det(states, 1e-4, false); 
  
  std::vector<RootInterval> intervals = find_root_intervals(det, -1e6);
  EXPECT_EQ(intervals.size(), 3);

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


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

