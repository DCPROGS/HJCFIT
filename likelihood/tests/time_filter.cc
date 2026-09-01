/***********************
    HJCFIT computes missed-events likelihood as described in
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

#include "HJCFITConfig.h"
#include <random>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include <time.h>
#include <cstdlib>
#include "../time_filter.h"
using namespace HJCFIT;

t_uint  const Nmax[2] = {6, 100};
t_real const tau = 1;
t_real const alpha = 10;

class TestTimeFilter : public ::testing::TestWithParam<t_int> { 
  public:
    //! Seeded from the environment, and fixed by default.
    //!
    //! This used to seed from std::random_device (or the clock), so every run
    //! drew different series and a failure could not be reproduced from the
    //! log -- the reported parameter index names the test case, not the data.
    //! It failed roughly one CI run in fifteen, on a different parameter each
    //! time, which made every red run ambiguous.
    //!
    //! Set HJCFIT_TEST_SEED to sweep other seeds deliberately; the series is
    //! printed on failure, so a known seed makes the failing case recoverable.
    TestTimeFilter() {
      char const * const env = std::getenv("HJCFIT_TEST_SEED");
      mersenne.seed(env ? static_cast<unsigned int>(std::strtoul(env, nullptr, 10)) : 42u);
    }
  protected:
    std::mt19937 mersenne;
};


//! Creates a fake time series with known number of sub-resolution steps.
//! \param[in] _N: Total number of times
//! \param[in] _n: Number of sub-resolution events
//! \param[in] _tau: Intervall below which two subsequent events cannot be detected
//! \param[in] _alpha: _alpha*_tau is the max interval between events
//! \param[in] _rng: random number generator engine.
template<class T>
  t_rvector fake_time_series(t_uint _N, t_uint _n, t_real _tau, t_real _alpha, T && _rng) {
  
    typedef std::uniform_real_distribution<t_real> t_rdist;
    t_rdist __supercrit(_tau*1.01, _tau*_alpha);
    t_rdist __subcrit(_tau*1e-4, _tau*0.999);


    std::vector<t_real> intervals(_N);
    std::generate(intervals.begin(), intervals.end()-_n, [&] { return __supercrit(_rng); });
    std::generate(intervals.end()-_n, intervals.end(), [&] { return __subcrit(_rng); });
    std::shuffle(intervals.begin(), intervals.end(), _rng);

    t_rvector result(_N+1);
    result(0) = 0e0;
    for(t_uint i(0); i < _N; ++i) result[i+1] = result[i] + intervals[i];
    return result;
  }

t_rvector::Index nbfiltered(t_rvector const &_vector, t_real _tau) {
 
  t_rvector intervals = _vector.tail(_vector.size()-1) - _vector.head(_vector.size()-1);
  // Find the first detectable *open* interval, stepping by two because the
  // series alternates open and shut -- exactly as interval_filter_impl does.
  //
  // Both of its early returns have to be modelled here. Only the first was,
  // and the second is reachable: when the first detectable open interval is
  // the last interval in the series, the filter yields nothing, because a
  // lone trailing opening is not a series. That happens in about one draw in
  // thirty thousand, which is why this disagreed with the implementation
  // rarely enough to look like flakiness rather than a defect in the test.
  t_rvector::Index start(0);
  for(; start < intervals.size() and intervals(start) < _tau; start += 2);
  if(start >= intervals.size()) return 0;
  if(start == intervals.size() - 1) return 0;
  t_rvector::Index result = (intervals.array() >= _tau).count();
  for(t_rvector::Index i(0); i < intervals.size(); ++i) 
    if(intervals(i) < _tau) {
      t_uint sub = 0;
      for(; i < intervals.size() and intervals(i) < _tau; ++i, ++sub);
      if(i != intervals.size() and sub % 2 == 1) --result;
    }  
  return result + 1;
}


TEST_P(TestTimeFilter, nbfiltered) {
  typedef std::uniform_int_distribution<t_uint> t_idist;
  t_rvector::Index const n = t_idist(Nmax[0], Nmax[1])(this->mersenne);
  t_rvector::Index const N = t_idist(Nmax[0], Nmax[1])(this->mersenne) + n;
  t_rvector const series = fake_time_series(N, n, tau, alpha, this->mersenne); 
  t_rvector const intervals = series.tail(N) - series.head(N);
  EXPECT_EQ(series.size(), N+1);
  EXPECT_EQ(series(0), 0e0);
  EXPECT_EQ((intervals.array() < tau).count(), n)
    << "Series of " << N << " has " 
    << (intervals.array() < tau).count() 
    << " sub-resolution intervals, rather than "
    << n << "."; 
  t_rvector const filtered = time_filter(series, tau);
  t_uint const nf = filtered.size();
  if(nf > 1) {
    EXPECT_TRUE(((filtered.tail(nf-1) - filtered.head(nf-1)).array() > tau).all());
  }
  EXPECT_EQ(intervals.size(), N);
  EXPECT_EQ(nbfiltered(series, tau), filtered.size()) << series.transpose();
}

INSTANTIATE_TEST_CASE_P(random, TestTimeFilter, ::testing::Range(t_int(0), t_int(300)));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

