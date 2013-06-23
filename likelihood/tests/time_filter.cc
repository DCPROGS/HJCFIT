#include <random>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include <time.h>
#include "../time_filter.h"
using namespace DCProgs;

t_int  const Nmax[2] = {6, 100};
t_real const tau = 1;
t_real const alpha = 10;

class TestTimeFilter : public ::testing::TestWithParam<t_int> { 
  public:
    TestTimeFilter() {
#   ifdef HAS_CXX11_RANDOM_DEVICE
      std::random_device rd;
      mersenne.seed(rd()); 
#   else 
      mersenne.seed(static_cast<unsigned int>(std::time(nullptr))); 
#   endif
    }
  protected:
    std::mt19937 mersenne;
};


//! Creates a fake time series with known number of critical steps.
//! \param[in] _N: Total number of times
//! \param[in] _n: Number of critical events
//! \param[in] _tau: Intervall below which two subsequent events cannot be detected
//! \param[in] _alpha: _alpha*_tau is the max interval between events
//! \param[in] _rng: random number generator engine.
template<class T>
  t_rvector fake_time_series(t_int _N, t_int _n, t_real _tau, t_real _alpha, T && _rng) {
  
    typedef std::uniform_real_distribution<t_real> t_rdist;
    t_rdist __supercrit(_tau*1.01, _tau*_alpha);
    t_rdist __subcrit(_tau*1e-4, _tau*0.999);


    std::vector<t_real> intervals(_N);
    std::generate(intervals.begin(), intervals.end()-_n, [&] { return __supercrit(_rng); });
    std::generate(intervals.end()-_n, intervals.end(), [&] { return __subcrit(_rng); });
    std::shuffle(intervals.begin(), intervals.end(), _rng);

    t_rvector result(_N+1);
    result(0) = 0e0;
    for(t_int i(0); i < _N; ++i) result[i+1] = result[i] + intervals[i];
    return result;
  }

t_int nbfiltered(t_rvector const &_vector, t_real _tau) {
 
  auto intervals = (_vector.tail(_vector.size()-1) - _vector.head(_vector.size()-1)).array();
  // Check special case where time series disappears.
  if(intervals(0) < _tau) { 
    t_int i(2);
    for(; i < intervals.size() and intervals(i) < _tau; i += 2);
    if(i >= intervals.size()) return 0;
  }
  t_int i = 0;
  t_int result = (intervals >= _tau).count();
  for(t_int i(0); i < intervals.size(); ++i) 
    if(intervals(i) < _tau) {
      t_int sub = 0;
      for(; i < intervals.size() and intervals(i) < _tau; ++i, ++sub);
      if(i != intervals.size() and sub % 2 == 1) --result;
    }  
  return result + 1;
}


TEST_P(TestTimeFilter, nbfiltered) {
  typedef std::uniform_int_distribution<t_int> t_idist;
  t_int const n = t_idist(Nmax[0], Nmax[1])(this->mersenne);
  t_int const N = t_idist(Nmax[0], Nmax[1])(this->mersenne) + n;
  t_rvector const series = fake_time_series(N, n, tau, alpha, this->mersenne); 
  auto const intervals = (series.tail(N) - series.head(N)).array();
// EXPECT_EQ(series.size(), N+1);
// EXPECT_EQ(series(0), 0e0);
// EXPECT_EQ((intervals < tau).count(), n)
//   << "Series of " << N << " has " 
//   << (intervals < tau).count() 
//   << " sub-critical intervals, rather than "
//   << n << "."; 
// t_rvector const filtered = time_filter(series, tau);
// t_int const nf = filtered.size();
// EXPECT_TRUE(((filtered.tail(nf-1) - filtered.head(nf-1)).array() > tau).all());
// EXPECT_EQ(intervals.size(), N);
// EXPECT_EQ(nbfiltered(series, tau), filtered.size()) << series.transpose();
}

INSTANTIATE_TEST_CASE_P(random, TestTimeFilter, ::testing::Range(0, 300));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

