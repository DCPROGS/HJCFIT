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

// The log-likelihood must be a function of the Q-matrix and the data alone.
//
// It is fed to a Simplex search, so any dependence on how many threads happened
// to be available steers the optimiser and makes fitted rate constants
// irreproducible between machines. Two people fitting the same mechanism to the
// same data would get different answers with nothing to warn them.
//
// The existing coherence check in likelihood.cc, TestLikelihood.vector_vs_real,
// looks like it covers this but does not: its bursts are one and three
// intervals long, so neither the high-level path (more than 100 bursts) nor the
// low-level one (more than 100 intervals in a burst) is ever entered. These
// cases are sized to enter both.

#include <HJCFITConfig.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#if defined(_OPENMP)
# include <omp.h>
#endif
#include "../likelihood.h"

using namespace HJCFIT;

namespace {

  QMatrix classic_qmatrix() {
    QMatrix qmatrix;
    qmatrix.matrix.resize(5, 5);
    qmatrix.matrix <<  -3050,        50,  3000,      0,    0,
                     2./3., -1502./3.,     0,    500,    0,
                        15,         0, -2065,     50, 2000,
                         0,     15000,  4000, -19000,    0,
                         0,         0,    10,      0,  -10;
    qmatrix.nopen = 2;
    return qmatrix;
  }

  //! Bursts with an odd number of intervals, as the likelihood requires.
  //! Interval lengths vary so the running product actually moves through the
  //! rescaling thresholds rather than sitting at one magnitude.
  t_Bursts make_bursts(std::size_t nbursts, std::size_t intervals_per_burst) {
    t_Bursts bursts;
    bursts.reserve(nbursts);
    for(std::size_t b = 0; b < nbursts; ++b) {
      t_Burst burst;
      burst.reserve(intervals_per_burst);
      for(std::size_t i = 0; i < intervals_per_burst; ++i)
        burst.push_back(2.0e-4 + 1.0e-4 * static_cast<t_real>((i + b) % 7));
      bursts.push_back(burst);
    }
    return bursts;
  }

  //! Log10Likelihood reads the thread count in its constructor, so it has to be
  //! rebuilt after omp_set_num_threads for the change to take effect.
  t_real evaluate(t_Bursts const &bursts, QMatrix const &qmatrix, int nthreads) {
#   if defined(_OPENMP)
      omp_set_num_threads(nthreads);
#   endif
    Log10Likelihood likelihood(bursts, qmatrix.nopen, 1e-4);
    return likelihood(qmatrix);
  }

  //! Distance in representable doubles between two values, so a report says how
  //! far apart the results were rather than only that they differed.
  long ulp_distance(t_real a, t_real b) {
    if(a == b) return 0;
    long n = 0;
    t_real lo = a < b ? a : b, hi = a < b ? b : a;
    while(lo < hi and n < 1000000) { lo = std::nextafter(lo, hi); ++n; }
    return n;
  }

  std::vector<int> thread_counts() {
#   if defined(_OPENMP)
      std::vector<int> const wanted{1, 2, 4, 8};
      std::vector<int> usable;
      int const cap = omp_get_max_threads();
      for(int n : wanted) if(n <= cap or n == 1) usable.push_back(n);
      return usable;
#   else
      return std::vector<int>{1};
#   endif
  }
}

//! Many short bursts: parallelised across bursts.
TEST(Determinism, across_threads_many_bursts) {
  QMatrix const qmatrix = classic_qmatrix();
  t_Bursts const bursts = make_bursts(150, 11);

  std::vector<int> const counts = thread_counts();
  t_real const reference = evaluate(bursts, qmatrix, counts.front());
  for(int n : counts) {
    t_real const value = evaluate(bursts, qmatrix, n);
    EXPECT_EQ(value, reference)
      << "log-likelihood changed with " << n << " threads (reference used "
      << counts.front() << "): " << std::setprecision(20) << value << " vs "
      << reference << ", " << ulp_distance(value, reference) << " ulp";
  }
}

//! Few long bursts: parallelised within each burst, which is the path that
//! splits the matrix product across threads.
TEST(Determinism, across_threads_long_bursts) {
  QMatrix const qmatrix = classic_qmatrix();
  t_Bursts const bursts = make_bursts(4, 501);

  std::vector<int> const counts = thread_counts();
  t_real const reference = evaluate(bursts, qmatrix, counts.front());
  for(int n : counts) {
    t_real const value = evaluate(bursts, qmatrix, n);
    EXPECT_EQ(value, reference)
      << "log-likelihood changed with " << n << " threads (reference used "
      << counts.front() << "): " << std::setprecision(20) << value << " vs "
      << reference << ", " << ulp_distance(value, reference) << " ulp";
  }
}

//! operator() is documented as the sum over bursts of what vector() returns.
//! Exercised here on burst counts that straddle the 100-burst switch between
//! the two code paths, which the existing coherence test never reaches.
TEST(Determinism, sum_of_vector_equals_scalar) {
  QMatrix const qmatrix = classic_qmatrix();

  for(std::size_t nbursts : {4u, 150u}) {
    for(std::size_t nintervals : {11u, 501u}) {
      t_Bursts const bursts = make_bursts(nbursts, nintervals);
      Log10Likelihood likelihood(bursts, qmatrix.nopen, 1e-4);

      t_rvector const per_burst = likelihood.vector(qmatrix);
      ASSERT_EQ(static_cast<std::size_t>(per_burst.size()), nbursts);

      t_real sum(0);
      for(t_int i = 0; i < per_burst.size(); ++i) sum += per_burst(i);

      t_real const scalar = likelihood(qmatrix);
      EXPECT_EQ(sum, scalar)
        << "sum(vector(Q)) != operator()(Q) for " << nbursts << " bursts of "
        << nintervals << " intervals: " << std::setprecision(20) << sum
        << " vs " << scalar << ", " << ulp_distance(sum, scalar) << " ulp";
    }
  }
}

int main(int argc, char **argv) {
#if defined(_OPENMP)
  std::cout << "[determinism] OpenMP enabled, omp_get_max_threads() = "
            << omp_get_max_threads() << std::endl;
#else
  std::cout << "[determinism] built without OpenMP: this test cannot "
               "demonstrate thread dependence" << std::endl;
#endif
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
