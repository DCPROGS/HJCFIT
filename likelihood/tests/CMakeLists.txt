########################
#   DCProgs computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

if(tests)

  add_gtest(qmatrix qmatrix.cc likelihood)
  add_gtest(idealG idealG.cc likelihood)
  add_gtest(idealG_exponentials idealg_exponentials.cc likelihood)
  add_gtest(idealG_laplacian idealG_laplacian.cc likelihood)
  add_gtest(occupancies occupancies.cc likelihood)
  add_gtest(time_filter time_filter.cc likelihood)
  add_gtest(laplace_survivor laplace_survivor.cc likelihood)
  add_gtest(determinant_equation determinant_equation.cc likelihood)
  add_gtest(recursion recursion_formula.cc likelihood)
  add_gtest(exact_survivor exact_survivor.cc likelihood)
  if(DCPROGS_USE_MPFR AND NOT DCPROGS_LONG_DOUBLE)
    add_gtest(mpfr_root_bounds mpfr_root_bounds.cc likelihood)
  endif(DCPROGS_USE_MPFR AND NOT DCPROGS_LONG_DOUBLE)
  add_gtest(root_finder_intervals
           "root_finder_intervals.cc;random_matrix.h;random_matrix.cc" 
           likelihood)
  add_gtest(root_finder_lower_bound
           "root_finder_lower_bound.cc;random_matrix.h;random_matrix.cc"
           likelihood)
  add_gtest(asymptotes asymptotes.cc likelihood)
  add_gtest(brentq
           "brentq.cc;random_matrix.h;random_matrix.cc"
           likelihood)
  add_gtest(likelihood likelihood.cc likelihood)

endif(tests)
