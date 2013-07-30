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

  cxx_test(qmatrix qmatrix.cc likelihood)
  cxx_test(idealG idealG.cc likelihood)
  cxx_test(idealG_exponentials idealg_exponentials.cc likelihood)
  cxx_test(idealG_laplacian idealG_laplacian.cc likelihood)
  cxx_test(occupancies occupancies.cc likelihood)
  cxx_test(time_filter time_filter.cc likelihood)
  cxx_test(laplace_survivor laplace_survivor.cc likelihood)
  cxx_test(determinant_equation determinant_equation.cc likelihood)
  cxx_test(recursion recursion_formula.cc likelihood)
  cxx_test(exact_survivor exact_survivor.cc likelihood)
  cxx_test(root_finder_intervals
           "root_finder_intervals.cc;random_matrix.h;random_matrix.cc" 
           likelihood)
  cxx_test(root_finder_lower_bound
           "root_finder_lower_bound.cc;random_matrix.h;random_matrix.cc"
           likelihood)
  cxx_test(asymptotes asymptotes.cc likelihood)
  cxx_test(brentq
           "brentq.cc;random_matrix.h;random_matrix.cc"
           likelihood)
  cxx_test(likelihood likelihood.cc likelihood)

endif(tests)
