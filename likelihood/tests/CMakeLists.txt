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

endif(tests)
