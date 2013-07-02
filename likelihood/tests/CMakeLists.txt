if(tests)

  cxx_test(state_matrix state_matrix.cc likelihood)
  cxx_test(idealG idealG.cc likelihood)
  cxx_test(idealG_exponentials idealg_exponentials.cc likelihood)
  cxx_test(idealG_laplacian idealG_laplacian.cc likelihood)
  cxx_test(equilibrium equilibrium.cc likelihood)
  cxx_test(time_filter time_filter.cc likelihood)
  cxx_test(determinant_equation determinant_equation.cc likelihood)
  cxx_test(recursion recursion_formula.cc likelihood)
  cxx_test(exactG exactG.cc likelihood)
  cxx_test(root_finder_intervals
           "root_finder_intervals.cc;random_matrix.h;random_matrix.cc" 
           likelihood)
  cxx_test(root_finder_lower_bound
           "root_finder_lower_bound.cc;random_matrix.h;random_matrix.cc"
           likelihood)

endif(tests)
