if(tests)

  cxx_test(state_matrix state_matrix.cc likelihood)
  cxx_test(idealG idealG.cc likelihood)
  cxx_test(idealG_exponentials idealg_exponentials.cc likelihood)
  cxx_test(idealG_laplacian idealG_laplacian.cc likelihood)
  cxx_test(equilibrium equilibrium.cc likelihood)
  cxx_test(time_filter time_filter.cc likelihood)
  cxx_test(asymptotes asymptotes.cc likelihood)

endif(tests)
