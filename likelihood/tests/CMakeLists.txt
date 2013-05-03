if(tests)
  cxx_test(idealG idealG.cc likelihood)
  cxx_test(state_matrix state_matrix.cc likelihood)
endif(tests)
