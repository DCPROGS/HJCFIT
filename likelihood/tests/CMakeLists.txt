if(tests)
  # on windows + visual studio, have to copy dll to output dir.
# if(MSVC) 
#   ExternalProject_Add_Step(googletest likelitest
#                            COMMAND 
#                              Copy.exe ${BINARY_DIR}/${CMAKE_CFG_INTDIR}/gtest.dll ${EXECUTABLE_OUTPUT_PATH}
#                            DEPENDEES build)
# endif(MSVC)
  cxx_test(state_matrix state_matrix.cc likelihood)
  cxx_test(idealG idealG.cc likelihood)
  cxx_test(idealG_exponentials idealg_exponentials.cc likelihood)
  cxx_test(idealG_laplacian idealG_laplacian.cc likelihood)
  cxx_test(equilibrium equilibrium.cc likelihood)
  cxx_test(time_filter time_filter.cc likelihood)

endif(tests)
