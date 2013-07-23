if(pythonBindings)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_SWIG_FLAGS "")

  set_source_files_properties(likelihood.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall")
  set(SWIG_MODULE_likelihood_EXTRA_DEPS determinant_equation.swg qmatrix.swg root_finder.swg
                                        asymptotes.swg exact_survivor.swg approx_survivor.swg
                                        missed_eventsG.swg idealg.swg chained.swg time_filter.swg)
  swig_add_module(likelihood python likelihood.i)
  set_target_properties(_likelihood PROPERTIES LINKER_LANGUAGE CXX)
  swig_link_libraries(likelihood ${PYTHON_LIBRARIES} likelihood)

  
  install(TARGETS _likelihood DESTINATION ${PYINSTALL_DIRECTORY}/dcprogs/likelihood)
  install(FILES  random.py _likelihood_methods.py optimization.py
                 ${CMAKE_CURRENT_BINARY_DIR}/likelihood.py __init__.py
          DESTINATION ${PYINSTALL_DIRECTORY}/dcprogs/likelihood)
endif(pythonBindings)

if(pythonBindings AND tests)

  feature_test(qmatrix features/qmatrix.feature)
  feature_test(determinant_equation features/determinant_equation.feature)
  feature_test(find_roots features/find_roots.feature)
  feature_test(asymptotes features/asymptotes.feature)
  feature_test(idealg features/idealg.feature)
  feature_test(missed_eventsG features/missed_eventsG.feature)
  feature_test(approx_survivor features/approx_survivor.feature)
  
endif(pythonBindings AND tests)
