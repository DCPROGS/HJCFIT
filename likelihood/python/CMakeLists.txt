if(pythonBindings)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_SWIG_FLAGS "")

  set_source_files_properties(likelihood.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall")
  set(SWIG_MODULE_likelihood_EXTRA_DEPS asymptotes.swg)
  swig_add_module(likelihood python likelihood.i asymptotes.swg)
  set_target_properties(_likelihood PROPERTIES LINKER_LANGUAGE CXX)
  swig_link_libraries(likelihood ${PYTHON_LIBRARIES} likelihood)
  install(TARGETS _likelihood DESTINATION ${CMAKE_PYINSTALL_PREFIX})
  install(FILES  __init__.py ${CMAKE_CURRENT_BINARY_DIR}/likelihood.py
          DESTINATION ${CMAKE_PYINSTALL_PREFIX})

endif(pythonBindings)

if(pythonBindings AND tests)

  feature_test(state_matrix features/state_matrix.feature)
  feature_test(asymptotes features/asymptotes.feature)
  
endif(pythonBindings AND tests)
