if(pythonBindings)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_SWIG_FLAGS "")

  set_source_files_properties(likelihood.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall")
  swig_add_module(likelihood python likelihood.i)
  set_target_properties(_likelihood PROPERTIES LINKER_LANGUAGE CXX)
  swig_link_libraries(likelihood ${PYTHON_LIBRARIES} likelihood)
  install(TARGETS _likelihood DESTINATION ${PYINSTALL_DIR}/dcprogs)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/likelihood.py DESTINATION ${PYINSTALL_DIR}/dcprogs)
  install(FILES __init__.py DESTINATION ${PYINSTALL_DIR}/dcprogs)

endif(pythonBindings)
